# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


from memoized_property import memoized_property
from serializable import Serializable

from .common import bio_seq_to_str


def _cryptic_motif_score(candidate):
    """Average of a cryptic-exon candidate's donor and acceptor motif
    scores. Stored in evidence as a motif score, not a probability.
    """
    donor = candidate.donor_score or 0.0
    acceptor = candidate.acceptor_score or 0.0
    return (donor + acceptor) / 2.0


class MutationEffect(Serializable):
    """
    Base class for mutation effects.
    """

    def __init__(self, variant):
        self.variant = variant

    def __str__(self):
        return "%s(variant=%s)" % (self.__class__.__name__, self.variant)

    def __repr__(self):
        return str(self)

    def __lt__(self, other):
        """
        Effects are ordered by their associated variants, which have
        comparison implement in terms of their chromosomal locations.
        """
        return self.variant < other.variant

    @property
    def short_description(self):
        """
        A short but human-readable description of the effect.
        Defaults to class name for most of the non-coding effects,
        but is more informative for coding ones.
        """
        return self.__class__.__name__.lower()

    transcript = None
    gene = None

    @property
    def original_protein_sequence(self):
        """Amino acid sequence of a coding transcript (without the nucleotide
        variant/mutation)
        """
        if self.transcript:
            return self.transcript.protein_sequence
        else:
            return None

    @property
    def gene_name(self):
        if self.gene:
            return self.gene.name
        else:
            return None

    @property
    def gene_id(self):
        if self.gene:
            return self.gene.id
        else:
            return None

    @property
    def transcript_name(self):
        if self.transcript:
            return self.transcript.name
        else:
            return None

    @property
    def transcript_id(self):
        if self.transcript:
            return self.transcript.id
        else:
            return None

    # It's convenient to have a property which tells us:
    # 1) is this a variant effect overlapping a transcript?
    # 2) does that transcript have a coding sequence?
    # 3) does the variant affect the coding sequence?
    modifies_coding_sequence = False

    # Additionally:
    # 4) Does the change to the coding sequence result in a change to the
    # protein sequence?
    modifies_protein_sequence = False

    mutant_protein_sequence = None

    # which residues in the mutant protein sequence are different from
    # the original?
    aa_mutation_start_offset = None
    aa_mutation_end_offset = None

    # Optional full :class:`~varcode.MutantTranscript` attached to
    # this effect (#269 / #296 / #335). SV effects and splice
    # exon-skip candidates populate this during construction. Other
    # effect kinds (Substitution, FrameShift, etc.) leave it None by
    # default; an Isovar-style :class:`PhaseResolver` can attach an
    # assembled-contig MutantTranscript post-hoc via
    # :func:`varcode.phasing.apply_phase_resolver_to_effects`.
    mutant_transcript = None


class MultiOutcomeEffect(MutationEffect):
    """Marker base class for effects that represent a set of plausible
    outcomes rather than a single deterministic effect.

    Subclasses must expose:

    * :attr:`candidates` — tuple of :class:`~varcode.effect_candidates.EffectCandidate`
      objects in producer order. Each entry pairs an
      inner :class:`MutationEffect` (concrete or placeholder) with
      its provenance — ``source`` (producer) and ``evidence`` dict.
      The :attr:`effects` helper unwraps to the
      inner Effects when callers don't need provenance.
    * :attr:`priority_class` — effect class whose priority this set
      adopts (read by :func:`varcode.effects.effect_priority`).

    Downstream consumers filter for multi-outcome results with
    ``isinstance(effect, MultiOutcomeEffect)``, so new wrappers (RNA
    evidence #259, germline-aware #268, SV-at-breakpoint) implement
    the same protocol uniformly (#382).

    External integrations (RNA evidence, SpliceAI scoring, etc.)
    attach extra candidates post-hoc via the ``_extra_candidates``
    slot — subclasses that override :attr:`candidates` must include
    those extras in their returned tuple. The
    :meth:`_combine_with_extra_candidates` helper does the right
    thing.

    Picking *the* candidate
    -----------------------

    Two orthogonal "best candidate" notions are available; pick the
    one that matches your question:

    * **Most likely**: the first candidate after producer ordering.
      Producers preserve their own deterministic order.
      :attr:`most_likely_candidate` returns the wrapped
      :class:`EffectCandidate` (provenance + inner effect);
      :attr:`most_likely_effect` returns just the inner
      :class:`MutationEffect`. Always equal to ``candidates[0]`` /
      ``effects[0]``.

    * **Highest priority**: top by varcode's effect-priority ordering
      (see :func:`~varcode.effects.effect_priority`) — the
      most protein-disruptive candidate regardless of producer order.
      :attr:`highest_priority_candidate` and
      :attr:`highest_priority_effect` are the analogous accessors.
      Use this for clinical / functional filtering ("flag if any
      candidate is at least a frameshift"), since a disruptive
      candidate sitting behind a less-disruptive primary candidate
      should still light up.

    The two coincide when producer order and priority ranking agree,
    which is common but not guaranteed. Pick consciously.
    """

    @property
    def effects(self):
        """Tuple of inner :class:`MutationEffect` objects, in
        :attr:`candidates` order. Convenience for callers that don't
        need per-candidate provenance — equivalent to
        ``tuple(c.effect for c in self.candidates)``.
        """
        return tuple(c.effect for c in self.candidates)

    @property
    def most_likely_candidate(self):
        """The first :class:`EffectCandidate` in producer order.
        Pairs the inner effect with its ``source`` / ``evidence``
        provenance.

        For just the inner :class:`MutationEffect`, use
        :attr:`most_likely_effect`. For the most protein-disruptive
        candidate (independent of producer order), use
        :attr:`highest_priority_candidate`.
        """
        return self.candidates[0]

    @property
    def most_likely_effect(self):
        """The :class:`MutationEffect` of :attr:`most_likely_candidate`.
        Equivalent to ``most_likely_candidate.effect`` /
        ``effects[0]`` — given here so callers that don't need
        provenance don't have to reach through the wrapper.
        """
        return self.candidates[0].effect

    @property
    def highest_priority_candidate(self):
        """The :class:`EffectCandidate` whose inner effect has the
        highest :func:`~varcode.effects.effect_priority` (most
        protein-disruptive). Pure priority ranking — producer order
        deliberately doesn't factor in, so a frameshift sitting
        behind a less-disruptive primary candidate still surfaces
        here.

        Ties on priority resolve to the first matching entry of
        :attr:`candidates`, preserving the subclass's candidate order.

        Behavior is deterministic.
        """
        from .effect_ordering import effect_priority
        return max(
            self.candidates,
            key=lambda c: effect_priority(c.effect),
        )

    @property
    def highest_priority_effect(self):
        """The inner :class:`MutationEffect` of
        :attr:`highest_priority_candidate`. Use when you want the
        worst-case effect for clinical / functional filtering and
        don't need provenance.
        """
        return self.highest_priority_candidate.effect

    def _combine_with_extra_candidates(self, base_candidates):
        """Append externally-attached candidates (#259) to a derived
        ``candidates`` tuple.

        Subclasses that compute :attr:`candidates` dynamically call
        this helper on the tuple they construct so that
        :func:`varcode.rna_evidence.apply_rna_evidence_to_effects`
        and other post-hoc evidence producers can attach observed candidates
        without subclassing or monkey-patching.
        """
        extra = getattr(self, "_extra_candidates", ())
        if not extra:
            return tuple(base_candidates)
        return tuple(base_candidates) + tuple(extra)


class Intergenic(MutationEffect):
    """Variant has unknown effect if it occurs between genes"""
    short_description = "intergenic"


class Intragenic(MutationEffect):
    """
    Variant within boundaries of a gene but does not overlap
    introns or exons of any transcript. This seems very peculiar but
    apparently does happen sometimes, maybe some genes have two distinct sets
    of exons which are never simultaneously expressed?
    """
    short_description = "intragenic"

    def __init__(self, variant, gene):
        MutationEffect.__init__(self, variant)
        self.gene = gene


class TranscriptMutationEffect(Intragenic):
    def __init__(self, variant, transcript):
        Intragenic.__init__(self, variant, gene=transcript.gene)
        self.transcript = transcript

    def __str__(self):
        return "%s(variant=%s, transcript_name=%s, transcript_id=%s)" % (
            self.__class__.__name__,
            self.variant,
            self.transcript.name,
            self.transcript.id)


class Failure(TranscriptMutationEffect):
    """
    Special placeholder effect for when we want to suppress errors but still
    need to create a non-empty list of effects for each variant.
    """
    pass


class NoncodingTranscript(TranscriptMutationEffect):
    """
    Any mutation to a transcript with a non-coding biotype
    """
    short_description = "non-coding-transcript"


class IncompleteTranscript(TranscriptMutationEffect):
    """
    Any mutation to an incompletely annotated transcript with a coding biotype
    """
    short_description = "incomplete"


class FivePrimeUTR(TranscriptMutationEffect):
    """
    Any mutation to the 5' untranslated region (before the start codon) of
    coding transcript.
    """
    short_description = "5' UTR"


class ThreePrimeUTR(TranscriptMutationEffect):
    """
    Any mutation to the 3' untranslated region (after the stop codon) of
    coding transcript.
    """
    short_description = "3' UTR"


class Intronic(TranscriptMutationEffect):
    """
    Mutation in an intronic region of a coding transcript
    """
    def __init__(self, variant, transcript, nearest_exon, distance_to_exon):
        TranscriptMutationEffect.__init__(self, variant, transcript)
        self.nearest_exon = nearest_exon
        self.distance_to_exon = distance_to_exon

    short_description = "intronic"


class SpliceSite(object):
    """
    Parent class for all splice site mutations.
    """
    pass


class IntronicSpliceSite(Intronic, SpliceSite):
    """
    Mutations near exon boundaries, excluding the first two and last two
    nucleotides in an intron, since those are  known to more confidently
    affect splicing and are given their own effect classes below.
    """
    def __init__(self, variant, transcript, nearest_exon, distance_to_exon):
        Intronic.__init__(
            self, variant, transcript, nearest_exon, distance_to_exon)

    short_description = "intronic-splice-site"


class SpliceDonor(IntronicSpliceSite):
    """
    Mutation in the first two intron residues.
    """
    def __init__(self, variant, transcript, nearest_exon, distance_to_exon):
        IntronicSpliceSite.__init__(
            self, variant, transcript, nearest_exon, distance_to_exon)

    short_description = "splice-donor"


class SpliceAcceptor(IntronicSpliceSite):
    """
    Mutation in the last two intron residues.
    """
    short_description = "splice-acceptor"


# =====================================================================
# Splice-mechanism effects (#382).
#
# What the splice machinery does when a disrupted splice signal forces
# a choice. Each subclass represents one mechanism and carries the
# protein-level consequence of that mechanism in its own vocab —
# ``aa_ref`` / ``aa_alt`` / ``mutant_protein_sequence`` /
# ``mutant_transcript``, populated when varcode can resolve the protein
# math, ``None`` otherwise. The "unresolved" state is "protein fields
# are None"; no parallel placeholder class hierarchy.
#
# Every subclass references back to the underlying splice-signal
# disruption via :attr:`splice_signal` — a SpliceDonor / SpliceAcceptor /
# IntronicSpliceSite / ExonicSpliceSite Effect describing *where* in
# the transcript the splice signal was hit. The mechanism Effect
# describes *what* happens downstream of that disruption.
#
# Used by :class:`~varcode.splice_outcomes.SpliceOutcomeSet` to
# fill :attr:`EffectCandidate.effect` for each plausible mechanism a
# splice disruption could produce. Class identity = mechanism — no
# parallel enum.
# =====================================================================


class SpliceMechanismEffect(TranscriptMutationEffect):
    """Base for effects describing what the splice machinery does
    downstream of a disrupted splice signal.

    Concrete subclasses: :class:`NormalSplicing`, :class:`ExonSkipping`,
    :class:`IntronRetention`, :class:`CrypticDonor`,
    :class:`CrypticAcceptor`. Class identity is the mechanism;
    instances of this base are not constructed directly.

    Protein vocab fields are optional — ``None`` means "predicted
    mechanism, exact protein not computable from cached cDNA alone."
    When resolved, the fields follow the same conventions as
    :class:`KnownAminoAcidChange` so consumers can read them
    uniformly.

    Parameters
    ----------
    variant : Variant
    transcript : pyensembl.Transcript
    splice_signal : MutationEffect
        The splice-signal-disruption Effect (``SpliceDonor`` /
        ``SpliceAcceptor`` / ``IntronicSpliceSite`` /
        ``ExonicSpliceSite``) describing *where* the disruption was
        in the transcript. Lets consumers reach
        ``effect.splice_signal.nearest_exon``,
        ``effect.splice_signal.distance_to_exon``, etc. without
        carrying those fields redundantly on every mechanism subclass.
    protein_effect : MutationEffect or None
        Classified protein-level consequence of this splice mechanism
        (e.g. ``Deletion``, ``FrameShift``, ``PrematureStop``) when the
        protein sequence was resolved. The mechanism object keeps its
        own class identity while delegating severity metadata to this
        classified consequence.
    mutant_transcript : MutantTranscript or None
    aa_ref, aa_alt : str or None
    aa_mutation_start_offset, aa_mutation_end_offset : int or None
    mutant_protein_sequence : str or None
    """

    def __init__(
            self, variant, transcript, splice_signal,
            protein_effect=None,
            mutant_transcript=None,
            aa_ref=None, aa_alt=None,
            aa_mutation_start_offset=None,
            aa_mutation_end_offset=None,
            mutant_protein_sequence=None):
        TranscriptMutationEffect.__init__(self, variant, transcript)
        self.splice_signal = splice_signal
        self.protein_effect = protein_effect
        self.mutant_transcript = mutant_transcript
        self.aa_ref = aa_ref
        self.aa_alt = aa_alt
        self.aa_mutation_start_offset = aa_mutation_start_offset
        self.aa_mutation_end_offset = aa_mutation_end_offset
        self.mutant_protein_sequence = mutant_protein_sequence

    @property
    def resolved(self) -> bool:
        """True when the protein math succeeded for this mechanism —
        i.e. :attr:`mutant_protein_sequence` is populated. Convenience
        for consumers that want to branch on "do we have a concrete
        protein for this mechanism?"
        """
        return self.mutant_protein_sequence is not None

    @property
    def priority_class(self):
        if self.protein_effect is None:
            return None
        return self.protein_effect.__class__

    @property
    def modifies_coding_sequence(self):
        return getattr(self.protein_effect, "modifies_coding_sequence", False)

    @property
    def modifies_protein_sequence(self):
        return getattr(self.protein_effect, "modifies_protein_sequence", False)


class NormalSplicing(SpliceMechanismEffect):
    """Canonical splicing proceeds despite the disruption.

    The splice signal was hit, but the spliceosome handles it; the
    protein consequence (if any) is whatever the underlying
    nucleotide change normally produces — a :class:`Substitution`
    for an :class:`ExonicSpliceSite`, no protein change for a
    purely intronic disruption. The underlying coding change is
    carried as :attr:`coding_effect`; when there isn't one (intronic
    variant), it's ``None``.

    Note: this class is special among :class:`SpliceMechanismEffect`
    subclasses — it doesn't itself transform the protein, so all the
    protein-vocab attributes (``aa_ref``, ``aa_alt``,
    ``mutant_protein_sequence``, ``mutant_transcript``,
    ``aa_mutation_*_offset``) are *delegating descriptors* that read
    through to ``coding_effect``. Construction therefore calls
    :class:`TranscriptMutationEffect`'s ``__init__`` directly rather
    than :class:`SpliceMechanismEffect`'s — the latter would set
    those names as plain instance attributes and shadow the
    descriptors. The other four mechanisms (skipping, retention,
    cryptic donor/acceptor) actually transform the protein and
    carry their own protein vocab as instance state.
    """

    def __init__(self, variant, transcript, splice_signal,
                 coding_effect=None):
        TranscriptMutationEffect.__init__(self, variant, transcript)
        self.splice_signal = splice_signal
        self.coding_effect = coding_effect

    # Delegating descriptors. Read-only by construction (no setters)
    # — NormalSplicing never owns these fields directly; they're a
    # view onto ``coding_effect``. Attempting to assign to one is a
    # programming error and rightly raises AttributeError.

    @property
    def aa_ref(self):
        return getattr(self.coding_effect, "aa_ref", None)

    @property
    def aa_alt(self):
        return getattr(self.coding_effect, "aa_alt", None)

    @property
    def mutant_protein_sequence(self):
        return getattr(self.coding_effect, "mutant_protein_sequence", None)

    @property
    def aa_mutation_start_offset(self):
        return getattr(self.coding_effect, "aa_mutation_start_offset", None)

    @property
    def aa_mutation_end_offset(self):
        return getattr(self.coding_effect, "aa_mutation_end_offset", None)

    @property
    def mutant_transcript(self):
        return getattr(self.coding_effect, "mutant_transcript", None)

    @property
    def priority_class(self):
        if self.coding_effect is None:
            return None
        return self.coding_effect.__class__

    @property
    def modifies_coding_sequence(self):
        return getattr(self.coding_effect, "modifies_coding_sequence", False)

    @property
    def modifies_protein_sequence(self):
        return getattr(self.coding_effect, "modifies_protein_sequence", False)

    @property
    def short_description(self):
        if self.coding_effect is None:
            return "normal-splicing"
        return "normal-splicing (%s)" % self.coding_effect.short_description


class ExonSkipping(SpliceMechanismEffect):
    """The affected exon is excluded from the mature transcript.

    When in-frame (the exon's length is divisible by 3), the result
    is a clean amino-acid deletion plus an optional boundary-codon
    reshape. When out-of-frame, a frameshift propagates from the
    new exon junction through the rest of the transcript.

    Parameters beyond the base class
    --------------------------------
    affected_exon : pyensembl.Exon
        The skipped exon.
    in_frame : bool
        ``True`` when ``len(affected_exon) % 3 == 0`` — the
        exon-skip is codon-aligned and the downstream sequence
        keeps its frame.

    Note on ``aa_ref`` / ``aa_alt`` semantics
    -----------------------------------------
    The AA-level fields come from
    :func:`~varcode.effects.classify.classify_from_protein_diff`,
    so what they mean depends on ``in_frame``:

    * **In-frame skip** (``in_frame=True``): ``aa_ref`` is the
      contiguous run of amino acids removed from the reference
      protein (including the boundary-codon reshape when the exon
      boundary lands mid-codon). ``aa_alt`` is the reshaped boundary
      codon's translation (or ``""`` for a clean codon-aligned skip).
    * **Out-of-frame skip** (``in_frame=False``): ``aa_ref`` is the
      reference AA at the divergence point (where the frame breaks
      after the new exon junction); ``aa_alt`` is the first
      frameshifted AA. The downstream protein continues in the new
      frame until a premature stop — read the full sequence from
      :attr:`mutant_protein_sequence`, not from ``aa_alt``.
    """

    def __init__(self, variant, transcript, splice_signal, affected_exon,
                 in_frame, mutant_transcript=None,
                 protein_effect=None,
                 aa_ref=None, aa_alt=None,
                 aa_mutation_start_offset=None,
                 aa_mutation_end_offset=None,
                 mutant_protein_sequence=None):
        SpliceMechanismEffect.__init__(
            self, variant, transcript, splice_signal,
            protein_effect=protein_effect,
            mutant_transcript=mutant_transcript,
            aa_ref=aa_ref, aa_alt=aa_alt,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_mutation_end_offset=aa_mutation_end_offset,
            mutant_protein_sequence=mutant_protein_sequence)
        self.affected_exon = affected_exon
        self.in_frame = in_frame

    @property
    def short_description(self):
        exon_id = getattr(self.affected_exon, "exon_id", "?")
        if not self.resolved:
            suffix = "" if self.in_frame else ", out-of-frame"
            return "exon-skip %s (predicted%s)" % (exon_id, suffix)
        if not self.aa_ref or self.aa_mutation_start_offset is None:
            # Resolved but no usable AA range — fall back to a bare label.
            return "exon-skip %s" % exon_id
        pos = self.aa_mutation_start_offset + 1
        if self.in_frame:
            # HGVS deletion: p.{first}{start}_{last}{end}del
            end = self.aa_mutation_start_offset + len(self.aa_ref)
            return "exon-skip %s (p.%s%d_%s%ddel)" % (
                exon_id, self.aa_ref[0], pos, self.aa_ref[-1], end)
        # HGVS frameshift: p.{aa}{pos}fs — distinct notation from del,
        # not concatenated with it.
        return "exon-skip %s (p.%s%dfs)" % (exon_id, self.aa_ref[0], pos)


class IntronRetention(SpliceMechanismEffect):
    """The intron stays in the mature transcript; translation
    typically hits a premature stop inside the retained intron.

    Parameters beyond the base class
    --------------------------------
    retained_intron_start, retained_intron_end : int or None
        Genomic coordinates of the retained intron (1-based
        inclusive). Populated when the adjacent-exon geometry
        identifies which intron is being retained; ``None`` when
        that can't be determined (e.g. transcript missing the
        adjacent exon, splice boundary ambiguous).
    side : str or None
        ``"donor"`` or ``"acceptor"`` — which splice boundary failed,
        causing the intron to be retained. ``None`` when the side
        can't be inferred from the splice classification.
    """

    def __init__(self, variant, transcript, splice_signal,
                 retained_intron_start=None, retained_intron_end=None,
                 side=None,
                 mutant_transcript=None,
                 protein_effect=None,
                 aa_ref=None, aa_alt=None,
                 aa_mutation_start_offset=None,
                 aa_mutation_end_offset=None,
                 mutant_protein_sequence=None):
        SpliceMechanismEffect.__init__(
            self, variant, transcript, splice_signal,
            protein_effect=protein_effect,
            mutant_transcript=mutant_transcript,
            aa_ref=aa_ref, aa_alt=aa_alt,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_mutation_end_offset=aa_mutation_end_offset,
            mutant_protein_sequence=mutant_protein_sequence)
        self.retained_intron_start = retained_intron_start
        self.retained_intron_end = retained_intron_end
        self.side = side

    @property
    def short_description(self):
        side = self.side or "?"
        if not self.resolved:
            return "intron-retention (%s side, predicted)" % side
        if self.aa_mutation_start_offset is None or not self.aa_ref:
            return "intron-retention (%s side)" % side
        return "intron-retention (%s side, p.%s%dfs*)" % (
            side, self.aa_ref[0], self.aa_mutation_start_offset + 1)


class CrypticSpliceSiteEffect(SpliceMechanismEffect):
    """Base for cryptic donor / acceptor effects. Carries the
    cryptic motif's position and score, plus the resulting exon
    length delta."""

    direction = None  # "donor" / "acceptor", set by subclass

    def __init__(self, variant, transcript, splice_signal, affected_exon,
                 cryptic_genomic_position=None, motif_score=None,
                 exon_length_delta=None,
                 mutant_transcript=None,
                 protein_effect=None,
                 aa_ref=None, aa_alt=None,
                 aa_mutation_start_offset=None,
                 aa_mutation_end_offset=None,
                 mutant_protein_sequence=None):
        SpliceMechanismEffect.__init__(
            self, variant, transcript, splice_signal,
            protein_effect=protein_effect,
            mutant_transcript=mutant_transcript,
            aa_ref=aa_ref, aa_alt=aa_alt,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_mutation_end_offset=aa_mutation_end_offset,
            mutant_protein_sequence=mutant_protein_sequence)
        self.affected_exon = affected_exon
        self.cryptic_genomic_position = cryptic_genomic_position
        self.motif_score = motif_score
        self.exon_length_delta = exon_length_delta

    @property
    def short_description(self):
        if not self.resolved:
            return "cryptic-%s (predicted)" % self.direction
        delta = self.exon_length_delta or 0
        kind = "ext" if delta > 0 else "trunc" if delta < 0 else "shift"
        return "cryptic-%s %s%dbp" % (
            self.direction, kind, abs(delta))


class CrypticDonor(CrypticSpliceSiteEffect):
    """A cryptic GT donor is used instead of the disrupted canonical
    donor; the exon is extended (cryptic site downstream) or
    truncated (cryptic site upstream)."""
    direction = "donor"


class CrypticAcceptor(CrypticSpliceSiteEffect):
    """A cryptic AG acceptor is used instead of the disrupted
    canonical acceptor; the exon is extended (cryptic site upstream)
    or truncated (cryptic site downstream)."""
    direction = "acceptor"


class Exonic(TranscriptMutationEffect):
    """
    Any mutation which affects the contents of an exon (coding region or UTRs)
    """
    pass


class ExonLoss(Exonic):
    """
    Deletion of one or more exons in a transcript.
    """
    def __init__(self, variant, transcript, exons):
        Exonic.__init__(self, variant, transcript)
        self.exons = exons

    def __str__(self):
        return "ExonLoss(%s, %s, %s)" % (
            self.variant,
            self.transcript,
            "+".join(str(exon) for exon in self.exons))

    short_description = "exon-loss"

    @property
    def modifies_protein_sequence(self):
        # TODO: distinguish between exon loss in the CDS and UTRs
        return True

    @property
    def modifies_coding_sequence(self):
        # TODO: distinguish between exon loss in the CDS and UTRs
        return True


class ExonicSpliceSite(Exonic, SpliceSite, MultiOutcomeEffect):
    """
    Mutation in the last three nucleotides before an intron
    or in the first nucleotide after an intron.

    Expresses the two plausible outcomes of a splice-adjacent exonic
    variant — splice-signal disruption vs. the underlying coding
    change if splicing proceeds normally — via the
    :class:`MultiOutcomeEffect` protocol (#299, #382).
    :attr:`candidates` yields a 2-tuple of
    :class:`~varcode.effect_candidates.EffectCandidate` objects;
    ``alternate_effect`` stays on the instance as a first-class field
    for back-compat with callers that depended on it.

    This is the **lightweight 2-outcome form**. Callers that want
    the richer exon-skipping / intron-retention / cryptic-splice
    candidate set opt into ``splice_outcomes=True`` and get a
    :class:`~varcode.splice_outcomes.SpliceOutcomeSet` instead.
    """
    def __init__(self, variant, transcript, exon, alternate_effect):
        Exonic.__init__(self, variant, transcript)
        self.exon = exon
        self.alternate_effect = alternate_effect

    def __str__(self):
        return "ExonicSpliceSite(exon=%s, alternate_effect=%s)" % (
            self.exon,
            self.alternate_effect)

    short_description = "exonic-splice-site"

    @property
    def mutant_protein_sequence(self):
        """
        TODO: determine when exonic splice variants cause exon skipping
        vs. translation of the underlying modified coding sequence.

        For now just pretending like there is no effect on splicing.
        """
        return self.alternate_effect.mutant_protein_sequence

    @property
    def modifies_protein_sequence(self):
        return self.alternate_effect.modifies_protein_sequence

    @property
    def modifies_coding_sequence(self):
        return self.alternate_effect.modifies_coding_sequence

    @property
    def candidates(self):
        """The two plausible outcomes wrapped as
        :class:`~varcode.effect_candidates.EffectCandidate`:

        Position [0] is the splice-disruption outcome (this effect
        itself); position [1] is the coding change if splicing
        proceeds. Returns a 1-tuple when ``alternate_effect`` is None.
        Extra candidates attached via :meth:`_combine_with_extra_candidates`
        (e.g. from RNA evidence) come after.
        """
        from ..effect_candidates import EffectCandidate
        base = [EffectCandidate(effect=self, source="varcode")]
        if self.alternate_effect is not None:
            base.append(EffectCandidate(
                effect=self.alternate_effect, source="varcode"))
        return self._combine_with_extra_candidates(tuple(base))

    @property
    def priority_class(self):
        """ExonicSpliceSite priority is used directly — no
        delegation, since this class IS the splice-adjacent effect
        rather than a wrapper around one.
        """
        return ExonicSpliceSite


class CodingMutation(Exonic):
    """
    Base class for all mutations which result in a modified coding sequence.
    """
    def __str__(self):
        fields = [
            ("variant", self.variant),
            ("transcript_name", self.transcript.name),
            ("transcript_id", self.transcript.id),
            ("effect_description", self.short_description)
        ]
        fields_str = ", ".join([
            "%s=%s" % (field_name, field_value)
            for (field_name, field_value)
            in fields
        ])
        return "%s(%s)" % (self.__class__.__name__, fields_str)

    @property
    def modifies_coding_sequence(self):
        return True


class Silent(CodingMutation):
    """Mutation to an exon of a coding region which doesn't change the
    amino acid sequence.
    """
    def __init__(
            self,
            variant,
            transcript,
            aa_pos,
            aa_ref):
        """
        Parameters
        ----------
        variant : Variant

        transcript : Transcript

        aa_pos : int
            Offset of first synonymous codon in protein sequence

        aa_ref : str or Bio.Seq
            Reference amino acid(s) at offset
        """
        CodingMutation.__init__(
            self,
            variant=variant,
            transcript=transcript)
        self.aa_pos = aa_pos
        self.aa_ref = bio_seq_to_str(aa_ref)

    @property
    def short_description(self):
        # HGVS p.= notation for synonymous changes: "p.{aa_ref}{1-indexed pos}="
        # (http://varnomen.hgvs.org/). If aa_ref is empty (the silent
        # change was fully captured by the shared prefix/suffix and we
        # have nothing to name), fall back to a positional-only form.
        if self.aa_ref:
            return "p.%s%d=" % (self.aa_ref, self.aa_pos + 1)
        return "p.%d=" % (self.aa_pos + 1)


class AlternateStartCodon(Silent):
    """Change to the start codon (e.g. ATG>CTG) but without changing the
    starting amino acid from methionine.
    """
    def __init__(
            self,
            variant,
            transcript,
            ref_codon,
            alt_codon):
        Silent.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_pos=0,
            aa_ref=transcript.protein_sequence[0])
        self.ref_codon = bio_seq_to_str(ref_codon)
        self.alt_codon = bio_seq_to_str(alt_codon)

    @property
    def short_description(self):
        return "alternate-start-codon (%s>%s)" % (
            self.ref_codon, self.alt_codon)


class NonsilentCodingMutation(CodingMutation):
    """
    All coding mutations other than silent codon substitutions
    """

    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_mutation_end_offset,
            aa_ref):
        """
        variant : Variant

        transcript : Transcript

        aa_mutation_start_offset : int
            Offset of first modified amino acid in protein (starting from 0)

        aa_mutation_end_offset : int
            Offset after last mutated amino acid (half-open coordinates)

        aa_ref : str
            Amino acid string of what used to be at aa_mutation_start_offset
            in the wildtype (unmutated) protein.
        """
        CodingMutation.__init__(
            self,
            variant=variant,
            transcript=transcript)
        self.aa_mutation_start_offset = aa_mutation_start_offset
        self.aa_mutation_end_offset = aa_mutation_end_offset
        self.aa_ref = bio_seq_to_str(aa_ref)

    @property
    def modifies_protein_sequence(self):
        return True


class StartLoss(NonsilentCodingMutation):
    """
    When a start codon is lost it's difficult to determine if there is
    an alternative Kozak consensus sequence (either before or after the
    original) from which an alternative start codon can be inferred.

    TODO:
        - look for downstream alternative start codon to predict
          new coding sequence (probably also requires matching
          pattern of preceding ~6nt)
        - If an alternative start codon is changed to ATG then
          we should make a StrongerStartCodon effect which is effectively
          silent
        - If ATG is changed to the two common alternative codons then
          we should make a WeakerStartCodon effect which is also
          effectively silent.
    """
    def __init__(
            self,
            variant,
            transcript):
        NonsilentCodingMutation.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=0,
            aa_mutation_end_offset=1,
            aa_ref=transcript.protein_sequence[0])

    @property
    def mutant_protein_sequence(self):
        # TODO: scan downstream in the cDNA sequence to predict an alternate
        # start codon using a Kozak sequence template
        return None

    @property
    def short_description(self):
        return "p.%s1? (start-loss)" % (self.transcript.protein_sequence[0],)


class KnownAminoAcidChange(NonsilentCodingMutation):
    """
    Coding mutations in which we can predict what the new/mutant protein
    sequence will be.
    """
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_ref,
            aa_alt):
        NonsilentCodingMutation.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_mutation_end_offset=aa_mutation_start_offset + len(aa_alt),
            aa_ref=aa_ref)
        self.aa_alt = bio_seq_to_str(aa_alt)

    @property
    def short_description(self):
        if len(self.aa_ref) == 0:
            return "p.%dins%s" % (
                self.aa_mutation_start_offset,
                self.aa_alt)
        elif len(self.aa_alt) == 0:
            return "p.%s%ddel" % (
                self.aa_ref,
                self.aa_mutation_start_offset)
        else:
            return "p.%s%d%s" % (
                self.aa_ref,
                self.aa_mutation_start_offset + 1,
                self.aa_alt)

    @memoized_property
    def mutant_protein_sequence(self):
        original = self.original_protein_sequence
        prefix = original[:self.aa_mutation_start_offset]
        suffix = original[self.aa_mutation_start_offset + len(self.aa_ref):]
        return prefix + self.aa_alt + suffix


class Substitution(KnownAminoAcidChange):
    """
    Single amino acid substitution, e.g. BRAF-001 V600E
    """
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_ref,
            aa_alt):
        if len(aa_ref) != 1:
            raise ValueError(
                "Simple substitution can't have aa_ref='%s'" % (aa_ref,))
        if len(aa_alt) != 1:
            raise ValueError(
                "Simple substitution can't have aa_alt='%s'" % (aa_alt,))
        KnownAminoAcidChange.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt=aa_alt)


class ComplexSubstitution(KnownAminoAcidChange):
    """
    In-frame substitution of multiple amino acids, e.g. TP53-002 p.391FY>QQQ
    Can change the length of the protein sequence but since it has
    non-empty ref and alt strings, is more complicated than an insertion or
    deletion alone.
    """
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_ref,
            aa_alt):
        if len(aa_ref) == 1 and len(aa_alt) == 1:
            raise ValueError(
                "ComplexSubstitution can't have aa_ref='%s' and aa_alt='%s'" % (
                    aa_ref, aa_alt))
        KnownAminoAcidChange.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt=aa_alt)


class Insertion(KnownAminoAcidChange):
    """
    In-frame insertion of one or more amino acids.
    """
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_alt):
        KnownAminoAcidChange.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref="",
            aa_alt=aa_alt)


class Deletion(KnownAminoAcidChange):
    """
    In-frame deletion of one or more amino acids.
    """

    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_ref):
        KnownAminoAcidChange.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt="")


class PrematureStop(KnownAminoAcidChange):
    """In-frame insertion of codons containing a stop codon. May also involve
    insertion/deletion/substitution of other amino acids preceding the stop."""
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_ref="",
            aa_alt=""):
        """
        Insertion of premature stop codon, possibly preceded by a substitution
        of `aa_ref` amino acids for `aa_alt` alternative residues.
        """
        if "*" in aa_ref:
            raise ValueError(
                ("Unexpected aa_ref = '%s', should only include amino acids "
                 "before the new stop codon.") % aa_ref)
        if "*" in aa_alt:
            raise ValueError(
                ("Unexpected aa_ref = '%s', should only include amino acids "
                 "before the new stop codon.") % aa_alt)
        KnownAminoAcidChange.__init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt=aa_alt)
        self.stop_codon_offset = aa_mutation_start_offset + len(aa_alt)

        if self.stop_codon_offset >= len(transcript.protein_sequence):
            raise ValueError(
                ("Premature stop codon cannot be at position %d"
                 " since the original protein of %s has length %d") % (
                    self.stop_codon_offset,
                    transcript,
                    len(transcript.protein_sequence)))

    @property
    def short_description(self):
        # When aa_ref is empty, the mutation is an insertion (no residues
        # were replaced) — use HGVS "ins" notation so the description is
        # unambiguous. Otherwise fall through to "{ref}{pos}{alt}*".
        if len(self.aa_ref) == 0:
            return "p.%dins%s*" % (
                self.aa_mutation_start_offset + 1,
                self.aa_alt)
        return "p.%s%d%s*" % (
            self.aa_ref,
            self.aa_mutation_start_offset + 1,
            self.aa_alt)

    @memoized_property
    def mutant_protein_sequence(self):
        prefix = self.original_protein_sequence[:self.aa_mutation_start_offset]
        return prefix + self.aa_alt


class StopLoss(KnownAminoAcidChange):
    def __init__(
            self,
            variant,
            transcript,
            aa_ref,
            aa_alt):
        # StopLoss assumes that we deleted some codons ending with a
        # stop codon
        if "*" in aa_ref:
            raise ValueError(
                "StopLoss aa_ref '%s' should not contain '*'" % (
                    aa_ref,))
        if len(aa_alt) == 0:
            raise ValueError(
                "If no amino acids added by StopLoss then it should be Silent")
        # subtract 1 for the stop codon
        n_ref_amino_acids = len(aa_ref)
        protein_length = len(transcript.protein_sequence)
        aa_mutation_start_offset = protein_length - n_ref_amino_acids
        KnownAminoAcidChange.__init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_alt=aa_alt,
            aa_ref=aa_ref)

    @property
    def extended_protein_sequence(self):
        """Deprecated name for aa_alt"""
        return self.aa_alt

    @property
    def short_description(self):
        return "p.%s*%d%s (stop-loss)" % (
            self.aa_ref,
            self.aa_mutation_start_offset + 1,
            self.extended_protein_sequence)


class FrameShift(KnownAminoAcidChange):
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            shifted_sequence):
        """Frameshift mutation preserves all the amino acids up to
        aa_mutation_start_offset and then replaces the rest of the protein with
        new (frameshifted) sequence. Unlike an insertion, where we denote with
        aa_ref as the chracter before the variant sequence, a frameshift starts
        at aa_ref.
        """
        aa_ref = transcript.protein_sequence[aa_mutation_start_offset:]
        KnownAminoAcidChange.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt=shifted_sequence)

    @property
    def shifted_sequence(self):
        return self.aa_alt

    @memoized_property
    def mutant_protein_sequence(self):
        original_aa_sequence = self.original_protein_sequence
        prefix = original_aa_sequence[:self.aa_mutation_start_offset]
        return prefix + self.shifted_sequence

    @property
    def short_description(self):
        return "p.%s%dfs" % (
            self.aa_ref[0],
            self.aa_mutation_start_offset + 1)


class FrameShiftTruncation(PrematureStop, FrameShift):
    """
    A frame-shift mutation which immediately introduces a stop codon.
    """
    def __init__(
            self,
            variant,
            transcript,
            stop_codon_offset,
            aa_ref=""):
        PrematureStop.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=stop_codon_offset,
            aa_ref=aa_ref,
            aa_alt="")

    @property
    def shifted_sequence(self):
        return self.aa_alt

    @property
    def short_description(self):
        return "p.%s%dfs*" % (
            self.aa_ref,
            self.aa_mutation_start_offset + 1)


# =====================================================================
# Structural-variant effects (#252 / PR 10)
#
# These describe consequences of SVs at the transcript level. A single
# SV commonly produces multiple plausible consequences (fusion vs exon
# loss vs readthrough), so each of these is a ``MultiOutcomeEffect``
# subclass exposing a :attr:`candidates` tuple per #299 / #382.
#
# The classes are deliberately thin wrappers — the important interface
# is :attr:`candidates`, which carries :class:`EffectCandidate` objects
# with source / evidence provenance so external tools (RNA evidence,
# SpliceAI, long-read assembly) can annotate them without subclassing.
# =====================================================================


class StructuralVariantEffect(TranscriptMutationEffect, MultiOutcomeEffect):
    """Base class for effects of a :class:`StructuralVariant` on a
    specific transcript. Subclasses set ``short_description``; the
    :attr:`candidates` tuple is what downstream callers read.

    :attr:`candidates` always returns a tuple of
    :class:`~varcode.effect_candidates.EffectCandidate` objects.
    Primary candidates carry ``source="varcode"``; cryptic-exon
    candidates attached via :meth:`_attach_cryptic_candidates` carry
    ``source="varcode_motif"``; splice-outcome candidates attached
    via :meth:`_attach_splice_outcomes` carry ``source="varcode_splice"``.
    """

    def __init__(
            self, variant, transcript,
            primary_effects=None, mutant_transcript=None):
        TranscriptMutationEffect.__init__(self, variant, transcript)
        # Primary classifications: either the subclass's own effect
        # (``self``) or a caller-supplied tuple of
        # :class:`MutationEffect` instances. Lifted into
        # :class:`EffectCandidate` shape on access.
        self._primary_effects = (
            tuple(primary_effects) if primary_effects is not None
            else (self,))
        self.mutant_transcript = mutant_transcript
        # Cryptic-exon candidates nominated by
        # :func:`varcode.cryptic_exons.enumerate_from_structural_variant`
        # (#337). Kept separate from primary effects because their
        # entries carry different source / evidence.
        self._cryptic_candidates = ()
        # Splice-outcome candidates attached by the SV annotator when
        # an SV breakpoint lands in a canonical splice window (#341).
        # Pre-constructed :class:`~varcode.EffectCandidate` tuples; the
        # annotator re-sources them as ``"varcode_splice"`` and
        # enriches evidence with ``sv_type`` before attaching.
        self._splice_candidates = ()

    @property
    def candidates(self):
        """Unified :class:`~varcode.EffectCandidate` view over primary
        SV classifications, attached cryptic candidates, and any
        splice-outcome candidates (#339, #337, #341, #382).

        Primary candidates carry ``source="varcode"``; cryptic-exon
        candidates carry ``source="varcode_motif"``; splice-outcome
        candidates carry ``source="varcode_splice"``. Each source
        marks provenance so external scorers (SpliceAI, Pangolin,
        RNA evidence) can filter before rescoring.
        """
        from ..effect_candidates import EffectCandidate
        sv_type = getattr(self.variant, "sv_type", None)
        base_evidence = {"sv_type": sv_type} if sv_type is not None else {}
        primary = tuple(
            EffectCandidate(
                effect=effect,
                source="varcode",
                evidence=base_evidence)
            for effect in self._primary_effects)
        cryptic = tuple(
            EffectCandidate(
                effect=c,
                source="varcode_motif",
                evidence={
                    **base_evidence,
                    "donor_score": c.donor_score,
                    "acceptor_score": c.acceptor_score,
                    "motif_score": _cryptic_motif_score(c),
                    "interval_start": c.interval_start,
                    "interval_end": c.interval_end,
                })
            for c in self._cryptic_candidates)
        return self._combine_with_extra_candidates(
            primary + cryptic + tuple(self._splice_candidates))

    def _attach_cryptic_candidates(self, cryptic_candidates):
        """Attach cryptic-exon candidates (#337). Called by the SV
        annotator after effect construction so the candidates appear
        as additional :class:`EffectCandidate` entries on
        :attr:`candidates` without changing the primary classification.
        """
        self._cryptic_candidates = tuple(cryptic_candidates)

    def _attach_splice_outcomes(self, splice_outcomes):
        """Attach splice-outcome candidates that the SV annotator
        generated by feeding a synthesized splice-disrupting effect
        into :func:`enumerate_splice_outcomes` (#341). Each entry is
        a pre-constructed :class:`~varcode.EffectCandidate`; the annotator
        has already re-sourced them as ``"varcode_splice"`` and
        enriched evidence with the SV type.
        """
        self._splice_candidates = tuple(splice_outcomes)


class LargeDeletion(StructuralVariantEffect):
    """A deletion (``<DEL>`` / ``<CN0>``) that removes one or more
    exons — or an entire gene. Carries the list of affected exons
    for downstream analysis; the single ``outcome`` is this effect
    itself (callers that want to add RNA evidence or SpliceAI
    scoring construct additional outcomes and wrap the result)."""

    short_description = "sv-deletion"

    def __init__(
            self, variant, transcript, affected_exons,
            primary_effects=None, mutant_transcript=None):
        StructuralVariantEffect.__init__(
            self, variant, transcript,
            primary_effects=primary_effects,
            mutant_transcript=mutant_transcript)
        self.affected_exons = tuple(affected_exons)

    def __str__(self):
        return "LargeDeletion(%s, %s, exons=%d)" % (
            self.variant.short_description,
            self.transcript.id,
            len(self.affected_exons))


class LargeDuplication(StructuralVariantEffect):
    """A tandem duplication (``<DUP>``) overlapping exons.
    Biologically may produce a copy-number increase, a fused
    reading frame if junctions land in-frame, or a regulatory
    effect. Varcode reports the affected exons and leaves scoring
    to downstream tools."""

    short_description = "sv-duplication"

    def __init__(
            self, variant, transcript, affected_exons,
            primary_effects=None, mutant_transcript=None):
        StructuralVariantEffect.__init__(
            self, variant, transcript,
            primary_effects=primary_effects,
            mutant_transcript=mutant_transcript)
        self.affected_exons = tuple(affected_exons)


class Inversion(StructuralVariantEffect):
    """An inversion (``<INV>``) that flips a stretch of a transcript.
    Depending on whether breakpoints fall in exons or introns, the
    consequence is very different (exonic: likely disruptive; purely
    intronic: may or may not affect splicing). Reported as a single
    outcome here; subclasses / callers can enrich with cryptic-
    splice candidates per PR 11."""

    short_description = "sv-inversion"


class GeneFusion(StructuralVariantEffect):
    """A breakend (``<BND>``) whose mate lies in another
    protein-coding gene — the canonical fusion shape.

    Carries the two partner transcripts (5' and 3') and, when the
    annotator has enough context, a :class:`MutantTranscript` built
    from :class:`ReferenceSegment` entries describing the fused
    allele. Predicting the exact fused-protein sequence requires
    knowing which exons are retained, which typically needs RNA
    evidence — outcomes beyond "this is a plausible fusion" are
    left to downstream tools that attach :class:`EffectCandidate`
    objects with their own producer ``source`` tag.
    """

    short_description = "sv-gene-fusion"

    def __init__(
            self, variant, transcript, partner_transcript,
            mutant_transcript=None, primary_effects=None):
        StructuralVariantEffect.__init__(
            self, variant, transcript,
            primary_effects=primary_effects,
            mutant_transcript=mutant_transcript)
        self.partner_transcript = partner_transcript


class TranslocationToIntergenic(StructuralVariantEffect):
    """A breakend whose mate lies in intergenic space. The
    downstream consequence depends on whether the intergenic region
    contains cryptic splice / ORF signals — reported as a single
    outcome here, with PR 11's cryptic-exon enumerator adding
    candidate outcomes when applicable."""

    short_description = "sv-translocation-intergenic"


class CrypticExonCandidate(MutationEffect):
    """A region where an SV has brought novel sequence into range
    of the transcript, and motif scoring flags a plausible new
    splice acceptor / donor pair. Produced by PR 11's cryptic-exon
    enumerator; attached as additional :class:`EffectCandidate` entries on
    SV effects rather than standalone.

    Not a :class:`TranscriptMutationEffect` because the candidate
    region may not overlap any existing transcript — it's a *new*
    exon hypothesis. Carries the contig / interval and the motif
    scores as plain fields; external predictors (SpliceAI,
    Pangolin) attach their own scores via the enclosing
    :class:`EffectCandidate.evidence` dict.
    """

    short_description = "sv-cryptic-exon-candidate"

    def __init__(
            self, variant, contig, interval_start, interval_end,
            donor_score=None, acceptor_score=None):
        MutationEffect.__init__(self, variant)
        self.contig = contig
        self.interval_start = interval_start
        self.interval_end = interval_end
        self.donor_score = donor_score
        self.acceptor_score = acceptor_score


# =====================================================================
# Haplotype effects — joint consequence of multiple cis variants on
# the same transcript (#269).
#
# When a phase resolver (VCF PS, Isovar assembly, long-read haplotype
# caller) confirms two or more variants sit on the same physical copy
# of a transcript, their combined effect on the protein isn't just the
# sum of per-variant effects — they can share a codon, one can rescue
# another's frameshift, etc. HaplotypeEffect carries the joint
# MutantTranscript (via apply_variants_to_transcript) and coexists
# with the per-variant effects on the same EffectCollection so
# downstream consumers can pick the granularity they need.
# =====================================================================


class HaplotypeEffect(TranscriptMutationEffect, MultiOutcomeEffect):
    """Joint effect of two or more cis variants on the same
    transcript (#269).

    Emitted by :meth:`VariantCollection.effects` when a
    ``phase_resolver`` (VCF ``PS`` tag, Isovar assembly, or any
    other :class:`PhaseResolver`) groups cis variants together and
    :func:`varcode.mutant_transcript.apply_variants_to_transcript`
    can build a combined mutant cDNA. The per-variant effects stay
    on the collection — this is additive, not a replacement.

    Downstream consumers that want haplotype-level context (peptide
    prediction, ASE, compound het interpretation) read the
    :attr:`mutant_transcript` here. Consumers that want per-variant
    HGVS or annotation granularity iterate the per-variant effects
    as before.
    """

    def __init__(
            self, variants, transcript, mutant_transcript,
            phase_source=None):
        # The HaplotypeEffect is scoped to a set of variants — the
        # :class:`MutationEffect.variant` slot takes the first one
        # as a back-compat anchor so existing single-variant
        # consumers don't blow up, but real consumers read
        # :attr:`variants`.
        assert len(variants) >= 2, (
            "HaplotypeEffect requires ≥ 2 cis variants; got %d" % len(variants))
        TranscriptMutationEffect.__init__(self, variants[0], transcript)
        self.variants = tuple(variants)
        self.mutant_transcript = mutant_transcript
        # Which resolver produced the phase grouping (e.g. "vcf_ps"
        # for DNA-PS-tag phasing, "read_phasing" for an RNA / long-
        # read source). Useful for downstream filtering and for
        # auditing whether the cis call came from DNA-only phasing
        # or RNA-assembly evidence.
        self.phase_source = phase_source

    @property
    def candidates(self):
        """Single-outcome wrapping: the haplotype IS the outcome,
        lifted into the uniform :class:`EffectCandidate` shape so
        consumers iterate :attr:`candidates` the same way across all
        :class:`MultiOutcomeEffect` subclasses (#382). Extra
        candidates attached via :meth:`_combine_with_extra_candidates`
        (e.g. from RNA evidence) come after.
        """
        from ..effect_candidates import EffectCandidate
        evidence = {}
        if self.phase_source is not None:
            evidence["phase_source"] = self.phase_source
        base = (EffectCandidate(
            effect=self, source="varcode", evidence=evidence),)
        return self._combine_with_extra_candidates(base)

    @property
    def mutant_protein_sequence(self):
        """Joint protein translated from the haplotype's cDNA, or
        ``None`` if the edits don't land after the CDS start."""
        return getattr(
            self.mutant_transcript, "mutant_protein_sequence", None)

    def __str__(self):
        return "HaplotypeEffect(%d variants, transcript=%s, source=%s)" % (
            len(self.variants),
            self.transcript.name if self.transcript else "?",
            self.phase_source or "unknown")

    @property
    def short_description(self):
        """HGVS haplotype notation: ``[c.1A>T;c.5G>C]`` for cis.

        Per-variant short descriptions are synthesized from each
        variant's genomic coords when the per-variant effect isn't
        available here — the transcript-level HGVS would require
        re-running the per-variant classifier, which
        ``HaplotypeEffect`` deliberately doesn't do (that work lives
        on the per-variant effects that coexist in the collection).
        """
        parts = [v.short_description for v in self.variants]
        return "[" + ";".join(parts) + "]"


# =====================================================================
# Phase candidate set — possibility set when somatic + germline
# share a window and phase between them is unknown (#268).
#
# Sibling of HaplotypeEffect: that one captures the *known-cis* case
# (multiple variants composed on the same haplotype), this one
# captures the *unknown-phase* case (multiple hypotheses, each with
# its own resulting effect). Both inherit the same MultiOutcomeEffect
# protocol so consumers iterate ``outcomes`` uniformly. This isn't a
# wrapper around a reference-relative effect — it's the primary
# effect class for the unknown-phase case, just like
# SpliceOutcomeSet is for splice ambiguity.
# =====================================================================


class PhaseCandidateSet(TranscriptMutationEffect, MultiOutcomeEffect):
    """Possibility set across phase hypotheses when a somatic variant
    and one or more germline variants share a window on a transcript
    and phase between them is unknown.

    The somatic effect at this locus depends on which haplotype the
    somatic landed on, and without phase data we can't say which.
    The honest output is the set of plausible effects, one per
    hypothesis. ``most_likely`` returns the highest-priority
    candidate; :attr:`candidates` exposes the full set with
    per-hypothesis evidence keys (``phase_state``, ``haplotype``,
    ``germline_variants``) so consumers — including RNA-evidence
    integrations from #259 — can align across axes.

    Emitted by :func:`varcode.germline.predict_germline_aware_effect`.
    Sibling of :class:`HaplotypeEffect`: that one captures the
    *known-cis* multi-variant case; this one captures the
    *unknown-phase* possibility set.
    """

    def __init__(
            self,
            variant,
            transcript,
            candidates,
            hypotheses,
            germline_variants):
        TranscriptMutationEffect.__init__(self, variant, transcript)
        if len(candidates) != len(hypotheses):
            raise ValueError(
                "PhaseCandidateSet needs one candidate per "
                "hypothesis; got %d candidates and %d hypotheses." % (
                    len(candidates), len(hypotheses)))
        if not candidates:
            raise ValueError(
                "PhaseCandidateSet requires at least one candidate.")
        self._candidates_raw = tuple(candidates)
        self._hypotheses = tuple(hypotheses)
        self.germline_variants = tuple(germline_variants)

    @property
    def candidates(self):
        """One :class:`EffectCandidate` per hypothesis, sorted by
        underlying effect priority (most-severe first), carrying the
        phase metadata needed to align with RNA-evidence outcomes
        (#259, #382).

        ``evidence`` keys: ``phase_state`` (``"phased"`` /
        ``"implicit"`` / ``"unknown"`` / ``"too_many_hypotheses"``),
        ``haplotype`` (opaque tag), ``germline_variants`` (tuple of
        the cis germline variants on that hypothesis's haplotype).
        """
        from ..effect_candidates import EffectCandidate
        from .effect_ordering import effect_priority
        paired = sorted(
            zip(self._candidates_raw, self._hypotheses),
            key=lambda pair: -effect_priority(pair[0]))
        base = tuple(
            EffectCandidate(
                effect=candidate,
                source="varcode_germline",
                evidence={
                    "phase_state": hypothesis.phase_state,
                    "haplotype": hypothesis.haplotype,
                    "germline_variants": tuple(hypothesis.cis),
                })
            for candidate, hypothesis in paired)
        return self._combine_with_extra_candidates(base)

    # Note on accessor semantics for this class: the base-class
    # "producer order" contract is replaced here by
    # "sorted highest-impact-first" via
    # :func:`effect_priority`. Consequently
    # :attr:`most_likely_candidate` (== ``candidates[0]``) and
    # :attr:`highest_priority_candidate` coincide for this class —
    # both return the most protein-disruptive hypothesis, which is
    # also the conservative single-effect representation downstream
    # consumers see in ``short_description`` etc.

    @property
    def short_description(self):
        """``"?<most_likely_effect>"`` — the leading ``?`` flags the
        ambiguity. Consumers wanting the full possibility set read
        :attr:`candidates`."""
        return "?" + self.most_likely_effect.short_description

    @property
    def mutant_protein_sequence(self):
        """Most-likely candidate's mutant protein. Consumers iterating
        over hypotheses pull per-candidate sequences from
        :attr:`candidates`."""
        return getattr(self.most_likely_effect, "mutant_protein_sequence", None)

    @property
    def modifies_protein_sequence(self):
        return getattr(
            self.most_likely_effect, "modifies_protein_sequence", False)

    @property
    def modifies_coding_sequence(self):
        return getattr(
            self.most_likely_effect, "modifies_coding_sequence", False)

    def __str__(self):
        return (
            "PhaseCandidateSet(variant=%s, transcript=%s, "
            "%d hypotheses, germline=%s)" % (
                self.variant,
                self.transcript.name if self.transcript else "?",
                len(self._hypotheses),
                ",".join(g.short_description for g in self.germline_variants)))
