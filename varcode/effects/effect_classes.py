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


class MultiOutcomeEffect(MutationEffect):
    """Marker base class for effects that represent a set of plausible
    outcomes rather than a single deterministic effect.

    Subclasses must expose:

    * ``candidates`` — sequence of :class:`MutationEffect` instances,
      sorted most-plausible-first. (Kept for back-compat with 2.x
      callers.)
    * ``most_likely`` — the top candidate (i.e. ``candidates[0]``).
    * ``priority_class`` — effect class whose priority this set adopts
      (read by :func:`varcode.effects.effect_priority`).

    **Harmonized interface (#299):** new code should read
    :attr:`outcomes` instead of ``candidates``. Each entry is an
    :class:`~varcode.outcomes.Outcome` carrying the effect plus
    provenance (probability, source, evidence dict). The default
    implementation wraps ``candidates`` with ``source="varcode"``
    and no probability — external scorers (SpliceAI, Pangolin),
    RNA-evidence callers (Isovar), and long-read assembly tools
    override to attach their own scores without subclassing.

    Downstream consumers filter for multi-outcome results with
    ``isinstance(effect, MultiOutcomeEffect)``, so new wrappers (RNA
    evidence #259, germline-aware #268, SV-at-breakpoint) implement
    the same protocol without downstream code churn.
    """

    @property
    def outcomes(self):
        """Tuple of :class:`~varcode.outcomes.Outcome` objects,
        most-plausible-first. Default implementation wraps
        :attr:`candidates` under ``source="varcode"``; subclasses
        (or external integrations) override to attach probabilities
        and evidence.
        """
        from ..outcomes import outcomes_from_candidates
        return outcomes_from_candidates(self.candidates)


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
# Predicted-but-uncomputed placeholder effects (#339).
#
# Used by :class:`~varcode.splice_outcomes.SpliceOutcomeSet` and, later,
# :class:`StructuralVariantEffect` to fill :attr:`Outcome.effect` with a
# real :class:`MutationEffect` when the protein-level outcome cannot be
# computed from cached transcript cDNA alone (intron retention requires
# intron sequence; cryptic splice requires flanking genomic sequence).
# These types carry no aa_ref / aa_alt / mutant_protein_sequence — they
# are markers so that ``outcome.effect`` still satisfies the
# MutationEffect interface and downstream iteration doesn't have to
# branch on None.
# =====================================================================


class PredictedIntronRetention(TranscriptMutationEffect):
    """Placeholder effect: intron retention predicted, exact protein
    not computable from cached transcript cDNA.

    Emitted as the :attr:`Outcome.effect` of a SpliceOutcomeSet's
    ``INTRON_RETENTION`` outcome. The biologically expected outcome is
    a premature stop codon inside the retained intron; consumers that
    need the exact protein sequence require intron genomic sequence
    (see #296).
    """

    short_description = "predicted-intron-retention"


class PredictedCrypticSpliceSite(TranscriptMutationEffect):
    """Placeholder effect: a cryptic donor or acceptor is expected to
    be used in place of the disrupted canonical signal.

    ``direction`` is ``"donor"`` or ``"acceptor"``. Exact protein
    consequence requires flanking genomic sequence; emitted as the
    :attr:`Outcome.effect` of a SpliceOutcomeSet's ``CRYPTIC_DONOR`` /
    ``CRYPTIC_ACCEPTOR`` outcome.
    """

    def __init__(self, variant, transcript, direction):
        TranscriptMutationEffect.__init__(self, variant, transcript)
        if direction not in ("donor", "acceptor"):
            raise ValueError(
                "PredictedCrypticSpliceSite direction must be "
                "'donor' or 'acceptor', got %r" % (direction,))
        self.direction = direction

    @property
    def short_description(self):
        return "predicted-cryptic-%s" % self.direction


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
    :class:`MultiOutcomeEffect` protocol (see #299). ``candidates``
    yields ``(self, alternate_effect)`` (both :class:`MutationEffect`
    instances), and ``alternate_effect`` stays on the instance as a
    first-class field for back-compat with callers that depended on
    the 2.x shape.

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
        """The two plausible outcomes, most-likely-first:

        ``(self, alternate_effect)`` — position [0] is the
        splice-disruption outcome (this effect itself), position [1]
        is the coding change if splicing proceeds. Returns a
        1-tuple ``(self,)`` when ``alternate_effect`` is None.

        Elements are :class:`MutationEffect` instances directly,
        not :class:`SpliceCandidate` objects — that's a different
        (richer) shape used by :class:`SpliceOutcomeSet`. Downstream
        code that iterates ``candidates`` here reads the effect
        class / ``short_description`` / ``aa_ref`` etc. on each
        element, same as on any ``MutationEffect``.
        """
        if self.alternate_effect is None:
            return (self,)
        return (self, self.alternate_effect)

    @property
    def most_likely(self):
        """The splice-disruption outcome is the primary classification
        (this effect itself). Follows the ``MultiOutcomeEffect``
        contract of ``most_likely == candidates[0]``.
        """
        return self

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
# subclass exposing an ``outcomes`` tuple per #299.
#
# The classes are deliberately thin wrappers — the important interface
# is ``outcomes``, which carries the ``Outcome`` objects with
# probability / source / evidence so external tools (RNA evidence,
# SpliceAI, long-read assembly) can score them without subclassing.
# =====================================================================


class StructuralVariantEffect(TranscriptMutationEffect, MultiOutcomeEffect):
    """Base class for effects of a :class:`StructuralVariant` on a
    specific transcript. Subclasses set ``short_description``; the
    ``outcomes`` property is what downstream callers read.

    The ``candidates`` tuple (back-compat 2.x shape) returns
    ``(self,)`` when no specific alternates have been nominated, or
    a richer tuple when the subclass carries explicit alternate
    effects. ``MultiOutcomeEffect.outcomes`` lifts that to the
    unified :class:`Outcome` shape automatically.
    """

    def __init__(
            self, variant, transcript, candidates=None, mutant_transcript=None):
        TranscriptMutationEffect.__init__(self, variant, transcript)
        self._candidates = (
            tuple(candidates) if candidates is not None else (self,))
        self.mutant_transcript = mutant_transcript

    @property
    def candidates(self):
        return self._candidates

    @property
    def most_likely(self):
        return self._candidates[0]

    @property
    def outcomes(self):
        """Unified :class:`~varcode.Outcome` view over
        :attr:`candidates` (#339).

        Each outcome's ``effect`` is always a :class:`MutationEffect`
        (SV candidates are MutationEffect instances by construction —
        the default ``(self,)`` is already correct). ``evidence``
        carries the SV-type tag so consumers can filter by
        ``DEL`` / ``DUP`` / ``INV`` / ``BND`` without re-reading the
        variant.
        """
        from ..outcomes import Outcome
        sv_type = getattr(self.variant, "sv_type", None)
        evidence = {"sv_type": sv_type} if sv_type is not None else {}
        return tuple(
            Outcome(
                effect=candidate,
                source="varcode",
                description=getattr(candidate, "short_description", None),
                evidence=evidence)
            for candidate in self._candidates)


class LargeDeletion(StructuralVariantEffect):
    """A deletion (``<DEL>`` / ``<CN0>``) that removes one or more
    exons — or an entire gene. Carries the list of affected exons
    for downstream analysis; the single ``outcome`` is this effect
    itself (callers that want to add RNA evidence or SpliceAI
    scoring construct additional outcomes and wrap the result)."""

    short_description = "sv-deletion"

    def __init__(
            self, variant, transcript, affected_exons,
            candidates=None, mutant_transcript=None):
        StructuralVariantEffect.__init__(
            self, variant, transcript,
            candidates=candidates, mutant_transcript=mutant_transcript)
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
            candidates=None, mutant_transcript=None):
        StructuralVariantEffect.__init__(
            self, variant, transcript,
            candidates=candidates, mutant_transcript=mutant_transcript)
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
    left to downstream tools that attach :class:`Outcome` objects
    with ``source="isovar"`` / ``"longread_assembly"`` / etc.
    """

    short_description = "sv-gene-fusion"

    def __init__(
            self, variant, transcript, partner_transcript,
            mutant_transcript=None, candidates=None):
        StructuralVariantEffect.__init__(
            self, variant, transcript,
            candidates=candidates, mutant_transcript=mutant_transcript)
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
    enumerator; attached as additional :class:`Outcome` entries on
    SV effects rather than standalone.

    Not a :class:`TranscriptMutationEffect` because the candidate
    region may not overlap any existing transcript — it's a *new*
    exon hypothesis. Carries the contig / interval and the motif
    scores as plain fields; external predictors (SpliceAI,
    Pangolin) attach their own scores via the enclosing
    :class:`Outcome.evidence` dict.
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
