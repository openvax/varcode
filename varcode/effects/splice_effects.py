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

"""Splice-signal disruption and splice-mechanism consequence classes."""

from __future__ import annotations

from typing import TYPE_CHECKING

from .base import Exonic, Intronic, MultiOutcomeEffect, TranscriptMutationEffect

if TYPE_CHECKING:
    from ..effect_candidates import EffectCandidate
    from ..mutant_transcript import MutantTranscript


class SpliceSite(object):
    """Marker base for the four splice-signal *disruption* effects:
    :class:`SpliceDonor`, :class:`SpliceAcceptor`,
    :class:`IntronicSpliceSite`, and :class:`ExonicSpliceSite`.

    A ``SpliceSite`` effect answers **where** a variant hit a splice
    signal — a position fact (which boundary, how close). It does
    *not* describe the protein consequence: from DNA alone you can't
    say how the spliceosome will respond. That downstream protein-level
    outcome is modeled separately by :class:`SpliceMechanismEffect`
    subclasses (exon skipping, intron retention, cryptic-site use,
    normal splicing).

    So the two halves of splice handling are:

    * ``SpliceSite``           → *where* the signal was disrupted.
    * ``SpliceMechanismEffect`` → *what* the spliceosome does, and the
      protein change that results.

    :func:`varcode.splice_outcomes.enumerate_splice_outcomes` bridges
    them: given a ``SpliceSite`` disruption it produces a
    :class:`~varcode.splice_outcomes.SpliceOutcomeSet` of plausible
    mechanisms.
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

    @property
    def short_description(self) -> str:
        return "intronic-splice-site"


class SpliceDonor(IntronicSpliceSite):
    """
    Mutation in the first two intron residues.
    """
    def __init__(self, variant, transcript, nearest_exon, distance_to_exon):
        IntronicSpliceSite.__init__(
            self, variant, transcript, nearest_exon, distance_to_exon)

    @property
    def short_description(self) -> str:
        return "splice-donor"


class SpliceAcceptor(IntronicSpliceSite):
    """
    Mutation in the last two intron residues.
    """

    @property
    def short_description(self) -> str:
        return "splice-acceptor"


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
# disruption via :attr:`splice_signal` — a :class:`SpliceSite` effect
# (SpliceDonor / SpliceAcceptor / IntronicSpliceSite / ExonicSpliceSite)
# describing *where* in the transcript the splice signal was hit. The
# mechanism Effect describes *what* happens downstream of that
# disruption.
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
    splice_signal : SpliceSite
        The splice-signal-disruption Effect *instance* describing
        *where* the disruption was — a :class:`SpliceSite` (one of
        ``SpliceDonor`` / ``SpliceAcceptor`` / ``IntronicSpliceSite`` /
        ``ExonicSpliceSite``). Lets consumers reach
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

    def _get_modifies_coding_sequence(self) -> bool | None:
        return getattr(self.protein_effect, "modifies_coding_sequence", False)

    def _get_modifies_protein_sequence(self) -> bool | None:
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
    def mutant_protein_sequence(self) -> str | None:
        return getattr(self.coding_effect, "mutant_protein_sequence", None)

    @mutant_protein_sequence.setter
    def mutant_protein_sequence(self, value: str | None) -> None:
        raise AttributeError("NormalSplicing.mutant_protein_sequence is read-only")

    @property
    def aa_mutation_start_offset(self) -> int | None:
        return getattr(self.coding_effect, "aa_mutation_start_offset", None)

    @aa_mutation_start_offset.setter
    def aa_mutation_start_offset(self, value: int | None) -> None:
        raise AttributeError("NormalSplicing.aa_mutation_start_offset is read-only")

    @property
    def aa_mutation_end_offset(self) -> int | None:
        return getattr(self.coding_effect, "aa_mutation_end_offset", None)

    @aa_mutation_end_offset.setter
    def aa_mutation_end_offset(self, value: int | None) -> None:
        raise AttributeError("NormalSplicing.aa_mutation_end_offset is read-only")

    @property
    def mutant_transcript(self) -> MutantTranscript | None:
        return getattr(self.coding_effect, "mutant_transcript", None)

    @mutant_transcript.setter
    def mutant_transcript(self, value: MutantTranscript | None) -> None:
        raise AttributeError("NormalSplicing.mutant_transcript is read-only")

    @property
    def priority_class(self):
        if self.coding_effect is None:
            return None
        return self.coding_effect.__class__

    def _get_modifies_coding_sequence(self) -> bool | None:
        return getattr(self.coding_effect, "modifies_coding_sequence", False)

    def _get_modifies_protein_sequence(self) -> bool | None:
        return getattr(self.coding_effect, "modifies_protein_sequence", False)

    @property
    def short_description(self):
        if self.coding_effect is None:
            return "normal-splicing"
        return "normal-splicing (%s)" % self.coding_effect.short_description

    def _format_short_description(self, deletion_notation):
        if self.coding_effect is None:
            return self.short_description
        return "normal-splicing (%s)" % self.coding_effect.format_short_description(
            deletion_notation=deletion_notation)


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

    direction: str | None = None  # "donor" / "acceptor", set by subclass

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

    @property
    def short_description(self) -> str:
        return "exon-loss"

    _default_modifies_protein_sequence = True
    _default_modifies_coding_sequence = True


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

    Note: as of varcode 6.0, ``ExonicSpliceSite`` no longer surfaces
    as the top-level effect at the user-facing API — splice-disrupting
    variants are always wrapped in
    :class:`~varcode.splice_outcomes.SpliceOutcomeSet`. ``ExonicSpliceSite``
    survives as ``splice_set.disrupted_signal_class`` (a type) and as
    the ``splice_signal`` reference on each candidate (an instance).
    """
    def __init__(self, variant, transcript, exon, alternate_effect):
        Exonic.__init__(self, variant, transcript)
        self.exon = exon
        self.alternate_effect = alternate_effect

    def __str__(self):
        return "ExonicSpliceSite(exon=%s, alternate_effect=%s)" % (
            self.exon,
            self.alternate_effect)

    @property
    def short_description(self) -> str:
        return "exonic-splice-site"

    def _get_mutant_protein_sequence(self) -> str | None:
        """
        TODO: determine when exonic splice variants cause exon skipping
        vs. translation of the underlying modified coding sequence.

        For now just pretending like there is no effect on splicing.
        """
        return self.alternate_effect.mutant_protein_sequence

    def _get_modifies_protein_sequence(self) -> bool | None:
        return self.alternate_effect.modifies_protein_sequence

    def _get_modifies_coding_sequence(self) -> bool | None:
        return self.alternate_effect.modifies_coding_sequence

    def _get_candidates(self) -> tuple[EffectCandidate, ...]:
        """The two plausible outcomes wrapped as
        :class:`~varcode.effect_candidates.EffectCandidate`:

        Position [0] is the splice-disruption outcome (this effect
        itself); position [1] is the coding change if splicing
        proceeds. Returns a 1-tuple when ``alternate_effect`` is None.
        Additional candidates appended via :meth:`add_candidates`
        (e.g. from RNA evidence) follow automatically.
        """
        from ..effect_candidates import EffectCandidate
        base = [EffectCandidate(effect=self, source="varcode")]
        if self.alternate_effect is not None:
            base.append(EffectCandidate(
                effect=self.alternate_effect, source="varcode"))
        return tuple(base)

    @property
    def priority_class(self):
        """ExonicSpliceSite priority is used directly — no
        delegation, since this class IS the splice-adjacent effect
        rather than a wrapper around one.
        """
        return ExonicSpliceSite



__all__ = [
    "SpliceSite",
    "IntronicSpliceSite",
    "SpliceDonor",
    "SpliceAcceptor",
    "SpliceMechanismEffect",
    "NormalSplicing",
    "ExonSkipping",
    "IntronRetention",
    "CrypticSpliceSiteEffect",
    "CrypticDonor",
    "CrypticAcceptor",
    "ExonLoss",
    "ExonicSpliceSite",
]
