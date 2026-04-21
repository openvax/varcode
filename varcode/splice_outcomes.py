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

"""
Multi-effect candidate model for splice-disrupting variants — see
openvax/varcode#262.

When a variant disrupts a splice signal, the downstream protein
outcome is not deterministic from DNA alone. Instead of forcing a
single Effect, this module wraps splice effects in a
:class:`SpliceOutcomeSet` carrying multiple :class:`SpliceCandidate`
outcomes with rough plausibility scores.

Limitations of this prototype (documented openly so callers know what
the scores mean):

* Plausibility values are **hand-tuned heuristics**, not real
  probabilities. They are derived from rough literature consensus
  about which outcomes dominate at canonical donor/acceptor positions
  vs exonic splice sites vs flanking signals. Real scoring requires
  models like SpliceAI/SpliceTransformer or RNA evidence.
* :class:`SpliceOutcome.EXON_SKIPPING` and :class:`SpliceOutcome.NORMAL_SPLICING`
  produce concrete mutant protein sequences using the cDNA already
  available from PyEnsembl.
* :class:`SpliceOutcome.INTRON_RETENTION` and the cryptic-splice outcomes
  produce candidates with a predicted-class label
  (e.g. "likely PrematureStop") but no exact mutant protein. They
  require intron / flanking genomic sequence that PyEnsembl does
  not cache by default. A future PR with genomic-FASTA support
  would fill these in.
* This is a **prototype**. Once the foundational
  :class:`MutantTranscript` abstraction (#271) lands, this code will
  be re-expressed as a thin layer over it.
"""

from dataclasses import dataclass
from enum import Enum
from typing import Optional

from serializable import DataclassSerializable

from .effects.classify import classify_from_protein_diff
from .effects.codon_tables import codon_table_for_transcript, translate_sequence
from .effects.effect_classes import (
    ExonLoss,
    ExonicSpliceSite,
    Intronic,
    IntronicSpliceSite,
    MultiOutcomeEffect,
    MutationEffect,
    PredictedCrypticSpliceSite,
    PredictedIntronRetention,
    SpliceAcceptor,
    SpliceDonor,
)
from .mutant_transcript import MutantTranscript, TranscriptEdit
from .outcomes import Outcome


class SpliceOutcome(Enum):
    """Possible biological outcome when a splice signal is disrupted."""

    NORMAL_SPLICING = "normal_splicing"
    """The disruption is partial or leaky; canonical splicing still
    occurs and the variant's effect is whatever the underlying
    nucleotide change produces (substitution, frameshift, etc.)."""

    EXON_SKIPPING = "exon_skipping"
    """The spliceosome can't recognize the disrupted donor/acceptor
    and skips the entire affected exon. The flanking exons join
    directly. If the skipped exon is divisible by 3, the result is
    an in-frame deletion of those amino acids; otherwise a frameshift."""

    INTRON_RETENTION = "intron_retention"
    """The spliceosome fails entirely; the intron stays in the
    transcript. Almost always produces a premature stop codon
    within the intronic sequence."""

    CRYPTIC_DONOR = "cryptic_donor"
    """A nearby weak GT donor site is used instead of the disrupted
    canonical one. Truncates or extends the affected exon."""

    CRYPTIC_ACCEPTOR = "cryptic_acceptor"
    """A nearby weak AG acceptor site is used instead of the
    disrupted canonical one. Truncates or extends the affected exon."""


@dataclass(frozen=True)
class SpliceCandidate(DataclassSerializable):
    """One possible outcome of a splice-disrupting variant.

    Plausibility is a rough hand-tuned score, not a probability.
    Use the relative ordering of candidates rather than the absolute
    values for downstream filtering.

    Serialization is inherited from
    :class:`~serializable.DataclassSerializable`; ``to_json`` /
    ``from_json`` round-trip the ``coding_effect`` polymorphically
    via the standard ``__class__`` / ``__module__`` stamping from
    ``serializable.helpers``.
    """

    outcome: SpliceOutcome
    plausibility: float
    description: str
    """Human-readable summary of the outcome."""

    coding_effect: Optional[MutationEffect] = None
    """Concrete coding effect (Substitution, Deletion, FrameShift,
    PrematureStop, etc.) when computable. None for candidates whose
    protein math requires intron/flanking sequence we don't have
    cached (intron retention, cryptic splice)."""

    predicted_class_name: Optional[str] = None
    """When ``coding_effect`` is None, the predicted Effect class name
    (e.g. ``'PrematureStop'``) inferred from the outcome and exon
    properties. Useful for filtering downstream."""

    @property
    def has_protein(self) -> bool:
        return self.coding_effect is not None and getattr(
            self.coding_effect, "mutant_protein_sequence", None
        ) is not None

    @property
    def short_description(self) -> str:
        if self.coding_effect is not None:
            return "%s (%s, p=%.2f)" % (
                self.coding_effect.short_description,
                self.outcome.value,
                self.plausibility,
            )
        return "%s (%s, p=%.2f)" % (
            self.predicted_class_name or "?",
            self.outcome.value,
            self.plausibility,
        )

    def to_dict(self):
        """Emit the dataclass fields with :class:`SpliceOutcome`
        converted to its string value. ``python-serializable`` has no
        built-in enum branch in
        :func:`~serializable.helpers.to_serializable_repr`, so we
        stringify here to keep the JSON round-trip stable. Any
        future dataclass field typed as an :class:`Enum` will need
        the same treatment until the upstream helper grows an
        enum branch (#343).
        """
        d = super().to_dict()
        d["outcome"] = self.outcome.value
        return d

    @classmethod
    def from_dict(cls, state_dict):
        """Inverse of :meth:`to_dict`: reparse the string
        ``outcome`` into the :class:`SpliceOutcome` enum before
        delegating to the dataclass constructor (#343)."""
        state_dict = dict(state_dict)
        outcome = state_dict.get("outcome")
        if isinstance(outcome, str):
            state_dict["outcome"] = SpliceOutcome(outcome)
        return super().from_dict(state_dict)


class SpliceOutcomeSet(MultiOutcomeEffect):
    """A splice-disrupting variant's effect, expressed as a set of
    plausible outcomes rather than a single Effect.

    Implements the :class:`MultiOutcomeEffect` protocol: downstream
    consumers can write ``isinstance(effect, MultiOutcomeEffect)`` to
    catch this and any future multi-outcome wrapper (RNA evidence,
    germline-aware, etc.) uniformly — see #299.

    For back-compat, :attr:`short_description` delegates to the most
    plausible candidate. Callers that want all candidates iterate
    :attr:`candidates`.

    Constructed by :func:`enumerate_splice_outcomes` when a caller
    passes ``splice_outcomes=True`` to ``Variant.effects()`` or
    ``VariantCollection.effects()``.
    """

    def __init__(self, variant, transcript, candidates, disrupted_signal_class=None):
        MultiOutcomeEffect.__init__(self, variant)
        self.transcript = transcript
        # Sort candidates by plausibility descending so .candidates[0]
        # is always the most likely.
        self.candidates = tuple(sorted(
            candidates, key=lambda c: c.plausibility, reverse=True))
        # Class of the original splice-disrupting effect this set
        # replaced (SpliceDonor, SpliceAcceptor, ExonicSpliceSite, or
        # IntronicSpliceSite). Used for priority lookup.
        self.disrupted_signal_class = disrupted_signal_class

    # Serialization is inherited from Serializable (the base of
    # MutationEffect). to_dict emits {variant, transcript, candidates,
    # disrupted_signal_class} via introspection; to_serializable_repr
    # recursively handles nested MutationEffects, SpliceCandidate
    # dataclasses, and the disrupted_signal_class type object via
    # serializable.helpers' standard __class__ / __module__ stamping
    # (#343).

    @property
    def priority_class(self):
        """Delegate priority lookup to the disrupted-signal class so a
        SpliceOutcomeSet sorts as if it were the original splice effect.

        Read by :func:`varcode.effects.effect_priority`.
        """
        return self.disrupted_signal_class

    @property
    def most_likely(self) -> SpliceCandidate:
        return self.candidates[0]

    @property
    def outcomes(self):
        """Unified :class:`~varcode.Outcome` view over :attr:`candidates`
        (#339).

        Each outcome's ``effect`` is guaranteed to be a real
        :class:`MutationEffect` — for candidates whose protein math
        requires genomic sequence we don't cache (intron retention,
        cryptic splice), the effect is a
        :class:`PredictedIntronRetention` or
        :class:`PredictedCrypticSpliceSite` placeholder with
        ``mutant_protein_sequence = None`` instead of ``None`` itself.

        ``evidence`` carries the :class:`SpliceOutcome` enum under key
        ``"splice_outcome"`` so splice-aware consumers can still
        dispatch on the biological outcome kind without reaching
        back into :attr:`candidates`.
        """
        return tuple(
            Outcome(
                effect=self._outcome_effect(candidate),
                probability=candidate.plausibility,
                source="varcode",
                description=candidate.description,
                evidence={"splice_outcome": candidate.outcome})
            for candidate in self.candidates)

    def _outcome_effect(self, candidate):
        """Return the :class:`MutationEffect` that backs ``candidate``
        for :attr:`outcomes`. Falls back to placeholder classes when
        the candidate has no computable ``coding_effect``.
        """
        if candidate.coding_effect is not None:
            return candidate.coding_effect
        transcript = self.transcript
        if candidate.outcome is SpliceOutcome.INTRON_RETENTION:
            return PredictedIntronRetention(self.variant, transcript)
        if candidate.outcome is SpliceOutcome.CRYPTIC_DONOR:
            return PredictedCrypticSpliceSite(
                self.variant, transcript, direction="donor")
        if candidate.outcome is SpliceOutcome.CRYPTIC_ACCEPTOR:
            return PredictedCrypticSpliceSite(
                self.variant, transcript, direction="acceptor")
        if candidate.outcome is SpliceOutcome.EXON_SKIPPING:
            # Exon-skipping without a computed protein — we know an
            # exon is affected but the specific exon / skipped-protein
            # math wasn't resolvable at build time. Emit ExonLoss with
            # an empty exons tuple so ``.effect`` still satisfies the
            # MutationEffect interface; consumers that need the
            # specific exon should reach into ``candidates`` or
            # construct a SpliceOutcomeSet with RNA evidence attached.
            return ExonLoss(self.variant, transcript, exons=())
        # NORMAL_SPLICING candidate with no coding_effect: the variant
        # is intronic and has no coding change if splicing proceeds.
        # Represent with Intronic so .effect still satisfies the
        # MutationEffect contract.
        return Intronic(
            self.variant,
            transcript,
            nearest_exon=None,
            distance_to_exon=None)

    @property
    def alternate_effect(self):
        """Back-compat shim matching :attr:`ExonicSpliceSite.alternate_effect`.

        Resolves to the ``coding_effect`` of the ``NORMAL_SPLICING``
        candidate — i.e. "the coding change that applies when
        splicing proceeds normally." Returns ``None`` when the
        ``NORMAL_SPLICING`` candidate has no coding_effect (e.g. the
        variant is intronic and has no underlying coding change) or
        when the outcome set doesn't include ``NORMAL_SPLICING``.
        """
        for candidate in self.candidates:
            if candidate.outcome is SpliceOutcome.NORMAL_SPLICING:
                return candidate.coding_effect
        return None

    @property
    def candidate_proteins(self):
        """Mapping from each :class:`SpliceOutcome` to its computed
        mutant protein sequence (empty string when the protein is not
        computable from cDNA alone — i.e. intron retention and
        cryptic splice).
        """
        return {
            outcome.evidence["splice_outcome"]: _protein_str(outcome.effect)
            for outcome in self.outcomes
        }

    @property
    def mutant_protein_sequences(self):
        """Set of distinct non-empty mutant protein sequences across
        all candidates.
        """
        return {
            protein for protein in (
                _protein_str(o.effect) for o in self.outcomes)
            if protein
        }

    @property
    def short_description(self) -> str:
        return "splice-set:%s" % self.most_likely.short_description

    def __str__(self) -> str:
        return "SpliceOutcomeSet(variant=%s, transcript=%s, candidates=[%s])" % (
            self.variant,
            getattr(self.transcript, "name", None),
            ", ".join(c.short_description for c in self.candidates),
        )

    def __repr__(self) -> str:
        return str(self)


# ---------------------------------------------------------------------
# Plausibility tables.
#
# Hand-tuned heuristics. The numbers don't need to be precise — what
# matters is the relative ordering. Documented openly so reviewers and
# callers know what to expect. See module docstring for the limitation.
# ---------------------------------------------------------------------


# Disrupted canonical donor (intronic +1/+2): exon skipping dominates,
# intron retention common, leaky normal splicing rare.
_PLAUSIBILITY_SPLICE_DONOR = {
    SpliceOutcome.EXON_SKIPPING: 0.50,
    SpliceOutcome.INTRON_RETENTION: 0.30,
    SpliceOutcome.CRYPTIC_DONOR: 0.10,
    SpliceOutcome.NORMAL_SPLICING: 0.10,
}

# Disrupted canonical acceptor (intronic -1/-2): same rough
# distribution as donor.
_PLAUSIBILITY_SPLICE_ACCEPTOR = {
    SpliceOutcome.EXON_SKIPPING: 0.50,
    SpliceOutcome.INTRON_RETENTION: 0.30,
    SpliceOutcome.CRYPTIC_ACCEPTOR: 0.10,
    SpliceOutcome.NORMAL_SPLICING: 0.10,
}

# Disrupted exonic splice site (last 3 of exon): often splicing still
# proceeds (the disruption is in the exon, where the spliceosome can
# tolerate more variation), so normal splicing is more competitive.
_PLAUSIBILITY_EXONIC_SPLICE_SITE = {
    SpliceOutcome.NORMAL_SPLICING: 0.50,
    SpliceOutcome.EXON_SKIPPING: 0.30,
    SpliceOutcome.CRYPTIC_DONOR: 0.15,
    SpliceOutcome.INTRON_RETENTION: 0.05,
}

# Intronic splice site (positions +3 to +6 or -3): less critical
# region; normal splicing dominates. Two variants — donor-side (+3 to +6)
# and acceptor-side (-3) — which differ only in the cryptic direction.
_PLAUSIBILITY_INTRONIC_SPLICE_SITE_DONOR = {
    SpliceOutcome.NORMAL_SPLICING: 0.70,
    SpliceOutcome.EXON_SKIPPING: 0.20,
    SpliceOutcome.INTRON_RETENTION: 0.05,
    SpliceOutcome.CRYPTIC_DONOR: 0.05,
}

_PLAUSIBILITY_INTRONIC_SPLICE_SITE_ACCEPTOR = {
    SpliceOutcome.NORMAL_SPLICING: 0.70,
    SpliceOutcome.EXON_SKIPPING: 0.20,
    SpliceOutcome.INTRON_RETENTION: 0.05,
    SpliceOutcome.CRYPTIC_ACCEPTOR: 0.05,
}

# Order matters for isinstance matching: more specific classes first.
# SpliceDonor and SpliceAcceptor are subclasses of IntronicSpliceSite
# (see effect_classes.py), so they must be checked before the
# IntronicSpliceSite fallback in _plausibility_table_for — otherwise
# a SpliceDonor would incorrectly fall into the intronic-window table.
# IntronicSpliceSite itself is deliberately held out of this tuple and
# dispatched separately so side detection can pick the right cryptic
# direction.
_PLAUSIBILITY_TABLES = (
    (SpliceDonor, _PLAUSIBILITY_SPLICE_DONOR),
    (SpliceAcceptor, _PLAUSIBILITY_SPLICE_ACCEPTOR),
    (ExonicSpliceSite, _PLAUSIBILITY_EXONIC_SPLICE_SITE),
)


def _intronic_splice_side_is_acceptor(splice_effect):
    """True if an IntronicSpliceSite effect is on the acceptor side of
    its nearest exon.

    The classifier emits IntronicSpliceSite for positions +3–6 (donor
    side, after exon in transcript order) and -3 (acceptor side,
    before exon in transcript order). The side determines whether a
    cryptic donor or cryptic acceptor is the relevant alternative.

    Caller contract: ``splice_effect`` must be an ``IntronicSpliceSite``
    (guaranteed by :func:`_plausibility_table_for`), so ``variant`` and
    ``nearest_exon`` are always present.
    """
    exon = splice_effect.nearest_exon
    variant = splice_effect.variant
    if exon.strand == "+":
        return variant.trimmed_base1_start < exon.start
    # Reverse strand: acceptor-side intronic variants sit past the
    # genomic end of the exon (which is the 5' end in transcript order).
    return variant.trimmed_base1_end > exon.end


def _plausibility_table_for(splice_effect):
    """Return the plausibility table that applies to this splice
    effect, or None for effects we don't wrap.

    Uses :func:`isinstance` rather than exact-class dispatch so
    subclasses are handled correctly.
    """
    for cls, table in _PLAUSIBILITY_TABLES:
        if isinstance(splice_effect, cls):
            return table
    if isinstance(splice_effect, IntronicSpliceSite):
        if _intronic_splice_side_is_acceptor(splice_effect):
            return _PLAUSIBILITY_INTRONIC_SPLICE_SITE_ACCEPTOR
        return _PLAUSIBILITY_INTRONIC_SPLICE_SITE_DONOR
    return None


# ---------------------------------------------------------------------
# Outcome construction
# ---------------------------------------------------------------------


def _build_normal_splicing_candidate(splice_effect, plausibility):
    """Normal splicing: the underlying coding effect.

    For ExonicSpliceSite this is the .alternate_effect. For other
    splice classes there is no underlying coding effect (the variant
    is intronic) — return a candidate with no protein but the same
    plausibility for completeness.
    """
    underlying = getattr(splice_effect, "alternate_effect", None)
    if underlying is not None:
        return SpliceCandidate(
            outcome=SpliceOutcome.NORMAL_SPLICING,
            plausibility=plausibility,
            description=(
                "Disruption is partial or leaky; canonical splicing "
                "occurs and the underlying coding change applies."),
            coding_effect=underlying,
            predicted_class_name=type(underlying).__name__,
        )
    return SpliceCandidate(
        outcome=SpliceOutcome.NORMAL_SPLICING,
        plausibility=plausibility,
        description=(
            "Disruption is partial or leaky; canonical splicing "
            "occurs. Variant is intronic and has no coding impact "
            "if splicing proceeds normally."),
        coding_effect=None,
        predicted_class_name="Intronic",
    )


def _build_exon_skipping_candidate(splice_effect, plausibility):
    """Exon skipping: the affected exon is excluded from the transcript.

    Builds a :class:`MutantTranscript` by deleting the exon's cDNA
    range from the transcript sequence, translates, and classifies
    the effect via :func:`classify_from_protein_diff`. This replaces
    the three ad-hoc helpers from the original prototype (#305):
    boundary-codon reconstruction, start-loss detection, and
    in-frame vs. out-of-frame distinction all fall out of the diff
    naturally.
    """
    transcript = getattr(splice_effect, "transcript", None)
    exon = _affected_exon(splice_effect)
    if transcript is None or exon is None:
        return SpliceCandidate(
            outcome=SpliceOutcome.EXON_SKIPPING,
            plausibility=plausibility,
            description="Affected exon is skipped from the transcript.",
            coding_effect=None,
            predicted_class_name="ExonLoss",
        )

    mt = _build_exon_skip_mutant_transcript(
        splice_effect.variant, transcript, exon)
    exon_length = exon.end - exon.start + 1

    if mt is None or mt.mutant_protein_sequence is None:
        # Couldn't compute protein (incomplete transcript, exon not
        # found in transcript, etc.).
        predicted = "FrameShift" if exon_length % 3 != 0 else "Deletion"
        return SpliceCandidate(
            outcome=SpliceOutcome.EXON_SKIPPING,
            plausibility=plausibility,
            description="Exon %s is skipped." % getattr(exon, "exon_id", "?"),
            coding_effect=None,
            predicted_class_name=predicted,
        )

    ref_protein = str(transcript.protein_sequence)
    mut_protein = mt.mutant_protein_sequence
    coding_effect = classify_from_protein_diff(
        variant=splice_effect.variant,
        transcript=transcript,
        ref_protein=ref_protein,
        mut_protein=mut_protein,
        length_delta=-exon_length,
        mutant_transcript=mt,
    )
    predicted_class_name = type(coding_effect).__name__

    if exon_length % 3 == 0:
        description = (
            "Exon %s is skipped (in-frame, %d aa removed)." % (
                getattr(exon, "exon_id", "?"), exon_length // 3))
    else:
        description = (
            "Exon %s is skipped (out of frame, frameshift in the "
            "joined transcript)." % getattr(exon, "exon_id", "?"))

    return SpliceCandidate(
        outcome=SpliceOutcome.EXON_SKIPPING,
        plausibility=plausibility,
        description=description,
        coding_effect=coding_effect,
        predicted_class_name=predicted_class_name,
    )


def _build_intron_retention_candidate(splice_effect, plausibility):
    """Intron retention: stub candidate without exact protein.

    Concrete protein computation requires intron genomic sequence
    that PyEnsembl does not cache by default. We label the predicted
    class as PrematureStop (the typical outcome) but leave
    coding_effect None.
    """
    return SpliceCandidate(
        outcome=SpliceOutcome.INTRON_RETENTION,
        plausibility=plausibility,
        description=(
            "Intron is retained in the mature transcript. Almost "
            "always produces a premature stop codon within the "
            "intronic sequence; exact mutant protein requires "
            "intron genomic sequence not cached by PyEnsembl."),
        coding_effect=None,
        predicted_class_name="PrematureStop",
    )


def _build_cryptic_splice_candidate(splice_effect, plausibility, outcome):
    """Cryptic donor/acceptor: stub candidate without exact protein.

    Detecting the cryptic site requires scanning flanking genomic
    sequence. Stub for now; full implementation requires genomic
    FASTA support.
    """
    direction = "donor" if outcome is SpliceOutcome.CRYPTIC_DONOR else "acceptor"
    return SpliceCandidate(
        outcome=outcome,
        plausibility=plausibility,
        description=(
            "A cryptic %s site nearby may be used in place of the "
            "disrupted canonical signal. Truncates or extends the "
            "affected exon; exact mutant protein requires flanking "
            "genomic sequence not cached by PyEnsembl." % direction),
        coding_effect=None,
        predicted_class_name="ComplexSubstitution",
    )


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------


def _build_exon_skip_mutant_transcript(variant, transcript, exon):
    """Build a :class:`MutantTranscript` representing an exon skip.

    Deletes the exon's cDNA range from the transcript sequence and
    translates the result. Adjusts the CDS start offset when the
    skipped exon precedes the coding region (5'UTR exon). When the
    skipped exon contains the start codon, returns a MutantTranscript
    with an empty protein so the classifier emits StartLoss.
    """
    if not transcript.complete:
        return None
    try:
        exon_start_in_tx = _exon_start_offset_in_transcript(
            transcript, exon)
    except (StopIteration, ValueError):
        return None

    exon_length = exon.end - exon.start + 1
    full_sequence = str(transcript.sequence)

    post_skip_cdna = (
        full_sequence[:exon_start_in_tx]
        + full_sequence[exon_start_in_tx + exon_length:]
    )

    cds_start = min(transcript.start_codon_spliced_offsets)

    # Adjust CDS start for the deletion.
    if exon_start_in_tx + exon_length <= cds_start:
        # Exon is entirely in the 5' UTR → shift CDS start left.
        new_cds_start = cds_start - exon_length
    elif exon_start_in_tx >= cds_start:
        # Exon is entirely after the CDS start → no shift.
        new_cds_start = cds_start
    else:
        # Exon overlaps the CDS start → start codon lost.
        new_cds_start = None

    mut_protein = ""
    if new_cds_start is not None and new_cds_start < len(post_skip_cdna):
        codon_table = codon_table_for_transcript(transcript)
        coding = post_skip_cdna[new_cds_start:]
        truncated = coding[:(len(coding) // 3) * 3]
        try:
            mut_protein = translate_sequence(
                truncated, codon_table=codon_table, to_stop=True)
        except ValueError:
            mut_protein = ""

    edit = TranscriptEdit(
        cdna_start=exon_start_in_tx,
        cdna_end=exon_start_in_tx + exon_length,
        alt_bases="",
        source_variant=variant,
    )

    return MutantTranscript(
        reference_transcript=transcript,
        edits=(edit,),
        cdna_sequence=post_skip_cdna,
        mutant_protein_sequence=mut_protein,
        annotator_name="splice_outcomes",
    )


def _affected_exon(splice_effect):
    """Best-effort extraction of the affected exon from a splice effect.

    ExonicSpliceSite carries it directly. SpliceDonor/SpliceAcceptor
    carry the nearest exon. IntronicSpliceSite likewise.
    """
    if hasattr(splice_effect, "exon"):
        return splice_effect.exon
    if hasattr(splice_effect, "nearest_exon"):
        return splice_effect.nearest_exon
    return None


def _protein_str(effect):
    """Return ``effect.mutant_protein_sequence`` as a str, or ``""``
    if the effect carries no protein (placeholder stubs, start-loss,
    etc.). Centralizes the None/empty/str coercion that
    :attr:`SpliceOutcomeSet.candidate_proteins` used to inline.
    """
    protein = getattr(effect, "mutant_protein_sequence", None)
    return str(protein) if protein else ""


def _exon_start_offset_in_transcript(transcript, exon):
    """Offset (in transcript coordinates) of the first base of the
    given exon.
    """
    offset = 0
    for ex in transcript.exons:
        if ex.exon_id == exon.exon_id:
            return offset
        offset += ex.end - ex.start + 1
    raise ValueError("Exon %s not found in transcript %s" % (
        exon.exon_id, transcript.transcript_id))


# ---------------------------------------------------------------------
# Public entry point: enumerate outcomes for a splice effect
# ---------------------------------------------------------------------


def enumerate_splice_outcomes(splice_effect):
    """Wrap a splice-disrupting Effect in a SpliceOutcomeSet.

    Recognized splice effect classes are SpliceDonor, SpliceAcceptor,
    ExonicSpliceSite, and IntronicSpliceSite. Unrecognized classes
    pass through unchanged.

    Parameters
    ----------
    splice_effect : MutationEffect
        Output of varcode's existing splice classification.

    Returns
    -------
    SpliceOutcomeSet or the original effect
        SpliceOutcomeSet wrapping the candidate outcomes when the
        input is a splice-disrupting effect; otherwise the input
        unchanged.
    """
    table = _plausibility_table_for(splice_effect)
    if table is None:
        return splice_effect

    candidates = []
    for outcome, plausibility in table.items():
        if outcome is SpliceOutcome.NORMAL_SPLICING:
            candidates.append(_build_normal_splicing_candidate(
                splice_effect, plausibility))
        elif outcome is SpliceOutcome.EXON_SKIPPING:
            candidates.append(_build_exon_skipping_candidate(
                splice_effect, plausibility))
        elif outcome is SpliceOutcome.INTRON_RETENTION:
            candidates.append(_build_intron_retention_candidate(
                splice_effect, plausibility))
        elif outcome in (
                SpliceOutcome.CRYPTIC_DONOR,
                SpliceOutcome.CRYPTIC_ACCEPTOR):
            candidates.append(_build_cryptic_splice_candidate(
                splice_effect, plausibility, outcome))

    return SpliceOutcomeSet(
        variant=splice_effect.variant,
        transcript=getattr(splice_effect, "transcript", None),
        candidates=candidates,
        disrupted_signal_class=type(splice_effect),
    )


def wrap_splice_effects_in_collection(effect_collection):
    """Apply :func:`enumerate_splice_outcomes` to every splice-related
    Effect in an EffectCollection. Non-splice effects pass through.
    """
    new_effects = [enumerate_splice_outcomes(e) for e in effect_collection]
    return effect_collection.clone_with_new_elements(new_effects)


# Priority integration happens via the `priority_class` attribute on
# SpliceOutcomeSet — see `effect_priority` in effects.effect_ordering.

# Serialization is inherited from Serializable / DataclassSerializable
# (see #343). serializable.helpers.to_serializable_repr stamps each
# nested object's {__module__, __name__} and from_serializable_dict
# resolves them via _lookup_value, so SpliceCandidate.coding_effect
# (polymorphic MutationEffect) and SpliceOutcomeSet.disrupted_signal_class
# (a class object) both round-trip without the hand-rolled registries
# this module used to carry.
