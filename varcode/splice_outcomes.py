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

from .effects.effect_classes import (
    Deletion,
    ExonicSpliceSite,
    IntronicSpliceSite,
    MutationEffect,
    SpliceAcceptor,
    SpliceDonor,
    StartLoss,
)


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
class SpliceCandidate:
    """One possible outcome of a splice-disrupting variant.

    Plausibility is a rough hand-tuned score, not a probability.
    Use the relative ordering of candidates rather than the absolute
    values for downstream filtering.
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


class SpliceOutcomeSet(MutationEffect):
    """A splice-disrupting variant's effect, expressed as a set of
    plausible outcomes rather than a single Effect.

    For back-compat, :attr:`short_description` delegates to the most
    plausible candidate. Callers that want all candidates iterate
    :attr:`candidates`.

    Constructed by :func:`enumerate_splice_outcomes` when a caller
    passes ``splice_outcomes=True`` to ``Variant.effects()`` or
    ``VariantCollection.effects()``.
    """

    def __init__(self, variant, transcript, candidates, disrupted_signal_class=None):
        MutationEffect.__init__(self, variant)
        self.transcript = transcript
        # Sort candidates by plausibility descending so .candidates[0]
        # is always the most likely.
        self.candidates = tuple(sorted(
            candidates, key=lambda c: c.plausibility, reverse=True))
        # Class of the original splice-disrupting effect this set
        # replaced (SpliceDonor, SpliceAcceptor, ExonicSpliceSite, or
        # IntronicSpliceSite). Used for priority lookup.
        self.disrupted_signal_class = disrupted_signal_class

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
    def candidate_proteins(self):
        """Mapping from each :class:`SpliceOutcome` to its computed
        mutant protein sequence (empty string when the protein is not
        computable from cDNA alone — i.e. intron retention and
        cryptic splice).

        Returns
        -------
        dict[SpliceOutcome, str]
            One entry per candidate. The string is empty for stubs.
        """
        result = {}
        for candidate in self.candidates:
            coding = candidate.coding_effect
            if coding is not None:
                protein = getattr(coding, "mutant_protein_sequence", "")
                result[candidate.outcome] = str(protein) if protein else ""
            else:
                result[candidate.outcome] = ""
        return result

    @property
    def mutant_protein_sequences(self):
        """Set of distinct non-empty mutant protein sequences across
        all candidates.

        Useful for downstream consumers (neoantigen prediction,
        isovar/vaxrank) that want to enumerate every protein this
        splice disruption might produce.
        """
        return {
            p for p in self.candidate_proteins.values() if p
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
# Today all four are leaf classes so ordering doesn't matter in
# practice, but the iteration pattern in _plausibility_table_for()
# is future-proof against subclassing.
_PLAUSIBILITY_TABLES = (
    (SpliceDonor, _PLAUSIBILITY_SPLICE_DONOR),
    (SpliceAcceptor, _PLAUSIBILITY_SPLICE_ACCEPTOR),
    (ExonicSpliceSite, _PLAUSIBILITY_EXONIC_SPLICE_SITE),
    # IntronicSpliceSite requires side detection to pick the right
    # cryptic direction; handled in _plausibility_table_for().
)


def _intronic_splice_side_is_acceptor(splice_effect):
    """True if an IntronicSpliceSite effect is on the acceptor side of
    its nearest exon.

    The classifier emits IntronicSpliceSite for positions +3–6 (donor
    side, after exon in transcript order) and -3 (acceptor side,
    before exon in transcript order). The side determines whether a
    cryptic donor or cryptic acceptor is the relevant alternative.
    """
    exon = getattr(splice_effect, "nearest_exon", None)
    variant = getattr(splice_effect, "variant", None)
    if exon is None or variant is None:
        return False
    strand = getattr(exon, "strand", "+")
    variant_start = variant.trimmed_base1_start
    variant_end = variant.trimmed_base1_end
    if strand == "+":
        return variant_start < exon.start
    # Reverse strand: acceptor-side intronic variants sit past the
    # genomic end of the exon (which is the 5' end in transcript order).
    return variant_end > exon.end


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

    Computes the resulting protein by removing the exon's amino acids
    (in-frame skip) or by triggering a frameshift (out-of-frame skip).
    Reports the change as a Deletion when in-frame and as a frameshift
    label otherwise.
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

    exon_contains_start_codon = _exon_contains_start_codon(transcript, exon)
    exon_length = exon.end - exon.start + 1
    if exon_contains_start_codon:
        # Losing the exon that contains the start codon means the
        # protein as annotated can't be produced; label this as
        # StartLoss. We can't easily construct the StartLoss effect
        # here without more context, so return a stub with the right
        # predicted class.
        coding_effect = _build_start_loss_effect(
            splice_effect.variant, transcript)
        predicted_class_name = "StartLoss"
        description = (
            "Exon %s is skipped and contains the start codon; "
            "annotated translation initiation is lost." % (
                getattr(exon, "exon_id", "?"),))
    elif exon_length % 3 == 0:
        coding_effect = _build_in_frame_exon_skip_effect(
            splice_effect.variant, transcript, exon)
        predicted_class_name = "Deletion"
        description = (
            "Exon %s is skipped (in-frame, %d aa removed)." % (
                getattr(exon, "exon_id", "?"), exon_length // 3))
    else:
        coding_effect = _build_out_of_frame_exon_skip_effect(
            splice_effect.variant, transcript, exon)
        predicted_class_name = "FrameShift"
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


def _exon_contains_start_codon(transcript, exon):
    """True if the given exon contains the annotated start codon."""
    if not transcript.complete:
        return False
    try:
        start_offsets = transcript.start_codon_spliced_offsets
    except Exception:
        return False
    if not start_offsets:
        return False
    exon_start = _exon_start_offset_in_transcript_or_none(transcript, exon)
    if exon_start is None:
        return False
    exon_length = exon.end - exon.start + 1
    exon_end = exon_start + exon_length - 1
    first_start = min(start_offsets)
    return exon_start <= first_start <= exon_end


def _exon_start_offset_in_transcript_or_none(transcript, exon):
    """Like :func:`_exon_start_offset_in_transcript` but returns None
    on failure instead of raising."""
    try:
        return _exon_start_offset_in_transcript(transcript, exon)
    except (StopIteration, ValueError):
        return None


def _build_start_loss_effect(variant, transcript):
    """Build a StartLoss effect stub when the start-codon exon is
    skipped. Returns None if the transcript can't accommodate it.
    """
    if not transcript.complete:
        return None
    try:
        return StartLoss(variant=variant, transcript=transcript)
    except Exception:
        return None


def _build_out_of_frame_exon_skip_effect(variant, transcript, skipped_exon):
    """Compute the mutant protein for an out-of-frame exon skip.

    Builds the post-skip cDNA by joining all exons except the skipped
    one, translating from the original start codon, and following the
    frameshift until the first stop. Returns a Deletion-style effect
    with aa_ref set to the joined skipped region; leaves the mutant
    protein accessible via a custom subclass-like shim.
    """
    if not transcript.complete:
        return None
    try:
        exon_start_in_tx = _exon_start_offset_in_transcript(
            transcript, skipped_exon)
    except (StopIteration, ValueError):
        return None
    cds_start_offset = min(transcript.start_codon_spliced_offsets)
    if exon_start_in_tx < cds_start_offset:
        return None
    exon_length = skipped_exon.end - skipped_exon.start + 1
    # Construct the post-skip cDNA by excising the exon's nucleotides
    # from the full transcript sequence.
    full_sequence = str(transcript.sequence)
    post_skip_cdna = (
        full_sequence[:exon_start_in_tx]
        + full_sequence[exon_start_in_tx + exon_length:]
    )
    # Translate from the start codon through the frameshift to first stop.
    coding_from_start = post_skip_cdna[cds_start_offset:]
    protein = _translate_to_first_stop(coding_from_start)
    if not protein:
        return None
    # Compute the aa position where the frameshift begins — this is
    # where the exon used to start, in aa coordinates.
    aa_frameshift_start = (exon_start_in_tx - cds_start_offset) // 3
    if aa_frameshift_start >= len(transcript.protein_sequence):
        return None
    return _ExonSkipFrameshiftEffect(
        variant=variant,
        transcript=transcript,
        aa_frameshift_start=aa_frameshift_start,
        protein=protein,
    )


def _translate_to_first_stop(cdna):
    """Translate a cDNA string to protein, stopping at the first stop
    codon. Returns the protein string without the stop symbol.
    """
    from Bio.Seq import Seq
    n_codons = len(cdna) // 3
    truncated = str(cdna[:n_codons * 3])
    protein = str(Seq(truncated).translate(to_stop=True))
    return protein


class _ExonSkipFrameshiftEffect(MutationEffect):
    """Internal helper effect representing an exon-skip-induced
    frameshift. Carries the computed mutant protein sequence but isn't
    intended as a public effect class — it's wrapped inside the
    SpliceCandidate.
    """
    def __init__(self, variant, transcript, aa_frameshift_start, protein):
        MutationEffect.__init__(self, variant)
        self.transcript = transcript
        self.aa_mutation_start_offset = aa_frameshift_start
        self.aa_ref = str(
            transcript.protein_sequence[aa_frameshift_start:])
        # aa_alt is the new amino acids after the frameshift point.
        self.aa_alt = protein[aa_frameshift_start:]
        self.mutant_protein_sequence = protein

    @property
    def short_description(self):
        return "p.%s%dfs*%d" % (
            self.aa_ref[:1] if self.aa_ref else "?",
            self.aa_mutation_start_offset + 1,
            len(self.aa_alt),
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


def _build_in_frame_exon_skip_effect(variant, transcript, skipped_exon):
    """Construct a Deletion effect representing an in-frame exon skip.

    Computes the amino acid range corresponding to the skipped exon
    and emits a Deletion. Falls back to None if the math doesn't
    work out (e.g. exon spans the start codon).
    """
    if not transcript.complete:
        return None
    try:
        exon_start_in_tx = _exon_start_offset_in_transcript(
            transcript, skipped_exon)
    except (StopIteration, ValueError):
        return None
    cds_start_offset = min(transcript.start_codon_spliced_offsets)
    if exon_start_in_tx < cds_start_offset:
        # Exon overlaps the 5' UTR or start codon; the simple
        # amino-acid math doesn't apply.
        return None
    aa_start = (exon_start_in_tx - cds_start_offset) // 3
    exon_length = skipped_exon.end - skipped_exon.start + 1
    n_aa_removed = exon_length // 3
    aa_end = aa_start + n_aa_removed
    if aa_end > len(transcript.protein_sequence):
        return None
    aa_ref = str(transcript.protein_sequence[aa_start:aa_end])
    if not aa_ref:
        return None
    return Deletion(
        variant=variant,
        transcript=transcript,
        aa_mutation_start_offset=aa_start,
        aa_ref=aa_ref,
    )


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
