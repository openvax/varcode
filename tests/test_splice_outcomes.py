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
Tests for the splice outcome possibility-set prototype
(openvax/varcode#262).

Coverage:
  - Default behavior unchanged (back-compat)
  - Opt-in wraps splice effects in SpliceOutcomeSet
  - Each canonical splice signal class produces the expected outcome set
  - Plausibility ordering is stable
  - Per-outcome candidate construction (normal splicing, exon
    skipping, intron retention stubs, cryptic splice stubs)
  - SpliceOutcomeSet integrates with EffectCollection
  - Multi-allelic and reverse-strand variants work too
"""

import pytest
from pyensembl import cached_release

import varcode
from varcode import (
    SpliceCandidate,
    SpliceOutcome,
    SpliceOutcomeSet,
    Variant,
)
from varcode.effects import (
    ExonicSpliceSite,
    IntronicSpliceSite,
    SpliceAcceptor,
    SpliceDonor,
)
from varcode.splice_outcomes import enumerate_splice_outcomes


ensembl_grch38 = cached_release(81)
CFTR_TRANSCRIPT_ID = "ENST00000003084"
BRCA1_TRANSCRIPT_ID = "ENST00000357654"


# --------------------------------------------------------------------
# Back-compat
# --------------------------------------------------------------------


def test_default_behavior_unchanged():
    # No kwarg -> same as today: SpliceDonor effect, no wrapping.
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effect = variant.effect_on_transcript(transcript)
    assert effect.__class__ is SpliceDonor
    assert not isinstance(effect, SpliceOutcomeSet)


def test_default_for_collection_unchanged():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    effects = variant.effects()  # default: splice_outcomes=False
    classes = {type(e) for e in effects}
    assert SpliceOutcomeSet not in classes


# --------------------------------------------------------------------
# Opt-in wraps splice effects
# --------------------------------------------------------------------


def test_opt_in_wraps_splice_donor():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    effects = variant.effects(splice_outcomes=True)
    cftr_effect = next(
        e for e in effects
        if getattr(e, "transcript", None) is not None
        and e.transcript.id == CFTR_TRANSCRIPT_ID
    )
    assert isinstance(cftr_effect, SpliceOutcomeSet)
    assert cftr_effect.disrupted_signal_class is SpliceDonor


def test_opt_in_wraps_splice_acceptor():
    # CFTR exon 4 acceptor -1 with canonical ref G.
    variant = Variant("7", 117530898, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effect = variant.effects(splice_outcomes=True)
    target = next(e for e in effect if e.transcript is transcript)
    assert isinstance(target, SpliceOutcomeSet)
    assert target.disrupted_signal_class is SpliceAcceptor


def test_opt_in_wraps_exonic_splice_site():
    # CFTR exon 4 ends with AAG. G->T at -1 disrupts the MAG signal.
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert isinstance(target, SpliceOutcomeSet)
    assert target.disrupted_signal_class is ExonicSpliceSite


def test_opt_in_wraps_intronic_splice_site():
    # CFTR exon 4 +1 with NON-canonical ref A is downgraded to
    # IntronicSpliceSite (post-2.0.0 sequence-aware classification).
    variant = Variant("7", 117531115, "A", "G", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert isinstance(target, SpliceOutcomeSet)
    assert target.disrupted_signal_class is IntronicSpliceSite


def test_opt_in_passes_through_non_splice_effects():
    # Pure substitution that doesn't touch any splice signal: should
    # not be wrapped.
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", ensembl_grch38)
    effects = variant.effects(splice_outcomes=True)
    # At least one effect should be a Substitution, not wrapped.
    classes = [type(e).__name__ for e in effects]
    assert "Substitution" in classes


# --------------------------------------------------------------------
# Plausibility ordering and candidate composition
# --------------------------------------------------------------------


def test_splice_donor_candidate_set_has_expected_outcomes():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effect = variant.effects(splice_outcomes=True)
    target = next(e for e in effect if e.transcript is transcript)
    outcomes = {c.outcome for c in target.candidates}
    assert outcomes == {
        SpliceOutcome.EXON_SKIPPING,
        SpliceOutcome.INTRON_RETENTION,
        SpliceOutcome.CRYPTIC_DONOR,
        SpliceOutcome.NORMAL_SPLICING,
    }


def test_splice_acceptor_candidate_set_uses_cryptic_acceptor():
    variant = Variant("7", 117530898, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    outcomes = {c.outcome for c in target.candidates}
    # SpliceAcceptor disruption uses CRYPTIC_ACCEPTOR not CRYPTIC_DONOR.
    assert SpliceOutcome.CRYPTIC_ACCEPTOR in outcomes
    assert SpliceOutcome.CRYPTIC_DONOR not in outcomes


def test_candidates_sorted_by_plausibility_descending():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    plaus = [c.plausibility for c in target.candidates]
    assert plaus == sorted(plaus, reverse=True), (
        "Candidates should be ordered most-plausible-first, got %r"
        % plaus)


def test_most_likely_for_splice_donor_is_exon_skipping():
    # Per the plausibility table, EXON_SKIPPING dominates SpliceDonor.
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert target.most_likely.outcome is SpliceOutcome.EXON_SKIPPING


def test_most_likely_for_exonic_splice_site_is_normal_splicing():
    # ExonicSpliceSite gets NORMAL_SPLICING as the most-likely
    # outcome (the disruption is on the exon side and often tolerated).
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert target.most_likely.outcome is SpliceOutcome.NORMAL_SPLICING


def test_normal_splicing_carries_underlying_coding_effect():
    # ExonicSpliceSite has an alternate_effect (the coding change if
    # splicing proceeds). NORMAL_SPLICING candidate exposes it.
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    normal = next(
        c for c in target.candidates
        if c.outcome is SpliceOutcome.NORMAL_SPLICING
    )
    assert normal.coding_effect is not None
    assert "p." in normal.coding_effect.short_description


# --------------------------------------------------------------------
# Per-outcome detail
# --------------------------------------------------------------------


def test_intron_retention_candidate_predicts_premature_stop():
    # Intron retention typically produces a PrematureStop. Stub
    # without exact protein since we don't have intronic genomic seq.
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    intron = next(
        c for c in target.candidates
        if c.outcome is SpliceOutcome.INTRON_RETENTION
    )
    assert intron.predicted_class_name == "PrematureStop"
    assert intron.coding_effect is None


def test_cryptic_donor_candidate_is_a_stub():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    cryptic = next(
        c for c in target.candidates
        if c.outcome is SpliceOutcome.CRYPTIC_DONOR
    )
    assert cryptic.coding_effect is None
    assert "cryptic" in cryptic.description.lower()


def test_exon_skipping_for_in_frame_exon_emits_deletion():
    # CFTR exon 4 is 216 nucleotides = 72 codons (216 % 3 == 0), so
    # skipping it is in-frame. The candidate should report Deletion
    # of the exon's amino acids.
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    skip = next(
        c for c in target.candidates
        if c.outcome is SpliceOutcome.EXON_SKIPPING
    )
    # Either a Deletion was constructed, or the candidate falls back
    # to None with predicted_class_name still set. Both are valid.
    if skip.coding_effect is not None:
        assert skip.predicted_class_name == "Deletion"
        assert skip.coding_effect.aa_ref  # non-empty AA range
    else:
        assert skip.predicted_class_name in ("Deletion", "FrameShift", "ExonLoss")


# --------------------------------------------------------------------
# EffectCollection integration
# --------------------------------------------------------------------


def test_collection_iteration_after_wrapping():
    # The wrapped collection should still be iterable, indexable, and
    # produce SpliceOutcomeSet objects in place of the original splice
    # effects.
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    effects = variant.effects(splice_outcomes=True)
    items = list(effects)
    assert len(items) > 0
    splice_set_count = sum(1 for e in items if isinstance(e, SpliceOutcomeSet))
    assert splice_set_count >= 1


def test_short_description_uses_most_likely():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    desc = target.short_description
    assert desc.startswith("splice-set:")
    assert target.most_likely.outcome.value in desc


# --------------------------------------------------------------------
# Reverse-strand
# --------------------------------------------------------------------


def test_opt_in_works_on_reverse_strand_donor():
    # BRCA1 exon 12 reverse-strand donor at 43082403 with canonical ref C.
    variant = Variant("17", 43082403, "C", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(BRCA1_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert isinstance(target, SpliceOutcomeSet)
    assert target.disrupted_signal_class is SpliceDonor


# --------------------------------------------------------------------
# Direct enumerate_splice_outcomes tests
# --------------------------------------------------------------------


def test_enumerate_passes_through_non_splice():
    # Non-splice effect should pass through unchanged.
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(BRCA1_TRANSCRIPT_ID)
    sub_effect = variant.effect_on_transcript(transcript)
    assert type(sub_effect).__name__ == "Substitution"
    wrapped = enumerate_splice_outcomes(sub_effect)
    assert wrapped is sub_effect


# --------------------------------------------------------------------
# SpliceCandidate dataclass ergonomics
# --------------------------------------------------------------------


def test_splice_candidate_is_frozen():
    c = SpliceCandidate(
        outcome=SpliceOutcome.EXON_SKIPPING,
        plausibility=0.5,
        description="test",
    )
    try:
        c.plausibility = 0.9  # type: ignore
    except Exception:
        pass
    else:
        raise AssertionError("SpliceCandidate should be frozen")


def test_splice_candidate_equality():
    a = SpliceCandidate(
        outcome=SpliceOutcome.EXON_SKIPPING,
        plausibility=0.5,
        description="d",
    )
    b = SpliceCandidate(
        outcome=SpliceOutcome.EXON_SKIPPING,
        plausibility=0.5,
        description="d",
    )
    assert a == b


def test_package_level_exports():
    assert varcode.SpliceCandidate is SpliceCandidate
    assert varcode.SpliceOutcome is SpliceOutcome
    assert varcode.SpliceOutcomeSet is SpliceOutcomeSet


# --------------------------------------------------------------------
# MultiOutcomeEffect protocol (see #299 for the planned generalization).
# --------------------------------------------------------------------


def test_splice_outcome_set_is_a_multi_outcome_effect():
    # Downstream consumers filter multi-outcome results with
    # isinstance(e, MultiOutcomeEffect) so future wrappers (RNA
    # evidence #259, germline-aware #268, etc.) don't force churn.
    from varcode import MultiOutcomeEffect
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert isinstance(target, MultiOutcomeEffect)
    # Protocol surface: candidates, most_likely, priority_class.
    assert hasattr(target, "candidates") and len(target.candidates) > 0
    assert target.most_likely is target.candidates[0]
    assert target.priority_class is target.disrupted_signal_class


def test_multi_outcome_effect_exported_at_package_level():
    from varcode import MultiOutcomeEffect
    from varcode import effects
    assert MultiOutcomeEffect is effects.MultiOutcomeEffect
    # Confirm SpliceOutcomeSet is a subclass, not a duck.
    assert issubclass(SpliceOutcomeSet, MultiOutcomeEffect)


# --------------------------------------------------------------------
# Unified Outcome.effect contract (#339).
#
# After #339, ``outcome.effect`` is always a MutationEffect — never a
# SpliceCandidate — and consumers can read
# ``outcome.effect.mutant_protein_sequence`` and
# ``outcome.effect.short_description`` uniformly across SV, splice,
# and point-variant outcomes. The splice-outcome enum is carried as
# ``evidence["splice_outcome"]``.
# --------------------------------------------------------------------


def test_outcomes_effect_is_always_a_mutation_effect():
    """Every outcome's .effect is a :class:`MutationEffect`, including
    the intron-retention and cryptic-splice placeholders (which used
    to be represented by ``coding_effect=None`` in the candidate
    shape). Consumers can call ``effect.short_description`` and
    ``effect.mutant_protein_sequence`` without None-checking."""
    from varcode.effects.effect_classes import MutationEffect
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    for outcome in target.outcomes:
        assert isinstance(outcome.effect, MutationEffect), (
            "Outcome.effect must be MutationEffect, got %s"
            % type(outcome.effect).__name__)
        # short_description is always present and a str.
        assert isinstance(outcome.effect.short_description, str)


def test_outcomes_carry_splice_outcome_in_evidence():
    """Evidence dict carries the :class:`SpliceOutcome` enum so splice-
    aware consumers can dispatch without re-reading
    :attr:`candidates`."""
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    tags = {o.evidence["splice_outcome"] for o in target.outcomes}
    assert SpliceOutcome.EXON_SKIPPING in tags
    assert SpliceOutcome.INTRON_RETENTION in tags


def test_outcomes_carry_plausibility_and_description():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    for outcome in target.outcomes:
        # plausibility moved from SpliceCandidate to Outcome.probability.
        assert outcome.probability is not None
        assert 0.0 <= outcome.probability <= 1.0
        # Description carries the human sentence.
        assert outcome.description is None or isinstance(
            outcome.description, str)


def test_intron_retention_outcome_effect_is_placeholder_class():
    """INTRON_RETENTION outcomes used to hold ``coding_effect=None``;
    after #339 they carry a :class:`PredictedIntronRetention` so
    ``outcome.effect`` satisfies the MutationEffect contract even
    though no protein is computable."""
    from varcode.effects.effect_classes import PredictedIntronRetention
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    ir = next(
        o for o in target.outcomes
        if o.evidence["splice_outcome"] is SpliceOutcome.INTRON_RETENTION)
    assert isinstance(ir.effect, PredictedIntronRetention)
    assert ir.effect.mutant_protein_sequence is None
    assert "intron-retention" in ir.effect.short_description


def test_cryptic_donor_outcome_effect_is_placeholder_class():
    from varcode.effects.effect_classes import PredictedCrypticSpliceSite
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    cryptic = next(
        (o for o in target.outcomes
         if o.evidence["splice_outcome"] is SpliceOutcome.CRYPTIC_DONOR),
        None)
    if cryptic is None:
        pytest.skip("Plausibility table does not include CRYPTIC_DONOR for "
                    "this splice-class combination.")
    assert isinstance(cryptic.effect, PredictedCrypticSpliceSite)
    assert cryptic.effect.direction == "donor"
    assert cryptic.effect.mutant_protein_sequence is None


# --------------------------------------------------------------------
# Boundary-codon reconstruction when in-frame exon skip starts
# mid-codon (see #298). CFTR exon 6 (cDNA 875-1000, length 126) is
# in-frame but starts at codon phase 2, so the preceding exon
# contributes 2 bases to the boundary codon. Skipping it reshapes
# that codon from two flanking-exon bases.
# --------------------------------------------------------------------


def test_mid_codon_in_frame_exon_skip_reshapes_boundary():
    # Donor +1 disruption at CFTR exon 6's 3' end.
    variant = Variant("7", 117536674, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    if not isinstance(bare, SpliceDonor):
        pytest.skip(
            "Canonical donor G not present at 117536674; "
            "classifier emitted %s." % type(bare).__name__)
    effects = variant.effects(splice_outcomes=True)
    splice_set = next(e for e in effects if e.transcript is transcript)
    skip_candidate = next(
        c for c in splice_set.candidates
        if c.outcome is SpliceOutcome.EXON_SKIPPING
    )
    # Pre-#298, this returned a plain Deletion with an aa_ref that
    # didn't account for the reshaped boundary codon. Post-#298, the
    # effect is a Deletion / Substitution / ComplexSubstitution
    # depending on what the reshaped codon translates to.
    assert skip_candidate.coding_effect is not None, (
        "Expected a concrete coding_effect for mid-codon in-frame skip, "
        "got stub with predicted_class_name=%s"
        % skip_candidate.predicted_class_name)
    effect_class = type(skip_candidate.coding_effect).__name__
    assert effect_class in ("Deletion", "Substitution", "ComplexSubstitution"), (
        "Expected Deletion/Substitution/ComplexSubstitution, got %s"
        % effect_class)
    # Whichever class it is, the predicted_class_name label matches.
    assert skip_candidate.predicted_class_name == effect_class


def test_codon_aligned_exon_skip_still_pure_deletion():
    # CFTR exon 3 (cDNA 405-620, length 216) starts AND ends at codon
    # boundaries (phase 0, in-frame). No boundary reshaping needed;
    # result should stay a pure Deletion.
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    splice_set = next(e for e in effects if e.transcript is transcript)
    skip_candidate = next(
        c for c in splice_set.candidates
        if c.outcome is SpliceOutcome.EXON_SKIPPING
    )
    # CFTR exon 3 is codon-aligned, so we expect a clean Deletion.
    assert type(skip_candidate.coding_effect).__name__ == "Deletion"


# --------------------------------------------------------------------
# Serialization round-trip (see #295)
# --------------------------------------------------------------------


def test_splice_candidate_round_trip_with_coding_effect():
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    normal = next(
        c for c in target.candidates
        if c.outcome is SpliceOutcome.NORMAL_SPLICING
    )
    d = normal.to_dict()
    rt = SpliceCandidate.from_dict(d)
    assert rt.outcome is SpliceOutcome.NORMAL_SPLICING
    assert rt.plausibility == normal.plausibility
    assert rt.description == normal.description
    assert type(rt.coding_effect).__name__ == type(normal.coding_effect).__name__
    assert rt.coding_effect.aa_ref == normal.coding_effect.aa_ref


def test_splice_candidate_round_trip_stub():
    # Intron retention candidate has coding_effect=None but carries a
    # predicted_class_name.
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    intron = next(
        c for c in target.candidates
        if c.outcome is SpliceOutcome.INTRON_RETENTION
    )
    d = intron.to_dict()
    rt = SpliceCandidate.from_dict(d)
    assert rt.outcome is SpliceOutcome.INTRON_RETENTION
    assert rt.coding_effect is None
    assert rt.predicted_class_name == "PrematureStop"


def test_splice_outcome_set_json_round_trip():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)

    payload = target.to_json()
    restored = SpliceOutcomeSet.from_json(payload)

    assert type(restored) is SpliceOutcomeSet
    assert len(restored.candidates) == len(target.candidates)
    assert restored.disrupted_signal_class is target.disrupted_signal_class
    assert restored.most_likely.outcome is target.most_likely.outcome
    # Pairwise candidate comparison — tuple equality via __eq__ on the
    # frozen dataclass + the underlying effect's structural equality.
    for rt_c, og_c in zip(restored.candidates, target.candidates):
        assert rt_c.outcome is og_c.outcome
        assert rt_c.plausibility == og_c.plausibility
        assert rt_c.predicted_class_name == og_c.predicted_class_name
        assert (rt_c.coding_effect is None) == (og_c.coding_effect is None)
        if rt_c.coding_effect is not None:
            assert type(rt_c.coding_effect) is type(og_c.coding_effect)


def test_splice_outcome_set_dict_round_trip_is_idempotent():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)

    d = target.to_dict()
    restored = SpliceOutcomeSet.from_dict(d)
    # Serializable's to_dict is stable across round-trip — the restored
    # object's to_dict should match the original's.
    assert restored.to_dict() == d


def test_splice_outcome_set_round_trip_preserves_priority_class():
    # priority_class is a property derived from disrupted_signal_class,
    # so round-trip correctness depends on the class being rehydrated
    # to a class object (not a bare string).
    from varcode.effects import effect_priority
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    restored = SpliceOutcomeSet.from_json(target.to_json())
    assert effect_priority(restored) == effect_priority(target)


def test_non_splice_effects_are_not_multi_outcome():
    # Guard against future class-hierarchy rearrangements that might
    # accidentally mark deterministic effects as multi-outcome.
    from varcode import MultiOutcomeEffect
    from varcode.effects import Substitution, Silent, Intronic, MutationEffect
    for cls in (Substitution, Silent, Intronic, MutationEffect):
        assert not issubclass(cls, MultiOutcomeEffect), (
            "%s should not be a MultiOutcomeEffect" % cls.__name__)


# --------------------------------------------------------------------
# Priority integration: SpliceOutcomeSet sorts as if it were the
# disrupted-signal class (review feedback on PR #292).
# --------------------------------------------------------------------


def test_splice_outcome_set_sorts_as_disrupted_signal_class():
    # When wrapped, a SpliceDonor-backed SpliceOutcomeSet should have
    # the same priority as a bare SpliceDonor — higher than Intronic,
    # lower than Substitution. If the priority delegation is broken,
    # SpliceOutcomeSet gets priority -1 and sorts to the bottom.
    from varcode.effects import effect_priority

    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare_effect = variant.effect_on_transcript(transcript)
    assert isinstance(bare_effect, SpliceDonor)
    bare_priority = effect_priority(bare_effect)

    wrapped_effects = variant.effects(splice_outcomes=True)
    wrapped = next(e for e in wrapped_effects if e.transcript is transcript)
    assert isinstance(wrapped, SpliceOutcomeSet)
    wrapped_priority = effect_priority(wrapped)

    assert wrapped_priority == bare_priority, (
        "SpliceOutcomeSet priority (%d) must match the disrupted-"
        "signal class priority (%d); otherwise sorting and "
        "top_priority_effect() behave wrongly." % (
            wrapped_priority, bare_priority))


def test_splice_outcome_set_top_priority_works():
    # top_priority_effect on a collection containing SpliceOutcomeSet
    # should not pick a lower-priority non-splice effect.
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    effects = variant.effects(splice_outcomes=True)
    top = effects.top_priority_effect()
    # The wrapped SpliceDonor (or one of the splice-set variants) is
    # higher priority than Intronic/NoncodingTranscript from other
    # overlapping transcripts, so the top should be a splice-related
    # effect.
    top_class_name = type(top).__name__
    assert top_class_name in ("SpliceOutcomeSet", "SpliceDonor"), (
        "Expected a splice-related effect at top priority, got %s"
        % top_class_name)


# --------------------------------------------------------------------
# Acceptor-side IntronicSpliceSite emits CRYPTIC_ACCEPTOR, not DONOR.
# --------------------------------------------------------------------


def test_acceptor_side_intronic_splice_site_uses_cryptic_acceptor():
    # CFTR exon 4 acceptor -3 (3bp before exon.start). A variant here
    # with NON-canonical ref (not A, the canonical MAG component) is
    # classified as IntronicSpliceSite. The splice set should include
    # CRYPTIC_ACCEPTOR (the relevant cryptic direction for the
    # acceptor side), not CRYPTIC_DONOR.
    from varcode.effects import IntronicSpliceSite
    # chr7:117530896 is -3 before CFTR exon 4 (forward strand).
    # Use a non-canonical ref for the -3 position so it's
    # IntronicSpliceSite (not SpliceAcceptor which covers -1/-2).
    # At distance -3, the position isn't required to be canonical
    # anyway — the classifier emits IntronicSpliceSite for this window.
    variant = Variant("7", 117530896, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    assert isinstance(bare, IntronicSpliceSite) and \
        not isinstance(bare, (SpliceDonor, SpliceAcceptor))
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    outcomes = {c.outcome for c in target.candidates}
    assert SpliceOutcome.CRYPTIC_ACCEPTOR in outcomes, \
        "Acceptor-side IntronicSpliceSite should use CRYPTIC_ACCEPTOR"
    assert SpliceOutcome.CRYPTIC_DONOR not in outcomes, \
        "Acceptor-side IntronicSpliceSite should not use CRYPTIC_DONOR"


def test_donor_side_intronic_splice_site_uses_cryptic_donor():
    # Mirror test for donor-side IntronicSpliceSite at +3 after CFTR
    # exon 4 end (117531117 = exon.end + 3).
    from varcode.effects import IntronicSpliceSite
    variant = Variant("7", 117531117, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    assert isinstance(bare, IntronicSpliceSite) and \
        not isinstance(bare, (SpliceDonor, SpliceAcceptor))
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    outcomes = {c.outcome for c in target.candidates}
    assert SpliceOutcome.CRYPTIC_DONOR in outcomes
    assert SpliceOutcome.CRYPTIC_ACCEPTOR not in outcomes


# --------------------------------------------------------------------
# Multi-protein surface: candidate_proteins and mutant_protein_sequences
# --------------------------------------------------------------------


def test_candidate_proteins_maps_each_outcome_to_a_protein():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    proteins = target.candidate_proteins
    # Every candidate outcome appears as a key.
    outcomes = {c.outcome for c in target.candidates}
    assert set(proteins.keys()) == outcomes
    # EXON_SKIPPING for an in-frame exon should have a non-empty
    # protein (reference minus the skipped AAs). CFTR exon 4 is 216
    # nucleotides = 72 codons = in-frame.
    assert proteins[SpliceOutcome.EXON_SKIPPING], \
        "Expected a concrete mutant protein for in-frame exon skipping"
    # INTRON_RETENTION and CRYPTIC are stubs → empty string for now.
    assert proteins[SpliceOutcome.INTRON_RETENTION] == ""
    assert proteins[SpliceOutcome.CRYPTIC_DONOR] == ""


def test_mutant_protein_sequences_collects_distinct_proteins():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    proteins = target.mutant_protein_sequences
    assert isinstance(proteins, set)
    assert len(proteins) >= 1
    # Reference protein should be in there or a proper subset of
    # reference (exon-skipped version is shorter).
    ref = str(transcript.protein_sequence)
    # The in-frame exon skip removes exon 4 AAs; resulting protein
    # should be shorter than reference.
    shortest = min(proteins, key=len)
    assert len(shortest) < len(ref)


# --------------------------------------------------------------------
# Out-of-frame exon skip now produces a real mutant protein
# --------------------------------------------------------------------


def test_out_of_frame_exon_skip_produces_mutant_protein():
    # CFTR exon 5 is 90 nucleotides = 30 codons, BUT exon 5 is not
    # out of frame — need a different exon. Use a variant known to
    # target an out-of-frame exon. We'll discover one empirically
    # by finding an exon whose length is not divisible by 3.
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    target_exon = None
    for exon in transcript.exons[2:]:
        length = exon.end - exon.start + 1
        if length % 3 != 0:
            target_exon = exon
            break
    if target_exon is None:
        pytest.skip("No out-of-frame exon found in CFTR beyond exon 2")

    # Construct a donor-side disrupting variant at this exon's end
    # (+1 position after the exon in + strand coords).
    donor_plus_1 = target_exon.end + 1
    # Use a canonical-ref SNV to ensure SpliceDonor classification.
    variant = Variant("7", donor_plus_1, "G", "A", ensembl_grch38)
    bare = variant.effect_on_transcript(transcript)
    if not isinstance(bare, SpliceDonor):
        pytest.skip(
            "Canonical donor G not present at %d; classifier emitted %s "
            "rather than SpliceDonor." % (donor_plus_1, type(bare).__name__))
    effects = variant.effects(splice_outcomes=True)
    splice_set = next(e for e in effects if e.transcript is transcript)
    skip_candidate = next(
        c for c in splice_set.candidates
        if c.outcome is SpliceOutcome.EXON_SKIPPING
    )
    # Out-of-frame skip should now carry a mutant protein.
    assert skip_candidate.coding_effect is not None, (
        "Out-of-frame exon skip should produce a concrete mutant "
        "protein, not a stub")
    protein = skip_candidate.coding_effect.mutant_protein_sequence
    assert isinstance(protein, str)
    assert len(protein) > 0
    # The frameshifted protein should differ from the reference
    # after the skip point.
    assert protein != str(transcript.protein_sequence)


# --------------------------------------------------------------------
# has_protein property on SpliceCandidate
# --------------------------------------------------------------------


def test_has_protein_is_true_for_candidates_with_coding_effect():
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    # NORMAL_SPLICING has a Substitution coding_effect with a protein.
    normal = next(
        c for c in target.candidates
        if c.outcome is SpliceOutcome.NORMAL_SPLICING
    )
    assert normal.has_protein is True


def test_has_protein_is_false_for_stub_candidates():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    intron = next(
        c for c in target.candidates
        if c.outcome is SpliceOutcome.INTRON_RETENTION
    )
    assert intron.has_protein is False


# ------------------------------------------------------------------
# SpliceCandidate.mutant_transcript exposure (#305).
# ------------------------------------------------------------------


def test_exon_skipping_candidate_exposes_mutant_transcript():
    """EXON_SKIPPING candidates now carry the :class:`MutantTranscript`
    they were built from. Consumers get the full cDNA + protein
    without re-deriving from coding_effect fields (#305)."""
    from varcode import MutantTranscript
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    exon_skip = next(
        c for c in target.candidates
        if c.outcome is SpliceOutcome.EXON_SKIPPING
    )
    assert isinstance(exon_skip.mutant_transcript, MutantTranscript)
    # The MutantTranscript has a materialized cDNA and the expected
    # edit shape (one deletion of the skipped exon).
    assert exon_skip.mutant_transcript.cdna_sequence is not None
    assert len(exon_skip.mutant_transcript.edits) == 1
    edit = exon_skip.mutant_transcript.edits[0]
    assert edit.alt_bases == ""
    # The skipped-exon cDNA delta matches the exon length on the
    # transcript.
    assert edit.cdna_end > edit.cdna_start
    assert exon_skip.mutant_transcript.annotator_name == "splice_outcomes"


def test_stub_candidates_have_no_mutant_transcript():
    """INTRON_RETENTION / CRYPTIC_* have ``mutant_transcript=None``
    when no ``genomic_sequence`` provider is passed. With a provider,
    INTRON_RETENTION populates the transcript (see the #296 tests
    below). CRYPTIC_* remain None until a future pass."""
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    intron = next(
        c for c in target.candidates
        if c.outcome is SpliceOutcome.INTRON_RETENTION
    )
    assert intron.mutant_transcript is None


# ------------------------------------------------------------------
# Intron-retention with a genomic_sequence provider (#296).
# ------------------------------------------------------------------


def _synthetic_intron_provider(stop_offset=30, pad_with="A"):
    """Return a callable that acts as a ``genomic_sequence`` provider.

    The returned sequence starts with ``pad_with * stop_offset`` then
    ``TAA`` (a stop codon) then further ``pad_with`` padding — so any
    ORF translating through the intron hits a premature stop at
    roughly ``stop_offset`` bases in.
    """
    def provider(contig, start, end):
        length = end - start + 1
        core = (pad_with * stop_offset) + "TAA" + (pad_with * length)
        return core[:length]
    return provider


def test_intron_retention_with_provider_populates_mutant_transcript():
    """When a ``genomic_sequence`` provider is passed to
    :func:`enumerate_splice_outcomes`, the INTRON_RETENTION candidate
    stops being a stub — it carries a full :class:`MutantTranscript`
    with the intron inserted and the protein truncated at the first
    in-intron stop codon (#296)."""
    from varcode import MutantTranscript
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    assert isinstance(bare, SpliceDonor)
    splice_set = enumerate_splice_outcomes(
        bare, genomic_sequence=_synthetic_intron_provider())
    intron = next(
        c for c in splice_set.candidates
        if c.outcome is SpliceOutcome.INTRON_RETENTION
    )
    assert isinstance(intron.mutant_transcript, MutantTranscript)
    mt = intron.mutant_transcript
    # cDNA is strictly longer than the reference (intron inserted).
    assert len(mt.cdna_sequence) > len(str(transcript.sequence))
    # Protein is truncated — shorter than the reference protein.
    assert mt.mutant_protein_sequence is not None
    assert len(mt.mutant_protein_sequence) < len(
        str(transcript.protein_sequence))
    # coding_effect is populated (was None in the stub case).
    assert intron.coding_effect is not None


def test_intron_retention_without_provider_still_stub():
    """Back-compat: no ``genomic_sequence`` kwarg → stub behavior
    preserved."""
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    splice_set = enumerate_splice_outcomes(bare)  # no provider
    intron = next(
        c for c in splice_set.candidates
        if c.outcome is SpliceOutcome.INTRON_RETENTION
    )
    assert intron.mutant_transcript is None
    assert intron.coding_effect is None
    assert intron.predicted_class_name == "PrematureStop"


def _synthetic_cryptic_donor_provider(motif_offset=65, motif="CAGGTAAGT"):
    """Provider that embeds a full-consensus donor motif at the given
    offset from scan-window start. Used to verify the cryptic donor
    scan picks up the implanted site."""
    def provider(contig, start, end):
        length = end - start + 1
        seq = bytearray(b"A" * length)
        for i, c in enumerate(motif):
            if motif_offset + i < length:
                seq[motif_offset + i] = ord(c)
        return seq.decode()
    return provider


def test_cryptic_donor_with_provider_picks_up_implanted_motif():
    """With a genomic_sequence provider returning a strong donor
    motif a fixed distance from the canonical boundary, the
    CRYPTIC_DONOR candidate finds it, computes a MutantTranscript
    with the new exon boundary, and classifies the resulting coding
    effect (#296)."""
    from varcode import MutantTranscript
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    # Place the motif 15 bp into the intron of the scan window
    # (canonical sits at offset 50; 65 = canonical + 15).
    splice_set = enumerate_splice_outcomes(
        bare, genomic_sequence=_synthetic_cryptic_donor_provider(
            motif_offset=65))
    cryptic = next(
        c for c in splice_set.candidates
        if c.outcome is SpliceOutcome.CRYPTIC_DONOR
    )
    assert isinstance(cryptic.mutant_transcript, MutantTranscript)
    mt = cryptic.mutant_transcript
    # cDNA is longer (exon extended into intron).
    delta = len(mt.cdna_sequence) - len(str(transcript.sequence))
    assert delta > 0, "Cryptic donor should extend the exon"
    assert cryptic.coding_effect is not None
    # Description names the motif score and genomic position.
    assert "motif score" in cryptic.description
    assert "Cryptic donor" in cryptic.description


def test_cryptic_acceptor_with_provider_picks_up_implanted_motif():
    """Same as the donor test but for the acceptor side. Disrupted
    canonical acceptor + embedded cryptic acceptor motif → the
    CRYPTIC_ACCEPTOR candidate materializes."""
    # CFTR exon 4 acceptor -1: variant at 117530898.
    variant = Variant("7", 117530898, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    if not isinstance(bare, SpliceAcceptor):
        pytest.skip("Variant didn't classify as SpliceAcceptor")

    # Acceptor consensus is YAG|G (4 bp). Put "CAGG" (full-match) at
    # offset 40 — 10 bp before the canonical boundary at 50.
    def provider(contig, start, end):
        length = end - start + 1
        seq = bytearray(b"A" * length)
        motif = b"CAGG"
        pos = 40
        for i, c in enumerate(motif):
            if pos + i < length:
                seq[pos + i] = c
        return seq.decode()

    splice_set = enumerate_splice_outcomes(bare, genomic_sequence=provider)
    cryptic = next(
        c for c in splice_set.candidates
        if c.outcome is SpliceOutcome.CRYPTIC_ACCEPTOR
    )
    assert cryptic.mutant_transcript is not None
    assert "Cryptic acceptor" in cryptic.description


def test_cryptic_without_provider_still_stub():
    """Back-compat: without a provider, cryptic candidates stay stubs."""
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    splice_set = enumerate_splice_outcomes(bare)
    cryptic = next(
        c for c in splice_set.candidates
        if c.outcome is SpliceOutcome.CRYPTIC_DONOR
    )
    assert cryptic.mutant_transcript is None
    assert cryptic.coding_effect is None


def test_cryptic_with_no_motif_in_range_is_stub():
    """If the scan window contains no above-threshold motif, the
    cryptic candidate falls back to the stub (nothing to materialize)."""
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)

    # Provider returns nothing resembling a donor/acceptor motif.
    def provider(contig, start, end):
        return "N" * (end - start + 1)

    splice_set = enumerate_splice_outcomes(bare, genomic_sequence=provider)
    cryptic = next(
        c for c in splice_set.candidates
        if c.outcome is SpliceOutcome.CRYPTIC_DONOR
    )
    assert cryptic.mutant_transcript is None
    assert cryptic.predicted_class_name == "ComplexSubstitution"


def test_intron_retention_provider_can_raise_without_crashing():
    """If the provider raises (unknown contig, missing FASTA, etc.),
    the builder falls back to the stub candidate rather than
    propagating the error."""
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)

    def failing_provider(contig, start, end):
        raise IOError("FASTA not available")
    splice_set = enumerate_splice_outcomes(
        bare, genomic_sequence=failing_provider)
    intron = next(
        c for c in splice_set.candidates
        if c.outcome is SpliceOutcome.INTRON_RETENTION
    )
    # Falls back to stub.
    assert intron.mutant_transcript is None
    assert intron.coding_effect is None
