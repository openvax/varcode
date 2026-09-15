"""Tests for graph-valid hypotheses and baseline-aware merging."""

from types import SimpleNamespace

import pytest

from varcode.effect_hypotheses import (
    ClassifiedOutcome,
    EffectHypothesis,
    RealizedTranscriptProduct,
    SpliceAxis,
    SpliceOption,
    SpliceSiteKey,
    enumerate_splice_plans,
    merge_classified_outcomes,
    order_realized_candidates,
)


def effect(name="Deletion", priority=5):
    cls = type(name, (), {})
    value = cls()
    value.aa_mutation_start_offset = 3
    value.aa_ref = "D"
    value.aa_alt = ""
    value.priority = priority
    return value


def product(protein, cdna="AAACCC", junctions=()):
    return RealizedTranscriptProduct(
        transcript=None,
        cdna_sequence=cdna,
        protein_sequence=protein,
        junction_signature=junctions)


def test_splice_enumeration_omits_site_removed_by_another_choice():
    exon4 = (4, 1)
    acceptor = SpliceAxis(
        key=SpliceSiteKey("acceptor", 4),
        anchor_runs=(exon4,),
        options=(
            SpliceOption("normal", 0, required_runs=(exon4,)),
            SpliceOption("exon_skip", 1, skipped_runs=(exon4,)),
        ))
    donor = SpliceAxis(
        key=SpliceSiteKey("donor", 4),
        anchor_runs=(exon4,),
        options=(
            SpliceOption("normal", 0, required_runs=(exon4,)),
            SpliceOption("exon_skip", 1, skipped_runs=(exon4,)),
        ))

    plans = enumerate_splice_plans((acceptor, donor), (exon4,))
    mechanisms = [
        tuple(option.mechanism for _, option in plan.choices)
        for plan in plans]

    assert mechanisms == [("normal", "normal"), ("exon_skip",),
                          ("exon_skip",)]
    assert all(len(plan.choices) != 2 or plan.kept_runs for plan in plans)


def test_splice_rule_rank_is_not_reported_as_probability():
    axis = SpliceAxis(
        key=SpliceSiteKey("donor", 4),
        anchor_runs=(4,),
        options=(SpliceOption("exon_skip", 0, skipped_runs=(4,)),))

    (plan,) = enumerate_splice_plans((axis,), (4,))

    assert plan.ordinal_key == (0,)
    assert plan.probability is None


def test_calibrated_splice_probabilities_are_multiplied():
    axes = tuple(
        SpliceAxis(
            key=SpliceSiteKey("donor", exon),
            anchor_runs=(exon,),
            options=(SpliceOption(
                "normal", 0, required_runs=(exon,), probability=p),))
        for exon, p in ((4, 0.8), (5, 0.5)))

    (plan,) = enumerate_splice_plans(axes, (4, 5))

    assert plan.probability == pytest.approx(0.4)


def test_same_mutant_different_patient_baseline_does_not_merge():
    mutant = product("ABCDE")
    outcome_a = ClassifiedOutcome(
        hypothesis=EffectHypothesis(phase=(("g", "trans"),)),
        effect=effect(),
        baseline=product("ABCDEFG"),
        mutant=mutant,
        change_aa_offset=3)
    outcome_b = ClassifiedOutcome(
        hypothesis=EffectHypothesis(phase=(("g", "cis"),)),
        effect=effect(),
        baseline=product("ABCXEFG"),
        mutant=mutant,
        change_aa_offset=3)

    candidates = merge_classified_outcomes((outcome_a, outcome_b))

    assert len(candidates) == 2


def test_identical_comparisons_merge_and_preserve_hypotheses():
    baseline = product("ABCDEFG")
    mutant = product("ABCDE")
    outcomes = tuple(
        ClassifiedOutcome(
            hypothesis=EffectHypothesis(
                phase=(("g", phase),), phase_probability=0.5,
                enumeration_index=index),
            effect=effect(),
            baseline=baseline,
            mutant=mutant,
            change_aa_offset=3)
        for index, phase in enumerate(("trans", "cis")))

    (candidate,) = merge_classified_outcomes(outcomes)

    assert len(candidate.hypotheses) == 2
    assert candidate.phase_mass == 1.0


def test_candidate_probability_combines_calibrated_splice_and_phase():
    baseline = product("ABCDEFG")
    mutant = product("ABCDE")
    plan = SimpleNamespace(ordinal_key=(0,), probability=0.4)
    outcomes = tuple(
        ClassifiedOutcome(
            EffectHypothesis(
                phase=(("g", phase),), phase_probability=0.5,
                splice_plan=plan),
            effect(), baseline, mutant, 3)
        for phase in ("trans", "cis"))

    (candidate,) = merge_classified_outcomes(outcomes)

    assert candidate.probability == pytest.approx(0.4)


def test_junction_signature_is_part_of_product_identity():
    baseline = product("ABCDEFG")
    outcomes = (
        ClassifiedOutcome(
            EffectHypothesis(), effect(), baseline,
            product("ABCDE", junctions=((4, 6),)), 3),
        ClassifiedOutcome(
            EffectHypothesis(), effect(), baseline,
            product("ABCDE", junctions=((4, 7),)), 3),
    )

    assert len(merge_classified_outcomes(outcomes)) == 2


def test_ordering_uses_ordinal_rank_when_scores_are_uncalibrated():
    baseline = product("ABCDEFG")
    mutant = product("ABCDE")
    low_priority = effect("Deletion", priority=5)
    high_priority = effect("FrameShift", priority=10)
    plans = [
        SimpleNamespace(ordinal_key=(0,), probability=None),
        SimpleNamespace(ordinal_key=(1,), probability=None),
    ]
    outcomes = tuple(
        ClassifiedOutcome(
            EffectHypothesis(splice_plan=plan, enumeration_index=index),
            current_effect, baseline, mutant, 3)
        for index, (plan, current_effect) in enumerate(
            zip(plans, (low_priority, high_priority))))
    candidates = merge_classified_outcomes(outcomes)

    ordered = order_realized_candidates(
        candidates, lambda value: value.priority)

    assert ordered[0].effect is low_priority
    assert ordered[0].probability is None
