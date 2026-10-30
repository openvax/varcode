"""Public traversal separates comparison identity from evidence paths (#653)."""

import pytest
from pyensembl import cached_release

from varcode import (
    EffectCandidate, EffectCollection, MultiOutcomeEffect, MutationEffect,
    Variant, classification_priority, get_effect_candidates, effect_priority,
    iter_candidate_effects, iter_candidate_paths,
)
from varcode.effect_hypotheses import (
    ClassifiedOutcome, EffectHypothesis, RealizedEffectCandidate,
    RealizedTranscriptProduct,
)
from varcode.effects import Deletion, Intronic
from varcode.effects.sequence_change import modification_status


class OutcomeSet(MultiOutcomeEffect):
    def __init__(self, *candidates):
        super().__init__(None)
        self.candidates = candidates


def realized(effect, baseline, phase):
    outcome = ClassifiedOutcome(
        hypothesis=EffectHypothesis(evidence={"phase": phase}),
        baseline=RealizedTranscriptProduct(None, None, baseline),
        mutant=RealizedTranscriptProduct(None, None, "MALT"),
        effect=effect)
    return RealizedEffectCandidate(effect, (outcome,), (), None, None, 0)


def test_bare_and_missing_effects_do_not_invent_provenance():
    effect = MutationEffect(None)
    assert get_effect_candidates(effect) == ()
    assert list(iter_candidate_paths(effect)) == [(effect, ())]
    assert list(iter_candidate_paths(None)) == [(None, ())]
    assert list(iter_candidate_effects(None)) == [None]


def test_ordinary_selected_comparison_survives_alternative_only_candidates():
    selected, alternate = MutationEffect(None), MutationEffect(None)
    wrapper = EffectCandidate(alternate, source="external")
    selected.candidates = (wrapper,)
    paths = list(iter_candidate_paths(selected))
    assert paths == [(selected, ()), (alternate, (wrapper,))]


def test_paths_retain_shared_effect_distinct_baselines_and_nested_self_wrapper():
    effect = MutationEffect(None)
    cis, trans = realized(effect, "MPATIENT", "cis"), realized(effect, "MREF", "trans")
    effect.realized_candidates = (cis, trans)
    nested = EffectCandidate(effect, source="rna", evidence={"read_count": 12})
    root = OutcomeSet(nested)
    paths = list(iter_candidate_paths(root))
    assert len(paths) == 2
    assert all(terminal is effect for terminal, _ in paths)
    assert paths[0][1][0] is nested
    assert paths[0][1][1] is cis
    assert paths[1][1][1] is trans
    assert paths[0][1][-1].outcomes[0].baseline.protein_sequence == "MPATIENT"
    assert list(iter_candidate_effects(root)) == [effect]
    assert nested.evidence == {"read_count": 12}
    assert cis.probability is None


def test_realized_children_and_extra_rna_candidates_are_not_double_counted():
    effect, other = MutationEffect(None), MutationEffect(None)
    primary = realized(effect, "MREF", "unknown")
    extra = EffectCandidate(other, source="rna")
    effect.realized_candidates = (primary,)
    effect.candidates = (primary, extra)
    effect._extra_candidates = (extra,)
    assert get_effect_candidates(effect) == (primary, extra)
    assert list(iter_candidate_paths(effect)) == [(effect, (primary,)), (other, (extra,))]


def test_empty_sets_and_longer_cycles_remain_unknown_but_self_is_concrete():
    empty, first, second = OutcomeSet(), OutcomeSet(), OutcomeSet()
    first.candidates = (EffectCandidate(second),)
    second.candidates = (EffectCandidate(first),)
    assert list(iter_candidate_paths(empty)) == [(None, ())]
    terminal, path = next(iter_candidate_paths(first))
    assert terminal is None
    assert len(path) == 2
    assert list(iter_candidate_effects(first)) == [None]
    first.candidates += (EffectCandidate(first, source="primary"),)
    assert list(iter_candidate_effects(first)) == [None, first]


def test_shared_leaf_keeps_each_path_but_effect_iteration_deduplicates_identity():
    leaf = MutationEffect(None)
    left, right = EffectCandidate(leaf, source="left"), EffectCandidate(leaf, source="right")
    root = OutcomeSet(left, EffectCandidate(OutcomeSet(right)))
    assert list(iter_candidate_paths(root)) == [
        (leaf, (left,)), (leaf, (root.candidates[1], right))]
    assert list(iter_candidate_effects(root)) == [leaf]


def test_deep_nesting_does_not_depend_on_python_recursion_limit():
    leaf = MutationEffect(None)
    root = leaf
    for _ in range(1100):
        root = OutcomeSet(EffectCandidate(root))
    terminal, path = next(iter_candidate_paths(root))
    assert terminal is leaf
    assert len(path) == 1100


@pytest.mark.parametrize("sequence", [None, "", "MALT"])
def test_traversal_does_not_fill_unavailable_or_empty_products(sequence):
    effect = MutationEffect(None)
    effect.mutant_protein_sequence = sequence
    root = OutcomeSet(EffectCandidate(effect))
    assert next(iter_candidate_effects(root)).mutant_protein_sequence is sequence


def test_real_cftr_preferred_intronic_keeps_coding_alternative_and_local_priority():
    genome = cached_release(81)
    variant = Variant("7", 117531117, "A", "C", genome)
    effect = variant.effect_on_transcript(genome.transcript_by_id("ENST00000003084"))
    assert isinstance(effect, Intronic)
    assert effect.modifies_coding_sequence is False
    paths = list(iter_candidate_paths(effect))
    assert len(paths) == len(effect.candidates)
    assert paths[0][0] is effect
    assert all(len(path) == 1 for _, path in paths)
    assert all(path[0] is candidate for (_, path), candidate in zip(paths, effect.candidates))
    alternate = next(e for e in iter_candidate_effects(effect) if isinstance(e, Deletion))
    assert classification_priority(alternate) > classification_priority(effect)
    assert effect_priority(effect) >= classification_priority(alternate)
    assert modification_status(effect, "modifies_coding_sequence") is True
    assert list(EffectCollection([effect]).drop_silent_and_noncoding()) == [effect]
