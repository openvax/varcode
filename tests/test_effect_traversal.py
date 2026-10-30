"""Public traversal preserves alternative products and candidate provenance."""

import pytest
from pyensembl import cached_release

from varcode import (
    EffectCandidate, EffectCollection, MutationEffect, Variant,
    candidate_children, effect_priority, iter_effect_paths, iter_effects, to_proteome,
)
from varcode.effect_hypotheses import merge_classified_outcomes
from varcode.effects import Deletion, Intronic, MultiOutcomeEffect, PrematureStop, Unresolved
from varcode.effects.sequence_change import modification_status
from varcode.orf import annotate_orf_variants

from . import test_proteome as fixtures
from .test_proteome import outcome, comparisons

effect = fixtures.effect
transcript = fixtures.transcript
reference_orf = fixtures.reference_orf


def test_default_intronic_result_retains_protein_changing_alternative():
    genome = cached_release(81)
    tx = genome.transcript_by_id("ENST00000003084")
    result = Variant("7", 117531117, "A", "C", genome).effect_on_transcript(tx)
    assert isinstance(result, Intronic) and not isinstance(result, MultiOutcomeEffect)
    effects = list(iter_effects(result))
    assert effects[0] is result
    assert len(effects) == len({id(e) for e in effects})
    assert result.modifies_protein_sequence is False
    deletions = [e for e in effects if isinstance(e, Deletion)]
    assert deletions and any(e.modifies_protein_sequence is True for e in deletions)
    assert any(e.baseline_transcript.cdna_sequence == tx.sequence for e in deletions)
    paths = list(iter_effect_paths(result))
    assert paths[0] == (result, ())
    assert any(e is result and path for e, path in paths)  # selected self-candidate
    for _, path in paths[1:]:
        for candidate in path:
            assert candidate.outcomes
            assert all(o.baseline.cdna_sequence == tx.sequence for o in candidate.outcomes)


def test_paths_keep_distinct_evidence_and_cycles_terminate(effect):
    root = MutationEffect(effect.variant)
    branch = MutationEffect(effect.variant)
    root.candidates = (EffectCandidate(branch, "first", {"read": "A"}),
                       EffectCandidate(branch, "second", {"read": "B"}))
    branch.candidates = (EffectCandidate(effect, "coding"), EffectCandidate(root, "back"))
    assert list(iter_effects(root)) == [root, branch, effect]
    paths = list(iter_effect_paths(root))
    coding = [path for current, path in paths if current is effect]
    assert [[c.source for c in path] for path in coding] == [["first", "coding"], ["second", "coding"]]
    assert [path[0].evidence["read"] for path in coding] == ["A", "B"]
    assert len(paths) == 7  # root, then branch/coding/cycle-close for each route
    assert [path[-1].source for current, path in paths if current is root and path] == ["back", "back"]


def test_equal_distinct_effects_are_not_merged(effect):
    other = type(effect).from_dict(effect.to_dict())
    root = MutationEffect(effect.variant)
    root.candidates = (EffectCandidate(effect), EffectCandidate(other))
    found = list(iter_effects(root))
    assert len(found) == 3 and found[1] is effect and found[2] is other


def test_nested_legacy_phase_and_splice_results():
    from varcode.germline import PhaseHypothesis
    from varcode.effects.effect_classes import PhaseCandidateSet
    genome = cached_release(81)
    tx = genome.transcript_by_id("ENST00000003084")
    somatic = Variant("7", 117531114, "G", "A", genome)
    splice = somatic.effect_on_transcript(tx, annotator="fast")
    root = PhaseCandidateSet(somatic, tx, [splice],
                             [PhaseHypothesis(cis=(), trans=(), phase_state="unknown", haplotype="A")], ())
    assert isinstance(root, MultiOutcomeEffect)
    paths = list(iter_effect_paths(root))
    assert any(len(path) > 1 for _, path in paths)
    assert any(c.source == "varcode_germline" for _, path in paths for c in path)
    assert len(list(iter_effects(root))) > len(candidate_children(root))


def test_shared_candidates_include_external_alternative_in_filter_rank_and_export(effect):
    baseline = Intronic(effect.variant, effect.transcript, effect.transcript.exons[0], 1)
    baseline.realized_candidates = merge_classified_outcomes([outcome(baseline, mutant="MAK")])
    baseline.candidates = baseline.realized_candidates
    stop = PrematureStop(effect.variant, effect.transcript, 1)
    external = EffectCandidate(stop, "external", {"observation": "stop"})
    baseline._extra_candidates = (external,)
    assert candidate_children(baseline)[-1] is external
    assert list(iter_effects(baseline)) == [baseline, stop]
    assert modification_status(baseline, "modifies_protein_sequence") is True
    assert effect_priority(baseline) == effect_priority(stop)
    assert list(EffectCollection([baseline]).drop_silent_and_noncoding()) == [baseline]
    exported = comparisons(to_proteome(baseline))
    assert len(exported) == 2
    assert any(c.to_dict()["payload"]["contexts"][-1].get("source") == "external" for c in exported)


def test_realized_paths_preserve_different_patient_baselines_and_unknowns(effect):
    unresolved = Unresolved(effect.variant, effect.transcript, "missing_sequence")
    root = effect
    outcomes = (outcome(effect, baseline="MAK", mutant="MTK"),
                outcome(effect, baseline="MVK", mutant="MTK", index=1),
                outcome(unresolved, baseline=None, mutant=None, index=2))
    root.realized_candidates = merge_classified_outcomes(outcomes)
    paths = [path for _, path in iter_effect_paths(root) if path]
    assert len(paths) == 3
    assert {c.outcomes[0].baseline.protein_sequence for path in paths for c in path} == {"MAK", "MVK", None}
    assert list(iter_effects(root)) == [root, unresolved]
    assert unresolved.mutant_protein_sequence is None
    assert unresolved.modifies_protein_sequence is None


def test_noncanonical_orf_uses_same_sequence_and_traversal_contract(effect, reference_orf):
    annotation = annotate_orf_variants([effect.variant], reference_orf)
    root = annotation.effect
    assert root.orf is not None and root.orf.orf_id == reference_orf.orf_id
    found = list(iter_effects(root))
    assert found and found[0] is root
    assert root.original_protein_sequence == "MAK"
    assert root.mutant_protein_sequence == "MTK"
    assert root.reference_transcript is effect.transcript
    assert root.modifies_protein_sequence is True


def test_lazy_errors_propagate_instead_of_dropping_candidates(effect):
    class Broken(MutationEffect):
        @property
        def candidates(self):
            raise ValueError("candidate sequence missing")

    for traversal in (iter_effects, iter_effect_paths):
        with pytest.raises(ValueError, match="candidate sequence missing"):
            list(traversal(Broken(effect.variant)))


def test_collection_delegates_retain_subclass_and_metadata(effect):
    from varcode import VariantCollection

    class Effects(EffectCollection):
        pass

    class Variants(VariantCollection):
        pass

    effects = Effects([effect], sources={"test"}, annotator="transcript_model", annotator_version="11")
    variants = Variants([effect.variant], source_to_metadata_dict={"test": {effect.variant: {"GT": "0/1"}}})
    for collection in (effects, variants):
        item, = collection
        assert collection[0] is item and collection[:1] == [item]
        for copied in (collection.filter(lambda _: True), collection.clone_with_new_elements([item]),
                       collection.groupby(lambda _: "key")["key"]):
            assert type(copied) is type(collection)
            assert copied.sources == collection.sources
            if isinstance(copied, EffectCollection):
                assert copied.annotator == "transcript_model" and copied.annotator_version == "11"
            else:
                assert copied.metadata == collection.metadata
