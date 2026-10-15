"""Effect archives preserve realized biology, uncertainty and shared links."""

import pytest
from pyensembl import cached_release

from varcode import EffectCollection, Variant, predict_transcript_model_effect
from varcode.effect_hypotheses import RealizedEffectCandidate, ClassifiedOutcome
from varcode.effects import Substitution
from varcode.effects.structural_serialization import _EffectReference


def assert_realized_annotation_equal(before, after):
    assert len(after.candidates) == len(before.candidates)
    assert after.candidates == after.realized_candidates
    assert after.candidates[0].effect is after
    assert after.most_likely_candidate is after.candidates[0]
    assert any(after.highest_priority_candidate is c for c in after.candidates)
    assert after.hypotheses == before.hypotheses
    assert after.mechanism_rank == before.mechanism_rank
    assert after.probability == before.probability
    assert after.splice_signal == before.splice_signal
    for original, restored in zip(before.candidates, after.candidates):
        assert isinstance(restored, RealizedEffectCandidate)
        assert restored.effect.short_description == original.effect.short_description
        assert restored.hypotheses == original.hypotheses
        assert restored.probability == original.probability
        assert restored.phase_mass == original.phase_mass
        assert len(restored.outcomes) == len(original.outcomes)
        for left, right in zip(original.outcomes, restored.outcomes):
            assert isinstance(right, ClassifiedOutcome)
            assert right.baseline == left.baseline
            assert right.mutant == left.mutant
            assert right.change_aa_offset == left.change_aa_offset
            assert right.change_cdna_offset == left.change_cdna_offset
    assert before.to_dict() == after.to_dict()


@pytest.fixture
def uncertain_effect():
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084")
    somatic = Variant("7", 117531115, "G", "A", genome)
    germline = Variant("7", 117531100, "T", "A", genome)
    return predict_transcript_model_effect((somatic,), transcript, germline_variants=(germline,))


def test_phase_and_splice_alternatives_keep_products_and_unknown_probabilities(uncertain_effect):
    effect = uncertain_effect
    assert len(effect.candidates) == 8
    assert {state for c in effect.candidates for h in c.hypotheses for _, state in h.phase} == {"cis", "trans"}
    assert effect.probability is None
    assert any(o.mutant.protein_sequence for c in effect.candidates for o in c.outcomes)
    restored = EffectCollection.from_json(EffectCollection([effect]).to_json())[0]
    assert_realized_annotation_equal(effect, restored)
    assert_realized_annotation_equal(effect, type(effect).from_dict(effect.to_dict()))


def test_single_variant_and_legacy_effect_json_remain_readable():
    variant = Variant("7", 117531100, "T", "A", genome=81)
    transcript = variant.genome.transcript_by_id("ENST00000003084")
    legacy = variant.effect_on_transcript(transcript, annotator="fast")
    assert "effect_graph_schema" not in legacy.to_dict()
    assert Substitution.from_json(legacy.to_json()).short_description == "p.L159M"
    effect = variant.effect_on_transcript(transcript, annotator="transcript_model")
    restored = Substitution.from_json(effect.to_json())
    assert_realized_annotation_equal(effect, restored)
    assert restored.candidates[0].outcomes[0].effect is restored


def test_graph_schema_and_references_are_validated(uncertain_effect):
    state = uncertain_effect.to_dict()
    cls = type(uncertain_effect)
    with pytest.raises(ValueError, match="Unsupported.*schema"):
        cls.from_dict({**state, "effect_graph_schema": 99})
    nodes = list(state["nodes"])
    root_cls, fields, attached = nodes[0]
    nodes[0] = (root_cls, fields, {**attached, "most_likely_candidate": _EffectReference(-1)})
    with pytest.raises(ValueError, match="Invalid.*reference"):
        cls.from_dict({**state, "nodes": nodes})
