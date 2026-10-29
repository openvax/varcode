"""RNA refines existing composed products without manufacturing phase or sequence."""

import pickle

import pytest
from pyensembl import cached_release

from varcode import (
    EffectCollection, GermlineContext, MutantTranscript, RNAEvidence, Variant,
    apply_rna_evidence_to_effects, make_rna_outcome,
)
from varcode.effects import Deletion, Intronic, Unresolved


@pytest.fixture
def inputs():
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084")
    return transcript, Variant("7", 117531115, "G", "A", genome)


def observation(effect, model_id="observed", read_count=7):
    return make_rna_outcome(effect, source="reads", transcript_model_id=model_id,
                            read_count=read_count)


def test_rna_support_selects_matching_normal_product_and_keeps_alternatives(inputs):
    transcript, variant = inputs
    original = variant.effect_on_transcript(transcript, annotator="transcript_model")
    normal = next(c.effect for c in original.candidates if isinstance(c.effect, Intronic))
    evidence = RNAEvidence([observation(normal)])
    results = variant.effects(annotator="transcript_model", rna_resolver=evidence,
                              transcript_filter=lambda t: t.id == transcript.id)
    result, = results
    assert isinstance(original, Deletion)
    assert isinstance(result, Intronic)
    assert result.mutant_transcript.cdna_sequence == transcript.sequence
    assert result.mutant_protein_sequence == transcript.protein_sequence
    assert len(result.candidates) == len(original.candidates)
    assert result.candidates[0].outcomes[0].mutant.evidence["rna_support"] == (
        {**normal.mutant_transcript.evidence,
         "source": "reads", "transcript_model_id": "observed", "read_count": 7,
         "variant": variant},)
    assert result.candidates[0].probability == original.candidates[-1].probability
    assert all("rna_support" not in o.mutant.evidence
               for c in original.candidates for o in c.outcomes)
    assert any(isinstance(c.effect, Deletion) for c in result.candidates)
    assert not hasattr(result, "_extra_candidates")


def test_same_rna_support_preserves_distinct_patient_baselines(inputs):
    transcript, variant = inputs
    germline = Variant("7", 117531113, "A", "C", transcript.genome)
    context = GermlineContext.from_variants([germline])
    original = variant.effect_on_transcript(
        transcript, annotator="transcript_model", germline=context)
    skipped = next(c.effect for c in original.candidates if isinstance(c.effect, Deletion))
    results = variant.effects(
        annotator="transcript_model", germline=context,
        rna_resolver=RNAEvidence([observation(skipped)]),
        transcript_filter=lambda t: t.id == transcript.id)
    supported = [o for c in results[0].candidates for o in c.outcomes
                 if o.mutant.evidence.get("rna_support")]
    assert len(supported) == 2
    assert {o.hypothesis.phase for o in supported} == {((germline, "cis"),), ((germline, "trans"),)}
    assert len({o.baseline.cdna_sequence for o in supported}) == 2
    assert len({o.mutant.cdna_sequence for o in supported}) == 1
    assert all(o.hypothesis.evidence["phase_state"] == "unknown" for o in supported)


@pytest.mark.parametrize("status", ["junction_fragment", "partial_start"])
def test_partial_rna_remains_an_observed_product_without_reference_completion(inputs, status):
    transcript, variant = inputs
    original = variant.effect_on_transcript(transcript, annotator="transcript_model")
    model = MutantTranscript.from_sequence(
        "ATGGCTTAA", reference_transcript=transcript, mutant_protein_sequence="MA",
        annotator_name="reads", evidence={"sequence_status": status})
    observed = Unresolved(variant, transcript, "observed_rna")
    observed.mutant_transcript = model
    evidence = RNAEvidence([observation(observed)])
    results = EffectCollection([original])
    apply_rna_evidence_to_effects(results, evidence)
    assert isinstance(results[0], Deletion)
    candidate = results[0].candidates[-1]
    assert isinstance(candidate.effect, Unresolved)
    assert candidate.effect.mechanism == "rna_transcript_comparison"
    assert candidate.effect.modifies_protein_sequence is None
    assert candidate.effect.original_protein_sequence is None
    assert candidate.effect.mutant_transcript is model
    assert candidate.outcomes[0].mutant.cdna_sequence == "ATGGCTTAA"
    assert candidate.outcomes[0].baseline.cdna_sequence is None
    assert candidate.hypotheses[0].phase == ()
    assert not hasattr(results[0], "observed_mutant_transcripts")
    apply_rna_evidence_to_effects(results, evidence)
    assert len(results[0].candidates) == len(original.candidates) + 1
    for restored in (EffectCollection.from_json(results.to_json()),
                     pickle.loads(pickle.dumps(results))):
        candidate = restored[0].candidates[-1]
        assert candidate.effect.mutant_transcript.cdna_sequence == "ATGGCTTAA"
        assert candidate.effect.mutant_transcript.evidence["sequence_status"] == status
        assert candidate.effect.original_protein_sequence is None
        assert candidate.hypotheses[0].evidence["rna_observation"]["read_count"] == 7


def test_reference_completion_is_not_full_observed_rna_support(inputs):
    transcript, variant = inputs
    original = variant.effect_on_transcript(transcript, annotator="transcript_model")
    normal = next(c.effect for c in original.candidates if isinstance(c.effect, Intronic))
    assumed = make_rna_outcome(normal, source="varcode_reference_completion",
                              extra_evidence={"reference_completion": {"assumed_spans": (0, 10)}})
    results = EffectCollection([original])
    apply_rna_evidence_to_effects(results, RNAEvidence([assumed]))
    assert isinstance(results[0], Deletion)
    assert not any(o.mutant.evidence.get("rna_support")
                   for c in results[0].candidates for o in c.outcomes)
    assert isinstance(results[0].candidates[-1].effect, Unresolved)


def test_joint_product_uses_rna_evidence_attached_to_nonprimary_member(inputs):
    from types import SimpleNamespace
    from varcode import VariantCollection

    transcript, _ = inputs
    variants = [Variant("7", 117530914, "G", "T", transcript.genome),
                Variant("7", 117530915, "T", "A", transcript.genome)]
    options = dict(annotator="transcript_model",
                   phase_resolver=SimpleNamespace(in_cis=lambda *args, **kwargs: True),
                   transcript_filter=lambda t: t.id == transcript.id)
    original, = [e for e in VariantCollection(variants).effects(**options)
                 if len(getattr(e, "variants", ())) == 2]
    model = original.mutant_transcript
    observed = Unresolved(variants[1], transcript, "observed_rna")
    observed.mutant_transcript = model
    result, = [e for e in VariantCollection(variants).effects(
        **options, rna_resolver=RNAEvidence([observation(observed)]))
        if len(getattr(e, "variants", ())) == 2]
    support = result.candidates[0].outcomes[0].mutant.evidence["rna_support"]
    assert support[0]["variant"] == variants[1]
    assert result.variants == tuple(variants)
    assert result.mutant_transcript.cdna_sequence == model.cdna_sequence
    assert result.mutant_protein_sequence == model.mutant_protein_sequence


def test_phase_resolver_partial_model_uses_same_refinement_as_rna_resolver(inputs):
    from types import SimpleNamespace

    transcript, variant = inputs
    model = MutantTranscript.from_sequence(
        "ATGGCTTAA", reference_transcript=transcript, mutant_protein_sequence="MA",
        annotator_name="reads", evidence={"sequence_status": "junction_fragment"})
    effects = variant.effects(
        annotator="transcript_model",
        phase_resolver=SimpleNamespace(mutant_transcript=lambda *args: model),
        transcript_filter=lambda t: t.id == transcript.id)
    result, = effects
    assert isinstance(result, Deletion)
    assert result.mutant_transcript.cdna_sequence != model.cdna_sequence
    candidate = result.candidates[-1]
    assert isinstance(candidate.effect, Unresolved)
    assert candidate.effect.mutant_transcript is model
    assert candidate.outcomes[0].baseline.cdna_sequence is None
