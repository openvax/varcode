"""Unified default routing and input-dependent experimental refusal."""

import pytest
from pyensembl import cached_release

from varcode import (
    EffectAnnotator,
    EffectCollection,
    FastEffectAnnotator,
    GermlineContext,
    StructuralVariant,
    Variant,
    VariantCollection,
    get_annotator,
    use_annotator,
)
from varcode.effects import Failure, LargeDeletion, Silent, Unresolved
from varcode.effects.effect_prediction import predict_variant_effect_on_transcript


@pytest.fixture
def transcript():
    return cached_release(81).transcript_by_id("ENST00000003084")


@pytest.fixture
def point_variant():
    return Variant("7", 117531115, "G", "A", cached_release(81))


def structural_variant(kind="DEL"):
    return StructuralVariant(
        "7", 117531100, sv_type=kind, end=117531200,
        genome=cached_release(81))


class PartialAnnotator:
    name = "test_partial"
    version = "1"

    def annotate_on_transcript(self, variant, transcript):
        if getattr(variant, "is_structural", False):
            return NotImplemented
        return Silent(variant, transcript, aa_pos=0, aa_ref="M")


def assert_unsupported(effect, name):
    assert isinstance(effect, Unresolved)
    assert effect.mechanism == "unsupported_annotation"
    assert name in effect.reason


@pytest.mark.parametrize("selection", [None, "fast", FastEffectAnnotator()])
def test_default_routes_sv_for_all_selection_forms(selection, transcript):
    variant = structural_variant()
    with use_annotator("fast"):
        direct = variant.effect_on_transcript(transcript, annotator=selection)
        effects = variant.effects(annotator=selection, raise_on_error=True)
    assert isinstance(direct, LargeDeletion)
    assert isinstance(effects[0], LargeDeletion)
    assert effects.annotator == "fast"
    assert isinstance(FastEffectAnnotator().annotate_on_transcript(
        variant, transcript), LargeDeletion)


def test_low_level_public_predictor_routes_sv(transcript):
    result = predict_variant_effect_on_transcript(
        structural_variant(), transcript, annotator="fast")
    assert isinstance(result, LargeDeletion)


def test_partial_plugin_needs_no_capability_metadata(point_variant, transcript):
    annotator = PartialAnnotator()
    assert isinstance(annotator, EffectAnnotator)
    with use_annotator(annotator):
        assert get_annotator(annotator.name) is annotator
        assert isinstance(point_variant.effect_on_transcript(transcript), Silent)
        assert_unsupported(
            structural_variant().effect_on_transcript(transcript), annotator.name)


@pytest.mark.parametrize("raise_on_error", [False, True])
def test_partial_results_keep_unknowns_and_provenance(point_variant, raise_on_error):
    variants = VariantCollection([point_variant, structural_variant()])
    effects = variants.effects(
        annotator=PartialAnnotator(), raise_on_error=raise_on_error)
    assert any(isinstance(effect, Silent) for effect in effects)
    unknowns = [effect for effect in effects if isinstance(effect, Unresolved)]
    assert unknowns
    for effect in unknowns:
        assert_unsupported(effect, "test_partial")
    assert effects.annotator == "test_partial"
    assert effects.annotator_version == "1"
    restored = EffectCollection.from_json(effects.to_json())
    assert restored.annotator == effects.annotator
    assert [effect.reason for effect in restored if isinstance(effect, Unresolved)] == [
        effect.reason for effect in unknowns]


@pytest.mark.parametrize("kind", ["DEL", "DUP", "INV", "INS", "CNV", "BND"])
def test_protein_diff_declines_structural_inputs(kind, transcript):
    variant = structural_variant(kind)
    annotator = get_annotator("protein_diff")
    assert annotator.annotate_on_transcript(variant, transcript) is NotImplemented
    assert_unsupported(variant.effect_on_transcript(
        transcript, annotator=annotator), "protein_diff")


def test_scoped_experiment_is_not_overridden_for_sv(transcript):
    variant = structural_variant()
    with use_annotator("protein_diff"):
        assert_unsupported(variant.effect_on_transcript(transcript), "protein_diff")
        effects = variant.effects(raise_on_error=True)
        assert effects.annotator == "protein_diff"
        assert all(isinstance(effect, Unresolved) for effect in effects)


def test_refusal_can_depend_on_transcript(point_variant, transcript):
    class TranscriptSubset(PartialAnnotator):
        def annotate_on_transcript(self, variant, selected_transcript):
            if selected_transcript.id == transcript.id:
                return NotImplemented
            return super().annotate_on_transcript(variant, selected_transcript)

    effects = point_variant.effects(annotator=TranscriptSubset(), raise_on_error=True)
    assert any(isinstance(effect, Silent) for effect in effects)
    refused = [effect for effect in effects if isinstance(effect, Unresolved)]
    assert [effect.transcript.id for effect in refused] == [transcript.id]


def test_missing_context_hook_does_not_ignore_germline(point_variant, transcript):
    annotator = PartialAnnotator()
    context = GermlineContext.from_variants([point_variant])
    result = point_variant.effect_on_transcript(
        transcript, annotator=annotator, germline=context)
    assert_unsupported(result, annotator.name)
    assert "germline context" in result.reason
    assert isinstance(point_variant.effect_on_transcript(
        transcript, annotator=annotator, germline=GermlineContext.empty()), Silent)


def test_context_hook_can_decline_without_fallback(point_variant, transcript):
    class ContextSubset(PartialAnnotator):
        def annotate_with_context(
                self, variant, selected_transcript, germline_ctx, phase_resolver=None):
            assert germline_ctx is context
            assert phase_resolver is resolver
            return NotImplemented

    context = GermlineContext.from_variants([point_variant])
    resolver = object()
    result = point_variant.effect_on_transcript(
        transcript, annotator=ContextSubset(), germline=context,
        phase_resolver=resolver)
    assert_unsupported(result, "test_partial")


def test_default_declines_unimplemented_sv_germline_composition(point_variant, transcript):
    result = structural_variant().effect_on_transcript(
        transcript, annotator="fast",
        germline=GermlineContext.from_variants([point_variant]))
    assert_unsupported(result, "fast")


def test_germline_internal_prediction_does_not_use_scoped_experiment(
        point_variant, transcript, monkeypatch):
    from varcode.germline import PhaseHypothesis, _classify_against_patient_baseline

    # Force the established default's joint-build fallback. It must remain
    # an internal default prediction even if an experiment is globally selected.
    monkeypatch.setattr(
        "varcode.mutant_transcript.apply_variants_to_transcript",
        lambda *args, **kwargs: None)
    hypothesis = PhaseHypothesis(cis=(), trans=(), phase_state="known")
    expected = point_variant.effect_on_transcript(transcript, annotator="fast")

    class Declining(PartialAnnotator):
        def annotate_on_transcript(self, variant, transcript):
            return NotImplemented

    with use_annotator(Declining()):
        actual = _classify_against_patient_baseline(
            point_variant, transcript, hypothesis)
    assert type(actual) is type(expected)
    assert actual.short_description == expected.short_description


def test_germline_helper_preserves_a_subclass_refusal(point_variant, transcript):
    from varcode import Completeness

    class Declining(FastEffectAnnotator):
        name = "test_declining"

        def annotate_on_transcript(self, variant, transcript):
            return NotImplemented

    context = GermlineContext.from_variants(
        [Variant("7", 117531050, "A", "G", cached_release(81))],
        completeness=Completeness.SPARSE)
    result = point_variant.effect_on_transcript(
        transcript, annotator=Declining(), germline=context)
    assert_unsupported(result, "test_declining")


@pytest.mark.parametrize("with_context", [False, True])
def test_transcript_model_declines_unresolved_cnv_direction(transcript, point_variant, with_context):
    context = GermlineContext.from_variants([point_variant]) if with_context else None
    result = structural_variant("CNV").effect_on_transcript(
        transcript, annotator="transcript_model", germline=context)
    assert_unsupported(result, "transcript_model")


def test_transcript_model_reports_missing_insertion_sequence(transcript):
    result = structural_variant("INS").effect_on_transcript(
        transcript, annotator="transcript_model")
    assert isinstance(result, Unresolved)
    assert result.mechanism == "sequence_unavailable"
    assert "alt_assembly" in result.reason


def test_errors_are_not_treated_as_unsupported(point_variant):
    class Broken(PartialAnnotator):
        def annotate_on_transcript(self, variant, transcript):
            raise ValueError("broken implementation")

    with pytest.raises(ValueError, match="broken implementation"):
        point_variant.effects(annotator=Broken(), raise_on_error=True)
    effects = point_variant.effects(annotator=Broken(), raise_on_error=False)
    assert effects
    assert all(isinstance(effect, Failure) for effect in effects)


def test_none_is_a_plugin_bug_not_an_unknown(point_variant):
    class Broken(PartialAnnotator):
        def annotate_on_transcript(self, variant, transcript):
            return None

    with pytest.raises(TypeError, match="MutationEffect or NotImplemented"):
        point_variant.effects(annotator=Broken())
