"""End-to-end tests for the experimental transcript model."""

import pickle
import subprocess
import sys

import pytest
from pyensembl import cached_release

from varcode import (
    Completeness,
    GermlineContext,
    StructuralVariant,
    TranscriptModelEffectAnnotator,
    Variant,
    get_default_annotator,
    get_annotator,
    predict_transcript_model_effect,
    use_annotator,
)
from varcode.effects.effect_classes import (
    Deletion,
    GeneFusion,
    Insertion,
    Substitution,
    Unresolved,
)

from .test_splice_graph import (
    LongIntronTranscript,
    Transcript,
    deletion,
    long_provider,
    point,
    provider,
)


ensembl_grch38 = cached_release(81)
CFTR_TRANSCRIPT_ID = "ENST00000003084"


def test_transcript_model_annotator_is_registered_but_not_default(request):
    annotator = get_annotator("transcript_model")

    assert annotator.name == "transcript_model"
    assert isinstance(annotator, TranscriptModelEffectAnnotator)
    configured = request.config.getoption("--annotator") or "fast"
    assert get_default_annotator() is get_annotator(configured)


def test_transcript_model_import_does_not_promote_experiment():
    subprocess.run([
        sys.executable, "-c",
        "import varcode.transcript_model; "
        "from varcode import get_default_annotator; "
        "assert get_default_annotator().name == 'fast'",
    ], check=True)


def test_realized_names_are_compatibility_aliases():
    from varcode import RealizedEffectAnnotator, predict_realized_effect
    from varcode.realized_effects import (
        RealizedEffectAnnotator as LegacyAnnotator,
        predict_realized_effect as legacy_predict,
    )

    assert get_annotator("realized") is get_annotator("transcript_model")
    assert RealizedEffectAnnotator is LegacyAnnotator is TranscriptModelEffectAnnotator
    assert predict_realized_effect is legacy_predict is predict_transcript_model_effect
    assert isinstance(pickle.loads(pickle.dumps(LegacyAnnotator())),
                      TranscriptModelEffectAnnotator)


@pytest.mark.parametrize("name", ["transcript_model", "realized"])
def test_transcript_model_selection_and_provenance(name):
    from varcode import EffectCollection

    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    variant = Variant("7", 117531100, "T", "A", ensembl_grch38)
    prior_default = get_default_annotator()
    with use_annotator(name):
        assert get_default_annotator().name == "transcript_model"
        effect = variant.effect_on_transcript(transcript)
        effects = variant.effects(raise_on_error=True)
    assert get_default_annotator() is prior_default
    assert type(effect) is Substitution
    assert effect.short_description == "p.L159M"
    assert effects.annotator == "transcript_model"
    assert EffectCollection.from_json(effects.to_json()).annotator == "transcript_model"
    assert variant.effects(annotator=name).annotator == "transcript_model"


def test_weak_exonic_variant_keeps_normal_coding_effect_first():
    transcript = Transcript()
    variant = point(9, "A", "G")

    result = predict_transcript_model_effect(
        (variant,), transcript, sequence_provider=provider())

    assert type(result) is Substitution
    assert result.aa_ref == "Q"
    assert result.aa_alt == "R"
    assert all(rank == 0 for rank in result.mechanism_rank)
    assert result.probability is None


def test_strong_donor_variant_returns_ordinary_deletion_as_top_effect():
    transcript = Transcript()
    variant = point(11, "G", "A")

    result = predict_transcript_model_effect(
        (variant,), transcript, sequence_provider=provider())

    assert type(result) is Deletion
    assert result.aa_ref == "Q"
    assert [
        type(candidate.effect).__name__ for candidate in result.candidates
    ] == ["Deletion", "FrameShift", "Unresolved", "Intronic"]
    unresolved = [
        candidate.effect for candidate in result.candidates
        if isinstance(candidate.effect, Unresolved)]
    assert [effect.mechanism for effect in unresolved] == ["cryptic"]


def test_same_mutant_with_different_phase_baselines_does_not_merge():
    transcript = Transcript()
    somatic = point(11, "G", "A")
    germline = point(9, "A", "G")

    result = predict_transcript_model_effect(
        (somatic,), transcript, germline_variants=(germline,),
        sequence_provider=provider())

    deletion_candidates = [
        candidate for candidate in result.candidates
        if type(candidate.effect) is Deletion]
    assert len(deletion_candidates) == 2
    assert {
        candidate.outcomes[0].baseline.protein_sequence
        for candidate in deletion_candidates
    } == {"MQA", "MRA"}
    assert {
        candidate.outcomes[0].mutant.protein_sequence
        for candidate in deletion_candidates
    } == {"MA"}


def test_sv_deletion_does_not_merge_away_patient_specific_baseline():
    transcript = Transcript()
    somatic = deletion(8, 10)
    germline = point(9, "A", "G")

    result = predict_transcript_model_effect(
        (somatic,), transcript, germline_variants=(germline,),
        sequence_provider=provider())

    assert type(result) is Deletion
    deletion_candidates = [
        candidate for candidate in result.candidates
        if type(candidate.effect) is Deletion]
    assert len(deletion_candidates) == 2
    assert {
        candidate.outcomes[0].baseline.protein_sequence
        for candidate in deletion_candidates
    } == {"MQA", "MRA"}


def test_multiple_somatic_variants_are_realized_on_one_product():
    transcript = Transcript()
    variants = (point(8, "C", "A"), point(9, "A", "G"))

    result = predict_transcript_model_effect(
        variants, transcript, sequence_provider=provider())

    assert type(result) is Substitution
    assert result.aa_ref == "Q"
    assert result.aa_alt == "R"
    assert result.variants == variants


def test_combined_phase_and_splice_hypotheses_obey_global_cap():
    transcript = Transcript()
    somatic = point(11, "G", "A")
    germline = point(9, "A", "G")

    with pytest.raises(ValueError, match="Combined phase/splice"):
        predict_transcript_model_effect(
            (somatic,), transcript, germline_variants=(germline,),
            sequence_provider=provider(), max_hypotheses=7)


def test_real_cftr_missense_matches_existing_varcode_result():
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    variant = Variant("7", 117531100, "T", "A", ensembl_grch38)

    result = predict_transcript_model_effect((variant,), transcript)

    assert type(result) is Substitution
    assert result.short_description == "p.L159M"
    assert len(result.candidates) == 1


def test_effects_dispatches_realized_annotator_with_germline_context():
    somatic = Variant("7", 117531100, "T", "A", ensembl_grch38)
    germline = Variant("7", 117531101, "T", "C", ensembl_grch38)
    context = GermlineContext.from_variants(
        [germline], reference_name="GRCh38")

    effects = somatic.effects(annotator="transcript_model", germline=context)
    result = next(
        effect for effect in effects
        if effect.transcript.id == CFTR_TRANSCRIPT_ID)

    assert type(result) is Substitution
    assert len(result.candidates) == 2
    assert {
        candidate.outcomes[0].baseline.protein_sequence[158]
        for candidate in result.candidates
    } == {"L", "S"}


def test_transcript_model_dispatch_preserves_sparse_germline_uncertainty_flag():
    somatic = Variant("7", 117531100, "T", "A", ensembl_grch38)
    distant = Variant("7", 117587799, "C", "T", ensembl_grch38)
    context = GermlineContext.from_variants(
        [distant], reference_name="GRCh38",
        completeness=Completeness.SPARSE)

    effects = somatic.effects(annotator="transcript_model", germline=context)
    result = next(
        effect for effect in effects
        if effect.transcript.id == CFTR_TRANSCRIPT_ID)

    assert result.germline_unknown is True


def test_real_cftr_donor_plus_one_enumerates_splice_mechanisms():
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)

    result = predict_transcript_model_effect((variant,), transcript)

    assert type(result) is Deletion
    assert [
        candidate.effect.mechanism
        for candidate in result.candidates
        if isinstance(candidate.effect, Unresolved)
    ] == ["intron_retention", "cryptic"]
    assert result.probability is None


@pytest.mark.parametrize(
    "sv_type,expected_class",
    [("DEL", Deletion), ("DUP", Insertion)])
def test_real_cftr_whole_exon_sv_is_classified_from_realized_product(
        sv_type, expected_class):
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    exon = transcript.exons[4]
    variant = StructuralVariant(
        contig="7",
        start=exon.start - 300,
        end=exon.end + 300,
        sv_type=sv_type,
        genome=ensembl_grch38)

    result = predict_transcript_model_effect((variant,), transcript)

    assert type(result) is expected_class
    assert len(result.candidates) == 1


@pytest.mark.parametrize(
    "reverse,position,ref,alt",
    [(False, 23, "G", "A"), (True, 178, "C", "T")])
def test_cryptic_donor_is_realized_on_mutated_haplotype_both_strands(
        reverse, position, ref, alt):
    transcript = LongIntronTranscript(reverse=reverse)
    variant = point(position, ref, alt)

    result = predict_transcript_model_effect(
        (variant,), transcript, sequence_provider=long_provider(reverse))

    cryptic = next(
        candidate for candidate in result.candidates
        if candidate.ordinal_key == (2,))
    assert not isinstance(cryptic.effect, Unresolved)
    evidence = cryptic.outcomes[0].mutant.evidence
    resolved = evidence["resolved_splice_choices"]
    assert len(resolved) == 1
    assert resolved[0][1] == "cryptic"
    assert resolved[0][3] > 0


def test_transcript_model_annotator_routes_bnd_through_verified_fusion_builder():
    genome = cached_release(95)
    otx1 = genome.transcript_by_id("ENST00000282549")
    variant = StructuralVariant(
        contig="2",
        start=63_053_516,
        sv_type="BND",
        alt="N]2:25955847]",
        mate_contig="2",
        mate_start=25_955_847,
        mate_orientation="]]",
        genome=genome)

    result = predict_transcript_model_effect((variant,), otx1)

    assert isinstance(result, GeneFusion)
    assert result.partner_transcript.id == "ENST00000264712"
    assert result.mutant_protein_sequence == (
        "MMSYLKQPPYGMNGLGLAGPAMDLLHPSVGYPETS")
