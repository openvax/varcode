"""Reproducible SV object archives with self links and attached evidence."""

import pytest
from pathlib import Path
from pyensembl import cached_release

from varcode import EffectCandidate, EffectCollection, StructuralVariant, load_vcf
from varcode.effects import (
    CrypticExonCandidate, GeneFusion, Inversion, LargeDeletion,
    LargeDuplication, StructuralVariantEffect, TranslocationToIntergenic,
)
from .test_sv_consequences import local_model  # isolate individual transcript models


@pytest.fixture
def transcript():
    return cached_release(81).transcript_by_id("ENST00000003084")


def assert_round_trip(effect):
    clone = type(effect).from_json(effect.to_json())
    assert clone.to_dict() == effect.to_dict()
    assert clone.variant == effect.variant
    assert clone.transcript == effect.transcript
    assert clone.mutant_protein_sequence == effect.mutant_protein_sequence
    assert clone.modifies_coding_sequence == effect.modifies_coding_sequence
    assert clone.modifies_protein_sequence == effect.modifies_protein_sequence
    collection = EffectCollection(
        [effect], distinct=True, annotator="fast", annotator_version="test",
        annotated_at="2026-09-23T00:00:00Z")
    restored = EffectCollection.from_json(collection.to_json())
    assert restored[0].to_dict() == effect.to_dict()
    assert restored.annotator == "fast"
    assert restored.annotator_version == "test"
    assert restored.annotated_at == collection.annotated_at
    return clone


@pytest.mark.parametrize("kind", ["DEL", "DUP", "INV", "CNV", "INS"])
@pytest.mark.parametrize("assembly", [None, "ATGAAATAG"])
@pytest.mark.usefixtures("local_model")
def test_predicted_sv_round_trip(transcript, kind, assembly):
    exon = transcript.exons[4]
    variant = StructuralVariant(
        "7", exon.start + 10, kind, end=exon.end - 10,
        genome=transcript.genome, alt_assembly=assembly)
    effect = variant.effect_on_transcript(transcript)
    restored = assert_round_trip(effect)
    assert restored.affected_exons == effect.affected_exons
    assert restored.mutant_transcript == effect.mutant_transcript
    assert [(c.source, c.evidence) for c in restored.candidates] == [
        (c.source, c.evidence) for c in effect.candidates]


@pytest.mark.usefixtures("local_model")
def test_direct_consequence_preserves_attached_model(transcript):
    exon = transcript.exons[4]
    variant = StructuralVariant("7", exon.start - 20, "DEL", end=exon.end + 20,
                                genome=transcript.genome)
    effect = variant.effect_on_transcript(transcript)
    assert not isinstance(effect, StructuralVariantEffect)
    restored = assert_round_trip(effect)
    assert restored.mutant_transcript == effect.mutant_transcript
    assert restored.affected_exons == effect.affected_exons


@pytest.mark.parametrize("cls", [StructuralVariantEffect, LargeDeletion,
                                  LargeDuplication, Inversion, TranslocationToIntergenic])
def test_legacy_self_candidate_round_trip(transcript, cls):
    variant = StructuralVariant("7", 117531100, "DEL", end=117531200, genome=transcript.genome)
    kwargs = {"affected_exons": transcript.exons[:1]} if cls in (LargeDeletion, LargeDuplication) else {}
    effect = cls(variant, transcript, **kwargs)
    restored = assert_round_trip(effect)
    assert restored.candidates[0].effect is restored


def test_fusion_all_candidate_sources_and_shared_self_links(transcript):
    partner = transcript.genome.transcript_by_id("ENST00000357654")
    variant = StructuralVariant("7", transcript.exons[4].end, "BND", genome=transcript.genome,
                                mate_contig="17", mate_start=partner.exons[4].start)
    effect = GeneFusion(variant, transcript, partner)
    other = TranslocationToIntergenic(variant, transcript)
    effect._attach_primary_effects((other,))
    effect._attach_cryptic_candidates((
        CrypticExonCandidate(variant, "7", 123, 145, donor_score=0.8, acceptor_score=0.9),))
    effect._attach_splice_outcomes((EffectCandidate(other, "varcode_splice", {"reads": 2}),))
    effect._extra_candidates = (EffectCandidate(effect, "RNA", {"witnesses": ["read1"]}),)
    effect.candidate_evidence = {"assumption": "reference_splicing"}
    restored = assert_round_trip(effect)
    assert restored.partner_transcript == partner
    assert restored.five_prime_transcript == transcript
    assert restored.three_prime_transcript == partner
    assert restored.candidates[0].effect is restored
    assert restored.candidates[-1].effect is restored
    assert restored.candidates[1].effect is restored.candidates[-2].effect
    assert [c.source for c in restored.candidates] == [
        "varcode", "varcode", "varcode_motif", "varcode_splice", "RNA"]


def test_unknown_schema_rejected(transcript):
    with pytest.raises(ValueError, match="Unsupported.*schema"):
        StructuralVariantEffect.from_dict({"structural_effect_schema": 999})


def test_esvee_predicted_fusions_and_collection_round_trip():
    variants = load_vcf(
        str(Path(__file__).parent / "data/osteosarc_esvee_somatic.vcf"),
        genome=81, parse_structural_variants=True)
    effects = variants.effects()
    assert any(isinstance(effect, GeneFusion) for effect in effects)
    restored = EffectCollection.from_json(effects.to_json())
    assert len(restored) == len(effects)
    for before, after in zip(effects, restored):
        assert type(before) is type(after)
        assert before.to_dict() == after.to_dict()
        assert before.mutant_protein_sequence == after.mutant_protein_sequence
