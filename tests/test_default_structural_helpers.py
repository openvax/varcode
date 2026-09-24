"""One structural implementation, without a second registered annotator."""

import importlib
import sys

import pytest
from pyensembl import cached_release

import varcode
from varcode import FastEffectAnnotator, StructuralVariant, get_annotator
from varcode.effects import (
    StructuralVariantEffect,
    TranslocationToIntergenic,
)
from varcode.effects.structural import predict_structural_variant_effect


@pytest.fixture
def transcript():
    return cached_release(81).transcript_by_id("ENST00000003084")


def structural_variant(kind, assembly=None):
    kwargs = {}
    if kind == "BND":
        kwargs = dict(
            alt="N]22:15500000]", mate_contig="22", mate_start=15_500_000,
            mate_orientation="]]")
    return StructuralVariant(
        "7", 117_531_100, sv_type=kind, end=117_531_200,
        genome=cached_release(81), alt_assembly=assembly, **kwargs)


@pytest.mark.parametrize("kind, expected", [
    ("DEL", StructuralVariantEffect),
    ("DUP", StructuralVariantEffect),
    ("INV", StructuralVariantEffect),
    ("CNV", StructuralVariantEffect),
    ("INS", StructuralVariantEffect),
    ("BND", TranslocationToIntergenic),
])
@pytest.mark.parametrize("assembly", [None, "ATG" + "GCT" * 20 + "TAA"])
def test_default_uses_structural_helpers(kind, expected, assembly, transcript, monkeypatch):
    monkeypatch.setitem(sys.modules, "varcode.annotators.structural_variant", None)
    variant = structural_variant(kind, assembly)
    effect = variant.effect_on_transcript(transcript, annotator="fast")
    direct = predict_structural_variant_effect(variant, transcript)
    assert isinstance(effect, expected)
    assert type(effect) is type(direct)
    assert effect.short_description == direct.short_description
    actual_mt = effect.mutant_transcript
    expected_mt = direct.mutant_transcript
    if actual_mt is None:
        assert kind in ("INS", "CNV") and assembly is None
        assert expected_mt is None
        assert effect.modifies_protein_sequence is None
        return
    assert actual_mt.annotator_name == expected_mt.annotator_name == "structural_variant"
    assert actual_mt.cdna_sequence == expected_mt.cdna_sequence
    assert actual_mt.mutant_protein_sequence == expected_mt.mutant_protein_sequence
    assert actual_mt.evidence == expected_mt.evidence
    assert [
        (s.source.sequence, s.start, s.end, s.strand, s.label)
        for s in actual_mt.reference_segments
    ] == [
        (s.source.sequence, s.start, s.end, s.strand, s.label)
        for s in expected_mt.reference_segments
    ]
    assert [
        (type(c.effect), c.short_description, c.source, c.evidence)
        for c in effect.candidates
    ] == [
        (type(c.effect), c.short_description, c.source, c.evidence)
        for c in direct.candidates
    ]
    if assembly:
        assert actual_mt.cdna_sequence == assembly


@pytest.mark.parametrize("name", ["StructuralVariantAnnotator", "UnsupportedVariantError"])
def test_retired_public_apis_are_removed(name):
    from varcode import annotators
    from varcode.annotators import registry

    for module in (varcode, annotators, registry):
        assert not hasattr(module, name)
        assert name not in getattr(module, "__all__", ())


def test_retired_structural_module_is_removed():
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("varcode.annotators.structural_variant")


@pytest.mark.parametrize("resolve", [get_annotator, varcode.resolve_annotator])
def test_structural_registry_name_is_not_an_alias(resolve):
    with pytest.raises(KeyError, match="structural_variant"):
        resolve("structural_variant")
    assert isinstance(get_annotator("fast"), FastEffectAnnotator)
