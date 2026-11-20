"""Verify deprecated 11.x aliases were removed in varcode 12.0 while saved
archive/provenance compatibility remains intact (#693, #735).
"""

import importlib
import warnings

import pandas as pd
import pytest
from pyensembl import cached_release

import varcode
from varcode import (
    ORFAnnotation,
    StructuralVariant,
    Variant,
    effect_serialization_classes,
    get_annotator,
    registered_annotators,
    resolve_annotator,
    set_default_annotator,
    to_proteome,
    use_annotator,
)
from varcode.csv_helpers import warn_on_version_drift
from varcode.effects import StopLoss
from varcode.orf import ORFTranscriptView
from varcode.util import random_variants
from varcode.vcf import dataframes_to_variant_collection


ensembl_grch37 = cached_release(75)
ensembl_grch38 = cached_release(81)


@pytest.mark.parametrize("removed_name", [
    "load_vcf_fast",
    "detect_loh",
    "RealizedEffectAnnotator",
    "predict_realized_effect",
    "realized_effects",
    "uorf",
    "UORFAnnotation",
    "ORFTranscript",
    "load_uorf_catalog",
    "annotate_uorf_variants",
    "annotate_uorf_product",
    "predict_uorf_effects",
])
def test_removed_package_exports_are_absent(removed_name):
    assert removed_name not in varcode.__all__
    assert not hasattr(varcode, removed_name)


@pytest.mark.parametrize("removed_module", [
    "varcode.realized_effects",
    "varcode.uorf",
])
def test_removed_legacy_modules_raise_module_not_found(removed_module):
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module(removed_module)


def test_ensembl_keyword_and_property_removed():
    with pytest.raises(TypeError, match="ensembl"):
        Variant("1", 10, "A", "G", ensembl=ensembl_grch38)

    variant = Variant("1", 10, "A", "G", genome=ensembl_grch38)
    assert not hasattr(variant, "ensembl")

    with pytest.raises(TypeError, match="ensembl"):
        StructuralVariant(
            "7", 117_480_000, "DEL", end=117_580_000, alt="<DEL>",
            ensembl=ensembl_grch38)

    with pytest.raises(TypeError, match="ensembl"):
        random_variants(1, random_seed=0, ensembl=ensembl_grch38)

    frame = pd.DataFrame(
        [("22", 51_179_700, "snv1", "C", "T", "100", "PASS", ".")],
        columns=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"],
    )
    with pytest.raises(TypeError, match="ensembl"):
        dataframes_to_variant_collection(
            [frame], source_path="legacy.vcf", info_parser=lambda _: {},
            variant_kwargs={"ensembl": "GRCh37"})


def test_stop_loss_extended_protein_sequence_removed():
    variant = Variant("1", 46501738, "G", "C", genome=ensembl_grch37)
    effect = variant.effect_on_transcript(
        ensembl_grch37.transcript_by_id("ENST00000361297"))
    assert isinstance(effect, StopLoss)
    assert effect.aa_alt
    assert not hasattr(effect, "extended_protein_sequence")


def test_realized_annotator_alias_removed_from_registry():
    assert "realized" not in registered_annotators()
    assert "realized" not in registered_annotators(include_deprecated=True)
    with pytest.raises(KeyError):
        get_annotator("realized")
    with pytest.raises(KeyError):
        resolve_annotator("realized")
    with pytest.raises(KeyError):
        set_default_annotator("realized")
    with pytest.raises(KeyError):
        with use_annotator("realized"):
            pass


def test_legacy_saved_json_and_provenance_still_load():
    # Old Variant JSON storing the genome under "ensembl" still loads silently.
    state = Variant("1", 10, "A", "G", genome=ensembl_grch38).to_dict()
    state["ensembl"] = state.pop("genome")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        restored = Variant.from_dict(state)
    assert restored.genome == ensembl_grch38
    assert not caught

    # Legacy serialized class keys for varcode.uorf remain decodeable via the
    # finite serialization allowlist.
    registry = effect_serialization_classes()
    assert registry[("varcode.uorf", "UORFAnnotation")] is ORFAnnotation
    assert registry[("varcode.uorf", "ORFTranscript")] is ORFTranscriptView

    # Saved provenance naming 'realized' still exports/validates cleanly.
    transcript = ensembl_grch38.transcript_by_id("ENST00000003084")
    variant = Variant("7", 117531100, "T", "A", genome=ensembl_grch38)
    effect = variant.effect_on_transcript(transcript)
    effect.annotator = "realized"
    with warnings.catch_warnings(record=True) as caught_prov:
        warnings.simplefilter("always")
        proteome = to_proteome([effect])
        warn_on_version_drift({"annotator": "realized"}, "12.0.0", "legacy.csv")
    assert proteome.proteins[0].producer.tool == "varcode"
    assert not caught_prov
