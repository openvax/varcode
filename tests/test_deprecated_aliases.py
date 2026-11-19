"""Deprecated names that varcode keeps until its next major release, 12.0 (#693).

Varcode itself no longer calls these. Each warns with ``DeprecationWarning``,
and pyproject.toml makes that warning an error anywhere outside
``pytest.warns``. Each test checks that the old name warns once, from the
caller's line, and still reaches the canonical one. Delete a test when its
name is removed.
"""

import warnings
from contextlib import contextmanager

import pandas as pd
import pytest
from pyensembl import cached_release

from varcode import StructuralVariant, Variant, load_vcf, load_vcf_fast
from varcode.effects import StopLoss
from varcode.util import random_variants
from varcode.vcf import dataframes_to_variant_collection

from .data import data_path

ensembl_grch37 = cached_release(75)
ensembl_grch38 = cached_release(81)


@contextmanager
def deprecated(match):
    """Expect one varcode 12.0 deprecation, matching ``match``, from this file."""
    with pytest.warns(DeprecationWarning, match=match) as record:
        yield
    assert [w.filename for w in record if "varcode 12.0" in str(w.message)] == [__file__]


def test_variant_ensembl_keyword_and_property_are_genome():
    with deprecated(r"Variant\(ensembl=\.\.\.\) is deprecated"):
        variant = Variant("1", 10, "A", "G", ensembl=ensembl_grch38)
    assert variant == Variant("1", 10, "A", "G", genome=ensembl_grch38)
    with deprecated(r"Variant\.ensembl is deprecated"):
        assert variant.ensembl is variant.genome


def test_variant_genome_wins_over_ensembl():
    with deprecated(r"Variant\(ensembl="):
        variant = Variant("1", 10, "A", "G", genome=ensembl_grch38, ensembl=ensembl_grch37)
    assert variant.genome == ensembl_grch38


def test_structural_variant_ensembl_keyword_is_genome():
    # One warning: StructuralVariant passes genome= on to Variant.
    with deprecated(r"StructuralVariant\(ensembl=\.\.\.\) is deprecated"):
        sv = StructuralVariant("7", 117_480_000, "DEL", end=117_580_000, alt="<DEL>",
                               ensembl=ensembl_grch38)
    assert sv.genome == ensembl_grch38


def test_variant_kwargs_ensembl_key_is_genome_for_every_row():
    # Structural and ordinary rows alike, with one warning per call.
    frame = pd.DataFrame([
        ("22", 51_179_178, "sv1", "A", "<DEL>", "100", "PASS", "SVTYPE=DEL;END=51179500"),
        ("22", 51_179_700, "snv1", "C", "T", "100", "PASS", "."),
    ], columns=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])

    def parse_info(value):
        return {} if value == "." else {"SVTYPE": "DEL", "END": 51_179_500}

    with deprecated(r"variant_kwargs\['ensembl'\] is deprecated"):
        variants = dataframes_to_variant_collection(
            [frame], source_path="legacy-ensembl.vcf", info_parser=parse_info,
            variant_kwargs={"ensembl": "GRCh37"}, parse_structural_variants=True)
    assert [type(v).__name__ for v in variants] == ["StructuralVariant", "Variant"]
    assert {v.reference_name for v in variants} == {"GRCh37"}


def test_old_json_ensembl_key_loads_without_a_warning():
    # Saved JSON can't be changed by its reader, so the old key is mapped silently.
    state = Variant("1", 10, "A", "G", genome=ensembl_grch38).to_dict()
    state["ensembl"] = state.pop("genome")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        variant = Variant.from_dict(state)
    assert variant.genome == ensembl_grch38
    assert not [w for w in caught if "varcode 12.0" in str(w.message)]


def test_random_variants_ensembl_keyword_is_genome():
    with deprecated(r"random_variants\(ensembl=\.\.\.\) is deprecated"):
        old = random_variants(2, random_seed=0, ensembl=ensembl_grch38)
    new = random_variants(2, random_seed=0, genome=ensembl_grch38)
    assert list(old) == list(new)
    with pytest.raises(ValueError, match="not both"):
        random_variants(2, ensembl=ensembl_grch38, genome=ensembl_grch38)


def test_stop_loss_extended_protein_sequence_is_aa_alt():
    variant = Variant("1", 46501738, "G", "C", genome=ensembl_grch37)
    effect = variant.effect_on_transcript(
        ensembl_grch37.transcript_by_id("ENST00000361297"))
    assert isinstance(effect, StopLoss)
    assert effect.aa_alt
    with deprecated(r"StopLoss\.extended_protein_sequence is deprecated"):
        assert effect.extended_protein_sequence == effect.aa_alt
    assert effect.short_description.endswith(effect.aa_alt + " (stop-loss)")


@pytest.mark.filterwarnings(
    "ignore:Skipped 1 reference-identical ALT allele:UserWarning")
def test_load_vcf_fast_warns_and_forwards_to_load_vcf():
    path = data_path("somatic_hg19_14muts.vcf")
    with pytest.warns(DeprecationWarning, match="renamed to load_vcf"):
        variants = load_vcf_fast(path)
    assert list(variants) == list(load_vcf(path))


def test_realized_annotator_and_module_aliases_warn_and_delegate():
    import pickle
    from varcode import (
        RealizedEffectAnnotator,
        TranscriptModelEffectAnnotator,
        get_annotator,
        get_default_annotator,
        predict_realized_effect,
        predict_transcript_model_effect,
        resolve_annotator,
        set_default_annotator,
        to_proteome,
        use_annotator,
    )
    from varcode.csv_helpers import warn_on_version_drift
    from varcode.realized_effects import (
        RealizedEffectAnnotator as LegacyAnnotator,
        predict_realized_effect as legacy_predict,
    )

    assert RealizedEffectAnnotator is LegacyAnnotator
    assert predict_realized_effect is legacy_predict

    with deprecated(r"RealizedEffectAnnotator is deprecated"):
        annotator = LegacyAnnotator()
    assert isinstance(annotator, TranscriptModelEffectAnnotator)
    assert isinstance(pickle.loads(pickle.dumps(annotator)), TranscriptModelEffectAnnotator)

    transcript = ensembl_grch38.transcript_by_id("ENST00000003084")
    variant = Variant("7", 117531100, "T", "A", genome=ensembl_grch38)
    with deprecated(r"predict_realized_effect is deprecated"):
        legacy_effect = predict_realized_effect([variant], transcript)
    assert legacy_effect == predict_transcript_model_effect([variant], transcript)

    with deprecated(r"Annotator 'realized' is deprecated"):
        assert get_annotator("realized") is get_annotator("transcript_model")
    with deprecated(r"Annotator 'realized' is deprecated"):
        assert resolve_annotator("realized") is get_annotator("transcript_model")

    prior_default = get_default_annotator().name
    try:
        with deprecated(r"Annotator 'realized' is deprecated"):
            set_default_annotator("realized")
        assert get_default_annotator().name == "transcript_model"
    finally:
        set_default_annotator(prior_default)

    with deprecated(r"Annotator 'realized' is deprecated"):
        with use_annotator("realized"):
            assert get_default_annotator().name == "transcript_model"
            assert variant.effects().annotator == "transcript_model"

    with deprecated(r"Annotator 'realized' is deprecated"):
        effects = variant.effects(annotator="realized")
    assert effects.annotator == "transcript_model"

    # Saved provenance naming 'realized' continues to load without warning.
    effect = variant.effect_on_transcript(transcript)
    effect.annotator = "realized"
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        proteome = to_proteome([effect])
        warn_on_version_drift({"annotator": "realized"}, "11.16.0", "legacy.csv")
    assert proteome.proteins[0].producer.tool == "varcode"
    assert not caught


def test_uorf_compatibility_module_warns_and_delegates():
    import orfdata
    from varcode import (
        ORFAnnotation,
        ORFTranscript,
        UORFAnnotation,
        annotate_orf_product,
        annotate_orf_variants,
        annotate_uorf_product,
        annotate_uorf_variants,
        load_uorf_catalog,
        predict_orf_effects,
        predict_uorf_effects,
    )
    from varcode.effect_hypotheses import RealizedTranscriptProduct

    with deprecated(r"load_uorf_catalog is deprecated"):
        catalog = load_uorf_catalog()
    record = catalog.get("c7riboseqorf56")
    transcript = cached_release(95).transcript_by_id(record.transcript_id)
    reference = record.bind(transcript)
    variant = Variant("7", 66205486, "C", "G", genome=transcript.genome)

    with deprecated(r"ORFTranscript is deprecated"):
        view = ORFTranscript(record, transcript, genome_build="GRCh38")
    assert view.reference_orf == reference

    with deprecated(r"annotate_uorf_variants is deprecated"):
        legacy = annotate_uorf_variants([variant], record, transcript, genome_build="GRCh38")
    canonical = annotate_orf_variants([variant], reference)
    assert legacy.effect == canonical.effect
    assert isinstance(canonical, UORFAnnotation)

    with deprecated(r"UORFAnnotation is deprecated"):
        constructed = UORFAnnotation(record=record, effect=canonical.effect, reference=reference)
    assert type(constructed) is ORFAnnotation

    product = RealizedTranscriptProduct(
        transcript, transcript.sequence, None,
        evidence={"base_origins": [None] * len(transcript.sequence)},
    )
    with deprecated(r"annotate_uorf_product is deprecated"):
        legacy_prod = annotate_uorf_product(record, transcript, product, genome_build="GRCh38")
    assert legacy_prod == annotate_orf_product(reference, product)

    subset = orfdata.ORFCatalog([record], catalog.metadata)
    with deprecated(r"predict_uorf_effects is deprecated"):
        legacy_predicted = predict_uorf_effects(variant, catalog=subset)
    assert [a.orf_id for a in legacy_predicted] == [
        a.orf_id for a in predict_orf_effects(variant, catalog=subset)
    ]
