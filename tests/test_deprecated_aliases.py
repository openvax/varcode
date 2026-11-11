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
