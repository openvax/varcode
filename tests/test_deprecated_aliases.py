"""Deprecated names that varcode keeps until its next major release (#693).

Varcode itself no longer calls these; each test checks that the old name still
reaches the canonical one. Delete a test when its name is removed.
"""

import pytest
from pyensembl import cached_release

from varcode import Variant, load_vcf, load_vcf_fast
from varcode.effects import StopLoss
from varcode.util import random_variants

from .data import data_path

ensembl_grch37 = cached_release(75)
ensembl_grch38 = cached_release(81)


def test_variant_ensembl_keyword_and_property_are_genome():
    variant = Variant("1", 10, "A", "G", ensembl=ensembl_grch38)
    assert variant == Variant("1", 10, "A", "G", genome=ensembl_grch38)
    assert variant.ensembl is variant.genome


def test_random_variants_ensembl_keyword_is_genome():
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
    assert effect.extended_protein_sequence == effect.aa_alt
    assert effect.short_description.endswith(effect.aa_alt + " (stop-loss)")


@pytest.mark.filterwarnings(
    "ignore:Skipped 1 reference-identical ALT allele:UserWarning")
def test_load_vcf_fast_warns_and_forwards_to_load_vcf():
    path = data_path("somatic_hg19_14muts.vcf")
    with pytest.warns(DeprecationWarning, match="renamed to load_vcf"):
        variants = load_vcf_fast(path)
    assert list(variants) == list(load_vcf(path))
