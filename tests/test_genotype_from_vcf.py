# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
End-to-end tests for VariantCollection genotype/zygosity access —
exercises the full path from VCF → Variant → Genotype (openvax/varcode#267).
"""

import os
import tempfile

import pytest

from varcode import (
    Genotype,
    SampleNotFoundError,
    Zygosity,
    load_vcf,
)


VCF_BODY = """##fileformat=VCFv4.1
##reference=GRCh38
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	tumor	normal
17	43082575	.	C	T	100	PASS	.	GT:AD:DP:GQ	0/1:10,5:15:99	0/0:20,0:20:99
7	117531114	.	G	T	100	PASS	.	GT:AD:DP:GQ	1/1:0,20:20:99	0/1:10,10:20:99
1	100	.	A	T,G	100	PASS	.	GT:AD:DP:GQ	1/2:5,5,5:15:99	0/1:10,5,0:15:99
17	43082576	.	C	A	100	PASS	.	GT:AD:DP:GQ:PS	0|1:8,7:15:99:100	./.:.:.:.:.
"""


@pytest.fixture
def multi_sample_vcf():
    fd, path = tempfile.mkstemp(suffix=".vcf")
    with os.fdopen(fd, "w") as f:
        f.write(VCF_BODY)
    yield path
    os.unlink(path)


# -------------------------------------------------------------------
# Sample discovery
# -------------------------------------------------------------------


def test_samples_property_lists_names(multi_sample_vcf):
    vc = load_vcf(multi_sample_vcf, genome="GRCh38")
    assert vc.samples == ["normal", "tumor"]


def test_has_sample_data_true_for_vcf_with_samples(multi_sample_vcf):
    vc = load_vcf(multi_sample_vcf, genome="GRCh38")
    assert vc.has_sample_data() is True


def test_has_sample_data_false_for_directly_constructed(multi_sample_vcf):
    # Directly-constructed variants have no sample_info metadata.
    from varcode import Variant, VariantCollection
    vc = VariantCollection(variants=[
        Variant("17", 43082575, "C", "T", "GRCh38"),
    ])
    assert vc.has_sample_data() is False
    assert vc.samples == []


# -------------------------------------------------------------------
# Per-variant genotype access
# -------------------------------------------------------------------


def test_genotype_for_heterozygous_sample(multi_sample_vcf):
    vc = load_vcf(multi_sample_vcf, genome="GRCh38")
    # chr17:43082575 C>T — tumor 0/1, normal 0/0.
    variant = next(
        v for v in vc if v.start == 43082575 and v.alt == "T"
    )
    tumor_gt = vc.genotype(variant, "tumor")
    assert tumor_gt.alleles == (0, 1)
    assert tumor_gt.is_called
    assert tumor_gt.zygosity_for_alt(1) is Zygosity.HETEROZYGOUS
    assert tumor_gt.allele_depths == (10, 5)
    assert tumor_gt.genotype_quality == 99

    normal_gt = vc.genotype(variant, "normal")
    assert normal_gt.alleles == (0, 0)
    assert normal_gt.zygosity_for_alt(1) is Zygosity.ABSENT


def test_genotype_for_homozygous_sample(multi_sample_vcf):
    vc = load_vcf(multi_sample_vcf, genome="GRCh38")
    variant = next(v for v in vc if v.start == 117531114)
    tumor_gt = vc.genotype(variant, "tumor")
    assert tumor_gt.alleles == (1, 1)
    assert tumor_gt.zygosity_for_alt(1) is Zygosity.HOMOZYGOUS


def test_genotype_phased_call_preserves_phase_info(multi_sample_vcf):
    vc = load_vcf(multi_sample_vcf, genome="GRCh38")
    variant = next(v for v in vc if v.start == 43082576)
    tumor_gt = vc.genotype(variant, "tumor")
    assert tumor_gt.phased is True
    assert tumor_gt.phase_set == 100


def test_genotype_nocall(multi_sample_vcf):
    vc = load_vcf(multi_sample_vcf, genome="GRCh38")
    variant = next(v for v in vc if v.start == 43082576)
    normal_gt = vc.genotype(variant, "normal")
    assert normal_gt.is_missing
    assert normal_gt.zygosity_for_alt(1) is Zygosity.MISSING


def test_genotype_unknown_sample_raises(multi_sample_vcf):
    vc = load_vcf(multi_sample_vcf, genome="GRCh38")
    variant = vc[0]
    with pytest.raises(SampleNotFoundError):
        vc.genotype(variant, "nonexistent_sample")


def test_sample_not_found_is_key_error_subclass():
    # Callers who only catch KeyError should still work.
    assert issubclass(SampleNotFoundError, KeyError)


def test_genotype_returns_none_for_variant_without_sample_info():
    # A variant that wasn't loaded from a multi-sample VCF.
    from varcode import Variant, VariantCollection
    v = Variant("17", 43082575, "C", "T", "GRCh38")
    vc = VariantCollection(variants=[v])
    assert vc.genotype(v, "anyone") is None


# -------------------------------------------------------------------
# Multi-allelic sites
# -------------------------------------------------------------------


def test_multiallelic_row_splits_into_separate_variants_with_own_genotype(
        multi_sample_vcf):
    vc = load_vcf(multi_sample_vcf, genome="GRCh38")
    # chr1:100 A>T,G — tumor has GT=1/2 (one T, one G).
    at = next(v for v in vc if v.start == 100 and v.alt == "T")
    ag = next(v for v in vc if v.start == 100 and v.alt == "G")

    # Relative to each alt, the tumor is heterozygous (carries one
    # copy of this alt and one copy of a different alt).
    assert vc.zygosity(at, "tumor") is Zygosity.HETEROZYGOUS
    assert vc.zygosity(ag, "tumor") is Zygosity.HETEROZYGOUS

    # Normal has GT=0/1 (one ref, one T). Relative to T: het.
    # Relative to G: absent (normal doesn't have the G alt).
    assert vc.zygosity(at, "normal") is Zygosity.HETEROZYGOUS
    assert vc.zygosity(ag, "normal") is Zygosity.ABSENT


# -------------------------------------------------------------------
# Convenience filters
# -------------------------------------------------------------------


def test_for_sample_filters_to_variants_carried_by_that_sample(
        multi_sample_vcf):
    vc = load_vcf(multi_sample_vcf, genome="GRCh38")

    tumor_variants = vc.for_sample("tumor")
    # Tumor carries all 4 variant-derived alts except normal-only ones.
    # Check concretely:
    # - 17:43082575 C>T: tumor 0/1 -> het, carried
    # - 7:117531114 G>T: tumor 1/1 -> hom, carried
    # - 1:100 A>T: tumor 1/2 -> carries T (alt #1)
    # - 1:100 A>G: tumor 1/2 -> carries G (alt #2)
    # - 17:43082576 C>A: tumor 0|1 -> het, carried
    assert len(tumor_variants) == 5

    normal_variants = vc.for_sample("normal")
    # - 17:43082575 C>T: normal 0/0 -> absent
    # - 7:117531114 G>T: normal 0/1 -> het, carried
    # - 1:100 A>T: normal 0/1 -> carries T
    # - 1:100 A>G: normal 0/1 -> absent (normal doesn't have G)
    # - 17:43082576 C>A: normal ./. -> missing, not carried
    assert len(normal_variants) == 2


def test_heterozygous_in_excludes_homozygous_calls(multi_sample_vcf):
    vc = load_vcf(multi_sample_vcf, genome="GRCh38")
    het_in_tumor = vc.heterozygous_in("tumor")
    starts = sorted(v.start for v in het_in_tumor)
    # Tumor is het at:
    # - 17:43082575 (0/1)
    # - 1:100 (1/2 — het for both T and G)
    # - 17:43082576 (0|1)
    # NOT at 7:117531114 (1/1 is hom, not het)
    assert starts == [100, 100, 43082575, 43082576]


def test_homozygous_alt_in(multi_sample_vcf):
    vc = load_vcf(multi_sample_vcf, genome="GRCh38")
    hom_in_tumor = vc.homozygous_alt_in("tumor")
    assert len(hom_in_tumor) == 1
    assert hom_in_tumor[0].start == 117531114


def test_for_sample_with_unknown_sample_raises(multi_sample_vcf):
    vc = load_vcf(multi_sample_vcf, genome="GRCh38")
    # Fail fast on typos rather than silently returning empty.
    with pytest.raises(SampleNotFoundError):
        vc.for_sample("nonexistent")


def test_filter_chain_composes(multi_sample_vcf):
    # Cross-sample queries fall out of set operations on the primitives.
    vc = load_vcf(multi_sample_vcf, genome="GRCh38")
    # "In tumor but not in normal" — somatic candidates.
    tumor_set = set(vc.for_sample("tumor"))
    normal_set = set(vc.for_sample("normal"))
    somatic = tumor_set - normal_set
    # Tumor carries: 17:43082575, 7:117531114, 1:100(T), 1:100(G), 17:43082576.
    # Normal carries: 7:117531114, 1:100(T).
    # Somatic = 17:43082575, 1:100(G), 17:43082576.
    assert len(somatic) == 3
    starts = sorted(v.start for v in somatic)
    assert starts == [100, 43082575, 43082576]


# -------------------------------------------------------------------
# Package-level exports
# -------------------------------------------------------------------


def test_package_level_exports():
    import varcode
    assert varcode.Genotype is Genotype
    assert varcode.Zygosity is Zygosity
    assert varcode.SampleNotFoundError is SampleNotFoundError
