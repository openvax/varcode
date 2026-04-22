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

"""Tests for :class:`VCFPhaseResolver` — DNA phase via VCF
``GT`` + ``PS`` FORMAT fields (#269).

Uses a hand-rolled VCF with a mix of phased, unphased, multi-allelic,
and cross-block cases. No Isovar or other external deps.
"""

import os
import tempfile

import pytest

from varcode import Variant, VCFPhaseResolver, load_vcf


# Two samples, two phase sets, mix of phased/unphased/homozygous/
# multi-allelic. Positions picked to land on CFTR (chr7) and BRCA1
# (chr17) so the variants are real enough to load.
VCF_BODY = """##fileformat=VCFv4.1
##reference=GRCh38
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor
7\t117531100\t.\tT\tA\t100\tPASS\t.\tGT:PS\t0|1:100
7\t117531114\t.\tG\tT\t100\tPASS\t.\tGT:PS\t0|1:100
7\t117531120\t.\tG\tA\t100\tPASS\t.\tGT:PS\t1|0:100
7\t117531130\t.\tC\tT\t100\tPASS\t.\tGT:PS\t0|1:200
7\t117531140\t.\tA\tG\t100\tPASS\t.\tGT\t0/1
7\t117531150\t.\tA\tT\t100\tPASS\t.\tGT\t1/1
17\t43082575\t.\tC\tT,G\t100\tPASS\t.\tGT:PS\t1|2:500
17\t43082580\t.\tA\tG\t100\tPASS\t.\tGT:PS\t0|1:500
"""


@pytest.fixture
def phased_vcf():
    fd, path = tempfile.mkstemp(suffix=".vcf")
    with os.fdopen(fd, "w") as f:
        f.write(VCF_BODY)
    yield path
    os.unlink(path)


# Convenience lookup helpers — keep tests terse.


def _lookup(vc, contig, start, ref, alt):
    """Find a specific Variant in the loaded collection."""
    for v in vc:
        if (str(v.contig) == contig and v.start == start
                and v.ref == ref and v.alt == alt):
            return v
    raise AssertionError(
        "Variant %s:%d %s>%s not in collection" % (contig, start, ref, alt))


# ------------------------------------------------------------------
# Phased-cis and phased-trans within the same phase set
# ------------------------------------------------------------------


def test_two_phased_variants_same_slot_are_cis(phased_vcf):
    vc = load_vcf(phased_vcf, genome="GRCh38")
    resolver = VCFPhaseResolver(vc, sample="tumor")
    v1 = _lookup(vc, "7", 117531100, "T", "A")
    v2 = _lookup(vc, "7", 117531114, "G", "T")
    # Both are GT=0|1 PS=100 → alt on haplotype slot 1 for each. Cis.
    assert resolver.in_cis(v1, v2) is True


def test_two_phased_variants_different_slots_are_trans(phased_vcf):
    vc = load_vcf(phased_vcf, genome="GRCh38")
    resolver = VCFPhaseResolver(vc, sample="tumor")
    v1 = _lookup(vc, "7", 117531100, "T", "A")    # GT=0|1 → alt at slot 1
    v3 = _lookup(vc, "7", 117531120, "G", "A")    # GT=1|0 → alt at slot 0
    # Same phase set (PS=100) but opposite haplotype slots. Trans.
    assert resolver.in_cis(v1, v3) is False


# ------------------------------------------------------------------
# Cross-block and unphased cases return None
# ------------------------------------------------------------------


def test_different_phase_sets_returns_none(phased_vcf):
    vc = load_vcf(phased_vcf, genome="GRCh38")
    resolver = VCFPhaseResolver(vc, sample="tumor")
    v_block_a = _lookup(vc, "7", 117531100, "T", "A")  # PS=100
    v_block_b = _lookup(vc, "7", 117531130, "C", "T")  # PS=200
    assert resolver.in_cis(v_block_a, v_block_b) is None


def test_unphased_variant_returns_none(phased_vcf):
    vc = load_vcf(phased_vcf, genome="GRCh38")
    resolver = VCFPhaseResolver(vc, sample="tumor")
    v_phased = _lookup(vc, "7", 117531100, "T", "A")   # GT=0|1 PS=100
    v_unphased = _lookup(vc, "7", 117531140, "A", "G")  # GT=0/1 no PS
    assert resolver.in_cis(v_phased, v_unphased) is None


# ------------------------------------------------------------------
# Homozygous-alt short-circuit — phase is deterministic
# ------------------------------------------------------------------


def test_homozygous_alt_is_cis_with_any_called_variant(phased_vcf):
    vc = load_vcf(phased_vcf, genome="GRCh38")
    resolver = VCFPhaseResolver(vc, sample="tumor")
    v_hom = _lookup(vc, "7", 117531150, "A", "T")      # GT=1/1
    v_het = _lookup(vc, "7", 117531100, "T", "A")      # GT=0|1
    # Every haplotype carries the hom-alt, so it's cis with any other
    # variant the sample carries on any haplotype.
    assert resolver.in_cis(v_hom, v_het) is True
    # Commutative.
    assert resolver.in_cis(v_het, v_hom) is True


def test_two_homozygous_alt_are_cis(phased_vcf):
    vc = load_vcf(phased_vcf, genome="GRCh38")
    resolver = VCFPhaseResolver(vc, sample="tumor")
    v_hom = _lookup(vc, "7", 117531150, "A", "T")      # GT=1/1
    # Use the same variant as a stand-in for "both hom".
    assert resolver.in_cis(v_hom, v_hom) is True


# ------------------------------------------------------------------
# Multi-allelic: split Variants map to the right haplotype slot
# ------------------------------------------------------------------


def test_multi_allelic_variants_map_to_correct_slots(phased_vcf):
    """GT=1|2 on a multi-allelic row: first ALT sits on slot 0,
    second ALT on slot 1. Each split Variant should resolve its own
    slot via alt_allele_index."""
    vc = load_vcf(phased_vcf, genome="GRCh38")
    resolver = VCFPhaseResolver(vc, sample="tumor")
    v_alt1 = _lookup(vc, "17", 43082575, "C", "T")     # first ALT, slot 0
    v_alt2 = _lookup(vc, "17", 43082575, "C", "G")     # second ALT, slot 1
    v_nearby = _lookup(vc, "17", 43082580, "A", "G")   # GT=0|1 PS=500, slot 1

    # v_alt2 is on slot 1, same as v_nearby → cis
    assert resolver.in_cis(v_alt2, v_nearby) is True
    # v_alt1 is on slot 0, v_nearby is on slot 1 → trans
    assert resolver.in_cis(v_alt1, v_nearby) is False


# ------------------------------------------------------------------
# phased_partners walks the collection
# ------------------------------------------------------------------


def test_phased_partners_returns_cis_variants_only(phased_vcf):
    vc = load_vcf(phased_vcf, genome="GRCh38")
    resolver = VCFPhaseResolver(vc, sample="tumor")
    v1 = _lookup(vc, "7", 117531100, "T", "A")
    partners = resolver.phased_partners(v1)
    # v2 at 117531114 (same slot, same PS) should be a partner.
    v2 = _lookup(vc, "7", 117531114, "G", "T")
    assert v2 in partners
    # v3 at 117531120 (trans) should NOT.
    v3 = _lookup(vc, "7", 117531120, "G", "A")
    assert v3 not in partners


# ------------------------------------------------------------------
# Provenance
# ------------------------------------------------------------------


def test_vcf_resolver_source_tag():
    from varcode import VariantCollection
    resolver = VCFPhaseResolver(VariantCollection([]), sample="x")
    assert resolver.source == "vcf_ps"


# ------------------------------------------------------------------
# Interface symmetry with IsovarPhaseResolver — same method set
# ------------------------------------------------------------------


def test_vcf_resolver_has_same_phase_api_as_isovar_resolver():
    from varcode import VariantCollection
    from varcode.phasing import IsovarPhaseResolver
    required = {"in_cis", "phased_partners", "source"}
    vcf_r = VCFPhaseResolver(VariantCollection([]), sample="x")
    # Can't instantiate IsovarPhaseResolver without a provider; check
    # class attrs instead.
    for attr in required:
        assert hasattr(vcf_r, attr)
        assert hasattr(IsovarPhaseResolver, attr)
