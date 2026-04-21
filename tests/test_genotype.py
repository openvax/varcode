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
Unit tests for the Genotype dataclass and GT parsing (openvax/varcode#267).

Tests that touch a real VCF end-to-end live in
``tests/test_genotype_from_vcf.py``.
"""

import pytest

from varcode import Genotype, Zygosity
from varcode.genotype import parse_gt_string


# ------------------------------------------------------------------
# parse_gt_string: one low-level function, many cases
# ------------------------------------------------------------------


def test_parse_gt_diploid_unphased():
    assert parse_gt_string("0/1") == ((0, 1), False)
    assert parse_gt_string("1/0") == ((1, 0), False)
    assert parse_gt_string("1/1") == ((1, 1), False)
    assert parse_gt_string("0/0") == ((0, 0), False)


def test_parse_gt_diploid_phased():
    assert parse_gt_string("0|1") == ((0, 1), True)
    assert parse_gt_string("1|0") == ((1, 0), True)


def test_parse_gt_multiallelic():
    # Multi-allelic: 1/2 means one copy of alt #1, one copy of alt #2.
    assert parse_gt_string("1/2") == ((1, 2), False)
    assert parse_gt_string("2|1") == ((2, 1), True)


def test_parse_gt_haploid():
    # Chromosomes X/Y in males, mitochondrial, etc. — single allele.
    assert parse_gt_string("1") == ((1,), False)
    assert parse_gt_string("0") == ((0,), False)


def test_parse_gt_missing():
    assert parse_gt_string("./.") == ((None, None), False)
    assert parse_gt_string(".") == ((None,), False)
    assert parse_gt_string("") == ((None,), False)
    assert parse_gt_string(None) == ((None,), False)


def test_parse_gt_partial_missing():
    # Half-called genotypes — one haplotype missing.
    assert parse_gt_string("./1") == ((None, 1), False)
    assert parse_gt_string("0/.") == ((0, None), False)
    assert parse_gt_string(".|1") == ((None, 1), True)


def test_parse_gt_polyploid():
    # Triploid (e.g. trisomy) or higher — tuple simply grows.
    alleles, phased = parse_gt_string("0/1/1")
    assert alleles == (0, 1, 1)
    assert phased is False


# ------------------------------------------------------------------
# Genotype construction from pyvcf-style sample info dicts
# ------------------------------------------------------------------


def test_genotype_from_sample_info_full():
    gt = Genotype.from_sample_info({
        "GT": "0/1",
        "AD": [10, 5],
        "DP": 15,
        "GQ": 99,
    })
    assert gt.raw_gt == "0/1"
    assert gt.alleles == (0, 1)
    assert gt.phased is False
    assert gt.allele_depths == (10, 5)
    assert gt.total_depth == 15
    assert gt.genotype_quality == 99
    assert gt.phase_set is None


def test_genotype_from_sample_info_phased_with_ps():
    gt = Genotype.from_sample_info({
        "GT": "0|1",
        "PS": 100,
        "AD": [8, 7],
        "DP": 15,
        "GQ": 99,
    })
    assert gt.phased is True
    assert gt.phase_set == 100


def test_genotype_from_sample_info_nocall_handles_none_values():
    # Pyvcf returns None for all fields when the call is ./.
    gt = Genotype.from_sample_info({
        "GT": "./.",
        "AD": None,
        "DP": None,
        "GQ": None,
    })
    assert gt.alleles == (None, None)
    assert gt.allele_depths is None
    assert gt.is_missing
    assert not gt.is_called


def test_genotype_from_sample_info_none_input():
    # Entirely missing sample_info dict (e.g. variant not constructed
    # from a VCF).
    gt = Genotype.from_sample_info(None)
    assert gt.is_missing
    assert gt.alleles == (None, None)


# ------------------------------------------------------------------
# General predicates (alt-agnostic)
# ------------------------------------------------------------------


def test_is_called_and_is_missing():
    assert Genotype.from_sample_info({"GT": "0/1"}).is_called
    assert not Genotype.from_sample_info({"GT": "./."}).is_called
    assert Genotype.from_sample_info({"GT": "./."}).is_missing
    # Partial-missing is still "called" because one allele is known.
    assert Genotype.from_sample_info({"GT": "./1"}).is_called


def test_ploidy():
    assert Genotype.from_sample_info({"GT": "0/1"}).ploidy == 2
    assert Genotype.from_sample_info({"GT": "1"}).ploidy == 1
    assert Genotype.from_sample_info({"GT": "0/1/1"}).ploidy == 3


def test_is_haploid():
    assert Genotype.from_sample_info({"GT": "1"}).is_haploid
    assert not Genotype.from_sample_info({"GT": "0/1"}).is_haploid


# ------------------------------------------------------------------
# Alt-relative zygosity (the business end of the API)
# ------------------------------------------------------------------


def test_zygosity_heterozygous_simple():
    gt = Genotype.from_sample_info({"GT": "0/1"})
    assert gt.zygosity_for_alt(1) is Zygosity.HETEROZYGOUS
    assert gt.carries_alt(1)
    assert gt.copies_of_alt(1) == 1


def test_zygosity_homozygous_alt_simple():
    gt = Genotype.from_sample_info({"GT": "1/1"})
    assert gt.zygosity_for_alt(1) is Zygosity.HOMOZYGOUS
    assert gt.carries_alt(1)
    assert gt.copies_of_alt(1) == 2


def test_zygosity_homozygous_ref_is_absent_for_any_alt():
    gt = Genotype.from_sample_info({"GT": "0/0"})
    # Relative to alt 1: sample doesn't have this alt.
    assert gt.zygosity_for_alt(1) is Zygosity.ABSENT
    assert not gt.carries_alt(1)
    assert gt.copies_of_alt(1) == 0


def test_zygosity_missing():
    gt = Genotype.from_sample_info({"GT": "./."})
    assert gt.zygosity_for_alt(1) is Zygosity.MISSING
    assert not gt.carries_alt(1)
    assert gt.copies_of_alt(1) == 0


def test_zygosity_multiallelic_querying_different_alts():
    # GT = 1/2 means one copy of alt #1, one copy of alt #2.
    # Querying alt 1: het (one copy of this alt, one of a different alt).
    # Querying alt 2: also het.
    # Querying alt 3 (not carried): absent.
    gt = Genotype.from_sample_info({"GT": "1/2"})
    assert gt.zygosity_for_alt(1) is Zygosity.HETEROZYGOUS
    assert gt.zygosity_for_alt(2) is Zygosity.HETEROZYGOUS
    assert gt.zygosity_for_alt(3) is Zygosity.ABSENT


def test_zygosity_multiallelic_homozygous_for_second_alt():
    # GT = 2/2 means both copies are alt #2.
    # Alt 2: hom. Alt 1: absent (sample doesn't have alt 1).
    gt = Genotype.from_sample_info({"GT": "2/2"})
    assert gt.zygosity_for_alt(2) is Zygosity.HOMOZYGOUS
    assert gt.zygosity_for_alt(1) is Zygosity.ABSENT


def test_zygosity_haploid_single_alt():
    # chrY in a male with an alt call: GT = 1.
    gt = Genotype.from_sample_info({"GT": "1"})
    # Single-allele calls with that allele equal to alt: all copies
    # are alt → classify as HOMOZYGOUS.
    assert gt.zygosity_for_alt(1) is Zygosity.HOMOZYGOUS
    assert gt.carries_alt(1)


def test_zygosity_partial_call():
    # GT = ./1: one allele missing, the other is alt #1.
    # Out of the called alleles (just one), all are alt #1 → HOMOZYGOUS.
    # This is the defensible read: we count called alleles only.
    gt = Genotype.from_sample_info({"GT": "./1"})
    assert gt.zygosity_for_alt(1) is Zygosity.HOMOZYGOUS


# ------------------------------------------------------------------
# Per-allele depth lookup
# ------------------------------------------------------------------


def test_depth_for_alt():
    gt = Genotype.from_sample_info({"GT": "0/1", "AD": [10, 5]})
    assert gt.depth_for_alt(0) == 10   # ref depth
    assert gt.depth_for_alt(1) == 5    # first alt depth
    assert gt.depth_for_alt(2) is None  # out of range


def test_depth_for_alt_returns_none_when_ad_missing():
    gt = Genotype.from_sample_info({"GT": "0/1"})
    assert gt.depth_for_alt(1) is None


def test_depth_for_alt_multiallelic():
    # AD has one entry per allele (ref + each alt).
    gt = Genotype.from_sample_info({"GT": "1/2", "AD": [3, 7, 5]})
    assert gt.depth_for_alt(0) == 3
    assert gt.depth_for_alt(1) == 7
    assert gt.depth_for_alt(2) == 5


# ------------------------------------------------------------------
# Dataclass-y ergonomics: hashable, equatable, frozen
# ------------------------------------------------------------------


def test_genotype_is_frozen():
    gt = Genotype.from_sample_info({"GT": "0/1"})
    with pytest.raises((AttributeError, Exception)):
        gt.alleles = (1, 1)  # type: ignore


def test_genotype_equality():
    a = Genotype.from_sample_info({"GT": "0/1", "AD": [10, 5], "DP": 15})
    b = Genotype.from_sample_info({"GT": "0/1", "AD": [10, 5], "DP": 15})
    assert a == b
    c = Genotype.from_sample_info({"GT": "1/1", "AD": [0, 15], "DP": 15})
    assert a != c


# ------------------------------------------------------------------
# DataclassSerializable round-trip (#272).
# ------------------------------------------------------------------


def test_genotype_json_round_trip():
    """``Genotype`` now inherits :class:`DataclassSerializable`, so
    ``to_json`` / ``from_json`` round-trip without any hand-rolled
    serialization code."""
    gt = Genotype.from_sample_info({
        "GT": "0/1",
        "AD": [10, 5],
        "DP": 15,
        "GQ": 99,
        "PS": 12345,
    })
    rt = Genotype.from_json(gt.to_json())
    assert rt == gt
    assert rt.alleles == (0, 1)
    assert rt.allele_depths == (10, 5)
    assert rt.total_depth == 15
    assert rt.genotype_quality == 99
    assert rt.phase_set == 12345
