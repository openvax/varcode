# Copyright (c) 2015. Mount Sinai School of Medicine
#
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

import os
from nose.tools import eq_
from varcode import load_vcf, Variant
from . import data_path

# Set to 1 to enable, 0 to disable.
# TODO: consider running in an in-process HTTP server instead for these tests.
RUN_TESTS_REQUIRING_INTERNET = bool(int(
    os.environ.get("RUN_TESTS_REQUIRING_INTERNET", 0)))

VCF_FILENAME = data_path("somatic_hg19_14muts.vcf")
VCF_EXTERNAL_URL = (
    "https://raw.githubusercontent.com/hammerlab/varcode/master/test/data/somatic_hg19_14muts.vcf")

def test_load_vcf_local():
    variants = load_vcf(VCF_FILENAME)
    assert variants.reference_names() == {"GRCh37"}
    assert len(variants) == 14
    
    variants = load_vcf(VCF_FILENAME + ".gz")
    assert variants.reference_names() == {"GRCh37"}
    assert len(variants) == 14
    
    variants = load_vcf("file://%s" % VCF_FILENAME)
    assert variants.reference_names() == {"GRCh37"}
    assert len(variants) == 14
    
    variants = load_vcf("file://%s.gz" % VCF_FILENAME)
    assert variants.reference_names() == {"GRCh37"}
    assert len(variants) == 14

if RUN_TESTS_REQUIRING_INTERNET:
    def test_load_vcf_external():
        variants = load_vcf(VCF_EXTERNAL_URL)
        assert variants.reference_names() == {"GRCh37"}
        assert len(variants) == 14
        
        variants = load_vcf(VCF_EXTERNAL_URL + ".gz")
        assert variants.reference_names() == {"GRCh37"}
        assert len(variants) == 14

def test_vcf_reference_name():
    variants = load_vcf(VCF_FILENAME)
    # after normalization, hg19 should be remapped to GRCh37
    assert variants.reference_names() == {"GRCh37"}


def test_reference_arg_to_load_vcf():
    variants = load_vcf(VCF_FILENAME)
    eq_(variants, load_vcf(VCF_FILENAME, ensembl_version=75))
    eq_(variants, load_vcf(VCF_FILENAME, reference_name="grch37"))
    eq_(variants, load_vcf(VCF_FILENAME, reference_name="GRCh37"))
    eq_(variants, load_vcf(VCF_FILENAME, reference_name="b37"))
    # TODO: actually make hg19 different from b37! They should use
    # different MT sequences
    eq_(variants, load_vcf(VCF_FILENAME, reference_name="hg19"))

def test_vcf_number_entries():
    # there are 14 mutations listed in the VCF, make sure they are all parsed
    variants = load_vcf(VCF_FILENAME)
    assert len(variants) == 14, \
        "Expected 14 mutations, got %d" % (len(variants),)

def _check_variant_gene_name(variant):
    expected_gene_names = variant.info['GE']
    assert variant.gene_names == expected_gene_names, \
        "Expected gene name %s for variant %s, got %s" % (
            expected_gene_names, variant, variant.gene_names)

def test_vcf_gene_names():
    variants = load_vcf(VCF_FILENAME)
    for variant in variants:
        yield (_check_variant_gene_name, variant)

def test_multiple_alleles_per_line():
    variants = load_vcf(data_path("multiallelic.vcf"))
    assert len(variants) == 2, "Expected 2 variants but got %s" % variants
    variant_list = list(variants)
    ensembl = variant_list[0].ensembl
    expected_variants = [
        Variant(1, 1431105, "A", "C", ensembl=ensembl),
        Variant(1, 1431105, "A", "G", ensembl=ensembl),
    ]
    eq_(set(variant_list), set(expected_variants))
