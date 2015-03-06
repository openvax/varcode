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

from nose.tools import eq_
from varcode import load_vcf, Variant

VCF_FILENAME = "data/somatic_hg19_14muts.vcf"

def test_vcf_reference_name():
    variants = load_vcf(VCF_FILENAME)
    # after normalization, hg19 should be remapped to GRCh37
    assert variants.reference_names() == {"GRCh37"}

def test_vcf_number_entries():
    # there are 14 mutations listed in the VCF, make sure they are all parsed
    variants = load_vcf(VCF_FILENAME)
    assert len(variants) == 14, \
        "Expected 14 mutations, got %d" % (len(variants),)

def _check_variant_gene_name(variant):
    expected_gene_names = variant.info['GE']
    assert variant.gene_names() == expected_gene_names, \
        "Expected gene name %s for variant %s, got %s" % (
            expected_gene_names, variant, variant.gene_names())

def test_vcf_gene_names():
    variants = load_vcf(VCF_FILENAME)
    for variant in variants:
        yield (_check_variant_gene_name, variant)

def test_multiple_alleles_per_line():
    variants = load_vcf("data/multiallelic.vcf")
    assert len(variants) == 2, "Expected 2 variants but got %s" % variants
    variant_list = list(variants)
    ensembl = variant_list[0].ensembl
    expected_variants = [
        Variant(1, 1431105, "A", "C", ensembl=ensembl),
        Variant(1, 1431105, "A", "G", ensembl=ensembl),
    ]
    eq_(variant_list, expected_variants)