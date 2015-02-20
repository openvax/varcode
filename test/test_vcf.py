# Copyright (c) 2014. Mount Sinai School of Medicine
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

from varcode import load_vcf

VCF_FILENAME = "data/somatic_hg19_14muts.vcf"

def test_vcf_reference_name():
    variants = load_vcf(VCF_FILENAME)
    # after normalization, hg19 should be remapped to GRCh37
    assert variants.reference_names() ==  { "GRCh37" }

def test_vcf_number_entries():
    # there are 14 mutations listed in the VCF, make sure they are all parsed
    variants = load_vcf(VCF_FILENAME)
    assert len(variants) == 14, \
        "Expected 14 mutations, got %d" % (len(variants),)

def _check_annotation_gene_name(annotation):
    variant = annotation.variant
    expected_gene_names = variant.info['GE']
    gene_names = [gene.name for gene in annotation.genes]
    assert expected_gene_names == gene_names, \
        "Expected gene name %s for variant %s, got %s" % (
            expected_gene_name, variant, gene_names)

def test_vcf_gene_names():
    variants = load_vcf(VCF_FILENAME)
    for annotation in variants.annotations():
        yield (_check_annotation_gene_name, annotation)
