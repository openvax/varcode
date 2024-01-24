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

import pytest 

from pyensembl import cached_release
from varcode import load_vcf, Variant

from .common import eq_
from .data import data_path



# Set to 1 to enable, 0 to disable.
# TODO: consider running in an in-process HTTP server instead for these tests.
RUN_TESTS_REQUIRING_INTERNET = bool(int(
    os.environ.get("RUN_TESTS_REQUIRING_INTERNET", 0)))

HG19_VCF_FILENAME = data_path("somatic_hg19_14muts.vcf")
HG19_VCF_EXTERNAL_URL = (
    "https://raw.githubusercontent.com/hammerlab/varcode/master/test/data/somatic_hg19_14muts.vcf")

# To load from the branch that introduced these changs:
# (needed before this gets merged to master, can be removed after)
# VCF_EXTERNAL_URL = (
#   "https://raw.githubusercontent.com/hammerlab/varcode/faster-vcf-parsing/test/data/somatic_hg19_14muts.vcf")

def test_load_vcf_local():
    variants = load_vcf(HG19_VCF_FILENAME)
    assert variants.reference_names() == {"GRCh37"}
    assert len(variants) == 14

    variants = load_vcf(HG19_VCF_FILENAME + ".gz")
    assert variants.reference_names() == {"GRCh37"}
    assert len(variants) == 14

    variants = load_vcf("file://%s" % HG19_VCF_FILENAME)
    assert variants.reference_names() == {"GRCh37"}
    assert len(variants) == 14

    variants = load_vcf("file://%s.gz" % HG19_VCF_FILENAME)
    assert variants.reference_names() == {"GRCh37"}
    assert len(variants) == 14

    # An extra slashe before an absolute path can confuse URL parsing.
    # Test that it can still be opened:
    variants = load_vcf("/%s" % HG19_VCF_FILENAME)
    assert variants.reference_names() == {"GRCh37"}
    assert len(variants) == 14

if RUN_TESTS_REQUIRING_INTERNET:
    def test_load_vcf_external():
        variants = load_vcf(HG19_VCF_FILENAME)
        eq_(variants.reference_names(), {"GRCh37"})
        eq_(variants.original_reference_names(), {"hg19"})
        eq_(len(variants), 14)

        variants = load_vcf(HG19_VCF_FILENAME + ".gz")
        eq_(variants.reference_names(), {"GRCh37"})
        eq_(len(variants), 14)

def test_vcf_reference_name():
    variants = load_vcf(HG19_VCF_FILENAME)

    # after normalization, hg19 should be remapped to GRCh37
    assert variants.reference_names() == {"GRCh37"}

def test_genome_arg_to_load_vcf_hg19():
    eq_(load_vcf(HG19_VCF_FILENAME),
        load_vcf(HG19_VCF_FILENAME, genome="hg19"))

def test_genome_arg_to_load_vcf_int_75():
    # if we use Ensembl 75 -- which is backed by GRCh37 -- then the two variant
    # collections will be the same as long as we also convert the contig names
    eq_(load_vcf(HG19_VCF_FILENAME),
        load_vcf(HG19_VCF_FILENAME, genome=75, convert_ucsc_contig_names=True))

    assert load_vcf(HG19_VCF_FILENAME) != load_vcf(
        HG19_VCF_FILENAME,
        genome=75,
        convert_ucsc_contig_names=False)

def test_genome_arg_to_load_vcf_cached_75():
    eq_(load_vcf(HG19_VCF_FILENAME),
        load_vcf(HG19_VCF_FILENAME,
                 genome=cached_release(75), convert_ucsc_contig_names=True))
    assert load_vcf(HG19_VCF_FILENAME) != load_vcf(
        HG19_VCF_FILENAME,
        genome=cached_release(75),
        convert_ucsc_contig_names=False)

def test_genome_arg_to_load_vcf_grch37():
    eq_(load_vcf(HG19_VCF_FILENAME),
        load_vcf(
            HG19_VCF_FILENAME,
            genome="grch37",
            convert_ucsc_contig_names=True))
    eq_(load_vcf(HG19_VCF_FILENAME), load_vcf(
        HG19_VCF_FILENAME,
        genome="GRCh37",
        convert_ucsc_contig_names=True))

    assert load_vcf(HG19_VCF_FILENAME) != load_vcf(
        HG19_VCF_FILENAME,
        genome="grch37",
        convert_ucsc_contig_names=False)

def test_genome_arg_to_load_vcf_b37():
    eq_(load_vcf(HG19_VCF_FILENAME),
        load_vcf(HG19_VCF_FILENAME, genome="b37", convert_ucsc_contig_names=True))

def test_vcf_number_entries():
    # there are 14 mutations listed in the VCF, make sure they are all parsed
    variants = load_vcf(HG19_VCF_FILENAME)
    assert len(variants) == 14, \
        "Expected 14 mutations, got %d" % (len(variants),)

def test_vcf_number_entries_duplicates():
    # There are 3 duplicated mutations listed in the VCF
    path_to_vcf_with_duplicates = data_path("duplicates.vcf")
    variants = load_vcf(
        path_to_vcf_with_duplicates,
        genome='hg38',
        distinct=True)
    assert len(variants) == 1
    variants = load_vcf(
        path_to_vcf_with_duplicates,
        genome='hg38',
        distinct=False)
    assert len(variants) == 3

def generate_vcf_gene_names():
    variants = load_vcf(HG19_VCF_FILENAME)
    for variant in variants:
        yield (variants, variant)

@pytest.mark.parametrize(['collection', 'variant'], generate_vcf_gene_names())
def test_vcf_gene_names(collection, variant):
    expected_gene_names = collection.metadata[variant]['info']['GE']
    assert variant.gene_names == expected_gene_names, \
        "Expected gene name %s for variant %s, got %s" % (
            expected_gene_names, variant, variant.gene_names)


def test_multiple_alleles_per_line():
    variants = load_vcf(data_path("multiallelic.vcf"))
    assert len(variants) == 2, "Expected 2 variants but got %s" % variants
    variant_list = list(variants)
    expected_variants = [
        Variant(1, 1431105, "A", "C", genome="GRCh37"),
        Variant(1, 1431105, "A", "G", genome="GRCh37"),
    ]
    eq_(set(variant_list), set(expected_variants))

def test_sample_info_genotype():
    variants = load_vcf(data_path("multiallelic.vcf"))
    assert len(variants) == 2, "Expected 2 variants but got %s" % variants
    eq_(variants.metadata[variants[0]]['sample_info']['metastasis']['GT'],
        '0/1')
    eq_(variants.metadata[variants[1]]['sample_info']['metastasis']['GT'],
        '0/1')
