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

from __future__ import print_function, division, absolute_import
from varcode import load_vcf, load_maf

from varcode.vcf_output import variants_to_vcf
from .data import data_path
import tempfile

TEST_FILENAMES_HUMAN = [
    'duplicates.maf',
    'multiallelic.vcf',
    'mutect-example.vcf',
    'ov.wustle.subset5.maf',
    'somatic_hg19_14muts.space_in_sample_name.vcf',
    'somatic_hg19_14muts.vcf',
    'strelka-example.vcf',
    'tcga_ov.head.maf',
    'tcga_ov.head.xychr.maf',
    # 'dbnsfp_validation_set.csv',      # csv
    # 'duplicates.vcf',                 # no ref genome header
    # 'mutect-example-headerless.vcf',  # no ref genome header
    # 'somatic_hg19_14muts.vcf.gz',     # gzip
]

TEST_FILENAMES_MOUSE = [
    'mouse_vcf_dbsnp_chr1_partial.vcf',
]

TEST_FILENAMES = TEST_FILENAMES_HUMAN + TEST_FILENAMES_MOUSE


def _merge_metadata_naive(variants):
    return {
        k: v
        for d in variants.source_to_metadata_dict.values()
        for k, v in d.items()
    }


def _do_roundtrip_test(filenames):

    def load_fn(filename):
        return {
            'vcf': load_vcf,
            'maf': load_maf
        }[filename.split('.')[-1]]

    def load_variants():
        variant_collections = []
        for filename in filenames:
            variant_collections.append(load_fn(filename)(data_path(filename)))
        return variant_collections[0].union(*variant_collections[1:])

    variants = load_variants()

    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        metadata = _merge_metadata_naive(variants)
        variants_to_vcf(variants, metadata, out=f)
        tmp_name = f.name
    reparsed_variants = load_vcf(tmp_name)

    # `==` checks the reference genome, which won't necessarily match.
    assert all(
        v1.contig == v2.contig and
        v1.start == v2.start and
        v1.ref == v2.ref and
        v1.start == v2.start
        for (v1, v2) in zip(variants, reparsed_variants))

    return (variants, reparsed_variants)

    # TODO: There is definitely more opportunity here to compare metadata
    # fields, with caveats.
    #
    # First, any variants from non-VCF sources (e.g., MAF files) will inevitably
    # lose some information through the change in representation (more importantly,
    # even if there is no loss in data, that data will be in a different format in
    # the new metadata dictionary). Thus, we should either ignore such variants
    # or only check certain fields.
    #
    # Second, without the original metadata headers in the VCF file, all metadata
    # information will be parsed as strings. Thus, for a simple comparison between
    # metadata (without the need to individually convert fields), we'd need to add
    # these headers to the output VCF file. See `vcf_output.py` for more info.


def test_single_file_roundtrip_conversion():
    for filename in TEST_FILENAMES:
        yield (_do_roundtrip_test, [filename])


def test_multiple_file_roundtrip_conversion():
    file_groups = (
        ['simple.1.vcf', 'simple.2.vcf'],  # basic multi-file test
        ['duplicates.maf', 'multiallelic.vcf'],  # dif. file formats
        ['duplicate-id.1.vcf', 'duplicate-id.2.vcf'],
        TEST_FILENAMES_HUMAN,
    )
    for file_group in file_groups:
        yield (_do_roundtrip_test, file_group)


def test_same_samples_produce_samples():
    """Ensures that, if a set of variants have the same samples, the reparsed
    collection will output these samples.
    """
    (variants, reparsed_variants) = _do_roundtrip_test(
        ['same-samples.1.vcf', 'same-samples.2.vcf'])

    original_metadata = _merge_metadata_naive(variants)
    reparsed_metadata = _merge_metadata_naive(reparsed_variants)

    sample_names = set(list(original_metadata.values())[0]['sample_info'].keys())
    assert all(
        set(d.get('sample_info', {}).keys()) == sample_names
        for d in reparsed_metadata.values())


def test_different_samples_produce_no_samples():
    """Ensures that, if a set of variants have different samples, the reparsed
    collection will not output any samples.

    See `vcf_output.py` for details as to why this is the way it's done for now.
    """
    (_, reparsed_variants) = _do_roundtrip_test(
        ['different-samples.1.vcf', 'different-samples.2.vcf'])

    metadata = _merge_metadata_naive(reparsed_variants)
    assert all(d.get('sample_info') is None for d in metadata.values())
