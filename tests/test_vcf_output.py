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

import tempfile

import pytest 

from varcode import load_vcf, load_maf
from varcode.vcf_output import variants_to_vcf

from .data import data_path


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



def _do_roundtrip_test(filenames, convert_ucsc_to_grch37=False):

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
    if convert_ucsc_to_grch37:
        variants = variants.clone_without_ucsc_data()

    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        metadata = _merge_metadata_naive(variants)
        variants_to_vcf(variants, metadata, out=f)
        tmp_name = f.name
    reparsed_variants = load_vcf(tmp_name)

    # `==` checks the reference genome, which won't necessarily match.
    for (v1, v2) in zip(variants, reparsed_variants):
        assert (
            v1.contig == v2.contig and
            v1.start == v2.start and
            v1.ref == v2.ref and
            v1.start == v2.start), (v1, v2)

    return (variants, reparsed_variants)

    # TODO:
    #   There is definitely more opportunity here to compare metadata
    #   fields, with caveats.
    #   ---
    #   First, any variants from non-VCF sources (e.g., MAF files) will inevitably
    #   lose some information through the change in representation (more importantly,
    #   even if there is no loss in data, that data will be in a different format in
    #   the new metadata dictionary). Thus, we should either ignore such variants
    #   or only check certain fields.
    #   ---
    #   Second, without the original metadata headers in the VCF file, all metadata
    #   information will be parsed as strings. Thus, for a simple comparison between
    #   metadata (without the need to individually convert fields), we'd need to add
    #   these headers to the output VCF file. See `vcf_output.py` for more info.


@pytest.mark.parametrize(['filename'], [(f,) for f in TEST_FILENAMES])
def test_roundtrip_serialization_single_file(filename):
    _do_roundtrip_test([filename])

FILENAME_PAIRS = (
    ['simple.1.vcf', 'simple.2.vcf'],  # basic multi-file VCF test
    ['duplicates.maf', 'ov.wustle.subset5.maf'],  # multiple MAF files
    ['duplicate-id.1.vcf', 'duplicate-id.2.vcf'],
)

@pytest.mark.parametrize(['file_group'], [(f,) for f in FILENAME_PAIRS])
def test_multiple_file_roundtrip_conversion(file_group):
    _do_roundtrip_test(file_group)

def test_multiple_file_roundtrip_conversion_mixed_references():
    # testing roundtrip serialization of hg19 VCF files
    # converted to GRCh37 combined with b37 MAFs
    _do_roundtrip_test(TEST_FILENAMES_HUMAN, convert_ucsc_to_grch37=True)

def test_same_samples_produce_samples():
    """test_same_samples_produce_samples

    Ensures that, if a set of variants have the same samples, the reparsed
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
    """test_different_samples_produce_no_samples

    Ensures that, if a set of variants have different samples, the reparsed
    collection will not output any samples.

    See `vcf_output.py` for details as to why this is the way it's done for now.
    """
    (_, reparsed_variants) = _do_roundtrip_test(
        ['different-samples.1.vcf', 'different-samples.2.vcf'])

    metadata = _merge_metadata_naive(reparsed_variants)
    assert all(d.get('sample_info') is None for d in metadata.values())
