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


import warnings

import pytest 
from pyensembl import cached_release

from varcode import Variant, load_vcf
from varcode.cli.effects_script import arg_parser
from varcode.cli.variant_args import variant_collection_from_args
from varcode.reference import (
    ensembl_reference_aliases, infer_genome, infer_reference_name, most_recent_assembly_name)
from .common import eq_

## test cases are given as 
## expected response: list of inputs
reference_test_cases = {
    'NCBI36': [
        'ncbi36p2.fasta', 
        'b36.fasta', 
        '##reference=file:///var/lib/cwl/ncbi36/homo_sapiens.d1.vd1.fa'],
    'GRCh38': [
        'grch38p2.fasta', 
        '##reference=file:///var/lib/cwl/job367935311_index_001zdr/GRCh38.d1.vd1.fa',
        '##reference=file:///var/lib/cwl/job367935311_index_001zdr/GRCh38.job36.d1.vd1.fa',
    ],
}

# These resolve with a warning: the first file name names no assembly, and
# in the second "job36" contains NCBI36's alias "b36".
AMBIGUOUS_REFERENCE_PATHS = {
    '##reference=file:///var/lib/cwl/ncbi36/homo_sapiens.d1.vd1.fa':
        "could not be matched against filename",
    '##reference=file:///var/lib/cwl/job367935311_index_001zdr/GRCh38.job36.d1.vd1.fa':
        "More than one reference",
}

def test_most_recent_assembly():
    eq_(most_recent_assembly_name(['ncbi36', 'grch38']), 'grch38')
    eq_(most_recent_assembly_name(['ncbi36', 'grch38', '37mm']), 'grch38')
    eq_(most_recent_assembly_name(['ncbi36']), 'ncbi36')
    eq_(most_recent_assembly_name(['ncbi36', '35']), 'ncbi36')
def generate_reference_name_aliases():
    for assembly_name, aliases in ensembl_reference_aliases.items():
        candidate_list = [assembly_name] + list(aliases)
        for candidate in candidate_list:
            yield (
                candidate,
                assembly_name
            )

@pytest.mark.parametrize(['candidate', 'assembly_name'], list(generate_reference_name_aliases()))
def test_infer_reference_name_aliases(candidate, assembly_name):
    eq_(infer_reference_name(candidate), assembly_name)
    
def generate_reference_name_fasta_filenames():
    for assembly_name, aliases in reference_test_cases.items():
        candidate_list = [assembly_name] + list(aliases)
        for candidate in candidate_list:
            yield (
                candidate,
                assembly_name
            )

@pytest.mark.parametrize(['candidate', 'assembly_name'], list(generate_reference_name_fasta_filenames()))
def test_reference_name_fasta_filenames(candidate, assembly_name):
    if candidate in AMBIGUOUS_REFERENCE_PATHS:
        # Bypass memoization so the warning fires whatever ran first.
        with pytest.warns(UserWarning, match=AMBIGUOUS_REFERENCE_PATHS[candidate]):
            eq_(infer_reference_name.__wrapped__(candidate), assembly_name)
    else:
        eq_(infer_reference_name(candidate), assembly_name)


@pytest.mark.parametrize(["name", "reference_name"], [
    ("CanFam3.1", "CanFam3.1"),
    ("/ref/CanFam3.1.dna.toplevel.fa", "CanFam3.1"),
    ("GRCh37_hg19.fa", "GRCh37"),
    ("GRCm39_mm39.fa", "GRCm39"),
])
def test_ensembl_and_ucsc_names_of_one_assembly_are_not_ambiguous(name, reference_name):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert infer_reference_name.__wrapped__(name) == reference_name



@pytest.mark.parametrize(["name", "release", "species", "was_ucsc"], [
    ("GRCh38:93", 93, "homo_sapiens", False),
    (" GRCh38 : 93 ", 93, "homo_sapiens", False),
    ("B37:75", 75, "homo_sapiens", False),
    ("hg19:75", 75, "homo_sapiens", True),
    ("GRCm38:95", 95, "mus_musculus", False),
    ("GRCh38.p13:93", 93, "homo_sapiens", False),
    # Ensembl's file names put the release after a dot; other separators work too.
    ("GRCh38.93", 93, "homo_sapiens", False),
    ("GRCh38_93", 93, "homo_sapiens", False),
    ("GRCh38-93", 93, "homo_sapiens", False),
    ("GRCh38 93", 93, "homo_sapiens", False),
    ("hg38.93", 93, "homo_sapiens", True),
    ("mm10.95", 95, "mus_musculus", True),
    ("mm39:110", 110, "mus_musculus", True),
    ("mm39.110", 110, "mus_musculus", True),
])
def test_reference_name_with_release_chooses_that_release(name, release, species, was_ucsc):
    genome, converted = infer_genome(name)
    assert (genome.release, genome.species.latin_name, converted) == (release, species, was_ucsc)
    assert genome is cached_release(release, species=species)


@pytest.mark.parametrize(["name", "message"], [
    ("GRCh38:75", "Assembly 'GRCh38' conflicts with selected GRCh37"),
    ("GRCh38.75", "Assembly 'GRCh38' conflicts with selected GRCh37"),
    ("GRCh37:93", "Assembly 'GRCh37' conflicts with selected GRCh38"),
    ("GRCh38:40", "No genome for homo_sapiens in Ensembl release 40"),
    # mm39 is GRCm39; a GRCm38 release is an error, not a substitute (#668).
    ("mm39:95", "Assembly 'GRCm39' conflicts with selected GRCm38"),
])
def test_release_that_does_not_provide_the_assembly_is_rejected(name, message):
    # Previously the release was dropped and the latest release used instead (#512).
    with pytest.raises(ValueError, match=message):
        infer_genome(name)


@pytest.mark.parametrize(["name", "ucsc_name", "assembly"], [
    ("mm39", "mm39", "GRCm39"),
    ("MM39", "mm39", "GRCm39"),
    ("/data/mm39.fa", "mm39", "GRCm39"),
    ("##reference=file:///ref/mm39.fa", "mm39", "GRCm39"),
    ("mm10", "mm10", "GRCm38"),
    # UCSC mm9 is NCBI Build 37; "GRCm37" is not an assembly (#690).
    ("mm9", "mm9", "NCBIM37"),
    ("/data/mm9.fa", "mm9", "NCBIM37"),
])
def test_ucsc_mouse_names_select_the_same_assembly(name, ucsc_name, assembly):
    # A UCSC mouse build is its Ensembl assembly, so no coordinates change
    # (#668).
    genome, converted = infer_genome(name)
    assert (genome.reference_name, genome.species.latin_name, converted) == (
        assembly, "mus_musculus", True)
    variant = Variant("chr7", 100, "A", "T", genome=name)
    assert (variant.reference_name, variant.original_reference_name,
            variant.contig) == (assembly, ucsc_name, "7")


@pytest.mark.parametrize(["name", "assembly"], [
    ("GCA_000001635.1", "NCBIM37"),
    ("GCF_000001635.18", "NCBIM37"),
    ("GCA_000001635.2", "GRCm38"),
    ("GCA_000001635.8", "GRCm38"),
    ("GCF_000001635.26", "GRCm38"),
    ("GCA_000001635.9", "GRCm39"),
    ("GCF_000001635.27", "GRCm39"),
    ("/ref/GCA_000001635.9_GRCm39_genomic.fna", "GRCm39"),
    ("/ref/GCF_000001635.26_GRCm38.p6_genomic.fna", "GRCm38"),
])
def test_mouse_accessions_name_their_own_build(name, assembly):
    # Versioned accessions were read as "GRCm38, Ensembl release N" (#690).
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert infer_reference_name.__wrapped__(name) == assembly
    genome, converted = infer_genome(name)
    assert (genome.reference_name, converted) == (assembly, False)


def test_unversioned_mouse_accession_names_no_single_build():
    # GCA_000001635 is the series of every mouse build, not GRCm38 (#690).
    with pytest.raises(ValueError, match="Failed to infer genome assembly name"):
        infer_genome("GCA_000001635")


@pytest.mark.parametrize(["path", "reference_name"], [
    # 1000 Genomes b37 builds, named 1..22, X, Y, MT (#712).
    ("file:///refs/hs37d5.fa", "GRCh37"),
    ("hs37d5.fa.gz", "GRCh37"),
    ("/refs/hs37.fa", "GRCh37"),
    ("/refs/human_g1k_v37.fasta", "GRCh37"),
    ("/refs/human_g1k_v37_decoy.fasta", "GRCh37"),
    # Broad/GATK b37. The Broad's own bucket files it under hg19/, but its
    # contigs are named 1 and MT, not chr1 and chrM.
    ("file:///refs/Homo_sapiens_assembly19.fasta", "GRCh37"),
    ("gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta",
     "GRCh37"),
    # GRCh38 builds named chr1..chr22, chrX, chrY, chrM, as hg38 is.
    ("file:///refs/Homo_sapiens_assembly38.fasta", "hg38"),
    ("/refs/hs38DH.fa", "hg38"),
    ("/refs/hs38a.fa", "hg38"),
    ("/refs/GRCh38/hs38.fa", "hg38"),
])
def test_common_reference_file_names_name_their_build(path, reference_name):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert infer_reference_name.__wrapped__(path) == reference_name


@pytest.mark.parametrize("name,release,assembly,converted", [
    ("hs37d5", 75, "GRCh37", False),
    ("human_g1k_v37", 75, "GRCh37", False),
    ("Homo_sapiens_assembly19", 75, "GRCh37", False),
    ("Homo_sapiens_assembly38", 93, "GRCh38", True),
    ("hs38", 93, "GRCh38", True),
    ("hs38a", 93, "GRCh38", True),
    ("hs38DH", 93, "GRCh38", True),
])
@pytest.mark.parametrize("separator", [":", ".", "_", "-", " "])
def test_reference_file_aliases_preserve_requested_release(
        name, release, assembly, converted, separator):
    genome, was_ucsc = infer_genome("%s%s%s" % (name, separator, release))
    assert (genome.reference_name, genome.release, was_ucsc) == (assembly, release, converted)
    conflicting_release = 93 if assembly == "GRCh37" else 75
    with pytest.raises(ValueError, match="conflicts with selected"):
        infer_genome("%s%s%s" % (name, separator, conflicting_release))


@pytest.mark.parametrize(["reference", "chrom", "assembly", "original"], [
    ("file:///refs/hs37d5.fa", "1", "GRCh37", "GRCh37"),
    ("file:///refs/Homo_sapiens_assembly38.fasta", "chr1", "GRCh38", "hg38"),
])
def test_load_vcf_infers_the_build_of_common_reference_files(
        tmp_path, reference, chrom, assembly, original):
    path = tmp_path / "variants.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##reference=%s\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "%s\t1000000\t.\tA\tG\t.\tPASS\t.\n" % (reference, chrom))
    variant, = load_vcf(str(path))
    assert (variant.contig, variant.reference_name, variant.original_reference_name) == (
        "1", assembly, original)


def test_mm39_variants_from_vcf_api_and_cli_share_one_identity(tmp_path):
    path = tmp_path / "variants.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##reference=file:///ref/mm39.fa\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr7\t100\t.\tA\tT\t.\tPASS\t.\n")
    from_vcf, = load_vcf(str(path))
    from_api = Variant("chr7", 100, "A", "T", genome="mm39")
    from_cli, = variant_collection_from_args(arg_parser.parse_args(
        ["--genome", "mm39", "--variant", "chr7", "100", "A", "T"]))
    from_json = Variant.from_json(from_api.to_json())
    for variant in (from_vcf, from_api, from_cli, from_json):
        assert (variant.genome.reference_name, variant.reference_name,
                variant.original_reference_name, variant.contig) == (
            "GRCm39", "GRCm39", "mm39", "7")
    assert from_vcf == from_api == from_cli == from_json


def test_variant_keeps_the_chosen_release_through_serialization():
    variant = Variant("7", 140753336, "A", "T", genome="GRCh38:93")
    assert variant.genome.release == 93 and variant.reference_name == "GRCh38"
    assert Variant.from_json(variant.to_json()).genome.release == 93


@pytest.mark.parametrize(["name", "reference_name"], [
    ("GRCh38.p13", "GRCh38"),
    ("grch38.d1.vd1", "GRCh38"),
    ("Felis_catus_9.0", "Felis_catus_9.0"),
    ("CanFam3.1", "CanFam3.1"),
])
def test_names_ending_in_digits_are_not_read_as_releases(name, reference_name):
    genome, _ = infer_genome(name)
    assert genome.reference_name == reference_name


def test_paths_are_not_read_as_releases():
    genome, _ = infer_genome("##reference=file:///var/lib/cwl/job367935311_index_001zdr/GRCh38.d1.vd1.fa")
    assert genome.reference_name == "GRCh38"
