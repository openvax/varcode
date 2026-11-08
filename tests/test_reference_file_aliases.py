"""Common reference filenames must work through both VCF entry points (#712)."""

import warnings

import pytest

from varcode import SampleSpec, load_vcf, load_vcf_samples
from varcode.reference import infer_genome, infer_reference_name


REFERENCES = [
    ("hs37d5", "GRCh37", "GRCh37", "1", "MT", 75),
    ("human_g1k_v37", "GRCh37", "GRCh37", "1", "MT", 75),
    ("Homo_sapiens_assembly19", "GRCh37", "GRCh37", "1", "MT", 75),
    ("Homo_sapiens_assembly38", "hg38", "GRCh38", "chr1", "chrM", 93),
]


@pytest.mark.parametrize("name,resolved,assembly,chrom,mt,release", REFERENCES)
@pytest.mark.parametrize("template", ["{}", "/refs/{}.fasta", "file:///refs/{}.fa.gz"])
def test_common_reference_names(name, resolved, assembly, chrom, mt, release, template):
    value = template.format(name)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert infer_reference_name.__wrapped__(value) == resolved
    genome, converted = infer_genome(value)
    assert (genome.reference_name, converted) == (assembly, resolved == "hg38")


@pytest.mark.parametrize("name,resolved,assembly,chrom,mt,release", REFERENCES)
def test_alias_release_and_filename_precedence(name, resolved, assembly, chrom, mt, release):
    genome, converted = infer_genome("%s:%s" % (name, release))
    assert (genome.reference_name, genome.release, converted) == (
        assembly, release, resolved == "hg38")
    other_build = "GRCh38" if assembly == "GRCh37" else "GRCh37"
    # The actual filename wins over a directory naming a different build.
    assert infer_reference_name("/%s/%s.fasta" % (other_build, name)) == resolved


@pytest.mark.parametrize("name,resolved,assembly,chrom,mt,release", REFERENCES)
def test_vcf_loaders_infer_common_reference_filenames(
        tmp_path, name, resolved, assembly, chrom, mt, release):
    path = tmp_path / "variants.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##reference=file:///refs/%s.fasta.gz\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor\n"
        "%s\t100\t.\tA\tG\t.\tPASS\t.\tGT:DP:GQ\t0/1:30:60\n"
        "%s\t100\t.\tA\tG\t.\tPASS\t.\tGT:DP:GQ\t0/1:30:60\n"
        % (name, chrom, mt))
    variants = load_vcf(str(path), convert_ucsc_contig_names=None)
    assert {v.contig for v in variants} == {"1", "MT"}
    assert {v.reference_name for v in variants} == {assembly}
    assert {v.original_reference_name for v in variants} == {resolved}
    sample, = load_vcf_samples(path, kind="somatic", samples=[SampleSpec("tumor", role="tumor")])
    assert sample.assembly == assembly
    assert set(sample.somatic) == {
        ("small", "1", 100, "A", "G"), ("small", "MT", 100, "A", "G")}
    other_build = "GRCh38" if assembly == "GRCh37" else "GRCh37"
    with pytest.raises(ValueError, match="contradicts"):
        load_vcf_samples(path, assembly=other_build)
