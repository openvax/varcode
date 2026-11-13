"""T2T-CHM13v2.0 is inferred from the names its files go by (#723), and the rCRS
chrM of its ``_rCRS`` analysis set is annotated with GRCh38's MT genes (#722).

ORFdata publishes the names, each file's contig names and its mitochondrion
(openvax/orfdata#204), and annotates the rCRS as contig NC_012920.1
(openvax/orfdata#206). These tests run offline: ``orfdata.load_genome`` is
replaced by the CHM13 fixture genome with ORFdata's packaged rCRS annotation
appended, built the way ORFdata builds the installed one.
"""

import json
import re
import warnings
from importlib.resources import files

import orfdata
import pytest

from pyensembl import Genome

from varcode import InvalidContigError, ReferenceMismatchError, Variant, VariantCollection, load_vcf
from varcode.pinned_references import load_pinned_genome
from varcode.reference import infer_genome, infer_reference_name

from .native_reference_fixture import DATA

RCRS = "NC_012920.1"


@pytest.fixture(scope="module")
def chm13_with_rcrs(tmp_path_factory):
    """The CHM13 fixture with ORFdata's packaged rCRS MT annotation appended."""
    from orfdata._pinned_genome import PinnedEnsemblGenome

    t2t, = orfdata.pinned_genomes(species="human")
    assets = t2t.mitochondria[RCRS]["assets"]
    packaged = {role: files("orfdata").joinpath(asset["resource"]) for role, asset in assets.items()}
    directory = tmp_path_factory.mktemp("chm13-rcrs")
    gtf = directory / "chm13+rcrs.gtf.gz"
    gtf.write_bytes((DATA / "chm13.gtf.gz").read_bytes() + packaged["gtf"].read_bytes())
    metadata = json.loads((DATA / "chm13.json").read_text())
    metadata.pop("sources")
    metadata["reference_metadata"]["contig_compartments"] = {RCRS: "mitochondrial"}
    genome = PinnedEnsemblGenome(
        **metadata, pinned_name="t2t-chm13",
        gtf_path_or_url=str(gtf),
        transcript_fasta_paths_or_urls=[str(DATA / "chm13.cdna.fa.gz"), str(packaged["transcripts"])],
        protein_fasta_paths_or_urls=[str(DATA / "chm13.pep.fa.gz"), str(packaged["proteins"])],
        cache_directory_path=str(directory))
    genome.index()
    yield genome
    genome.close()


@pytest.fixture(autouse=True)
def installed(chm13_with_rcrs, monkeypatch):
    """Every T2T name opens the fixture, as an installed t2t-chm13 would."""
    def load_genome(name, **kwargs):
        assert orfdata.genomes._pinned(name).name == "t2t-chm13"
        return chm13_with_rcrs
    monkeypatch.setattr(orfdata, "load_genome", load_genome)
    return chm13_with_rcrs


def resolve(name):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        return infer_reference_name(name)


@pytest.mark.parametrize("path, name", [
    ("T2T-CHM13v2.0", "T2T-CHM13v2.0"),
    ("t2t-chm13v2.0", "T2T-CHM13v2.0"),
    ("hs1", "hs1"),
    ("hs1.fa", "hs1"),
    ("/refs/UCSC/hs1.fa.gz", "hs1"),
    ("##reference=file:///refs/hs1.fa", "hs1"),
    ("chm13v2.0.fa", "chm13v2.0"),
    ("chm13v2.0_maskedY.fa.gz", "chm13v2.0_maskedY"),
    ("chm13v2.0_noY.fa", "chm13v2.0_noY"),
    # The longest name wins: this is the rCRS analysis set, not chm13v2.0.
    ("chm13v2.0_maskedY_rCRS.fa", "chm13v2.0_maskedY_rCRS"),
    ("/data/chm13v2.0_maskedY.rCRS.fa.gz", "chm13v2.0_maskedY.rCRS"),
    ("##reference=file:///refs/GRCh38/chm13v2.0_maskedY_rCRS.fa", "chm13v2.0_maskedY_rCRS"),
    ("GCA_009914755.4_T2T-CHM13v2.0_genomic.fna", "GCA_009914755.4"),
    ("GCF_009914755.1_T2T-CHM13v2.0_genomic.fna", "GCF_009914755.1"),
    ("Homo_sapiens-GCA_009914755.4-unmasked.fa.gz", "GCA_009914755.4"),
    ("chm13.draft_v1.1.fasta.gz", "chm13.draft_v1.1"),
    ("T2T-CHM13v1.1.fa", "T2T-CHM13v1.1"),
    # A version may be set off by a separator.
    ("T2T-CHM13-v2.0.fa", "T2T-CHM13v2.0"),
    ("T2T-CHM13 v2.0", "T2T-CHM13v2.0"),
    ("chm13.draft.v1.1.fasta", "chm13.draft_v1.1"),
])
def test_t2t_names_are_found_in_names_paths_and_headers(path, name):
    assert resolve(path) == name


def test_a_t2t_name_only_in_a_directory_is_used_with_the_usual_warning():
    with pytest.warns(UserWarning, match="could not be matched against filename"):
        assert infer_reference_name("/refs/hs1/genome.fa") == "hs1"


@pytest.mark.parametrize("path", ["ths1.fa", "hs10.fa", "phs1x.fasta", "chm13.fa", "CHM13.fasta"])
def test_hs1_only_as_a_whole_word_and_bare_chm13_names_nothing(path):
    with pytest.raises(ValueError, match="Failed to infer genome assembly name"):
        infer_reference_name(path)


@pytest.mark.parametrize("path", [
    "chm13.draft_v1.0.fasta.gz", "GCA_009914755.2_T2T-CHM13v1.0_genomic.fna", "T2T-CHM13v1.0",
    "##reference=file:///refs/chm13.draft_v0.9.fasta",
    # Separated from the version, these used to match ORFdata's genome name, t2t-chm13.
    "T2T-CHM13-v1.0.fa", "T2T-CHM13_v1.0.fasta", "T2T-CHM13 v1.0", "t2t-chm13.v1.0.fa"])
def test_older_t2t_versions_are_refused_because_their_coordinates_differ(path):
    with pytest.raises(ValueError, match="is not a supported reference: .*differ from T2T-CHM13v2.0"):
        infer_reference_name(path)
    with pytest.raises(ValueError, match="not a supported reference"):
        infer_genome(path)


@pytest.mark.parametrize("path, name", [
    ("GRCh38.fa", "GRCh38"),
    ("/data/hs1_project/GRCh38.fa", "GRCh38"),
    ("Homo_sapiens_assembly38.fasta", "hg38"),
    ("hs38DH.fa", "hg38"),
    ("##reference=file:///refs/chm13/GRCh38_full_analysis_set_plus_decoy_hla.fa", "GRCh38"),
])
def test_grch38_names_do_not_resolve_to_t2t(path, name):
    assert resolve(path) == name


def test_a_file_name_naming_t2t_and_another_assembly_is_ambiguous():
    with pytest.raises(ValueError, match="names more than one reference"):
        infer_reference_name("GRCh38_and_chm13v2.0.fa")


def test_t2t_names_open_orfdatas_pinned_genome(installed):
    for name in ("hs1", "chm13v2.0_maskedY_rCRS.fa", "T2T-CHM13v2.0:2022_07", "GCF_009914755.1"):
        genome, converted = infer_genome(name)
        assert genome is installed and converted
    with pytest.raises(ValueError, match="one pinned annotation, ensembl-rapid 2022_07.*2023_03"):
        infer_genome("T2T-CHM13v2.0:2023_03")


def test_a_genome_that_is_not_installed_names_the_command(monkeypatch, tmp_path):
    from orfdata.genomes import load_genome as real

    monkeypatch.setattr(orfdata, "load_genome", real)
    monkeypatch.setenv("ORFDATA_CACHE", str(tmp_path))
    with pytest.raises(orfdata.NotInstalledError) as error:
        infer_genome("/refs/not-installed/chm13v2.0_noY.fa")
    message = str(error.value)
    assert message.startswith("chm13v2.0_noY (from '/refs/not-installed/chm13v2.0_noY.fa') is T2T-CHM13v2.0")
    assert "'orfdata genomes install t2t-chm13'" in message
    # Not through infer_genome, which memoizes: other tests open T2T-CHM13v2.0 by name.
    with pytest.raises(orfdata.NotInstalledError, match="^The annotation of T2T-CHM13v2.0 comes from ORFdata"):
        load_pinned_genome("T2T-CHM13v2.0", source="T2T-CHM13v2.0")
    assert isinstance(error.value, ValueError) and isinstance(error.value, FileNotFoundError)


@pytest.mark.parametrize("contig", ["chr1", "1", "CP068277.2", "NC_060925.1"])
def test_contigs_are_named_by_orfdatas_table(contig):
    variant = Variant(contig, 48121278, "G", "T", genome="hs1")
    assert (variant.contig, variant.reference_name, variant.original_reference_name) == (
        "1", "T2T-CHM13v2.0", "hs1")
    assert "ENST05220073379" in variant.transcript_ids


# Classic pathogenic variants, in rCRS coordinates, with Ensembl GRCh38 115's
# consequence for each on MT: the same molecule, so the same answer.
RCRS_VARIANTS = [
    (8993, "T", "G", "MT-ATP6", "p.L156R", "Substitution"),  # NARP
    (6340, "C", "T", "MT-CO1", "p.T146I", "Substitution"),
    (14484, "T", "C", "MT-ND6", "p.M64V", "Substitution"),  # LHON; MT-ND6 is on the minus strand
    (3243, "A", "G", "MT-TL1", "non-coding-transcript", "NoncodingTranscript"),  # MELAS
    # LHON. Polyadenylation completes MT-ND1's stop codon, so Ensembl's GTF
    # has none and varcode reads the transcript as incomplete, on GRCh38 too.
    (3460, "G", "A", "MT-ND1", "incomplete", "IncompleteTranscript"),
]


@pytest.mark.parametrize("position, ref, alt, gene, change, kind", RCRS_VARIANTS)
def test_the_rcrs_chrm_of_an_rcrs_file_is_annotated_with_grch38s_mt_genes(
        position, ref, alt, gene, change, kind):
    variant = Variant("chrM", position, ref, alt, genome="chm13v2.0_maskedY_rCRS.fa")
    assert (variant.contig, variant.chrm_is_rcrs) == (RCRS, True)
    effect = variant.effects().top_priority_effect()
    assert (effect.gene_name, effect.short_description, type(effect).__name__) == (gene, change, kind)
    assert effect.transcript.translation_table == 2


def test_chm13s_own_chrm_is_unannotated_and_the_warning_says_why():
    variants = [Variant("chrM", 3460, "G", "A", genome=name) for name in ("hs1", "chm13v2.0.fa")]
    assert {(v.contig, v.chrm_is_rcrs) for v in variants} == {("MT", False)}
    with pytest.raises(Exception, match="own mitochondrion has no annotation"):
        variants[0].effects()
    from varcode import VariantCollection

    with pytest.warns(UserWarning, match=r"T2T-CHM13v2.0:MT; retained as failures\. T2T-CHM13v2.0's own "
                                         r"mitochondrion has no annotation.*chrm_is_rcrs=True"):
        VariantCollection(variants).effects()
    # Saying explicitly that this chrM is the rCRS names it NC_012920.1.
    assert Variant("chrM", 3460, "G", "A", genome="hs1", chrm_is_rcrs=True).contig == RCRS


VCF = """##fileformat=VCFv4.2
##reference={reference}
##contig=<ID=chr1,length={chr1}>
##contig=<ID=chrM,length=16569>
##contig=<ID=chrEBV,length=171823>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t48121278\t.\tG\tT\t.\tPASS\t.
chrM\t3460\t.\tG\tA\t.\tPASS\t.
"""


def write_vcf(tmp_path, reference, chr1=248387328, key="reference"):
    """A VCF with T2T's ##contig lengths, GRCh38's (chr1=248956422) or,
    with chr1=None, none, as many VCFs have."""
    text = VCF.format(reference=reference, chr1=chr1).replace("##reference=", "##%s=" % key)
    if chr1 is None:
        text = "".join(line for line in text.splitlines(True) if not line.startswith("##contig"))
    path = tmp_path / "t2t.vcf"
    path.write_text(text)
    return path


def test_load_vcf_infers_t2t_and_its_mitochondrion_from_the_reference(tmp_path, installed):
    path = write_vcf(tmp_path, "file:///refs/chm13v2.0_maskedY_rCRS.fa")
    variants = load_vcf(path)
    assert [(v.contig, v.chrm_is_rcrs, v.original_reference_name) for v in variants] == [
        ("1", True, "chm13v2.0_maskedY_rCRS"), (RCRS, True, "chm13v2.0_maskedY_rCRS")]
    assert {v.genome for v in variants} == {installed}
    effects = variants.effects()
    assert {e.gene_name for e in effects} >= {"SLC5A9", "MT-ND1"}
    assert "MT-ND1" in variants[1].gene_names
    # The written VCF names the same file, so it reads back the same way.
    from varcode.vcf_output import variants_to_vcf

    written = tmp_path / "written.vcf"
    with open(written, "w") as out:
        variants_to_vcf(variants, variants.source_to_metadata_dict[str(path)], out=out)
    assert "##reference=chm13v2.0_maskedY_rCRS\n" in written.read_text()
    assert [v.contig for v in load_vcf(written)] == ["1", RCRS]


def test_load_vcf_reads_chrm_from_the_header_file_when_given_the_genome(tmp_path, installed):
    rcrs = load_vcf(write_vcf(tmp_path, "file:///refs/chm13v2.0_maskedY_rCRS.fa"), genome=installed)
    own = load_vcf(write_vcf(tmp_path, "file:///refs/hs1.fa"), genome=installed)
    assert [v.contig for v in rcrs] == ["1", RCRS] and [v.contig for v in own] == ["1", "MT"]


def test_load_vcf_refuses_a_grch38_vcf_read_as_t2t(tmp_path):
    path = write_vcf(tmp_path, "file:///refs/chm13v2.0.fa", chr1=248956422)
    with pytest.raises(ValueError, match="is not on T2T-CHM13v2.0: .*chr1 is 248,956,422 bp here but "
                                         "248,387,328 bp in T2T-CHM13v2.0"):
        load_vcf(path)


def test_vcf_contig_lengths_are_checked_against_the_table_only_for_t2t(tmp_path):
    # chrEBV isn't in the table, and chrM is the rCRS's 16,569 bp either way.
    assert len(load_vcf(write_vcf(tmp_path, "file:///refs/hs1.fa"))) == 2


# Negative controls (#669): another assembly's VCF, alleles and contigs are
# refused rather than annotated on T2T's transcripts.

@pytest.mark.parametrize("reference, named", [
    ("GRCh38", "GRCh38"),
    ("file:///refs/Homo_sapiens_assembly38.fasta", "GRCh38"),
    ("hg38.fa", "GRCh38"),
    ("file:///refs/chm13.draft_v1.0.fasta", "chm13.draft_v1.0"),  # T2T-CHM13v1.0: other coordinates
    ("file:///refs/T2T-CHM13-v1.0.fa", "T2T-CHM13v1.0"),
    # The file name decides, not a directory.
    ("file:///data/t2t-chm13_compare/Homo_sapiens_assembly38.fasta", "GRCh38"),
    ("file:///refs/hs1/GRCh38.fa", "GRCh38"),
])
@pytest.mark.parametrize("by_name", [False, True])
def test_load_vcf_refuses_a_header_naming_another_assembly_on_t2t(
        tmp_path, installed, reference, named, by_name):
    path = write_vcf(tmp_path, reference, chr1=None)
    genome = "T2T-CHM13v2.0" if by_name else installed
    with pytest.raises(ValueError, match="its ##reference, .*, names %s, but genome= is T2T-CHM13v2.0"
                                         % re.escape(named)):
        load_vcf(path, genome=genome)


@pytest.mark.parametrize("chr1", [None, 248387328])  # no ##contig lengths, or T2T's
def test_load_vcf_refuses_a_t2t_header_on_grch38(tmp_path, chr1):
    with pytest.raises(ValueError, match="names T2T-CHM13v2.0, but genome= is GRCh38"):
        load_vcf(write_vcf(tmp_path, "file:///refs/hs1.fa", chr1=chr1), genome=81)


@pytest.mark.parametrize("reference", [
    "file:///refs/hs1.fa",
    "file:///refs/chm13v2.0_maskedY_rCRS.fa",
    "file:///refs/T2T-CHM13-v2.0.fa",
    "T2T-CHM13v1.1",  # v2.0 without chrY: the same coordinates
    "file:///refs/reference.fa",  # names no assembly, so isn't checked
    # Directories don't count, nor does an assembly name inside another word.
    "file:///refs/GRCh38_to_T2T/chm13v2.0.fa",
    "file:///home/rob38/genome.fa",
    "file:///work/a1/b38c9f0e/reference.fa",
    "file:///refs/GRCh38_and_chm13v2.0.fa",  # names two assemblies
])
@pytest.mark.parametrize("chr1", [None, 248387328])
def test_load_vcf_accepts_a_header_that_agrees_or_names_nothing(tmp_path, installed, reference, chr1):
    assert len(load_vcf(write_vcf(tmp_path, reference, chr1=chr1), genome=installed)) == 2


def test_contig_lengths_outweigh_a_stale_reference(tmp_path, installed):
    # bcftools +liftover writes the target's ##contig lines but keeps the
    # source's ##reference; Picard writes the target's and CrossMap drops it.
    lifted = write_vcf(tmp_path, "file:///refs/Homo_sapiens_assembly38.fasta")
    assert [v.contig for v in load_vcf(lifted, genome=installed)] == ["1", "MT"]
    lifted = write_vcf(tmp_path, "file:///refs/hs1.fa", chr1=248956422)
    assert [v.contig for v in load_vcf(lifted, genome=81)] == ["1", "MT"]


def test_reference_vcf_key_none_ignores_a_wrong_header(tmp_path, installed):
    variants = load_vcf(write_vcf(tmp_path, "GRCh38", chr1=None), genome=installed, reference_vcf_key=None)
    assert [v.contig for v in variants] == ["1", "MT"]


def test_another_reference_key_is_checked_too(tmp_path, installed):
    # Keys other than reference hold a list of values.
    path = write_vcf(tmp_path, "GRCh38", chr1=None, key="assembly")
    with pytest.raises(ValueError, match="its ##assembly, 'GRCh38', names GRCh38"):
        load_vcf(path, genome=installed, reference_vcf_key="assembly")


@pytest.fixture(scope="module")
def custom_genome(tmp_path_factory):
    """A T2T annotation a user built, named as they chose."""
    genomes = []

    def build(reference_name):
        genome = Genome(
            reference_name=reference_name, annotation_name="custom",
            gtf_path_or_url=str(DATA / "chm13.gtf.gz"),
            transcript_fasta_paths_or_urls=[str(DATA / "chm13.cdna.fa.gz")],
            protein_fasta_paths_or_urls=[str(DATA / "chm13.pep.fa.gz")],
            cache_directory_path=str(tmp_path_factory.mktemp(reference_name)))
        genome.index()
        genomes.append(genome)
        return genome

    yield build
    for genome in genomes:
        genome.close()


def test_a_genome_named_by_a_t2t_alias_is_t2t(tmp_path, custom_genome):
    genome = custom_genome("hs1")
    for reference in ("file:///refs/chm13v2.0.fa", "file:///refs/hs1.fa"):
        assert len(load_vcf(write_vcf(tmp_path, reference, chr1=None), genome=genome)) == 2
    with pytest.raises(ValueError, match="names GRCh38, but genome= is T2T-CHM13v2.0"):
        load_vcf(write_vcf(tmp_path, "GRCh38", chr1=None), genome=genome)


@pytest.mark.parametrize("name", ["T2T-CHM13v2.0", "hs1", "chm13v2.0", "HS1"])
def test_every_name_of_a_t2t_genome_gets_the_assembly_rules(tmp_path, custom_genome, name):
    # The rules belong to the assembly, not to what the user called the annotation (#730).
    genome = custom_genome(name)
    for reference in ("reference.fa", "file:///refs/hs1.fa"):
        with pytest.raises(ValueError, match="##contig lengths disagree"):
            load_vcf(write_vcf(tmp_path, reference, chr1=248956422), genome=genome)  # GRCh38's
    assert [v.contig for v in load_vcf(write_vcf(tmp_path, "reference.fa", chr1=248387328), genome=genome)] == [
        "1", "MT"]
    variant = Variant("CP068277.2", 48121278, "G", "T", genome=genome, convert_ucsc_contig_names=True)
    assert (variant.contig, variant.reference_name) == ("1", name)
    assert Variant("chrM", 1, "A", "T", genome=genome, convert_ucsc_contig_names=True).contig == "MT"
    from varcode.pinned_references import mitochondrion_note

    assert mitochondrion_note("MT", name) is not None


@pytest.mark.parametrize("name", ["T2T-CHM13v1.0", "chm13v1.0"])
def test_another_version_of_t2t_stays_unpinned(custom_genome, name):
    from varcode.pinned_references import is_pinned_assembly, pinned_assembly

    assert not is_pinned_assembly(name) and pinned_assembly(name) is None
    variant = Variant("CP068277.2", 1000, "G", "T", genome=custom_genome(name),
                      convert_ucsc_contig_names=True)
    assert variant.contig == "CP068277.2"


def test_two_assemblies_neither_pinned_are_not_compared(tmp_path, custom_genome):
    genome = custom_genome("T2T-CHM13v1.0")
    path = write_vcf(tmp_path, "file:///refs/chm13.draft_v1.0.fasta", chr1=None)
    assert len(load_vcf(path, genome=genome)) == 2


def test_an_explicit_genome_leaves_header_inference_alone(tmp_path):
    # infer_reference_name memoizes, so a check that ran it with warnings
    # silenced would silence them for the rest of the process.
    reference = "file:///data/GRCh38/x/reference_669.fa"
    load_vcf(write_vcf(tmp_path, reference, chr1=None), genome=81)
    with pytest.warns(UserWarning, match="could not be matched against filename"):
        infer_reference_name(reference)


@pytest.mark.parametrize("contig, position, ref, alt, strand", [
    ("1", 48121278, "A", "T", "+"),  # SLC5A9; CHM13 has G
    ("17", 70030393, "C", "A", "-"),  # ABCA10; CHM13 has G
])
def test_another_assemblys_allele_is_refused_on_either_strand(installed, contig, position, ref, alt, strand):
    with pytest.raises(ReferenceMismatchError, match=re.escape("strand %s" % strand)):
        Variant(contig, position, ref, alt, installed).effects()


@pytest.mark.parametrize("contig", ["chrUn_KI270302v1", "KI270302.1", "chrEBV"])
def test_grch38_only_contigs_are_not_placed_on_t2t(installed, contig):
    # GRCh38's unplaced scaffolds and EBV decoy have no T2T counterpart.
    variant = Variant(contig, 100, "A", "T", installed)
    with pytest.raises(InvalidContigError):
        variant.effects()
    with pytest.warns(UserWarning, match="absent from the annotation.*retained as failures"):
        effects = VariantCollection([variant]).effects()
    assert [type(effect).__name__ for effect in effects] == ["Failure"]


def test_effects_record_the_pinned_annotation(tmp_path, installed):
    effects = VariantCollection([Variant("1", 48121278, "G", "T", installed)]).effects()
    path = tmp_path / "effects.csv"
    effects.to_csv(str(path))
    header = {line for line in path.read_text().splitlines() if line.startswith("# ")}
    assert {"# reference_name=T2T-CHM13v2.0", "# annotation_name=ensembl-rapid",
            "# annotation_version=2022_07", "# assembly_accession=GCA_009914755.4"} <= header


def test_pinned_contig_table_has_every_spelling():
    from varcode.pinned_references import pinned_contig_lengths, pinned_contig_name

    for chrm_is_rcrs, mitochondrion in ((False, "MT"), (True, RCRS)):
        assert pinned_contig_name("chrM", "T2T-CHM13v2.0", chrm_is_rcrs) == mitochondrion
        assert pinned_contig_lengths("T2T-CHM13v2.0", chrm_is_rcrs)[mitochondrion] == 16569
    assert pinned_contig_name("chrx", "T2T-CHM13v2.0") == "X"
    assert pinned_contig_name("CP086569.2", "T2T-CHM13v2.0") == "Y"
    assert pinned_contig_name("chr1", "GRCh38") is None


class NamedReader(dict):
    """A FASTA reader of whole sequences, under a file name."""

    def __init__(self, filename, sequences):
        super().__init__(sequences)
        self.filename = filename


def test_an_attached_fasta_is_read_by_orfdatas_contig_names(installed):
    from varcode import Genome

    chr1, mitochondrion = "ACGT" * 13000, "GATCACAGGT" + "A" * 16559

    def read(filename, sequences, contig, end):
        genome = Genome(installed, fasta=NamedReader(filename, sequences), verify=False)
        return genome.sequence(contig, 1, end)

    # NCBI's GenBank FASTA names chromosomes by accession.
    assert read("/refs/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna", {"CP068277.2": chr1}, "1", 8) == "ACGTACGT"
    # The _rCRS analysis set's chrM is the rCRS, which the annotation calls NC_012920.1.
    rcrs = {"chr1": chr1, "chrM": mitochondrion}
    assert read("/refs/chm13v2.0_maskedY_rCRS.fa", rcrs, RCRS, 10) == "GATCACAGGT"
    # hs1's chrM is CHM13's own, and an unnamed file's could be either: neither serves the rCRS.
    for filename in ("/refs/hs1.fa", "/refs/genome.fa", None):
        assert read(filename, rcrs, RCRS, 10) == "" and read(filename, rcrs, "1", 4) == "ACGT"
