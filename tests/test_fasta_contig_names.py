"""A chromosome FASTA is read by the annotation's contig names (#721).

``varcode.Genome(fasta=...)`` used to read the FASTA by the annotation's own
spelling, so a FASTA naming its contigs ``chr22`` gave an annotation that
uses ``22`` no DNA at all, and no warning.
"""

import gzip
import json
import warnings

import pytest
from pyensembl import Genome as AnnotationGenome

from varcode import Genome

from .native_reference_fixture import DATA, ReferenceWindows


class ListedWindows(ReferenceWindows):
    """CHM13 fixture DNA under ``prefix`` + its contig names, listed by ``keys()``."""

    def __init__(self, windows, contig=None, prefix=""):
        super().__init__([dict(w, contig=prefix + w["contig"]) for w in windows], contig)

    def keys(self):
        return sorted({window["contig"] for window in self.windows})


@pytest.fixture(scope="module")
def chm13(tmp_path_factory):
    from orfdata._pinned_genome import PinnedEnsemblGenome

    metadata = json.loads((DATA / "chm13.json").read_text())
    metadata.pop("sources")
    genome = PinnedEnsemblGenome(
        **metadata,
        gtf_path_or_url=str(DATA / "chm13.gtf.gz"),
        transcript_fasta_paths_or_urls=[str(DATA / "chm13.cdna.fa.gz")],
        protein_fasta_paths_or_urls=[str(DATA / "chm13.pep.fa.gz")],
        cache_directory_path=str(tmp_path_factory.mktemp("chm13-fasta-names")),
    )
    genome.index()
    windows = [window for window in json.loads(gzip.decompress(
        (DATA / "dna-windows.json.gz").read_bytes())) if window["genome"] == "chm13"]
    yield genome, windows
    genome.close()


@pytest.mark.parametrize("prefix", ["", "chr"])
def test_chm13_fasta_is_read_whichever_way_it_spells_contigs(chm13, prefix):
    # The issue's table: chr-named contigs gave '' and plain ones the DNA.
    annotation, windows = chm13
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # The transcript spot-check agrees.
        genome = Genome(annotation, fasta=ListedWindows(windows, prefix=prefix))
    assert genome.sequence("22", 21796720, 21796729) == "GTCCCTATGC"
    assert genome.reference_range("22", 21796720, 21796729) == "GTCCCTATGC"
    assert genome.reference_base("22", 21796720) == "G"


def test_reader_that_lists_no_contigs_is_read_by_exact_name(chm13):
    annotation, windows = chm13
    renamed = [dict(window, contig="chr" + window["contig"]) for window in windows]
    genome = Genome(annotation, fasta=ReferenceWindows(renamed), verify=False)
    assert genome.sequence("22", 21796720, 21796729) == ""
    assert genome.sequence("chr22", 21796720, 21796729) == "GTCCCTATGC"


def _annotation(tmp_path, reference_name, contigs):
    """A one-gene-per-contig annotation, so contigs() lists ``contigs``."""
    attributes = ('gene_id "g{i}"; gene_name "G{i}"; transcript_id "t{i}"; '
                  'transcript_name "T{i}"; transcript_biotype "lncRNA";')
    lines = []
    for i, contig in enumerate(contigs):
        for feature, extra in [("gene", ""), ("transcript", ""),
                               ("exon", ' exon_id "e{i}"; exon_number "1";')]:
            lines.append("%s\ttest\t%s\t1\t8\t.\t+\t.\t%s\n" % (
                contig, feature, (attributes + extra).format(i=i)))
    directory = tmp_path / ("annotation-%s-%s" % (reference_name, "-".join(contigs)))
    directory.mkdir()
    gtf = directory / "genes.gtf"
    gtf.write_text("".join(lines))
    genome = AnnotationGenome(
        reference_name=reference_name, annotation_name="fasta_names_test",
        gtf_path_or_url=str(gtf), cache_directory_path=str(directory / "cache"))
    genome.index()
    return genome


def _fasta(tmp_path, name, contigs):
    path = tmp_path / name
    path.write_text("".join(">%s\n%s\n" % item for item in contigs.items()))
    return path


def _mitochondrion(length):
    return ("GATCACAGGT" * (length // 10 + 1))[:length]


@pytest.mark.parametrize("annotation_contig, fasta_contig", [
    ("1", "chr1"),  # hg38.fa with Ensembl GRCh38
    ("chr1", "1"),  # the other way round
])
def test_pyfaidx_fasta_spelled_differently(tmp_path, annotation_contig, fasta_contig):
    annotation = _annotation(tmp_path, "GRCh38", [annotation_contig, "2"])
    path = _fasta(tmp_path, "reference.fa", {fasta_contig: "ACGTTGCA", "chr2": "GGGGCCCC"})
    genome = Genome(annotation, fasta=path)
    assert genome.sequence(annotation_contig, 2, 5) == "CGTT"
    assert genome.sequence("2", 1, 4) == "GGGG"
    assert genome.reference_base(annotation_contig, 1) == "A"
    # The reader keeps its own API: names, path, ownership.
    assert sorted(genome.fasta.keys()) == sorted([fasta_contig, "chr2"])
    assert annotation_contig in genome.fasta
    assert genome.to_dict()["fasta"] == str(path)
    reader = genome.fasta
    genome.close()
    assert reader.faidx.file.closed


def test_exact_name_wins_over_alias(tmp_path):
    annotation = _annotation(tmp_path, "GRCh38", ["1"])
    path = _fasta(tmp_path, "both.fa", {"chr1": "AAAAAAAA", "1": "CCCCCCCC"})
    genome = Genome(annotation, fasta=path, verify=False)
    assert genome.sequence("1", 1, 4) == "CCCC"
    assert genome.sequence("chr1", 1, 4) == "AAAA"


@pytest.mark.parametrize("as_path", [True, False])
def test_fasta_sharing_no_contig_raises(tmp_path, as_path):
    annotation = _annotation(tmp_path, "GRCh38", ["1", "2"])
    if as_path:
        fasta = _fasta(tmp_path, "other.fa", {"NC_045512.2": "ACGT"})
    else:
        fasta = {"NC_045512.2": "ACGT"}
    with pytest.raises(ValueError, match="shares no contig.*'NC_045512.2'.*'1', '2'"):
        Genome(annotation, fasta=fasta, verify=False)


def test_reader_with_listed_names_is_probed_by_a_shared_contig(tmp_path):
    # The probe used to try the annotation's first contig by its own name,
    # so a pre-opened chr-named reader raised TypeError.
    annotation = _annotation(tmp_path, "GRCh38", ["1", "2"])
    genome = Genome(annotation, fasta={"chr2": "TTTTGGGG"}, verify=False)
    assert genome.sequence("2", 5, 8) == "GGGG"
    assert genome.sequence("1", 1, 4) == ""


def test_assigning_a_reader_resolves_its_names(tmp_path):
    annotation = _annotation(tmp_path, "GRCh38", ["1"])
    genome = Genome(annotation)
    genome.fasta = {"chr1": "ACGTACGT"}
    assert genome.sequence("1", 1, 4) == "ACGT"
    genome.fasta = genome.fasta  # Reassigning the same reader is a no-op.
    assert genome.sequence("1", 5, 8) == "ACGT"
    with pytest.raises(ValueError, match="shares no contig"):
        genome.fasta = {"chrUn": "ACGT"}
    assert genome.sequence("1", 1, 4) == "ACGT"  # The old reader stays.


@pytest.mark.parametrize("length, readable", [(16571, False), (16569, True)])
def test_grch37_reads_chrm_as_mt_only_at_the_rcrs_length(tmp_path, length, readable):
    # hg19's chrM (NC_001807, 16,571 bp) is offset from the rCRS that
    # Ensembl's GRCh37 names MT (#707); GENCODE's GRCh37 chrM is the rCRS.
    annotation = _annotation(tmp_path, "GRCh37", ["1", "MT"])
    path = _fasta(tmp_path, "hg19.fa", {"chr1": "ACGTACGT", "chrM": _mitochondrion(length)})
    genome = Genome(annotation, fasta=path, verify=False)
    assert genome.sequence("1", 1, 4) == "ACGT"
    assert genome.sequence("MT", 1, 4) == ("GATC" if readable else "")


def test_grch37_never_reads_mt_as_an_annotation_chrm(tmp_path):
    annotation = _annotation(tmp_path, "GRCh37", ["1", "chrM"])
    path = _fasta(tmp_path, "b37.fa", {"1": "ACGTACGT", "MT": _mitochondrion(16569)})
    genome = Genome(annotation, fasta=path, verify=False)
    assert genome.sequence("chrM", 1, 4) == ""


@pytest.mark.parametrize("annotation_contig, fasta_contig", [("MT", "chrM"), ("chrM", "MT")])
def test_t2t_never_pairs_chrm_with_mt(tmp_path, annotation_contig, fasta_contig):
    # A T2T FASTA's chrM is CHM13's own mitochondrion in hs1 but the rCRS in
    # the _rCRS analysis sets; both are 16,569 bp.
    annotation = _annotation(tmp_path, "T2T-CHM13v2.0", ["1", annotation_contig])
    path = _fasta(tmp_path, "chm13v2.0_maskedY_rCRS.fa",
                  {"chr1": "ACGTACGT", fasta_contig: _mitochondrion(16569)})
    genome = Genome(annotation, fasta=path, verify=False)
    assert genome.sequence("1", 1, 4) == "ACGT"
    assert genome.sequence(annotation_contig, 1, 4) == ""


def test_t2t_reads_an_exactly_named_chrm(tmp_path):
    annotation = _annotation(tmp_path, "T2T-CHM13v2.0", ["1", "chrM"])
    path = _fasta(tmp_path, "rcrs.fa", {"chr1": "ACGTACGT", "chrM": _mitochondrion(16569)})
    genome = Genome(annotation, fasta=path, verify=False)
    assert genome.sequence("chrM", 1, 4) == "GATC"


def test_other_assemblies_read_chrm_as_mt(tmp_path):
    annotation = _annotation(tmp_path, "GRCm39", ["1", "MT"])
    path = _fasta(tmp_path, "mm39.fa", {"chr1": "ACGTACGT", "chrM": _mitochondrion(16299)})
    genome = Genome(annotation, fasta=path, verify=False)
    assert genome.sequence("MT", 1, 4) == "GATC"
