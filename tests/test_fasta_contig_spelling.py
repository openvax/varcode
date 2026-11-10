"""varcode.Genome reads a FASTA by its own contig spelling (#721).

A FASTA that names contigs ``chr22`` used to give an annotation that names
them ``22`` no DNA at all, and no warning.
"""

import gzip
import json
import warnings

import pytest
from pyensembl import cached_release

from varcode import Genome
from varcode.genome_sequence import reference_base, reference_range

from .native_reference_fixture import DATA

# rCRS positions 1-9 (NC_012920.1).
RCRS_START = "GATCACAGG"


class _Record:
    """One contig, ``sequence`` placed at 1-based ``offset``, with pyfaidx's
    slice and length protocol."""

    def __init__(self, sequence, offset=1, length=None):
        self.sequence = sequence
        self.offset = offset
        self.length = offset - 1 + len(sequence) if length is None else length

    def __getitem__(self, s):
        lo = max(0, s.start - self.offset + 1)
        hi = max(0, s.stop - self.offset + 1)
        return self.sequence[lo:hi]

    def __len__(self):
        return self.length


class _Fasta:
    """A reader that lists its names, as ``pyfaidx.Fasta`` does."""

    def __init__(self, records):
        self.records = records

    def keys(self):
        return self.records.keys()

    def __getitem__(self, name):
        return self.records[name]


@pytest.fixture(scope="module")
def grch38():
    return cached_release(81)


@pytest.fixture(scope="module")
def cftr_base(grch38):
    """(position, base) inside CFTR's first exon on chromosome 7."""
    transcript = grch38.transcript_by_id("ENST00000003084")
    return transcript.exons[0].start + 5, transcript.sequence[5].upper()


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
        cache_directory_path=str(tmp_path_factory.mktemp("chm13")))
    genome.index()
    yield genome
    genome.close()


def _other(base):
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[base]


def test_ucsc_named_fasta_serves_an_ensembl_annotation(grch38, cftr_base):
    position, base = cftr_base
    # A base the cDNA doesn't have shows the FASTA was read.
    window = "A" + _other(base) + "C"
    fasta = _Fasta({"chr7": _Record(window, offset=position - 1)})
    genome = Genome(grch38, fasta=fasta, verify=False)
    assert genome.fasta is fasta
    assert genome.sequence("7", position - 1, position + 1) == window
    assert genome.reference_base("7", position) == _other(base)
    assert reference_base(genome, "7", position) == _other(base)
    assert reference_range(genome, "7", position - 1, position + 1) == window
    assert Genome(genome).sequence("7", position, position) == _other(base)


def test_a_name_the_fasta_lists_wins(grch38, cftr_base):
    position, _ = cftr_base
    fasta = _Fasta({"7": _Record("A", offset=position),
                    "chr7": _Record("C", offset=position)})
    genome = Genome(grch38, fasta=fasta, verify=False)
    assert genome.sequence("7", position, position) == "A"
    assert genome.sequence("chr7", position, position) == "C"


def test_ambiguous_spellings_are_not_guessed(grch38, cftr_base):
    position, _ = cftr_base
    fasta = _Fasta({"chr7": _Record("A", offset=position),
                    "Chr7": _Record("C", offset=position)})
    genome = Genome(grch38, fasta=fasta, verify=False)
    assert genome.sequence("7", position, position) == ""


@pytest.mark.parametrize("release,records,contig,expected", [
    # hg38's chrM is the rCRS that Ensembl's GRCh38 names MT.
    (81, {"chrM": _Record(RCRS_START, length=16569)}, "MT", "GAT"),
    # Under GRCh37, a chrM of the rCRS length is the rCRS, as in GENCODE's
    # GRCh37 (#710); hg19's own chrM is 16,571 bp and offset from it (#707).
    (75, {"chrM": _Record(RCRS_START, length=16569)}, "MT", "GAT"),
    (75, {"chrM": _Record(RCRS_START, length=16571)}, "MT", ""),
    # hg19's chrM, which variants keep as chrM, never reads the rCRS MT.
    (75, {"MT": _Record(RCRS_START, length=16569)}, "chrM", ""),
])
def test_a_respelled_mitochondrion_is_read_only_when_it_is_the_annotations(
        release, records, contig, expected):
    genome = Genome(cached_release(release), fasta=_Fasta(records), verify=False)
    assert genome.sequence(contig, 1, 3) == expected


def test_t2t_mitochondrion_is_never_chosen_by_spelling(chm13):
    # chm13v2.0.fa's chrM is CHM13's own (CP068254.1); the _rCRS set's is the
    # rCRS. Both are 16,569 bp, with positions 576 to 577 bp apart (#722).
    fasta = _Fasta({"chr22": _Record("ACGT"), "chrM": _Record(RCRS_START, length=16569)})
    genome = Genome(chm13, fasta=fasta, verify=False)
    assert genome.sequence("MT", 1, 3) == ""
    assert genome.sequence("22", 1, 4) == "ACGT"


def test_t2t_consortium_fasta_serves_ensembl_annotation(chm13, tmp_path):
    """The case in #721: a chr-named FASTA file, read through pyfaidx."""
    windows = json.loads(gzip.decompress((DATA / "dna-windows.json.gz").read_bytes()))
    window = next(w for w in windows if w["genome"] == "chm13" and w["contig"] == "22")
    sequence = window["sequence"].upper()
    path = tmp_path / "chm13v2.0.fa"
    with open(path, "w") as handle:
        handle.write(">chr22\n%s%s\n" % ("N" * (window["start"] - 1), window["sequence"]))
    start = window["start"]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        genome = Genome(chm13, fasta=str(path))
    assert genome.sequence("22", start, start + 9) == sequence[:10]
    end = start + len(sequence) - 1
    assert genome.reference_range("22", end - 9, end) == sequence[-10:]
    genome.close()


def test_a_fasta_naming_none_of_the_contigs_raises(grch38, cftr_base):
    position, base = cftr_base
    # NCBI's FASTAs name chromosomes by accession (CP068277.2 is T2T's chr1).
    fasta = _Fasta({"CP068277.2": _Record(base, offset=position)})
    with pytest.raises(ValueError, match="names none of the .* contigs in the GRCh38"):
        Genome(grch38, fasta=fasta)
    assert Genome(grch38, fasta=fasta, verify=False).sequence("7", position, position) == ""


def test_verify_compares_a_respelled_fasta(grch38):
    class Wrong(_Record):
        def __getitem__(self, s):
            return "X" * (s.stop - s.start)

    names = ["chr%s" % c for c in list(range(1, 23)) + ["X", "Y", "M"]]
    fasta = _Fasta({name: Wrong("", length=10 ** 9) for name in names})
    with pytest.warns(UserWarning, match="disagree"):
        Genome(grch38, fasta=fasta)
