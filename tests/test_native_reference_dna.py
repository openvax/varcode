"""Native PyEnsembl DNA stays available through varcode.Genome (#488)."""

import inspect

import pytest
from pyensembl import Genome as AnnotationGenome

from varcode import Genome, Variant
from varcode.genome_sequence import reference_base, reference_range
from varcode.reference import infer_genome


@pytest.fixture
def native_genome(tmp_path):
    if "genome_fasta_path_or_url" not in inspect.signature(AnnotationGenome).parameters:
        pytest.skip("Native reference DNA requires PyEnsembl >= 2.11")
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">1\naacgttcaggtaccgattcgaacgt\n")
    genome = AnnotationGenome(
        reference_name="GRCh38", annotation_name="native_dna_test",
        genome_fasta_path_or_url=str(fasta),
        cache_directory_path=str(tmp_path / "cache"))
    yield genome
    genome.close()


def test_native_intronic_and_intergenic_queries(native_genome, tmp_path):
    gtf = tmp_path / "tiny.gtf"
    attributes = (
        'gene_id "g1"; gene_name "G1"; transcript_id "t1"; '
        'transcript_name "T1"; transcript_biotype "lncRNA";')
    gtf.write_text("".join(
        '1\ttest\t%s\t%d\t%d\t.\t-\t.\t%s\n' % (
            feature, start, end, attributes + extra)
        for feature, start, end, extra in [
            ("gene", 2, 13, ""), ("transcript", 2, 13, ""),
            ("exon", 2, 5, ' exon_id "e1"; exon_number "2";'),
            ("exon", 10, 13, ' exon_id "e2"; exon_number "1";')]))
    native = AnnotationGenome(
        reference_name="GRCh38", annotation_name="native_dna_annotation",
        gtf_path_or_url=str(gtf),
        genome_fasta_path_or_url=str(tmp_path / "reference.fa"),
        cache_directory_path=str(tmp_path / "annotated_cache"))
    try:
        native.index()
        wrapped = Genome(native)
        # This intron belongs to a minus-strand transcript; genomic reads
        # must still return plus-strand DNA, with inclusive endpoints.
        transcript = native.transcript_by_id("t1")
        assert transcript.contains("1", 6, 9)
        assert all(not exon.overlaps("1", 6, 9) for exon in transcript.exons)
        assert not native.transcripts_at_locus("1", 16, 20)
        assert wrapped.sequence("1", 6, 9) == "TCAG"
        assert wrapped.reference_range("1", 6, 9) == "TCAG"
        assert reference_base(wrapped, "1", 7) == "C"
        assert reference_range(wrapped, "1", 16, 20) == "ATTCG"
        assert wrapped.fasta is native.fasta
    finally:
        native.close()


def test_construction_inference_and_repr_do_not_open_native_dna(
        native_genome, monkeypatch):
    def unexpected_open(self):
        pytest.fail("Native FASTA opened before a sequence request")

    monkeypatch.setattr(AnnotationGenome, "fasta", property(unexpected_open))
    wrapped = Genome(native_genome)
    rewrapped = Genome(wrapped)
    assert infer_genome(native_genome) == (native_genome, False)
    assert infer_genome(wrapped) == (wrapped, False)
    assert Variant("1", 7, "C", "A", genome=rewrapped).genome is rewrapped
    assert "native FASTA configured" in repr(wrapped)
    wrapped.close()


def test_rewrap_tracks_native_reader_refresh(native_genome):
    wrapped = Genome(native_genome)
    rewrapped = Genome(wrapped)
    original = wrapped.fasta
    assert original is not None
    native_genome.clear_cache()
    assert wrapped.fasta is not original
    assert rewrapped.fasta is wrapped.fasta is native_genome.fasta
    assert str(original["1"][5:9]).upper() == "TCAG"
    original.close()


@pytest.mark.parametrize("rewrap", [False, True])
def test_closing_borrowing_wrapper_keeps_native_reader_usable(native_genome, rewrap):
    first = Genome(native_genome)
    second = Genome(first if rewrap else native_genome)
    reader = native_genome.fasta
    first.close()
    first.close()
    assert str(reader["1"][5:9]).upper() == "TCAG"
    assert second.sequence("1", 6, 9) == "TCAG"
    assert native_genome.fasta is reader


def test_explicit_override_and_rewrap_precede_native_dna(native_genome, tmp_path):
    override = tmp_path / "override.fa"
    override.write_text(">1\n" + "a" * 24 + "\n")
    wrapped = Genome(native_genome, fasta=override, verify=False)
    rewrapped = Genome(wrapped)
    assert rewrapped.fasta is wrapped.fasta
    assert rewrapped.sequence("1", 6, 9) == "AAAA"
    assert rewrapped.reference_range("1", 6, 9) == "AAAA"
    assert native_genome.sequence("1", 6, 9) == "TCAG"
    reader = wrapped.fasta
    rewrapped.close()
    assert not reader.faidx.file.closed
    assert wrapped.reference_base("1", 6) == "A"
    wrapped.close()
    assert reader.faidx.file.closed


def test_caller_provided_reader_is_borrowed(native_genome, tmp_path):
    from pyfaidx import Fasta

    path = tmp_path / "borrowed.fa"
    path.write_text(">1\ncccc\n")
    with Fasta(str(path)) as reader:
        wrapped = Genome(native_genome, fasta=reader, verify=False)
        wrapped.close()
        assert wrapped.sequence("1", 1, 4) == "CCCC"
        assert not reader.faidx.file.closed


def test_assigning_fasta_releases_owned_reader(native_genome, tmp_path):
    path = tmp_path / "override.fa"
    path.write_text(">1\ncccc\n")
    wrapped = Genome(native_genome, fasta=path, verify=False)
    reader = wrapped.fasta
    wrapped.fasta = None
    assert reader.faidx.file.closed
    assert wrapped.sequence("1", 6, 9) == "TCAG"


@pytest.mark.parametrize("remote", [False, True])
def test_uninstalled_native_dna_never_downloads(native_genome, tmp_path, monkeypatch, remote):
    from pyensembl.genome_fasta import GenomeFasta, MissingGenomeFastaError

    def unexpected_download(*args, **kwargs):
        pytest.fail("Reference lookup attempted a download")

    monkeypatch.setattr(GenomeFasta, "_download", unexpected_download)
    source = ("https://example.invalid/reference.fa" if remote
              else str(tmp_path / "missing.fa"))
    native = AnnotationGenome(
        reference_name="GRCh38", annotation_name="missing_native_dna",
        genome_fasta_path_or_url=source,
        cache_directory_path=str(tmp_path / "missing_cache"))
    monkeypatch.setattr(native, "transcripts_at_locus", lambda *args, **kwargs: [])
    wrapped = Genome(native)
    assert wrapped.fasta is None
    assert wrapped.reference_base("1", 7) == ""
    assert wrapped.reference_range("1", 6, 9) == ""
    with pytest.raises(MissingGenomeFastaError):
        wrapped.sequence("1", 6, 9)
    wrapped.close()


def test_equal_native_genomes_keep_distinct_dna(native_genome, tmp_path):
    path = tmp_path / "other.fa"
    path.write_text(">1\n" + "c" * 24 + "\n")
    other = AnnotationGenome(
        reference_name=native_genome.reference_name,
        annotation_name=native_genome.annotation_name,
        genome_fasta_path_or_url=str(path),
        cache_directory_path=str(tmp_path / "other_cache"))
    try:
        assert native_genome == other
        assert infer_genome(native_genome)[0] is native_genome
        assert infer_genome(other)[0] is other
        assert Genome(native_genome).sequence("1", 6, 9) == "TCAG"
        assert Genome(other).sequence("1", 6, 9) == "CCCC"
        assert Variant("1", 7, "C", "A", genome=other).genome is other
    finally:
        other.close()


@pytest.mark.parametrize("start,end", [(0, 3), (3, 2), (20, 30)])
def test_native_sequence_preserves_coordinate_errors(native_genome, start, end):
    with pytest.raises(ValueError):
        Genome(native_genome).sequence("1", start, end)


def test_native_sequence_preserves_missing_contig_error(native_genome):
    with pytest.raises(ValueError, match="Contig"):
        Genome(native_genome).sequence("missing", 1, 4)


def test_old_pyensembl_without_native_api(monkeypatch, tmp_path):
    # Simulate supported PyEnsembl versions before native DNA was added.
    monkeypatch.delattr(AnnotationGenome, "fasta", raising=False)
    monkeypatch.delattr(AnnotationGenome, "requires_genome_fasta", raising=False)
    native = AnnotationGenome(
        reference_name="GRCh38", annotation_name="old_api",
        cache_directory_path=str(tmp_path / "old_cache"))
    monkeypatch.setattr(native, "transcripts_at_locus", lambda *args, **kwargs: [])
    wrapped = Genome(native)
    assert wrapped.fasta is None
    assert wrapped.sequence("1", 1, 4) == ""
    assert wrapped.reference_base("1", 1) == ""
    wrapped.close()
