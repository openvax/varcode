"""Reference reads must reject uncovered SV flanks without per-base work (#407)."""

import logging

import pytest

from varcode.genome_sequence import reference_range


class Transcript:
    def __init__(self, intervals, sequence, strand="+"):
        self.exon_intervals = intervals
        self.sequence = sequence
        self.strand = strand
        positions = sorted(
            p for lo, hi in intervals for p in range(lo, hi + 1))
        if strand == "-":
            positions.reverse()
        self.offsets = {p: i for i, p in enumerate(positions)}
        self.offset_calls = 0

    def spliced_offset(self, position):
        self.offset_calls += 1
        if position not in self.offsets:
            raise ValueError("Not exonic")
        return self.offsets[position]


class Genome:
    def __init__(self, *transcripts):
        self.transcripts = transcripts
        self.queries = 0

    def transcripts_at_locus(self, contig, start, end):
        self.queries += 1
        return self.transcripts


@pytest.mark.parametrize("strand", ["+", "-"])
def test_uncovered_flank_does_not_map_every_base(strand):
    # A long exonic prefix, then a one-base gap shared by many isoforms.
    # Work must depend on the annotation size, not prefix length × isoforms.
    transcripts = [
        Transcript([(1, 1000), (1002, 1100)], "A" * 1099, strand)
        for _ in range(20)]
    genome = Genome(*transcripts)
    assert reference_range(genome, "1", 1, 1100) == ""
    assert genome.queries == 1
    assert sum(t.offset_calls for t in transcripts) <= 2 * len(transcripts)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_complementary_exons_keep_complete_range(strand):
    first = Transcript([(1, 3), (7, 9)],
                       "ACGCTA" if strand == "+" else "TAGCGT", strand)
    second = Transcript([(4, 6)], "TGA" if strand == "+" else "TCA", strand)
    # The union of both transcripts supplies all nine genomic positions.
    assert reference_range(Genome(first, second), "1", 1, 9) == "ACGTGACTA"
    assert first.offset_calls + second.offset_calls <= 10


def test_overlapping_exons_preserve_first_answer_and_diagnostics(caplog):
    first = Transcript([(1, 4)], "ACGT")
    second = Transcript([(3, 6)], "TTAA")
    with caplog.at_level(logging.DEBUG, logger="varcode.genome_sequence"):
        assert reference_range(Genome(first, second), "1", 1, 6) == "ACGTAA"
    assert "disagree" in caplog.text


@pytest.mark.parametrize("sequence", [None, ""])
def test_annotated_coverage_does_not_invent_missing_sequence(sequence):
    first = Transcript([(1, 3)], "ACG")
    second = Transcript([(4, 6)], sequence)
    assert reference_range(Genome(first, second), "1", 1, 6) == ""


def test_short_sequence_leaves_range_unresolved():
    first = Transcript([(1, 3)], "AC")
    second = Transcript([(4, 6)], "TGA")
    assert reference_range(Genome(first, second), "1", 1, 6) == ""


def test_short_sequence_can_still_supply_part_of_a_range():
    first = Transcript([(1, 3)], "AC")
    second = Transcript([(3, 6)], "GTGA")
    assert reference_range(Genome(first, second), "1", 1, 6) == "ACGTGA"


def test_adapter_without_exon_intervals_keeps_sequence_fallback():
    first = Transcript([(1, 3)], "ACG")
    second = Transcript([(4, 6)], "TGA")
    del second.exon_intervals
    assert reference_range(Genome(first, second), "1", 1, 6) == "ACGTGA"


def test_intronic_range_does_not_load_cdna():
    class UnloadedTranscript(Transcript):
        @property
        def sequence(self):
            pytest.fail("A discontinuous range must not load transcript FASTA")

        @sequence.setter
        def sequence(self, value):
            pass

    transcript = UnloadedTranscript([(1, 3), (5, 7)], None)
    assert reference_range(Genome(transcript), "1", 1, 7) == ""
