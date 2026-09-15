"""Tests for splice-site detection on realized forward/reverse layouts."""

from types import SimpleNamespace

import pytest

from varcode.genomic_layout import GenomicLayout
from varcode.nucleotides import reverse_complement
from varcode.splice_graph import (
    disrupted_splice_sites,
    enumerate_layout_splice_plans,
)
from varcode.transcript_layout import build_exon_runs, realize_exon_path


FORWARD = "ATGGTAGCAAGTAGGCTGTAGTAG"


def _long_forward_sequence():
    bases = list("A" * 200)
    for start, sequence in (
            (1, "ATG"), (4, "GT"),
            (18, "AG"), (20, "CAA"), (23, "GT"),
            (29, "CAGGTAAGT"),
            (118, "AG"), (120, "GCT"), (123, "GT"),
            (138, "AG"), (140, "TAG")):
        bases[start - 1:start - 1 + len(sequence)] = sequence
    return "".join(bases)


LONG_FORWARD = _long_forward_sequence()


class Exon:
    def __init__(self, start, end):
        self.start = start
        self.end = end


class Transcript:
    contig = "x"
    start = 1
    end = 24
    is_protein_coding = True
    complete = True
    protein_sequence = "MQA"
    stop_codon_spliced_offsets = (9, 10, 11)
    id = "synthetic-splice"
    name = "SYN-001"

    def __init__(self, reverse=False):
        self.gene = SimpleNamespace(id="synthetic-gene", name="SYN")
        self.on_backward_strand = reverse
        intervals = ((1, 3), (8, 10), (15, 17), (22, 24))
        if reverse:
            intervals = tuple(
                (25 - end, 25 - start)
                for start, end in reversed(intervals))
            self.start_codon_positions = (22, 23, 24)
        else:
            self.start_codon_positions = (1, 2, 3)
        self.start_codon_spliced_offsets = (0, 1, 2)
        self.exons = tuple(Exon(*interval) for interval in intervals)

    def spliced_offset(self, position):
        offset = 0
        exons = sorted(
            self.exons, key=lambda exon: exon.start,
            reverse=self.on_backward_strand)
        for exon in exons:
            if exon.start <= position <= exon.end:
                if self.on_backward_strand:
                    return offset + exon.end - position
                return offset + position - exon.start
            offset += exon.end - exon.start + 1
        raise ValueError(position)


class LongIntronTranscript(Transcript):
    end = 142
    id = "synthetic-long-intron"

    def __init__(self, reverse=False):
        self.gene = SimpleNamespace(id="synthetic-gene", name="SYN")
        self.on_backward_strand = reverse
        intervals = ((1, 3), (20, 22), (120, 122), (140, 142))
        if reverse:
            intervals = tuple(
                (201 - end, 201 - start)
                for start, end in reversed(intervals))
            self.start_codon_positions = (198, 199, 200)
        else:
            self.start_codon_positions = (1, 2, 3)
        self.start = min(start for start, _ in intervals)
        self.end = max(end for _, end in intervals)
        self.start_codon_spliced_offsets = (0, 1, 2)
        self.exons = tuple(Exon(*interval) for interval in intervals)


def long_provider(reverse=False):
    chromosome = reverse_complement(LONG_FORWARD) if reverse else LONG_FORWARD

    def get_sequence(contig, start, end):
        assert contig == "x"
        sequence = chromosome[start - 1:end]
        return sequence + "A" * (end - start + 1 - len(sequence))

    return get_sequence


def provider(reverse=False):
    chromosome = reverse_complement(FORWARD) if reverse else FORWARD

    def get_sequence(contig, start, end):
        assert contig == "x"
        return chromosome[start - 1:end]

    return get_sequence


def point(start, ref, alt):
    return SimpleNamespace(
        contig="x", trimmed_base1_start=start,
        trimmed_base1_end=start + len(ref) - 1,
        trimmed_ref=ref, trimmed_alt=alt, is_structural=False)


def deletion(start, end):
    return SimpleNamespace(
        contig="x", start=start, end=end,
        affected_start=start, affected_end=end,
        sv_type="DEL", alt_assembly=None, is_structural=True)


def statuses_for(variant=None, reverse=False, with_sequence=True):
    transcript = Transcript(reverse)
    layout = GenomicLayout.from_transcript(
        transcript,
        flank=0,
        sequence_provider=provider(reverse) if with_sequence else None)
    if variant is not None:
        layout = layout.apply_point_variant(
            variant, validate_reference=with_sequence)
    runs = build_exon_runs(transcript, layout)
    return disrupted_splice_sites(transcript, layout, runs)


@pytest.mark.parametrize("reverse", [False, True])
def test_reference_splice_sites_are_intact(reverse):
    assert statuses_for(reverse=reverse) == ()


def test_donor_canonical_and_weak_windows_are_distinct():
    strong = statuses_for(point(11, "G", "A"))
    weak_intronic = statuses_for(point(14, "G", "A"))
    weak_exonic = statuses_for(point(10, "A", "C"))

    assert [(s.key.side, s.key.exon_number, s.strength) for s in strong] == [
        ("donor", 2, "lost")]
    donor = next(
        status for status in weak_intronic
        if status.key.side == "donor" and status.key.exon_number == 2)
    assert donor.strength == "weak_intronic"
    assert weak_exonic[0].strength == "weak_exonic"


def test_reverse_strand_donor_uses_transcript_direction():
    # Forward position 11 maps to reverse-complement position 14.
    statuses = statuses_for(point(14, "C", "T"), reverse=True)

    assert [(s.key.side, s.key.exon_number, s.strength) for s in statuses] == [
        ("donor", 2, "lost")]


def test_tier_zero_detects_explicit_canonical_site_edit_without_fasta():
    statuses = statuses_for(
        point(11, "G", "A"), with_sequence=False)

    assert statuses[0].strength == "lost"


def test_exon_truncation_is_a_lost_splice_site():
    transcript = Transcript()
    layout = GenomicLayout.from_transcript(
        transcript, flank=0, sequence_provider=provider())
    layout = layout.apply_structural_variant(deletion(10, 14))

    statuses = disrupted_splice_sites(
        transcript, layout, build_exon_runs(transcript, layout))

    assert any(
        status.key.side == "donor"
        and status.key.exon_number == 2
        and status.strength == "lost"
        for status in statuses)


def test_strong_site_options_are_ordinal_not_probabilities():
    transcript = Transcript()
    layout = GenomicLayout.from_transcript(
        transcript, flank=0, sequence_provider=provider())
    layout = layout.apply_point_variant(point(11, "G", "A"))

    plans = enumerate_layout_splice_plans(transcript, layout)
    mechanisms = [
        tuple(option.mechanism for _, option in plan.choices)
        for plan in plans]

    assert mechanisms == [
        ("exon_skip",),
        ("intron_retention",),
        ("cryptic",),
        ("normal",),
    ]
    assert all(plan.probability is None for plan in plans)
    skipped = realize_exon_path(
        transcript, layout, kept_run_keys=plans[0].kept_runs)
    assert skipped.cdna_sequence == "ATGGCTTAG"


def test_two_sites_on_skipped_exon_do_not_form_impossible_cartesian_pair():
    transcript = Transcript()
    layout = GenomicLayout.from_transcript(
        transcript, flank=0, sequence_provider=provider())
    # Destroy exon 2's acceptor and donor.
    layout = layout.apply_variants((
        point(6, "A", "C"),
        point(11, "G", "A"),
    ))

    plans = enumerate_layout_splice_plans(transcript, layout)

    assert plans
    assert all(not (
        len(plan.choices) == 2
        and all(option.mechanism == "exon_skip"
                for _, option in plan.choices))
        for plan in plans)


def test_retained_intron_cannot_coexist_with_skipped_downstream_exon():
    transcript = Transcript()
    layout = GenomicLayout.from_transcript(
        transcript, flank=0, sequence_provider=provider())
    # Destroy the donors of exon 2 and exon 3. Retaining intron 2 requires
    # exon 3, so it cannot coexist with exon 3 skipping.
    layout = layout.apply_variants((
        point(11, "G", "A"),
        point(18, "G", "A"),
    ))

    plans = enumerate_layout_splice_plans(transcript, layout)

    assert all(not (
        dict((key.exon_number, option.mechanism) for key, option in plan.choices)
        == {2: "intron_retention", 3: "exon_skip"})
        for plan in plans)
