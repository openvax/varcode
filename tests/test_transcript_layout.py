"""Tests for exon projection and realized transcript classification."""

from types import SimpleNamespace

import pytest

from varcode.effects.effect_classes import Deletion, Insertion, Silent
from varcode.genomic_layout import GenomicLayout
from varcode.nucleotides import reverse_complement
from varcode.transcript_layout import (
    build_exon_runs,
    classify_products,
    realize_exon_path,
)


FORWARD = "ATGGGCAAGGGCTGGTAG"


class Exon:
    def __init__(self, start, end):
        self.start = start
        self.end = end


class Transcript:
    contig = "x"
    start = 1
    end = 18
    is_protein_coding = True
    complete = True
    protein_sequence = "MQA"
    stop_codon_spliced_offsets = (9, 10, 11)
    id = "synthetic"

    def __init__(self, reverse=False):
        self.gene = SimpleNamespace(id="synthetic-gene", name="SYN")
        self.on_backward_strand = reverse
        intervals = ((1, 3), (6, 8), (11, 13), (16, 18))
        if reverse:
            intervals = tuple(reversed(intervals))
            self.start_codon_positions = (16, 17, 18)
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


def provider(reverse=False):
    chromosome = reverse_complement(FORWARD) if reverse else FORWARD

    def get_sequence(contig, start, end):
        assert contig == "x"
        return chromosome[start - 1:end]

    return get_sequence


def structural(kind, start, end):
    return SimpleNamespace(
        contig="x", start=start, end=end,
        affected_start=start, affected_end=end,
        sv_type=kind, alt_assembly=None, is_structural=True)


@pytest.mark.parametrize("reverse", [False, True])
def test_reference_layout_reproduces_transcript_on_both_strands(reverse):
    transcript = Transcript(reverse)
    layout = GenomicLayout.from_transcript(
        transcript, flank=0, sequence_provider=provider(reverse))

    product = realize_exon_path(transcript, layout)

    assert product.cdna_sequence == "ATGCAAGCTTAG"
    assert product.protein_sequence == "MQA"
    assert product.junction_signature == ()
    assert product.start_codon_present is True


@pytest.mark.parametrize("reverse", [False, True])
def test_whole_exon_deletion_is_classified_from_realized_protein(reverse):
    transcript = Transcript(reverse)
    layout = GenomicLayout.from_transcript(
        transcript, flank=0, sequence_provider=provider(reverse))
    baseline = realize_exon_path(transcript, layout)
    # Forward exon 2 is 6..8; the corresponding transcript-order exon 2 on
    # the reverse fixture is 11..13.
    span = (11, 13) if reverse else (6, 8)
    variant = structural("DEL", *span)
    mutant_layout = layout.apply_structural_variant(variant)

    mutant = realize_exon_path(transcript, mutant_layout)
    consequence = classify_products(variant, transcript, baseline, mutant)

    assert mutant.cdna_sequence == "ATGGCTTAG"
    assert mutant.protein_sequence == "MA"
    assert type(consequence) is Deletion
    assert consequence.aa_ref == "Q"
    assert len(mutant.junction_signature) == 1


@pytest.mark.parametrize("reverse", [False, True])
def test_whole_exon_duplication_is_an_inserted_amino_acid(reverse):
    transcript = Transcript(reverse)
    layout = GenomicLayout.from_transcript(
        transcript, flank=0, sequence_provider=provider(reverse))
    baseline = realize_exon_path(transcript, layout)
    span = (11, 13) if reverse else (6, 8)
    variant = structural("DUP", *span)

    mutant = realize_exon_path(
        transcript, layout.apply_structural_variant(variant))
    consequence = classify_products(variant, transcript, baseline, mutant)

    assert mutant.cdna_sequence == "ATGCAACAAGCTTAG"
    assert mutant.protein_sequence == "MQQA"
    assert type(consequence) is Insertion
    assert consequence.aa_alt == "Q"


@pytest.mark.parametrize("reverse", [False, True])
def test_inverted_exon_is_not_treated_as_sense_exon(reverse):
    transcript = Transcript(reverse)
    layout = GenomicLayout.from_transcript(
        transcript, flank=0, sequence_provider=provider(reverse))
    span = (11, 13) if reverse else (6, 8)

    runs = build_exon_runs(
        transcript, layout.apply_structural_variant(
            structural("INV", *span)))

    assert [run.exon_number for run in runs] == [1, 3, 4]


def test_exonic_variant_removed_by_splicing_is_silent_not_intronic():
    transcript = Transcript()
    layout = GenomicLayout.from_transcript(
        transcript, flank=0, sequence_provider=provider())
    baseline = realize_exon_path(transcript, layout)
    # Model exon 2 skipping while classifying a coding SNV inside exon 2.
    mutant = realize_exon_path(
        transcript, layout, kept_run_keys=((1, 1), (3, 1), (4, 1)))
    variant = SimpleNamespace(
        contig="x", trimmed_base1_start=7, trimmed_base1_end=7,
        is_structural=False)
    # Compare the somatic variant against a baseline taking the same splice
    # path: its allele is excluded, so protein is unchanged on this molecule.
    consequence = classify_products(variant, transcript, mutant, mutant)

    assert type(consequence) is Silent
    assert consequence.excluded_from_mrna is True
