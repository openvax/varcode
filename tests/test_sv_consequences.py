"""Structural consequences use the spliced coding edit, not event geometry."""

import pytest
from pyensembl import cached_release

from varcode import StructuralVariant
from varcode.effects import (
    Deletion, FivePrimeUTR, FrameShift, Insertion, StartLoss,
    StructuralVariantEffect, ThreePrimeUTR, Unresolved, effect_priority,
)
from varcode.effects import structural


@pytest.fixture
def cftr():
    return cached_release(81).transcript_by_id("ENST00000003084")


@pytest.fixture
def local_model(monkeypatch):
    monkeypatch.setattr(structural, "_fusion_partners", lambda *args: ())
    monkeypatch.setattr(structural, "_enumerate_and_attach_cryptics", lambda *args: None)
    monkeypatch.setattr(structural, "_enumerate_and_attach_splice_outcomes", lambda *args: None)


def _annotate(tx, kind, start, end, **kwargs):
    variant = StructuralVariant(tx.contig, start, kind, end=end, genome=tx.genome, **kwargs)
    return variant.effect_on_transcript(tx)


@pytest.mark.parametrize("first,last,cls,protein_length", [
    (0, 2, StartLoss, None), (4, 4, Deletion, 1450), (4, 5, FrameShift, 171),
])
def test_issue_420_cftr_deletion_consequences(cftr, local_model, first, last, cls, protein_length):
    # Fully removed exons with breakpoints outside the splice windows.
    effect = _annotate(cftr, "DEL", cftr.exons[first].start - 20,
                       cftr.exons[last].end + 20)
    assert isinstance(effect, cls)
    protein = effect.mutant_transcript.mutant_protein_sequence
    assert (len(protein) if protein is not None else None) == protein_length
    assert effect.modifies_protein_sequence is True
    assert effect.modifies_coding_sequence is True
    assert effect.affected_exons == tuple(cftr.exons[first:last + 1])
    assert effect.mutant_transcript.evidence["sv_type"] == "DEL"
    if cls is not StartLoss:
        assert effect.mutant_protein_sequence == effect.mutant_transcript.mutant_protein_sequence


def test_partial_exon_keeps_conditional_consequence(cftr, local_model):
    exon = cftr.exons[4]
    effect = _annotate(cftr, "DEL", exon.start + 10, exon.start + 13)
    assert isinstance(effect, StructuralVariantEffect)
    candidate = effect.candidates[0]
    assert isinstance(candidate.effect, FrameShift)
    assert candidate.evidence["splice_ambiguous"] is True
    assert candidate.evidence["assumption"] == "reference_splicing"
    assert candidate.evidence["sv_type"] == "DEL"
    assert effect_priority(effect) == effect_priority(candidate.effect)
    assert effect.short_description == candidate.effect.short_description
    assert effect.mutant_transcript is candidate.effect.mutant_transcript


@pytest.mark.parametrize("tx_id", ["ENST00000003084", "ENST00000357654"])
@pytest.mark.parametrize("region,cls", [("5utr", FivePrimeUTR), ("3utr", ThreePrimeUTR)])
def test_utr_deletion_preserves_mapped_protein(tx_id, region, cls, local_model):
    tx = cached_release(81).transcript_by_id(tx_id)
    exon = tx.exons[0] if region == "5utr" else tx.exons[-1]
    low = (region == "5utr") == (tx.strand == "+")
    start, end = (exon.start, exon.start + 9) if low else (exon.end - 9, exon.end)
    result = _annotate(tx, "DEL", start, end)
    effect = result.most_likely_effect
    assert isinstance(effect, cls)
    assert effect.mutant_transcript.mutant_protein_sequence == tx.protein_sequence
    assert result.modifies_protein_sequence is False
    assert result.modifies_coding_sequence is False


@pytest.mark.parametrize("tx_id", ["ENST00000003084", "ENST00000357654"])
def test_start_loss_maps_all_three_bases_on_both_strands(tx_id, local_model):
    tx = cached_release(81).transcript_by_id(tx_id)
    positions = tx.start_codon_positions
    result = _annotate(tx, "DEL", min(positions), min(positions))
    assert isinstance(result.most_likely_effect, StartLoss)
    assert result.mutant_transcript.mutant_protein_sequence is None


@pytest.mark.parametrize("last,cls", [(4, Insertion), (5, FrameShift)])
def test_tandem_duplication_translates_exonic_body(cftr, local_model, last, cls):
    effect = _annotate(cftr, "DUP", cftr.exons[4].start - 20, cftr.exons[last].end + 20)
    assert isinstance(effect, cls)
    assert effect.mutant_protein_sequence == effect.mutant_transcript.mutant_protein_sequence
    assert effect.mutant_transcript.evidence["structure_assumption"] == "tandem_duplication"


@pytest.mark.parametrize("kind", ["INV", "INS", "CNV"])
def test_unknown_alleles_do_not_invent_a_protein(cftr, local_model, kind):
    effect = _annotate(cftr, kind, cftr.exons[4].start - 20, cftr.exons[4].end + 20)
    assert isinstance(effect.most_likely_effect, Unresolved)
    assert effect.modifies_protein_sequence is None
    assert effect.modifies_coding_sequence is None
    assert effect.mutant_protein_sequence is None
    if kind != "INV":
        assert effect.mutant_transcript is None


def test_unmapped_assembly_does_not_imply_start_loss(cftr, local_model):
    effect = _annotate(cftr, "DEL", cftr.exons[4].start - 20, cftr.exons[4].end + 20,
                       alt_assembly="ATGGCTTAA")
    assert isinstance(effect.most_likely_effect, Unresolved)
    assert effect.mutant_transcript.cdna_sequence == "ATGGCTTAA"
    assert effect.modifies_protein_sequence is None


@pytest.mark.parametrize("tx_id", ["ENST00000523403", "ENST00000264498"])
def test_non_atg_initiator_matches_consequence_protein(tx_id, local_model):
    tx = cached_release(81).transcript_by_id(tx_id)
    exon = tx.exons[-1]
    position = exon.start + 20 if tx.strand == "+" else exon.end - 20
    result = _annotate(tx, "DEL", position, position + 2)
    effect = result.most_likely_effect
    assert isinstance(effect, Deletion)
    assert effect.mutant_protein_sequence == effect.mutant_transcript.mutant_protein_sequence
    assert effect.mutant_protein_sequence.startswith("M")


@pytest.mark.parametrize("insert", ["A", "AA"])
def test_local_paired_deletion_does_not_drop_inserted_bases(cftr, local_model, tmp_path, insert):
    from .test_osteosarc_fusions import _load_esvee_records
    from varcode.transforms import pair_breakends

    low, high = cftr.exons[4].start + 10, cftr.exons[4].start + 32
    variants = _load_esvee_records(tmp_path, [
        ("chr7", str(low), "a", "C", "C" + insert + "[chr7:%d[" % high,
         "60", "PASS", "SVTYPE=DEL;MATEID=b"),
        ("chr7", str(high), "b", "A", "]chr7:%d]" % low + insert + "A",
         "60", "PASS", "SVTYPE=DEL;MATEID=a"),
    ])
    variant, = pair_breakends(variants)
    effect = variant.effect_on_transcript(cftr)
    assert isinstance(effect.most_likely_effect, Unresolved)
    assert effect.mutant_transcript.cdna_sequence is None
    assert effect.mutant_transcript.mutant_protein_sequence is None
    assert effect.mutant_transcript.evidence["sequence_status"] == "unresolved_junction_insertion"


def test_incomplete_coding_transcript_retains_exon_loss(local_model):
    from varcode.effects import ExonLoss

    tx = cached_release(81).transcript_by_id("ENST00000446805")
    exon = tx.exons[0]
    assert tx.is_protein_coding and not tx.complete
    effect = _annotate(tx, "DEL", exon.start - 20, exon.end + 20)
    assert isinstance(effect, ExonLoss)
    assert effect.exons == (exon,)
    assert effect.mutant_transcript.mutant_protein_sequence is None
