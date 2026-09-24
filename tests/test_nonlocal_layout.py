"""Nonlocal rearrangements must not be classified from clipped alleles."""

import pytest
from pyensembl import cached_release

from varcode import StructuralVariant
from varcode.effects import Unresolved
from varcode.genomic_layout import GenomicLayout, NonlocalStructuralEdit
from varcode.transcript_model import predict_transcript_model_effect
from .test_genomic_layout import structural, provider


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("kind", ["DUP", "INV"])
@pytest.mark.parametrize("start,end", [(1, 6), (6, 12), (1, 12)])
def test_layout_rejects_crossing_or_enclosing_event(strand, kind, start, end):
    layout = GenomicLayout.from_interval("a", 3, 10, strand=strand, sequence_provider=provider)
    with pytest.raises(NonlocalStructuralEdit, match="complete rearrangement junctions"):
        layout.apply_structural_variant(structural(kind, start, end))
    assert layout.length == 8


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("kind", ["DUP", "INV"])
def test_fully_represented_and_disjoint_events(strand, kind):
    layout = GenomicLayout.from_interval("a", 3, 10, strand=strand, sequence_provider=provider)
    local = layout.apply_structural_variant(structural(kind, 3, 10))
    assert local.length == (16 if kind == "DUP" else 8)
    assert layout.apply_structural_variant(structural(kind, 20, 30)) is layout


@pytest.mark.parametrize("tx_id", ["ENST00000003084", "ENST00000357654"])
@pytest.mark.parametrize("kind", ["DUP", "INV"])
@pytest.mark.parametrize("side", ["left", "right", "both"])
def test_transcript_model_nonlocal_events_are_unresolved(tx_id, kind, side):
    tx = cached_release(81).transcript_by_id(tx_id)
    start = tx.start - 100 if side in ("left", "both") else tx.start + 50
    end = tx.end + 100 if side in ("right", "both") else tx.end - 50
    variant = StructuralVariant(tx.contig, start, kind, end=end, genome=tx.genome)
    for effect in (predict_transcript_model_effect((variant,), tx),
                   variant.effect_on_transcript(tx, annotator="transcript_model")):
        assert isinstance(effect, Unresolved)
        assert effect.mechanism == "nonlocal_structural_variant"
        assert effect.modifies_coding_sequence is None
        assert effect.modifies_protein_sequence is None
        assert effect.mutant_protein_sequence is None


@pytest.mark.parametrize("kind", ["DUP", "INV"])
@pytest.mark.parametrize("tx_id", ["ENST00000003084", "ENST00000357654"])
def test_nonlocal_supplied_assembly_is_preserved(kind, tx_id):
    tx = cached_release(81).transcript_by_id(tx_id)
    variant = StructuralVariant(tx.contig, tx.start - 100, kind,
                                end=tx.end + 100, genome=tx.genome,
                                alt_assembly="ATGAAATAG")
    effect = variant.effect_on_transcript(tx, annotator="transcript_model")
    assert effect.mutant_transcript.cdna_sequence == "ATGAAATAG"
    assert effect.modifies_protein_sequence is None


def test_each_isoform_checks_its_own_layout():
    genome = cached_release(81)
    gene = genome.transcript_by_id("ENST00000357654").gene
    variant = StructuralVariant(gene.contig, gene.start - 100, "DUP", end=gene.end + 100, genome=genome)
    effects = variant.effects(annotator="transcript_model")
    coding = [e for e in effects if e.transcript and e.transcript.complete
              and e.transcript.is_protein_coding and e.gene_name == "BRCA1"]
    assert len(coding) > 1
    assert all(isinstance(e, Unresolved) for e in coding)
    assert all(e.modifies_protein_sequence is None for e in coding)
