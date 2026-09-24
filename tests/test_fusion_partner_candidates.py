"""Fusion proteins must not depend on discarding all but the first isoform."""

from types import SimpleNamespace

import pytest
from pyensembl import cached_release

from varcode import EffectCandidate, StructuralVariant
from varcode.effects import GeneFusion, StructuralVariantEffect, TranslocationToIntergenic
from varcode.effects import structural


@pytest.fixture
def fusion_case(monkeypatch):
    genome = cached_release(95)
    transcript = genome.transcript_by_id("ENST00000538197")  # CPEB2, forward
    variant = StructuralVariant(
        "4", 15012987, "BND", alt="N[4:2664478[",
        mate_contig="4", mate_start=2664478, genome=genome)
    # These tests isolate junction/isoform enumeration. Dedicated existing
    # tests cover the splice and cryptic layers added after it.
    monkeypatch.setattr(structural, "_enumerate_and_attach_cryptics", lambda *args: None)
    monkeypatch.setattr(structural, "_enumerate_and_attach_splice_outcomes", lambda *args: None)
    return genome, transcript, variant


def _fusions(effect):
    return [candidate.effect for candidate in effect.candidates
            if isinstance(candidate.effect, GeneFusion)]


def test_audited_cpeb2_fam193a_proteins_are_both_candidates(fusion_case):
    genome, transcript, variant = fusion_case
    effect = variant.effect_on_transcript(transcript)
    candidates = _fusions(effect)
    expected = [t.id for t in genome.transcripts_at_locus("4", 2664478, 2664478)
                if t.is_protein_coding and t.gene_name == "FAM193A" and t.strand == "+"]
    assert [c.partner_transcript.id for c in candidates] == expected
    assert len(candidates) == 7
    assert candidates[0] is effect  # existing primary remains compatible
    proteins = {c.partner_transcript.id: c.mutant_protein_sequence for c in candidates}
    assert len(proteins["ENST00000324666"]) == 1541
    assert len(proteins["ENST00000637812"]) == 1500
    assert proteins["ENST00000324666"] != proteins["ENST00000637812"]
    partner = genome.transcript_by_id("ENST00000637812")
    assert proteins[partner.id] == transcript.protein_sequence[:678] + partner.protein_sequence[693:]
    for candidate in candidates:
        assert candidate.transcript is transcript
        assert candidate.five_prime_transcript is transcript
        assert candidate.three_prime_transcript is candidate.partner_transcript
        model = candidate.mutant_transcript
        assert model.cdna_sequence == "".join(
            segment.source.sequence[segment.start:segment.end]
            for segment in model.reference_segments)


def test_same_transcript_pair_has_same_protein_from_either_end(fusion_case):
    genome, transcript, variant = fusion_case
    partner = genome.transcript_by_id("ENST00000637812")
    mate = StructuralVariant(
        "4", 2664478, "BND", alt="]4:15012987]N",
        mate_contig="4", mate_start=15012987, genome=genome)
    from_five = next(c for c in _fusions(variant.effect_on_transcript(transcript))
                     if c.partner_transcript.id == partner.id)
    from_three = next(c for c in _fusions(mate.effect_on_transcript(partner))
                      if c.partner_transcript.id == transcript.id)
    assert from_five.mutant_protein_sequence == from_three.mutant_protein_sequence
    assert from_three.five_prime_transcript.id == transcript.id
    assert from_three.three_prime_transcript.id == partner.id


class _Isoform:
    """Reuse a reference transcript while varying only ID or available sequence."""
    def __init__(self, transcript, name, sequence):
        self._transcript = transcript
        self.id = name
        self.sequence = sequence

    def __getattr__(self, name):
        return getattr(self._transcript, name)


def test_large_isoform_set_is_not_capped_or_merged_by_protein(fusion_case, monkeypatch):
    genome, transcript, variant = fusion_case
    partner = genome.transcript_by_id("ENST00000637812")
    partners = [_Isoform(partner, "isoform_%03d" % i, partner.sequence) for i in range(257)]
    monkeypatch.setattr(structural, "_coding_transcripts_at", lambda *args: partners)
    candidates = _fusions(variant.effect_on_transcript(transcript))
    assert [c.partner_transcript.id for c in candidates] == [p.id for p in partners]
    assert len({c.mutant_protein_sequence for c in candidates}) == 1


def test_duplicate_junction_and_partner_records_are_not_repeated(fusion_case, monkeypatch):
    genome, transcript, variant = fusion_case
    partner = genome.transcript_by_id("ENST00000637812")
    monkeypatch.setattr(structural, "_coding_transcripts_at", lambda *args: [partner, partner])
    duplicate = SimpleNamespace(junctions=variant.junctions * 2)
    effect = structural._fusion_across_junction(duplicate, transcript, None)
    assert len(_fusions(effect)) == 1


@pytest.mark.parametrize("sequence", [None, ""])
def test_missing_partner_sequence_keeps_an_unresolved_candidate(fusion_case, monkeypatch, sequence):
    genome, transcript, variant = fusion_case
    partner = genome.transcript_by_id("ENST00000637812")
    missing = _Isoform(partner, "missing_sequence", sequence)
    monkeypatch.setattr(structural, "_coding_transcripts_at", lambda *args: [partner, missing])
    candidates = _fusions(variant.effect_on_transcript(transcript))
    assert len(candidates) == 2
    assert candidates[0].mutant_protein_sequence is not None
    assert candidates[1].partner_transcript is missing
    assert candidates[1].mutant_transcript is None
    assert candidates[1].mutant_protein_sequence is None


def test_supplied_assembly_and_external_candidates_are_preserved(fusion_case):
    genome, transcript, variant = fusion_case
    variant.alt_assembly = "ATG" + "GCC" * 20 + "TAA"
    effect = variant.effect_on_transcript(transcript)
    candidates = _fusions(effect)
    assert len(candidates) == 7
    assert all(c.mutant_transcript.cdna_sequence == variant.alt_assembly for c in candidates)
    assert all(c.mutant_transcript.evidence == {"source": "alt_assembly"} for c in candidates)
    observed = EffectCandidate(candidates[-1], source="rna", evidence={"read_count": 4})
    effect._extra_candidates = (observed,)
    assert effect.candidates[-1] is observed
    assert len(effect.candidates) == 8


def test_opposite_strand_partner_isoforms_are_all_kept(monkeypatch):
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084")  # CFTR, +
    variant = StructuralVariant(
        "7", 117485000, "BND", alt="N]17:43120000]",
        mate_contig="17", mate_start=43120000, genome=genome)
    monkeypatch.setattr(structural, "_enumerate_and_attach_cryptics", lambda *args: None)
    effect = variant.effect_on_transcript(transcript)
    expected = {t.id for t in genome.transcripts_at_locus("17", 43120000, 43120000)
                if t.is_protein_coding and t.strand == "-"}
    assert {c.partner_transcript.id for c in _fusions(effect)} == expected
    assert len(expected) > 1


def test_same_role_join_still_has_no_fusion_candidates(fusion_case):
    genome, transcript, _ = fusion_case
    variant = StructuralVariant(
        "4", 15012987, "BND", alt="N]4:2664478]",
        mate_contig="4", mate_start=2664478, genome=genome)
    effect = variant.effect_on_transcript(transcript)
    assert isinstance(effect, TranslocationToIntergenic)
    assert not _fusions(effect)


def test_span_effect_is_retained_after_all_fusion_partners(fusion_case):
    genome, _, _ = fusion_case
    tmprss2 = max((t for t in genome.genes_by_name("TMPRSS2")[0].transcripts if t.is_protein_coding),
                 key=lambda t: len(t.sequence or ""))
    erg = max((t for t in genome.genes_by_name("ERG")[0].transcripts if t.is_protein_coding),
              key=lambda t: len(t.sequence or ""))
    def midpoint(transcript, exon):
        left, right = transcript.exons[exon - 1:exon + 1]
        return (min(left.end, right.end) + max(left.start, right.start)) // 2
    variant = StructuralVariant("21", midpoint(erg, 3), "DEL",
                                end=midpoint(tmprss2, 2) - 1, genome=genome)
    effect = variant.effect_on_transcript(tmprss2)
    assert len(_fusions(effect)) > 1
    assert isinstance(effect.candidates[-1].effect, StructuralVariantEffect)
    assert effect.candidates[-1].effect.affected_exons


def test_inversion_keeps_both_junction_directions_with_local_five_prime_first(fusion_case):
    genome, _, _ = fusion_case
    transcript = genome.transcript_by_id("ENST00000003084")  # CFTR, +
    partner_gene = genome.genes_by_name("CTTNBP2")[0]  # downstream, -
    far = (max(partner_gene.start, transcript.gene.end + 1) + partner_gene.end) // 2
    variant = StructuralVariant("7", 117485000, "INV", end=far, genome=genome)
    effect = variant.effect_on_transcript(transcript)
    candidates = _fusions(effect)
    assert effect.five_prime_transcript.id == transcript.id
    directions = {(c.five_prime_transcript.gene_name, c.three_prime_transcript.gene_name)
                  for c in candidates}
    assert directions == {("CFTR", "CTTNBP2"), ("CTTNBP2", "CFTR")}
    keys = [(c.five_prime_transcript.id, c.three_prime_transcript.id) for c in candidates]
    assert len(keys) == len(set(keys))
