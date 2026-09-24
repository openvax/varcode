"""SV change flags read annotated selenocysteine UGA codons as the reference does.

Ensembl marks Sec as U in the reference protein. UGA encodes Sec only with a
SECIS element in the same mRNA's 3' UTR, which isn't annotated, so decoding
is known only where a model keeps the transcript through its 3' end.
"""

from dataclasses import replace

import pytest
from pyensembl import cached_release

from varcode import (
    MutantTranscript, ReferenceSegment, StructuralVariant, make_fusion_outcome,
)
from varcode.effects import (
    EffectCollection, GeneFusion, LargeDeletion, StructuralVariantEffect,
)
from varcode.effects import structural

# GPX1 (-), GPX4 (+), SEPP1 (-, ten Sec), TXNRD1 (+, Sec is the penultimate residue).
SELENOPROTEINS = ["ENST00000419783", "ENST00000354171", "ENST00000514985", "ENST00000525566"]


@pytest.fixture
def isolate_span(monkeypatch):
    monkeypatch.setattr(structural, "_fusion_partners", lambda *args: ())
    monkeypatch.setattr(structural, "_enumerate_and_attach_cryptics", lambda *args: None)
    monkeypatch.setattr(structural, "_enumerate_and_attach_splice_outcomes", lambda *args: None)


def _transcript(transcript_id):
    return cached_release(81).transcript_by_id(transcript_id)


def _effect(transcript, cdna, segments=None, protein=None, evidence=None):
    variant = StructuralVariant(transcript.contig, transcript.start, "BND",
                                genome=transcript.genome)
    model = MutantTranscript(
        reference_transcript=transcript, reference_segments=segments,
        cdna_sequence=cdna, mutant_protein_sequence=protein, evidence=evidence)
    return StructuralVariantEffect(variant, transcript, mutant_transcript=model)


def _sec_offset(transcript, nth=0):
    start = min(transcript.start_codon_spliced_offsets)
    return start + 3 * [i for i, aa in enumerate(transcript.protein_sequence) if aa == "U"][nth]


def _with_codon(sequence, pos, codon):
    return sequence[:pos] + codon + sequence[pos + 3:]


@pytest.mark.parametrize("transcript_id", SELENOPROTEINS)
def test_unchanged_selenoprotein_is_unchanged(transcript_id):
    t = _transcript(transcript_id)
    effect = _effect(t, t.sequence, (ReferenceSegment(t, 0, len(t.sequence)),))
    assert effect.modifies_coding_sequence is False
    assert effect.modifies_protein_sequence is False
    assert len(EffectCollection([effect]).drop_silent_and_noncoding()) == 0


@pytest.mark.parametrize("transcript_id", SELENOPROTEINS)
@pytest.mark.parametrize("change,coding,protein", [
    ("sec_to_trp", True, True), ("sec_to_stop", True, True),
    ("missense_before", True, True), ("synonymous_before", True, False),
])
def test_changes_at_and_around_selenocysteine(transcript_id, change, coding, protein):
    t = _transcript(transcript_id)
    sec = _sec_offset(t)
    start = min(t.start_codon_spliced_offsets)
    cdna = t.sequence
    if change == "sec_to_trp":
        cdna = _with_codon(cdna, sec, "TGG")
    elif change == "sec_to_stop":
        cdna = _with_codon(cdna, sec, "TAA")
    elif change == "missense_before":
        cdna = _with_codon(cdna, start + 3, "TGG" if cdna[start + 3:start + 6] != "TGG" else "GCT")
    else:
        # Leucine CTN codons are four-fold degenerate at the third base.
        pos = next(i for i in range(start + 3, sec, 3) if cdna[i:i + 2] == "CT")
        cdna = _with_codon(cdna, pos, "CT" + ("A" if cdna[pos + 2] != "A" else "G"))
    effect = _effect(t, cdna, (ReferenceSegment(t, 0, len(t.sequence)),))
    assert effect.modifies_coding_sequence is coding
    assert effect.modifies_protein_sequence is protein


def test_change_after_a_decoded_selenocysteine_is_found():
    t = _transcript("ENST00000354171")
    pos = _sec_offset(t) + 6
    codon = "TGG" if t.sequence[pos:pos + 3] != "TGG" else "GCT"
    effect = _effect(t, _with_codon(t.sequence, pos, codon),
                     (ReferenceSegment(t, 0, len(t.sequence)),))
    assert effect.modifies_coding_sequence is True
    assert effect.modifies_protein_sequence is True


def test_split_segments_still_retain_the_three_prime_end():
    t = _transcript("ENST00000354171")
    cut = _sec_offset(t) + 1
    segments = (ReferenceSegment(t, 0, cut), ReferenceSegment(t, cut, len(t.sequence)))
    effect = _effect(t, t.sequence, segments)
    assert effect.modifies_protein_sequence is False


@pytest.mark.parametrize("transcript_id", SELENOPROTEINS)
@pytest.mark.parametrize("region,coding,protein", [
    ("5utr", False, False), ("3utr", False, None)])
def test_utr_deletions_respect_secis_uncertainty(transcript_id, region, coding, protein,
                                                 isolate_span):
    # A 3' UTR deletion may remove the SECIS element, so Sec decoding is unknown.
    t = _transcript(transcript_id)
    exon = t.exons[0] if region == "5utr" else t.exons[-1]
    at_low_end = (region == "5utr") == (t.strand == "+")
    start, end = (exon.start, exon.start + 9) if at_low_end else (exon.end - 9, exon.end)
    variant = StructuralVariant(t.contig, start, "DEL", end=end, genome=t.genome)
    effect = variant.effect_on_transcript(t)
    assert effect.affected_exons
    assert effect.modifies_coding_sequence is coding
    assert effect.modifies_protein_sequence is protein
    assert len(EffectCollection([effect]).drop_silent_and_noncoding()) == int(protein is None)


@pytest.mark.parametrize("breakpoint,coding,protein", [
    ("three_prime_utr", False, None), ("after_sec", True, True)])
def test_fusion_of_a_selenoprotein_five_prime_partner(breakpoint, coding, protein):
    t = _transcript("ENST00000354171")
    partner = _transcript("ENST00000003084")
    if breakpoint == "three_prime_utr":
        position = t.end - 10
    else:
        # GPX4 is on the forward strand: walk exons to a CDS base past the Sec codon.
        offset = _sec_offset(t) + 30
        for exon in t.exons:
            if offset < exon.end - exon.start + 1:
                position = exon.start + offset
                break
            offset -= exon.end - exon.start + 1
    model = structural._build_fusion_mutant_transcript(t, t, position, partner, partner.end - 10)
    variant = StructuralVariant(t.contig, position, "BND", genome=t.genome)
    effect = GeneFusion(variant, t, partner, mutant_transcript=model)
    assert effect.modifies_coding_sequence is coding
    assert effect.modifies_protein_sequence is protein


def test_three_prime_utr_duplication_leaves_cds_unchanged(isolate_span):
    t = _transcript("ENST00000354171")
    start = t.exons[-1].end - 20
    effect = StructuralVariant(t.contig, start, "DUP", end=start + 5,
                               genome=t.genome).effect_on_transcript(t)
    assert effect.modifies_coding_sequence is False
    assert effect.modifies_protein_sequence is None


def test_losing_the_whole_three_prime_utr_truncates_at_selenocysteine(isolate_span):
    # No 3' UTR, so no SECIS: the CDS is intact, but UGA now terminates.
    t = _transcript("ENST00000354171")
    variant = StructuralVariant(t.contig, max(t.stop_codon_positions) + 1, "DEL",
                                end=t.end, genome=t.genome)
    effect = variant.effect_on_transcript(t)
    assert effect.modifies_coding_sequence is False
    assert effect.modifies_protein_sequence is True


@pytest.mark.parametrize("transcript_id", SELENOPROTEINS)
def test_unmapped_import_of_selenoprotein_is_unknown_not_changed(transcript_id):
    # Without coordinates, only the reference Sec index marks a possible Sec UGA.
    t = _transcript(transcript_id)
    variant = StructuralVariant(t.contig, t.start, "BND", genome=t.genome)
    candidate = make_fusion_outcome(variant, t, sequence=t.sequence,
                                    cds_start=min(t.start_codon_spliced_offsets),
                                    transcript_model_id="observed-model")
    assert candidate.effect.modifies_coding_sequence is False
    assert candidate.effect.modifies_protein_sequence is None
    start = min(t.start_codon_spliced_offsets)
    missense = _with_codon(t.sequence, start + 3,
                           "TGG" if t.sequence[start + 3:start + 6] != "TGG" else "GCT")
    changed = make_fusion_outcome(variant, t, sequence=missense, cds_start=start,
                                  transcript_model_id="observed-model")
    assert changed.effect.modifies_protein_sequence is True


def test_peptide_ending_at_selenocysteine_is_not_a_truncation():
    # Exacto translates UGA as a stop, so its peptide ends at the Sec residue.
    t = _transcript("ENST00000354171")
    start, sec = min(t.start_codon_spliced_offsets), _sec_offset(t)
    effect = _effect(t, t.sequence, protein=t.protein_sequence[:(sec - start) // 3],
                     evidence={"protein_completeness": "start_to_stop",
                               "cds_start": start, "cds_end": sec + 3})
    assert effect.modifies_coding_sequence is False
    assert effect.modifies_protein_sequence is None


@pytest.mark.parametrize("reaches_three_prime_end", [True, False])
def test_partial_observation_reads_through_only_decoded_selenocysteine(reaches_three_prime_end):
    t = _transcript("ENST00000354171")
    sec = _sec_offset(t)
    left = sec - 30
    right = len(t.sequence) if reaches_three_prime_end else max(t.stop_codon_spliced_offsets) + 1
    cdna = t.sequence[left:right]
    pos = sec - left + 6
    cdna = _with_codon(cdna, pos, "TGG" if cdna[pos:pos + 3] != "TGG" else "GCT")
    effect = _effect(t, cdna, (ReferenceSegment(t, left, right),), evidence={
        "protein_completeness": "partial_start", "cds_start": 0,
        "cds_end": max(t.stop_codon_spliced_offsets) + 1 - left})
    assert effect.modifies_coding_sequence is (True if reaches_three_prime_end else None)
    assert effect.modifies_protein_sequence is (True if reaches_three_prime_end else None)
    sec_to_stop = _with_codon(effect.mutant_transcript.cdna_sequence, sec - left, "TAA")
    effect.mutant_transcript = replace(effect.mutant_transcript, cdna_sequence=sec_to_stop)
    assert effect.modifies_coding_sequence is True
    assert effect.modifies_protein_sequence is True


# TXNRD3 (CTG start, Sec) and FGF2 (CTG start): Ensembl writes the initiator as L.
@pytest.mark.parametrize("transcript_id", ["ENST00000523403", "ENST00000264498"])
def test_non_atg_initiator_is_not_a_protein_change(transcript_id):
    t = _transcript(transcript_id)
    assert t.coding_sequence[:3] == "CTG" and t.protein_sequence[0] == "L"
    effect = _effect(t, t.sequence, (ReferenceSegment(t, 0, len(t.sequence)),))
    assert effect.modifies_coding_sequence is False
    assert effect.modifies_protein_sequence is False
    effect.mutant_transcript = replace(
        effect.mutant_transcript, mutant_protein_sequence="M" + t.protein_sequence[1:])
    assert effect.modifies_protein_sequence is False
