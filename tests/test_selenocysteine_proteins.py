"""Stored mutant proteins read annotated selenocysteine UGA codons (#470).

Ensembl marks Sec as U in the reference protein. Translators read those codons
as Sec unless the mutant keeps no selenoprotein 3' UTR, where the SECIS element
lies. Also covers the CDS start shifted by upstream 5' UTR edits (#471).
"""

import pytest
from pyensembl import cached_release

from varcode import StructuralVariant, TranscriptEdit, Variant
from varcode.effects import GeneFusion
from varcode.effects import structural
from varcode.effects.codon_tables import translate_sequence
from varcode.effects.selenocysteine import edited_selenocysteine, reference_selenocysteine
from varcode.effects.sequence_change import structural_sequence_changes
from varcode.mutant_transcript import apply_variant_to_transcript, apply_variants_to_transcript
from varcode.splice_outcomes import _build_exon_skip_mutant_transcript

# GPX1 (-), GPX4 (+), SEPP1 (-, ten Sec), TXNRD1 (+, Sec is the penultimate residue).
SELENOPROTEINS = ["ENST00000419783", "ENST00000354171", "ENST00000514985", "ENST00000525566"]
COMPLEMENT = dict(zip("ACGT", "TGCA"))


def _transcript(transcript_id):
    return cached_release(81).transcript_by_id(transcript_id)


def _genomic(transcript, offset):
    for exon in transcript.exons:
        length = exon.end - exon.start + 1
        if offset < length:
            return exon.end - offset if transcript.strand == "-" else exon.start + offset
        offset -= length


def _variant(transcript, offset, ref, alt):
    """A variant given in transcript orientation at a cDNA offset."""
    if transcript.strand == "-":
        ref = "".join(COMPLEMENT[b] for b in reversed(ref))
        alt = "".join(COMPLEMENT[b] for b in reversed(alt))
        position = _genomic(transcript, offset + len(ref) - 1)
    else:
        position = _genomic(transcript, offset)
    return Variant(transcript.contig, position, ref, alt, transcript.genome)


def _snv(transcript, residue, base, alt):
    offset = min(transcript.start_codon_spliced_offsets) + 3 * residue + base
    return _variant(transcript, offset, transcript.sequence[offset], alt)


def _cases(transcript):
    """Point variants at, before, and after the first Sec codon."""
    u = transcript.protein_sequence.index("U")
    start = min(transcript.start_codon_spliced_offsets)

    def codon(i):
        return transcript.sequence[start + 3 * i:start + 3 * i + 3]

    def other(base):
        return "A" if base != "A" else "C"

    silent = next(i for i in range(u + 1, len(transcript.protein_sequence))
                  if codon(i)[:2] in ("CT", "GT", "TC", "CC", "AC", "GC", "CG", "GG"))
    return {
        "sec_to_trp": _snv(transcript, u, 2, "G"),
        "before_sec": _snv(transcript, 10, 1, other(codon(10)[1])),
        "after_sec": _snv(transcript, u + 1, 1, other(codon(u + 1)[1])),
        "silent_after_sec": _snv(transcript, silent, 2, other(codon(silent)[2])),
    }


def test_translate_sequence_reads_listed_tga_as_selenocysteine():
    assert translate_sequence("ATGTGAGGCTAA") == "M"
    assert translate_sequence("ATGTGAGGCTAA", selenocysteine={3}) == "MUG"
    # Only a TGA at a listed codon offset is read as Sec.
    assert translate_sequence("ATGTAAGGC", selenocysteine={3}) == "M"


@pytest.mark.parametrize("transcript_id", SELENOPROTEINS)
@pytest.mark.parametrize("case", ["sec_to_trp", "before_sec", "after_sec", "silent_after_sec"])
@pytest.mark.parametrize("annotator", ["protein_diff", "transcript_model"])
def test_annotators_agree_with_fast_on_selenoproteins(transcript_id, case, annotator):
    t = _transcript(transcript_id)
    variant = _cases(t)[case]
    expected = variant.effect_on_transcript(t, annotator="fast")
    effect = variant.effect_on_transcript(t, annotator=annotator)
    assert (type(effect), effect.short_description) == (type(expected), expected.short_description)


def test_gpx1_calls_are_pinned():
    t = _transcript("ENST00000419783")
    for annotator in ("fast", "protein_diff", "transcript_model"):
        assert _snv(t, 48, 2, "G").effect_on_transcript(
            t, annotator=annotator).short_description == "p.U49W"
        assert _snv(t, 53, 0, "T").effect_on_transcript(
            t, annotator=annotator).short_description == "p.R54W"


@pytest.mark.parametrize("transcript_id", SELENOPROTEINS)
def test_stored_protein_continues_past_selenocysteine(transcript_id):
    t = _transcript(transcript_id)
    mt = apply_variant_to_transcript(_cases(t)["after_sec"], t)
    protein = mt.mutant_protein_sequence
    assert len(protein) == len(t.protein_sequence)
    assert [i for i, aa in enumerate(protein) if aa == "U"] == [
        i for i, aa in enumerate(t.protein_sequence) if aa == "U"]
    assert sum(a != b for a, b in zip(protein, t.protein_sequence)) == 1


def test_edited_selenocysteine_codon_is_translated_literally():
    t = _transcript("ENST00000354171")
    sec, = reference_selenocysteine(t)
    assert edited_selenocysteine(t, ()) == {sec}
    assert edited_selenocysteine(t, (TranscriptEdit(sec + 2, sec + 3, "G"),)) == set()
    assert edited_selenocysteine(t, (TranscriptEdit(sec + 1, sec + 1, "A"),)) == set()
    # A 5' UTR deletion shifts the codon; an insertion at its first base precedes it.
    assert edited_selenocysteine(t, (TranscriptEdit(10, 12, ""),)) == {sec - 2}
    assert edited_selenocysteine(t, (TranscriptEdit(sec, sec, "AAA"),)) == {sec + 3}


def test_deleting_the_whole_three_prime_utr_leaves_uga_a_stop():
    t = _transcript("ENST00000354171")
    utr = max(t.stop_codon_spliced_offsets) + 1
    assert edited_selenocysteine(t, (TranscriptEdit(utr, len(t.sequence), ""),)) == set()
    assert edited_selenocysteine(t, (TranscriptEdit(utr, len(t.sequence) - 1, ""),))


@pytest.mark.parametrize("transcript_id", ["ENST00000003084", "ENST00000354171"])
def test_joint_translation_follows_a_five_prime_utr_indel(transcript_id):
    # #471: a 5' UTR deletion shifts the CDS start in the mutant cDNA.
    t = _transcript(transcript_id)
    start = min(t.start_codon_spliced_offsets)
    utr = start // 2
    deletion = _variant(t, utr, t.sequence[utr:utr + 2], t.sequence[utr])
    missense = _snv(t, 1, 1, "A" if t.sequence[start + 4] != "A" else "C")
    alone = apply_variant_to_transcript(missense, t).mutant_protein_sequence
    joint = apply_variants_to_transcript([deletion, missense], t)
    assert joint.mutant_protein_sequence == alone
    assert len(alone) == len(t.protein_sequence)


@pytest.mark.parametrize("transcript_id,index", [
    ("ENST00000514985", 2),   # SEPP1 (-): between its first and second Sec
    ("ENST00000525566", 15),  # TXNRD1 (+): just upstream of its Sec
])
def test_exon_skipping_keeps_selenocysteine(transcript_id, index):
    t = _transcript(transcript_id)
    exon = t.exons[index]
    variant = StructuralVariant(t.contig, exon.start, "DEL", end=exon.end, genome=t.genome)
    mt = _build_exon_skip_mutant_transcript(variant, t, exon)
    skipped = (exon.end - exon.start + 1) // 3
    assert mt.mutant_protein_sequence.count("U") == t.protein_sequence.count("U")
    assert len(mt.mutant_protein_sequence) == len(t.protein_sequence) - skipped


def _fusion(five_prime, five_prime_position, three_prime, three_prime_position):
    model = structural._build_fusion_mutant_transcript(
        five_prime, five_prime, five_prime_position, three_prime, three_prime_position)
    variant = StructuralVariant(five_prime.contig, five_prime_position, "BND",
                                genome=five_prime.genome)
    return model, GeneFusion(variant, five_prime, three_prime, mutant_transcript=model)


def test_fused_selenoprotein_three_prime_partner_keeps_selenocysteine():
    # The 3' partner brings its own 3' UTR, and so its SECIS.
    cftr, gpx4 = _transcript("ENST00000003084"), _transcript("ENST00000354171")
    model, _ = _fusion(cftr, cftr.exons[4].end, gpx4, gpx4.exons[1].start)
    protein = model.mutant_protein_sequence
    assert "U" in protein
    assert protein.endswith(gpx4.protein_sequence[-20:])


def test_fusion_in_five_prime_partner_utr_keeps_reference_protein_but_flags_unknown():
    gpx4, cftr = _transcript("ENST00000354171"), _transcript("ENST00000003084")
    model, effect = _fusion(gpx4, gpx4.end - 10, cftr, cftr.end - 10)
    assert model.mutant_protein_sequence == gpx4.protein_sequence
    # Whether the kept part of the 3' UTR still holds the SECIS is unknown.
    assert structural_sequence_changes(effect) == (False, None)


def test_fusion_that_drops_the_selenoprotein_utr_stops_at_selenocysteine():
    gpx4, cftr = _transcript("ENST00000354171"), _transcript("ENST00000003084")
    sec, = reference_selenocysteine(gpx4)
    breakpoint = _genomic(gpx4, sec + 30)
    model, effect = _fusion(gpx4, breakpoint, cftr, cftr.end - 10)
    assert model.mutant_protein_sequence == gpx4.protein_sequence[:gpx4.protein_sequence.index("U")]
    assert effect.modifies_protein_sequence is True

