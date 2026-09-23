"""VCF 4.5 junction inserts must survive fusion assembly (#485)."""

from types import SimpleNamespace

import pytest
from pyensembl import cached_release

from varcode import StructuralVariant
from varcode.effects.codon_tables import translate_sequence
from varcode.effects import GeneFusion
from varcode.effects.structural import _build_fusion_mutant_transcript
from varcode.nucleotides import reverse_complement
from varcode.sv_allele_parser import breakend_inserted_sequence

from .test_osteosarc_fusions import _load_esvee_records
from varcode.transforms import pair_breakends


@pytest.mark.parametrize("first_side,second_side", [
    ("left", "right"), ("left", "left"), ("right", "left"), ("right", "right")])
@pytest.mark.parametrize("insert", ["", "A", "AC", "ACG"])
def test_reciprocal_junction_sequence_in_all_orientations(first_side, second_side, insert):
    """The same traversal must give identical bases from either VCF record."""
    def record(contig, pos, mate, mate_pos, side, mate_side, bases):
        local = bases if side == "left" else reverse_complement(bases)
        bracket = "[" if mate_side == "right" else "]"
        locus = "%s%s:%d%s" % (bracket, mate, mate_pos, bracket)
        alt = "C" + local + locus if side == "left" else locus + local + "C"
        return StructuralVariant(contig, pos, "BND", ref="C", alt=alt,
                                 mate_contig=mate, mate_start=mate_pos)

    a = record("1", 100, "2", 200, first_side, second_side, insert)
    b = record("2", 200, "1", 100, second_side, first_side, reverse_complement(insert))
    junction = a.junctions[0]
    assert a.junction_inserted_sequence(junction) == insert
    assert b.junction_inserted_sequence(junction) == insert
    assert a.junction_inserted_sequence(junction[::-1]) == reverse_complement(insert)
    a.source_variants = (a, b)
    assert a.junction_inserted_sequence(junction) == insert


@pytest.mark.parametrize("ref,alt", [
    ("C", "A[2:200["), ("C", "]2:200]A"), ("CC", "CCA[2:200["),
    ("C", "C[2:200]"), ("C", "C[2:200[C"), ("C", "<BND>")])
def test_uninterpretable_replacement_is_not_an_empty_insert(ref, alt):
    assert breakend_inserted_sequence(ref, alt) is None


@pytest.mark.parametrize("insert,expected_length", [("", 1184), ("A", 1237), ("AA", 1188)])
@pytest.mark.parametrize("svtype", ["BND", "DEL"])
def test_znf236_galr1_from_raw_reciprocal_and_paired_records(tmp_path, insert, expected_length, svtype):
    """The real exonic Osteosarc junction, with zero/one/two inserted bases."""
    vc = _load_esvee_records(tmp_path, [
        ("chr18", "76920050", "a", "C", "C" + insert + "[chr18:77271549[",
         "60", "PASS", "SVTYPE=" + svtype + ";MATEID=b"),
        ("chr18", "77271549", "b", "A", "]chr18:76920050]" + insert + "A",
         "60", "PASS", "SVTYPE=" + svtype + ";MATEID=a"),
    ])
    genome = cached_release(95)
    five = genome.transcript_by_id("ENST00000253159")
    three = genome.transcript_by_id("ENST00000299727")
    cdna = five.sequence[:five.spliced_offset(76920050) + 1] + insert + three.sequence[three.spliced_offset(77271549):]
    coding = cdna[min(five.start_codon_spliced_offsets):]
    protein = translate_sequence(coding[:len(coding) // 3 * 3], to_stop=True)
    assert len(protein) == expected_length
    for variant in list(vc) + list(pair_breakends(vc)):
        for transcript in (five, three):
            fusion = next(c.effect for c in variant.effect_on_transcript(transcript).candidates
                          if isinstance(c.effect, GeneFusion)
                          and c.effect.five_prime_transcript.id == five.id
                          and c.effect.three_prime_transcript.id == three.id)
            model = fusion.mutant_transcript
            assert model.cdna_sequence == cdna
            assert model.mutant_protein_sequence == protein
            assert "".join(s.source.sequence[s.start:s.end]
                           for s in model.reference_segments) == cdna
            if insert:
                assert model.reference_segments[1].label == "junction_insertion"
                assert model.evidence["junction_insertion_status"] == "retained"


def _transcript(reverse=False):
    return SimpleNamespace(sequence="ATGAAACCCGGGTAA", complete=True,
                           start_codon_spliced_offsets=[0, 1, 2],
                           contig="1", strand="-" if reverse else "+",
                           on_backward_strand=reverse,
                           exons=[SimpleNamespace(start=100, end=114)])


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("insert,protein", [("", "MKPG"), ("A", "MKTRV"), ("AC", "MKTPG"), ("ACG", "MKTPG")])
def test_exonic_insert_translation_on_both_strands(reverse, insert, protein):
    five, three = _transcript(reverse), _transcript(reverse)
    five_pos, three_pos = (109, 108) if reverse else (105, 106)
    model = _build_fusion_mutant_transcript(five, five, five_pos, three, three_pos, insert)
    assert model.cdna_sequence == "ATGAAA" + insert + "CCCGGGTAA"
    assert model.mutant_protein_sequence == protein


@pytest.mark.parametrize("five_pos,three_pos,status", [
    (115, 99, "excluded_by_reference_splicing"),
    (105, 99, "unresolved_insertion_retention"),
    (115, 106, "unresolved_insertion_retention")])
def test_insertion_retention_requires_splicing_evidence(five_pos, three_pos, status):
    five, three = _transcript(), _transcript()
    model = _build_fusion_mutant_transcript(five, five, five_pos, three, three_pos, "AC")
    assert model.evidence["junction_inserted_sequence"] == "AC"
    if status.startswith("unresolved"):
        assert model.evidence["sequence_status"] == status
        assert model.cdna_sequence is None
        assert model.mutant_protein_sequence is None
    else:
        assert model.evidence["junction_insertion_status"] == status
        assert model.cdna_sequence == five.sequence + three.sequence


def test_conflicting_mates_do_not_invent_a_junction_sequence():
    a = StructuralVariant("1", 100, "BND", ref="C", alt="CA[2:200[", mate_contig="2", mate_start=200)
    b = StructuralVariant("2", 200, "BND", ref="C", alt="]1:100]TC", mate_contig="1", mate_start=100)
    a.source_variants = (a, b)
    assert a.junction_inserted_sequence(a.junctions[0]) is None
    five, three = _transcript(), _transcript()
    model = _build_fusion_mutant_transcript(five, five, 105, three, 106, None)
    assert model.cdna_sequence is None
    assert model.mutant_protein_sequence is None
    assert model.evidence["sequence_status"] == "unresolved_junction_insertion"
