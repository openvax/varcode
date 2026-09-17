"""A BND fallback describes only its retained reference fragment (#447)."""

import pytest
from pyensembl import cached_release

from varcode import StructuralVariant, parse_symbolic_alt
from varcode.effects import TranslocationToIntergenic
from varcode.effects.structural import _build_translocation_mutant_transcript


@pytest.fixture(params=["ENST00000003084", "ENST00000357654"])
def transcript(request):
    return cached_release(81).transcript_by_id(request.param)


def expected_offsets(transcript, position, side):
    """Independent per-base oracle in transcript order, including the anchor."""
    positions = [p for exon in transcript.exons
                 for p in (range(exon.start, exon.end + 1) if transcript.strand == "+"
                           else range(exon.end, exon.start - 1, -1))]
    kept = [i for i, p in enumerate(positions)
            if (p <= position if side == "left" else p >= position)]
    return min(kept), max(kept) + 1


@pytest.mark.parametrize("alt,side", [
    ("N[22:15500000[", "left"), ("N]22:15500000]", "left"),
    ("[22:15500000[N", "right"), ("]22:15500000]N", "right"),
    ("ACGT.", "left"), (".ACGT", "right"),
])
@pytest.mark.parametrize("location", ["exonic", "exon_start", "exon_end", "intronic"])
def test_both_strands_and_local_sides(transcript, alt, side, location):
    # An interior exon also permits meaningful suffix/prefix comparisons.
    exon = transcript.exons[2]
    if location == "exonic":
        position = (exon.start + exon.end) // 2
    elif location == "exon_start":
        position = exon.start
    elif location == "exon_end":
        position = exon.end
    else:
        other = transcript.exons[3]
        position = ((exon.end + other.start) // 2 if transcript.strand == "+"
                    else (other.end + exon.start) // 2)
    sv = parse_symbolic_alt(transcript.contig, position, "N", alt,
                            genome=transcript.genome)
    effect = sv.effect_on_transcript(transcript)
    assert isinstance(effect, TranslocationToIntergenic)
    mt = effect.mutant_transcript
    assert mt is not None
    (segment,) = mt.reference_segments
    assert segment.source is transcript
    assert (segment.start, segment.end) == expected_offsets(transcript, position, side)
    assert segment.strand == "+"  # source is already oriented transcript cDNA
    five_prime = (side == "left") == (transcript.strand == "+")
    assert segment.label == ("translocation_5p" if five_prime else "translocation_3p")
    assert segment.length < len(transcript.sequence)
    assert mt.cdna_sequence is None
    assert mt.mutant_protein_sequence is None
    assert mt.evidence["sequence_status"] == "retained_reference_fragment"


def test_cftr_regression_keeps_51_bases_not_entire_transcript():
    tx = cached_release(81).transcript_by_id("ENST00000003084")
    sv = parse_symbolic_alt("7", tx.exons[0].start + 50, "N", "N]22:15500000]",
                            genome=tx.genome)
    mt = sv.effect_on_transcript(tx).mutant_transcript
    (segment,) = mt.reference_segments
    assert (segment.start, segment.end) == (0, 51)
    assert segment.source.sequence[segment.start:segment.end] == tx.sequence[:51]


@pytest.mark.parametrize("alt", [None, "<BND>", "N[?:?[", "N[22:15500000]", ".", "N.N"])
def test_unknown_local_orientation_is_unresolved(transcript, alt):
    sv = StructuralVariant(transcript.contig, transcript.exons[0].start + 5,
                           "BND", alt=alt, genome=transcript.genome)
    assert sv.effect_on_transcript(transcript).mutant_transcript is None


@pytest.mark.parametrize("orientation", ["[[", "]]"])
def test_mate_orientation_does_not_determine_local_side(transcript, orientation):
    sv = StructuralVariant(transcript.contig, transcript.exons[0].start + 5,
                           "BND", mate_contig="22", mate_start=15_500_000,
                           mate_orientation=orientation, genome=transcript.genome)
    with pytest.warns(UserWarning, match="can't read the breakend orientation"):
        effect = sv.effect_on_transcript(transcript)
    assert isinstance(effect, TranslocationToIntergenic)
    assert effect.mutant_transcript is None


@pytest.mark.parametrize("location", ["before", "after", "wrong_contig"])
def test_unrelated_transcript_does_not_supply_a_segment(transcript, location):
    position = {"before": transcript.start - 1, "after": transcript.end + 1,
                "wrong_contig": transcript.start + 5}[location]
    contig = "1" if location == "wrong_contig" else transcript.contig
    sv = StructuralVariant(contig, position, "BND", alt="N.", genome=transcript.genome)
    assert _build_translocation_mutant_transcript(sv, transcript) is None


@pytest.mark.parametrize("alt", ["N]22:15500000]", "]22:15500000]N", "N.", ".N", "<BND>"])
def test_supplied_assembly_is_unchanged(transcript, alt):
    assembly = "ATGGCTTAA"
    sv = StructuralVariant(transcript.contig, transcript.exons[0].start + 5,
                           "BND", alt=alt, alt_assembly=assembly, genome=transcript.genome)
    mt = sv.effect_on_transcript(transcript).mutant_transcript
    assert mt.cdna_sequence == assembly
    assert mt.mutant_protein_sequence is None
    assert mt.evidence == {"source": "alt_assembly"}
    (segment,) = mt.reference_segments
    assert segment.source.sequence == assembly
    assert (segment.start, segment.end, segment.label) == (0, len(assembly), "alt_assembly")
