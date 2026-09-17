"""Symbolic POS is retained; only POS+1..END changes (#404, VCF 4.3)."""

from types import SimpleNamespace

import pytest
from pyensembl import cached_release

from varcode import StructuralVariant, parse_symbolic_alt, load_vcf
from varcode.effects.structural import (
    _build_deletion_mutant_transcript, _build_duplication_mutant_transcript,
    _build_inversion_mutant_transcript,
)
from varcode.nucleotides import reverse_complement


@pytest.mark.parametrize("kind", ["DEL", "DUP", "INV", "CNV", "CN0", "CN3"])
def test_parsed_span_excludes_padding_but_keeps_record_coordinates(kind):
    sv = parse_symbolic_alt("1", 100, "T", "<%s>" % kind, info={"END": 105})
    assert (sv.start, sv.end) == (100, 105)
    assert (sv.affected_start, sv.affected_end, sv.length) == (101, 105, 5)
    assert StructuralVariant.from_json(sv.to_json()) == sv


@pytest.mark.parametrize("kind", ["DEL", "DUP", "INV", "CNV"])
@pytest.mark.parametrize("info", [{}, {"END": 100}, {"END": 99}])
def test_no_nonempty_span_cannot_delete_or_copy_the_anchor(kind, info):
    with pytest.raises(ValueError, match="END greater than POS"):
        parse_symbolic_alt("1", 100, "T", "<%s>" % kind, info=info)


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("kind,builder", [
    ("DEL", _build_deletion_mutant_transcript),
    ("DUP", _build_duplication_mutant_transcript),
    ("INV", _build_inversion_mutant_transcript),
])
def test_exact_mutant_cdna_on_both_strands(kind, builder, reverse):
    sequence = "ACGTTGCAAGCTTAGGCTAC"
    tx = SimpleNamespace(exons=[SimpleNamespace(start=100, end=119)],
                         contig="1", start=100, end=119, sequence=sequence,
                         on_backward_strand=reverse)
    sv = parse_symbolic_alt("1", 109, "G", "<%s>" % kind, info={"END": 114})
    a, b = (5, 10) if reverse else (10, 15)
    if kind == "DEL":
        expected = sequence[:a] + sequence[b:]
    elif kind == "DUP":
        expected = sequence[:b] + sequence[a:b] + sequence[b:]
    else:
        expected = sequence[:a] + reverse_complement(sequence[a:b]) + sequence[b:]
    model = builder(sv, tx)
    if kind == "INV":
        # This builder preserves the layout but intentionally defers sequence.
        observed = "".join(reverse_complement(sequence[s.start:s.end])
                           if s.strand == "-" else sequence[s.start:s.end]
                           for s in model.reference_segments)
        assert observed == expected
    else:
        assert model.cdna_sequence == expected


def test_vcf_deletion_does_not_remove_cftr_exon_terminal_anchor(tmp_path):
    genome = cached_release(81)
    tx = genome.transcript_by_id("ENST00000003084")
    first, second = tx.exons[:2]
    path = tmp_path / "anchor.vcf"
    path.write_text(
        "##fileformat=VCFv4.3\n"
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "7\t%d\t.\tT\t<DEL>\t.\tPASS\tEND=%d\n" % (first.end, second.end))
    (sv,) = load_vcf(str(path), genome=genome, parse_structural_variants=True)
    effect = sv.effect_on_transcript(tx)
    assert first not in effect.affected_exons
    assert second in effect.affected_exons
    a, b = first.end - first.start + 1, second.end - second.start + 1
    assert effect.mutant_transcript.cdna_sequence == tx.sequence[:a] + tx.sequence[a+b:]


def test_direct_constructor_span_defaults_remain_explicit():
    sv = StructuralVariant("1", 100, "DEL", end=105)
    assert (sv.affected_start, sv.affected_end, sv.length) == (100, 105, 6)
