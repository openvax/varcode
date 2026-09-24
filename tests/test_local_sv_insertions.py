"""Local DEL/DUP junction reconstruction from reciprocal VCF breakends (#491).

ALT orientation and retained anchors follow VCF 4.5 sections 5.4–5.4.1:
https://samtools.github.io/hts-specs/VCFv4.5.pdf
"""

import pytest
from pyensembl import cached_release

from varcode import StructuralVariant
from varcode.effects import FivePrimeUTR, FrameShift, StartLoss, ThreePrimeUTR, Unresolved
from varcode.effects import structural
from varcode.effects.codon_tables import translate_sequence
from varcode.nucleotides import reverse_complement

from .test_structural_variant_annotator import _pair_direct_breakends


@pytest.fixture(autouse=True)
def local_model(monkeypatch):
    monkeypatch.setattr(structural, "_fusion_partners", lambda *args: ())
    monkeypatch.setattr(structural, "_enumerate_and_attach_cryptics", lambda *args: None)
    monkeypatch.setattr(structural, "_enumerate_and_attach_splice_outcomes", lambda *args: None)


@pytest.fixture(params=["ENST00000003084", "ENST00000357654"])
def tx(request):
    return cached_release(81).transcript_by_id(request.param)


def _paired(tx, kind, low, high, insert="AC", mate_insert=None):
    """Insert is in forward-genomic orientation in both reciprocal ALTs."""
    if mate_insert is None:
        mate_insert = insert
    if kind == "DEL":
        low_alt = "N" + insert + "[%s:%d[" % (tx.contig, high)
        high_alt = "]%s:%d]" % (tx.contig, low) + mate_insert + "N"
    else:
        low_alt = "]%s:%d]" % (tx.contig, high) + insert + "N"
        high_alt = "N" + mate_insert + "[%s:%d[" % (tx.contig, low)
    a = StructuralVariant(tx.contig, low, "BND", ref="N", alt=low_alt,
                          mate_contig=tx.contig, mate_start=high, genome=tx.genome,
                          info={"SVTYPE": kind, "MATEID": "b"})
    b = StructuralVariant(tx.contig, high, "BND", ref="N", alt=high_alt,
                          mate_contig=tx.contig, mate_start=low, genome=tx.genome,
                          info={"SVTYPE": kind, "MATEID": "a"})
    variant = _pair_direct_breakends(a, b)
    assert variant.sv_type == kind
    return variant


def _expected_cdna(tx, kind, low, high, insert):
    start, end = sorted((tx.spliced_offset(low), tx.spliced_offset(high)))
    if kind == "DEL":
        return tx.sequence[:start + 1] + insert + tx.sequence[end:]
    return tx.sequence[:end + 1] + insert + tx.sequence[start:end + 1] + tx.sequence[end + 1:]


def _assert_sequence(effect, expected, cds_start):
    model = effect.mutant_transcript
    assert model.cdna_sequence == expected
    assert "".join(s.source.sequence[s.start:s.end] for s in model.reference_segments) == expected
    coding = expected[cds_start:]
    protein = translate_sequence(coding[:len(coding) // 3 * 3], to_stop=True)
    assert model.mutant_protein_sequence == protein
    primary = getattr(effect, "most_likely_effect", effect)
    if not isinstance(primary, (FivePrimeUTR, ThreePrimeUTR)):
        assert primary.mutant_protein_sequence == protein


@pytest.mark.parametrize("kind", ["DEL", "DUP"])
@pytest.mark.parametrize("insert", ["", "A", "AC", "ACG"])
@pytest.mark.parametrize("distance", [22, 23])
def test_exonic_insert_sequence_and_frame(tx, kind, insert, distance):
    low, high = tx.exons[4].start + 10, tx.exons[4].start + 10 + distance
    variant = _paired(tx, kind, low, high, insert)
    effect = variant.effect_on_transcript(tx)
    oriented = insert if tx.strand == "+" else reverse_complement(insert)
    expected = _expected_cdna(tx, kind, low, high, oriented)
    _assert_sequence(effect, expected, min(tx.start_codon_spliced_offsets))
    delta = len(expected) - len(tx.sequence)
    assert isinstance(effect.most_likely_effect, FrameShift) == (delta % 3 != 0)
    assert effect.modifies_coding_sequence is True
    assert effect.modifies_protein_sequence is True
    model = effect.mutant_transcript
    assert model.evidence["assumption"] == "reference_splicing"
    assert model.evidence["splice_ambiguous"] is True
    if insert:
        assert model.evidence["junction_inserted_sequence"] == oriented
        assert model.evidence["junction_insertion_status"] == "retained"
        assert model.reference_segments[1].label == "junction_insertion"
    else:
        assert "junction_insertion_status" not in model.evidence


@pytest.mark.parametrize("kind", ["DEL", "DUP"])
def test_intronic_insert_is_excluded_with_reference_splicing(tx, kind):
    exon = tx.exons[4]
    variant = _paired(tx, kind, exon.start - 20, exon.end + 20)
    effect = variant.effect_on_transcript(tx)
    control = _paired(tx, kind, exon.start - 20, exon.end + 20, "").effect_on_transcript(tx)
    assert type(getattr(effect, "most_likely_effect", effect)) is type(
        getattr(control, "most_likely_effect", control))
    model = effect.mutant_transcript
    assert model.cdna_sequence == control.mutant_transcript.cdna_sequence
    assert model.mutant_protein_sequence == control.mutant_transcript.mutant_protein_sequence
    assert model.evidence["junction_insertion_status"] == "excluded_by_reference_splicing"
    assert all(s.label != "junction_insertion" for s in model.reference_segments)


@pytest.mark.parametrize("kind", ["DEL", "DUP"])
@pytest.mark.parametrize("exonic_side", ["low", "high"])
def test_mixed_retention_is_unresolved(tx, kind, exonic_side):
    exon = tx.exons[4]
    low = exon.start + 10 if exonic_side == "low" else exon.start - 20
    high = exon.end - 10 if exonic_side == "high" else exon.end + 20
    effect = _paired(tx, kind, low, high).effect_on_transcript(tx)
    assert isinstance(effect.most_likely_effect, Unresolved)
    assert effect.mutant_transcript.cdna_sequence is None
    assert effect.mutant_transcript.mutant_protein_sequence is None
    assert effect.mutant_transcript.evidence["sequence_status"] == "unresolved_insertion_retention"


@pytest.mark.parametrize("kind", ["DEL", "DUP"])
@pytest.mark.parametrize("problem", ["conflict", "unreadable"])
def test_unknown_junction_sequence_stays_unresolved(tx, kind, problem):
    exon = tx.exons[4]
    variant = _paired(tx, kind, exon.start + 10, exon.start + 32, mate_insert="AG")
    if problem == "unreadable":
        # The ALT replacement no longer contains its REF anchor.
        variant.source_variants[1]._sv_alt = variant.source_variants[1].symbolic_alt.replace("N", "T")
    effect = variant.effect_on_transcript(tx)
    assert isinstance(effect.most_likely_effect, Unresolved)
    assert effect.mutant_transcript.cdna_sequence is None
    assert effect.mutant_transcript.mutant_protein_sequence is None
    assert effect.mutant_transcript.evidence["sequence_status"] == "unresolved_junction_insertion"


@pytest.mark.parametrize("kind", ["DEL", "DUP"])
def test_exon_boundary_anchors_retain_insert(tx, kind):
    first, second = sorted(tx.exons[4:6], key=lambda e: e.start)
    # Both retained anchors are exonic even though the DEL's affected span
    # starts and ends in introns. Retention must use the junction anchors.
    low, high = first.end, second.start
    if kind == "DEL":
        first, second = sorted((tx.exons[3], tx.exons[5]), key=lambda e: e.start)
        low, high = first.end, second.start
    effect = _paired(tx, kind, low, high).effect_on_transcript(tx)
    insert = "AC" if tx.strand == "+" else "GT"
    expected = _expected_cdna(tx, kind, low, high, insert)
    _assert_sequence(effect, expected, min(tx.start_codon_spliced_offsets))
    assert effect.mutant_transcript.evidence["junction_insertion_status"] == "retained"


def test_deleted_start_remains_start_loss(tx):
    low, high = min(tx.start_codon_positions) - 1, max(tx.start_codon_positions) + 1
    effect = _paired(tx, "DEL", low, high, "ATG").effect_on_transcript(tx)
    assert isinstance(effect.most_likely_effect, StartLoss)
    assert effect.mutant_transcript.cdna_sequence is not None
    assert effect.mutant_transcript.mutant_protein_sequence is None


def test_intron_deletion_can_insert_exonic_bases(tx):
    first, second = sorted(tx.exons[4:6], key=lambda e: e.start)
    low, high = first.end, second.start
    effect = _paired(tx, "DEL", low, high).effect_on_transcript(tx)
    insert = "AC" if tx.strand == "+" else "GT"
    expected = _expected_cdna(tx, "DEL", low, high, insert)
    _assert_sequence(effect, expected, min(tx.start_codon_spliced_offsets))
    assert isinstance(getattr(effect, "most_likely_effect", effect), FrameShift)


@pytest.mark.parametrize("kind", ["DEL", "DUP"])
@pytest.mark.parametrize("region,cls", [("5utr", FivePrimeUTR), ("3utr", ThreePrimeUTR)])
def test_utr_insert_preserves_cds_mapping(tx, kind, region, cls):
    exon = tx.exons[0] if region == "5utr" else tx.exons[-1]
    at_low = (region == "5utr") == (tx.strand == "+")
    low, high = ((exon.start + 1, exon.start + 8) if at_low
                 else (exon.end - 8, exon.end - 1))
    effect = _paired(tx, kind, low, high).effect_on_transcript(tx)
    insert = "AC" if tx.strand == "+" else "GT"
    expected = _expected_cdna(tx, kind, low, high, insert)
    cds_start = min(tx.start_codon_spliced_offsets)
    if region == "5utr":
        cds_start += len(expected) - len(tx.sequence)
    _assert_sequence(effect, expected, cds_start)
    assert isinstance(effect.most_likely_effect, cls)
    assert effect.modifies_coding_sequence is False
    assert effect.modifies_protein_sequence is False


@pytest.mark.parametrize("kind", ["DEL", "DUP"])
def test_ambiguous_insert_is_preserved_without_guessing_protein(tx, kind):
    exon = tx.exons[4]
    low, high = exon.start + 10, exon.start + 32
    effect = _paired(tx, kind, low, high, "N").effect_on_transcript(tx)
    assert effect.mutant_transcript.cdna_sequence == _expected_cdna(tx, kind, low, high, "N")
    assert isinstance(effect.most_likely_effect, Unresolved)
    assert effect.mutant_transcript.mutant_protein_sequence is None


def test_nonlocal_insertion_retention_is_not_assumed(tx):
    variant = _paired(tx, "DEL", tx.start - 20, tx.end + 20)
    effect = variant.effect_on_transcript(tx)
    assert isinstance(effect.most_likely_effect, Unresolved)
    assert effect.mutant_transcript.cdna_sequence is None


def test_explicit_assembly_is_preserved(tx):
    exon = tx.exons[4]
    variant = _paired(tx, "DEL", exon.start + 10, exon.start + 32, mate_insert="AG")
    variant.alt_assembly = "ATGACGTAA"
    effect = variant.effect_on_transcript(tx)
    assert isinstance(effect.most_likely_effect, Unresolved)
    assert effect.mutant_transcript.cdna_sequence == "ATGACGTAA"
    assert effect.mutant_transcript.evidence["source"] == "alt_assembly"
