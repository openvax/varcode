"""A partial overlap cannot establish a local DUP/INV transcript (#405)."""

import pytest
from pyensembl import cached_release

from varcode import StructuralVariant, parse_symbolic_alt
from varcode.effects import Inversion, LargeDuplication
from varcode.effects import structural


@pytest.fixture(params=["ENST00000003084", "ENST00000357654"])
def transcript(request):
    return cached_release(81).transcript_by_id(request.param)


@pytest.fixture(autouse=True)
def no_optional_candidates(monkeypatch):
    # Isolate the span fallback; fusion preservation also has real esvee tests.
    monkeypatch.setattr(structural, "_fusion_partners", lambda *args: ())
    monkeypatch.setattr(structural, "_enumerate_and_attach_cryptics", lambda *args: None)
    monkeypatch.setattr(structural, "_enumerate_and_attach_splice_outcomes", lambda *args: None)


def make_sv(transcript, kind, outside, symbolic=False, assembly=None):
    start, end = transcript.start + 50, transcript.end - 50
    if outside in ("left", "both"):
        start = transcript.start - 100
    if outside in ("right", "both"):
        end = transcript.end + 100
    if symbolic:
        return parse_symbolic_alt(
            transcript.contig, start, "N", "<%s>" % kind,
            info={"END": end}, genome=transcript.genome)
    return StructuralVariant(
        transcript.contig, start, kind, end=end,
        genome=transcript.genome, alt_assembly=assembly)


@pytest.mark.parametrize("kind,cls", [("DUP", LargeDuplication), ("INV", Inversion)])
@pytest.mark.parametrize("outside", ["left", "right", "both"])
@pytest.mark.parametrize("symbolic", [False, True])
def test_nonlocal_event_keeps_type_without_inventing_transcript(
        transcript, kind, cls, outside, symbolic):
    sv = make_sv(transcript, kind, outside, symbolic)
    effect = sv.effect_on_transcript(transcript)
    assert isinstance(effect, cls)
    assert effect.mutant_transcript is None
    assert all(c.effect.mutant_transcript is None for c in effect.candidates)


@pytest.mark.parametrize("kind", ["DUP", "INV"])
@pytest.mark.parametrize("outside", ["left", "right", "both"])
def test_nonlocal_event_preserves_supplied_assembly(transcript, kind, outside):
    assembly = "ATG" + "GCT" * 10 + "TAA"
    sv = make_sv(transcript, kind, outside, assembly=assembly)
    effect = sv.effect_on_transcript(transcript)
    assert effect.mutant_transcript.cdna_sequence == assembly


@pytest.mark.parametrize("kind", ["DUP", "INV"])
def test_wholly_internal_event_retains_local_model(transcript, kind):
    sv = make_sv(transcript, kind, outside=None)
    effect = sv.effect_on_transcript(transcript)
    assert effect.mutant_transcript is not None
    if kind == "DUP":
        assert len(effect.mutant_transcript.cdna_sequence) > len(transcript.sequence)
    else:
        assert any(s.strand == "-" for s in effect.mutant_transcript.reference_segments)


@pytest.mark.parametrize("transcript_id", ["ENST00000439541", "ENST00000550411"])
@pytest.mark.parametrize("kind", ["DUP", "INV"])
def test_gene_containment_does_not_establish_short_isoform_structure(transcript_id, kind):
    # CBX5's long isoform spans both ends, but these short isoforms do not.
    genome = cached_release(95)
    tx = genome.transcript_by_id(transcript_id)
    sv = StructuralVariant("12", 54_259_010, kind, end=54_262_307, genome=genome)
    ends = [end for junction in sv.junctions for end in junction]
    assert all(structural._contains(tx.gene, end) for end in ends)
    assert any(not structural._contains(tx, end) for end in ends)
    effect = sv.effect_on_transcript(tx)
    assert isinstance(effect, (LargeDuplication, Inversion))
    assert effect.mutant_transcript is None
