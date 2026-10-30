"""Unavailable patient comparisons must not become reference-only predictions."""

import pickle
from dataclasses import replace
from types import SimpleNamespace

import pytest
from pyensembl import cached_release

from varcode import EffectCollection, GermlineContext, Variant, effect_priority
from varcode.effects import Unresolved
from varcode.effects.effect_classes import PhaseCandidateSet
from varcode.nucleotides import reverse_complement
from varcode.splice_outcomes import SpliceOutcomeSet

# Legacy wrapper contract; the default composition path has separate regression tests.
pytestmark = pytest.mark.usefixtures("legacy_annotator")


@pytest.fixture
def transcript():
    return cached_release(81).transcript_by_id("ENST00000003084")


def pair(transcript, case):
    genome = transcript.genome
    if case == "joint_splice":
        return (Variant("7", 117531115, "G", "A", genome),
                Variant("7", 117531113, "A", "C", genome))
    if case == "intronic_baseline":
        return (Variant("7", 117531115, "G", "A", genome),
                Variant("7", 117531116, "T", "A", genome))
    if case == "conflicting_edits":
        return (Variant("7", 117531100, "T", "A", genome),
                Variant("7", 117531100, "T", "C", genome))
    if case == "baseline_ref_mismatch":
        return (Variant("7", 117531115, "G", "A", genome),
                Variant("7", 117531113, "T", "C", genome))
    raise AssertionError(case)


@pytest.mark.parametrize("case", [
    "joint_splice", "intronic_baseline", "conflicting_edits", "baseline_ref_mismatch"])
@pytest.mark.parametrize("annotator", ["fast", "protein_diff"])
def test_cis_composition_failure_is_explicit(transcript, case, annotator):
    somatic, germline = pair(transcript, case)
    result = somatic.effect_on_transcript(
        transcript, annotator=annotator,
        germline=GermlineContext.from_variants([germline]),
        phase_resolver=SimpleNamespace(in_cis=lambda *args: True))
    assert isinstance(result, Unresolved)
    assert result.mechanism == "germline_composition"
    assert result.reason
    assert result.modifies_coding_sequence is None
    assert result.modifies_protein_sequence is None
    assert result.mutant_protein_sequence is None
    assert result.germline_phase_state == "phased"
    assert result.germline_variants_in_window == (germline,)
    reference = somatic.effect_on_transcript(transcript, annotator="fast")
    assert type(result.reference_effect) is type(reference)
    assert result.reference_effect.short_description == reference.short_description
    assert effect_priority(result) == effect_priority(reference)


@pytest.mark.parametrize("case", ["joint_splice", "intronic_baseline", "conflicting_edits"])
def test_known_trans_retains_reference_annotation(transcript, case):
    somatic, germline = pair(transcript, case)
    result = somatic.effect_on_transcript(
        transcript, germline=GermlineContext.from_variants([germline]),
        phase_resolver=SimpleNamespace(in_cis=lambda *args: False))
    reference = somatic.effect_on_transcript(transcript)
    assert type(result) is type(reference)
    assert result.short_description == reference.short_description
    assert not isinstance(result, Unresolved)


@pytest.mark.parametrize("case", ["joint_splice", "intronic_baseline", "conflicting_edits"])
def test_unknown_phase_keeps_unresolved_cis_and_reference_trans(transcript, case):
    somatic, germline = pair(transcript, case)
    result = somatic.effect_on_transcript(
        transcript, germline=GermlineContext.from_variants([germline]))
    assert isinstance(result, PhaseCandidateSet)
    assert len(result.candidates) == 2
    cis = next(c for c in result.candidates if c.evidence["germline_variants"])
    trans = next(c for c in result.candidates if not c.evidence["germline_variants"])
    assert isinstance(cis.effect, Unresolved)
    assert cis.effect.modifies_protein_sequence is None
    assert cis.evidence["germline_variants"] == (germline,)
    assert type(trans.effect) is type(cis.effect.reference_effect)
    assert trans.short_description == cis.effect.reference_effect.short_description


@pytest.mark.parametrize("known_phase", [True, False])
def test_unresolved_comparison_survives_native_roundtrip(transcript, known_phase):
    somatic, germline = pair(transcript, "joint_splice")
    resolver = SimpleNamespace(in_cis=lambda *args: True) if known_phase else None
    effects = somatic.effects(
        germline=GermlineContext.from_variants([germline]), phase_resolver=resolver,
        transcript_filter=lambda t: t.id == transcript.id, raise_on_error=True)
    for restored in (EffectCollection.from_json(effects.to_json()),
                     pickle.loads(pickle.dumps(effects))):
        assert len(restored) == 1
        result = restored[0]
        if not known_phase:
            candidate = next(c for c in result.candidates if c.evidence["germline_variants"])
            assert candidate.evidence["germline_variants"] == (germline,)
            assert candidate.evidence["phase_state"] == "unknown"
            result = candidate.effect
        else:
            assert result.germline_phase_state == "phased"
            assert result.germline_variants_in_window == (germline,)
        assert isinstance(result, Unresolved)
        assert isinstance(result.reference_effect, SpliceOutcomeSet)
        assert result.reference_effect.variant == somatic
        assert result.modifies_protein_sequence is None
        assert result.mutant_protein_sequence is None


def test_legacy_unresolved_payload_remains_readable(transcript):
    somatic, _ = pair(transcript, "joint_splice")
    result = Unresolved.from_dict(dict(
        variant=somatic, transcript=transcript,
        mechanism="sequence_unavailable", reason="No sequence"))
    assert result.reference_effect is None
    assert result.modifies_protein_sequence is None


@pytest.mark.parametrize("transcript_id", ["ENST00000003084", "ENST00000357654"])
def test_conflicting_start_edits_on_both_strands(transcript_id):
    transcript = cached_release(81).transcript_by_id(transcript_id)
    position = sorted(transcript.start_codon_positions,
                      reverse=transcript.on_backward_strand)[0]
    bases = ["A", "C", "T"]
    if transcript.on_backward_strand:
        bases = [reverse_complement(base) for base in bases]
    ref, inherited, tumor = bases
    germline = Variant(transcript.contig, position, ref, inherited, transcript.genome)
    somatic = Variant(transcript.contig, position, ref, tumor, transcript.genome)
    result = somatic.effect_on_transcript(
        transcript, germline=GermlineContext.from_variants([germline]),
        phase_resolver=SimpleNamespace(in_cis=lambda *args: True))
    assert isinstance(result, Unresolved)
    assert "applied together" in result.reason
    assert result.modifies_protein_sequence is None
    assert result.reference_effect.variant == somatic


@pytest.mark.parametrize("missing", ["baseline", "joint"])
def test_unavailable_protein_is_not_replaced_by_reference_prediction(transcript, monkeypatch, missing):
    from varcode import mutant_transcript
    original = mutant_transcript.apply_variants_to_transcript
    somatic = Variant("7", 117531100, "T", "A", transcript.genome)
    germline = Variant("7", 117531101, "T", "C", transcript.genome)

    def build(variants, reference):
        model = original(variants, reference)
        assert model is not None
        if (missing == "baseline") == (len(variants) == 1):
            return replace(model, mutant_protein_sequence=None)
        return model

    monkeypatch.setattr(mutant_transcript, "apply_variants_to_transcript", build)
    result = somatic.effect_on_transcript(
        transcript, germline=GermlineContext.from_variants([germline]),
        phase_resolver=SimpleNamespace(in_cis=lambda *args: True))
    assert isinstance(result, Unresolved)
    assert "Protein sequence is unavailable" in result.reason
    assert missing in result.reason
    assert result.reference_effect.short_description == "p.L159M"
    assert result.mutant_protein_sequence is None


def test_supported_phase_alternatives_survive_native_roundtrip(transcript):
    somatic = Variant("7", 117531100, "T", "A", transcript.genome)
    germline = Variant("7", 117531101, "T", "C", transcript.genome)
    effects = EffectCollection([somatic.effect_on_transcript(
        transcript, germline=GermlineContext.from_variants([germline]))])
    for restored in (EffectCollection.from_json(effects.to_json()),
                     pickle.loads(pickle.dumps(effects))):
        candidates = restored[0].candidates
        assert {(c.short_description, c.evidence["germline_variants"])
                for c in candidates} == {("p.L159M", ()), ("p.S159T", (germline,))}
        assert all(c.evidence["phase_state"] == "unknown" for c in candidates)
