"""Patient-baseline initiation is distinct from internal codon decoding (#634)."""

from types import SimpleNamespace

import pytest
from pyensembl import cached_release

from varcode import GermlineContext, TranscriptEdit, Variant
from varcode.effects import EffectCollection, StartLoss
from varcode.effects.effect_classes import PhaseCandidateSet
from varcode.mutant_transcript import _translate_from_cds, apply_variants_to_transcript
from varcode.nucleotides import reverse_complement


def start_variant(transcript, offset, ref, alt):
    positions = sorted(transcript.start_codon_positions,
                       reverse=transcript.on_backward_strand)
    if transcript.on_backward_strand:
        ref, alt = reverse_complement(ref), reverse_complement(alt)
    return Variant(transcript.contig, positions[offset], ref, alt, transcript.genome)


@pytest.mark.parametrize("transcript_id,start,alternate", [
    ("ENST00000003084", "CTG", "CCG"),
    ("ENST00000003084", "TTG", "TCG"),
    ("ENST00000357654", "CTG", "CCG"),
    ("ENST00000357654", "TTG", "TCG"),
    ("ENST00000361624", "GTG", "GCG"),
])
@pytest.mark.parametrize("annotator", ["fast", "protein_diff"])
def test_germline_alternate_start_then_somatic_start_loss(transcript_id, start, alternate, annotator):
    transcript = cached_release(81).transcript_by_id(transcript_id)
    germline = start_variant(transcript, 0, "A", start[0])
    somatic = start_variant(transcript, 1, "T", "C")
    baseline = apply_variants_to_transcript([germline], transcript)
    assert baseline.mutant_protein_sequence == str(transcript.protein_sequence)

    effect = somatic.effect_on_transcript(
        transcript, annotator=annotator,
        germline=GermlineContext.from_variants([germline]),
        phase_resolver=SimpleNamespace(in_cis=lambda *args: True))
    assert type(effect) is StartLoss
    assert (effect.ref_codon, effect.alt_codon) == (start, alternate)
    assert effect.short_description == "p.M1? (start-loss)"
    assert effect.germline_phase_state == ("implicit" if transcript.contig == "MT" else "phased")
    assert effect.codon_change_status == "mapped"
    assert effect.mutant_protein_sequence is None
    model = effect.mutant_transcript
    assert model.mutant_protein_sequence is None
    assert {edit.source_variant for edit in model.edits} == {germline, somatic}
    assert model.cdna_sequence[effect.alt_codon_offset:effect.alt_codon_offset + 3] == alternate
    restored = EffectCollection.from_json(EffectCollection([effect]).to_json())[0]
    assert (restored.ref_codon, restored.alt_codon) == (start, alternate)
    assert restored.mutant_transcript.cdna_sequence == model.cdna_sequence
    assert restored.mutant_protein_sequence is None


@pytest.mark.parametrize("annotator", ["fast", "transcript_model"])
def test_unknown_phase_keeps_distinct_patient_start_codons(annotator):
    transcript = cached_release(81).transcript_by_id("ENST00000003084")
    germline = start_variant(transcript, 0, "A", "C")
    somatic = start_variant(transcript, 1, "T", "C")
    effect = somatic.effect_on_transcript(
        transcript, germline=GermlineContext.from_variants([germline]), annotator=annotator)
    assert type(effect) is (PhaseCandidateSet if annotator == "fast" else StartLoss)
    assert len(effect.candidates) == 2
    assert {(outcome.effect.ref_codon, outcome.effect.alt_codon)
            for outcome in effect.candidates} == {("CTG", "CCG"), ("ATG", "ACG")}
    assert all(type(outcome.effect) is StartLoss for outcome in effect.candidates)
    if annotator == "fast":
        assert all(c.evidence["phase_state"] == "unknown" for c in effect.candidates)
    else:
        assert all(h.evidence["phase_state"] == "unknown"
                   for c in effect.candidates for h in c.hypotheses)

    without = somatic.effect_on_transcript(transcript)
    trans = somatic.effect_on_transcript(
        transcript, germline=GermlineContext.from_variants([germline]),
        phase_resolver=SimpleNamespace(in_cis=lambda *args: False))
    assert type(trans) is type(without) is StartLoss
    assert (trans.ref_codon, trans.alt_codon) == (without.ref_codon, without.alt_codon)
    assert (without.ref_codon, without.alt_codon) == ("ATG", "ACG")


@pytest.mark.parametrize("code,start,declared,protein", [
    (1, "CTG", None, "ML"), (1, "TTG", None, "ML"),
    (1, "GTG", None, "VL"), (2, "GTG", None, "ML"),
    (2, "ATA", None, "ML"), (2, "CTG", None, "LL"),
    (3, "ATA", None, "MT"), (3, "CTG", None, "TT"),
    (3, "CTG", "CTG", "MT"), (1, "GTG", "GTG", "ML"),
])
def test_initiation_obeys_genetic_code_without_changing_internal_codons(code, start, declared, protein):
    transcript = SimpleNamespace(start_codon_spliced_offsets=(3, 4, 5),
                                 translation_table=code, start_codon=declared)
    sequence = "CCC" + start + "CTGTAA"
    assert _translate_from_cds(sequence, transcript) == protein
    upstream = TranscriptEdit(1, 1, "AA")
    assert _translate_from_cds(sequence[:1] + "AA" + sequence[1:], transcript, (upstream,)) == protein


def test_deleted_start_does_not_promote_downstream_alternate_initiator():
    transcript = SimpleNamespace(start_codon_spliced_offsets=(3, 4, 5), translation_table=1)
    deletion = TranscriptEdit(3, 6, "")
    # Removing ATG from CCCATGCTGTAA exposes an internal CTG, not a retained start.
    assert _translate_from_cds("CCCCTGTAA", transcript, (deletion,)) == "L"
