"""Start-loss nucleotide context is available without asserting translation (#206)."""

import pytest
from pyensembl import cached_release

from varcode import Variant
from varcode.effects import EffectCollection, StartLoss
from varcode.nucleotides import reverse_complement


@pytest.mark.parametrize("annotator", ["fast", "protein_diff", "transcript_model"])
@pytest.mark.parametrize("transcript_id", ["ENST00000003084", "ENST00000357654"])
@pytest.mark.parametrize("alternate", ["ACG", "TAC"])
def test_start_loss_retains_oriented_codons_and_sequence(annotator, transcript_id, alternate):
    genome = cached_release(81)
    transcript = genome.transcript_by_id(transcript_id)
    ref, alt = "ATG", alternate
    if transcript.on_backward_strand:
        ref, alt = reverse_complement(ref), reverse_complement(alt)
    variant = Variant(transcript.contig, min(transcript.start_codon_positions), ref, alt, genome)
    effect = variant.effect_on_transcript(transcript, annotator=annotator)
    assert type(effect) is StartLoss
    assert effect.ref_codon == "ATG"
    assert effect.alt_codon == alternate
    assert effect.codon_change_status == "mapped"
    assert effect.ref_codon_offset == min(transcript.start_codon_spliced_offsets)
    assert effect.mutant_transcript.cdna_sequence[
        effect.alt_codon_offset:effect.alt_codon_offset + 3] == alternate
    assert effect.alt_codon_amino_acid == {"ACG": "T", "TAC": "Y"}[alternate]
    assert effect.mutant_protein_sequence is None
    assert effect.mutant_transcript.mutant_protein_sequence is None
    assert not effect.mutant_transcript.is_identical_to_reference
    edit, = effect.mutant_transcript.edits
    assert edit.source_variant == variant
    assert effect.short_description == "p.M1? (start-loss)"
    # Native effects and collection archives must retain available nucleotide context.
    restored = EffectCollection.from_json(EffectCollection([effect]).to_json())[0]
    assert restored.ref_codon == effect.ref_codon
    assert restored.alt_codon == effect.alt_codon
    assert restored.codon_change_status == "mapped"
    assert restored.mutant_transcript.cdna_sequence == effect.mutant_transcript.cdna_sequence
    assert restored.mutant_protein_sequence is None


@pytest.mark.parametrize("annotator", ["fast", "protein_diff", "transcript_model"])
def test_deleted_start_keeps_edit_without_reading_the_next_triplet_as_its_replacement(annotator):
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084")
    variant = Variant("7", 117480094, "CATG", "C", genome)
    effect = variant.effect_on_transcript(transcript, annotator=annotator)
    assert type(effect) is StartLoss
    assert effect.ref_codon == "ATG"
    assert effect.alt_codon is None
    assert effect.alt_codon_offset is None
    assert effect.codon_change_status == "start_removed"
    assert effect.alt_codon_amino_acid is None
    assert effect.mutant_transcript.cdna_sequence is not None
    assert effect.mutant_protein_sequence is None
    edit, = effect.mutant_transcript.edits
    assert edit.cdna_start == effect.ref_codon_offset
    assert edit.cdna_end == effect.ref_codon_offset + 3
    assert edit.alt_bases == ""


def test_unresolvable_reference_does_not_invent_alternate_codon():
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084")
    effect = StartLoss(Variant("7", 117480096, "G", "C", genome), transcript)
    assert effect.ref_codon == "ATG"
    assert effect.alt_codon is None
    assert effect.codon_change_status == "edit_unresolved"
    assert effect.mutant_transcript is None


def test_original_mouse_example_has_a_pinned_reference():
    genome = cached_release(95, species="mouse")
    transcript = genome.transcript_by_id("ENSMUST00000204003")
    effect = Variant("10", 58224929, "A", "G", genome).effect_on_transcript(transcript)
    assert transcript.on_backward_strand
    assert type(effect) is StartLoss
    assert (effect.ref_codon, effect.alt_codon) == ("ATG", "ACG")
    assert effect.alt_codon_amino_acid == "T"
    assert effect.mutant_protein_sequence is None


@pytest.mark.parametrize("annotator", ["fast", "protein_diff", "transcript_model"])
@pytest.mark.parametrize("ref,alt,codon", [("A", "AA", "AAT"), ("AT", "A", "AGC")])
def test_start_overlapping_indels_retain_the_mapped_triplet(annotator, ref, alt, codon):
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084")
    variant = Variant("7", 117480095, ref, alt, genome)
    effect = variant.effect_on_transcript(transcript, annotator=annotator)
    assert type(effect) is StartLoss
    assert (effect.ref_codon, effect.alt_codon) == ("ATG", codon)
    assert effect.codon_change_status == "mapped"
    assert abs(effect.mutant_transcript.total_length_delta) == 1
    assert effect.mutant_transcript.mutant_protein_sequence is None


@pytest.mark.parametrize("transcript_id,alternate,amino_acid", [
    ("ENST00000003084", "TGA", "*"),
    ("ENST00000361624", "TGA", "W"),
    ("ENST00000361624", "CTG", "L"),
])
def test_literal_codon_decoding_uses_the_transcripts_genetic_code(
        transcript_id, alternate, amino_acid):
    genome = cached_release(81)
    transcript = genome.transcript_by_id(transcript_id)
    variant = Variant(transcript.contig, min(transcript.start_codon_positions),
                      "ATG", alternate, genome)
    effect = variant.effect_on_transcript(transcript, annotator="transcript_model")
    assert type(effect) is StartLoss
    assert effect.alt_codon == alternate
    assert effect.alt_codon_amino_acid == amino_acid
    assert effect.mutant_protein_sequence is None


@pytest.mark.parametrize("annotator", ["fast", "transcript_model"])
def test_csv_retains_saved_start_loss_without_reannotation(annotator, tmp_path, monkeypatch):
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084")
    effect = Variant("7", 117480096, "T", "C", genome).effect_on_transcript(
        transcript, annotator=annotator)
    collection = EffectCollection([effect], annotator=annotator)
    path = tmp_path / "start-loss.csv"
    collection.to_csv(path)

    def unavailable(*args, **kwargs):
        raise AssertionError("CSV attempted reannotation")

    monkeypatch.setattr(Variant, "effect_on_transcript", unavailable)
    restored = EffectCollection.from_csv(path)
    assert restored.annotator == annotator
    assert restored[0].to_dict() == effect.to_dict()
    assert restored.to_dataframe().iloc[0]["alt_codon"] == "ACG"


def test_unmapped_and_unavailable_sequences_do_not_invent_codons():
    from varcode.mutant_transcript import MutantTranscript
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084")
    for model, status in [
        (MutantTranscript.from_sequence("ACG", reference_transcript=transcript), "edit_unresolved"),
        (MutantTranscript(reference_transcript=transcript), "sequence_unavailable"),
    ]:
        effect = StartLoss(None, transcript, mutant_transcript=model)
        assert effect.ref_codon == "ATG"
        assert effect.alt_codon is None
        assert effect.codon_change_status == status
        assert effect.mutant_transcript is not None


def test_patient_baseline_and_joint_edits_are_preserved():
    from varcode.effects.classify import classify_from_protein_diff
    from varcode.mutant_transcript import apply_variants_to_transcript
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084")
    germline = Variant("7", 117480095, "A", "C", genome)
    somatic = Variant("7", 117480096, "T", "C", genome)
    baseline = apply_variants_to_transcript([germline], transcript)
    joint = apply_variants_to_transcript([germline, somatic], transcript)
    effect = classify_from_protein_diff(
        somatic, transcript, "M" + baseline.mutant_protein_sequence[1:],
        joint.mutant_protein_sequence, mutant_transcript=joint,
        baseline_transcript=baseline)
    assert type(effect) is StartLoss
    assert (effect.ref_codon, effect.alt_codon) == ("CTG", "CCG")
    assert len(effect.mutant_transcript.edits) == 2
    assert effect.mutant_transcript.cdna_sequence == joint.cdna_sequence
    assert effect.mutant_transcript.mutant_protein_sequence is None


def test_source_declared_orf_retains_start_loss_context_and_provenance():
    pytest.importorskip("orfdata")
    from .test_non_aug_orfs import declared_orf
    reference = declared_orf("GTG")
    variant = reference.transcript.variant(3, "GTG", "GCG")
    effect = variant.effects(extra_orfs=[reference]).orf_annotations[0].effect
    assert type(effect) is StartLoss
    assert (effect.ref_codon, effect.alt_codon) == ("GTG", "GCG")
    assert effect.alt_codon_amino_acid == "A"
    restored = StartLoss.from_json(effect.to_json())
    assert restored.orf_annotation.record.source == reference.record.source
    assert restored.mutant_transcript.cdna_sequence == effect.mutant_transcript.cdna_sequence
    assert restored.mutant_protein_sequence is None
