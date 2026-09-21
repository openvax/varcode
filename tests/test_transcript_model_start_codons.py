"""Initiation codons in realized products (openvax/varcode#476)."""

import pytest
from pyensembl import cached_release

from varcode import Variant
from varcode.effect_hypotheses import SplicePlan
from varcode.effects import AlternateStartCodon, StartLoss, Substitution
from varcode.genomic_layout import GenomicLayout
from varcode.nucleotides import reverse_complement
from varcode.transcript_layout import (
    build_exon_runs,
    classify_products,
    realize_exon_path,
    realize_splice_plan,
)

from .test_splice_graph import point
from .test_transcript_layout import Transcript, provider, structural


@pytest.mark.parametrize("annotator", ["fast", "protein_diff", "transcript_model"])
@pytest.mark.parametrize("transcript_id,alt_codon", [
    ("ENST00000003084", "CTG"),
    ("ENST00000003084", "TTG"),
    ("ENST00000357654", "CTG"),
    ("ENST00000357654", "TTG"),
    ("ENST00000361624", "ATT"),
    ("ENST00000361624", "ATC"),
    ("ENST00000361624", "ATA"),
    ("ENST00000361624", "GTG"),
])
def test_recognized_alternate_start_preserves_protein(
        annotator, transcript_id, alt_codon):
    genome = cached_release(81)
    transcript = genome.transcript_by_id(transcript_id)
    ref, alt = "ATG", alt_codon
    if transcript.on_backward_strand:
        ref, alt = reverse_complement(ref), reverse_complement(alt)
    variant = Variant(
        transcript.contig, min(transcript.start_codon_positions), ref, alt,
        genome=genome)

    effect = variant.effect_on_transcript(transcript, annotator=annotator)

    assert type(effect) is AlternateStartCodon
    assert effect.ref_codon == "ATG"
    assert effect.alt_codon == alt_codon
    assert effect.modifies_coding_sequence is True
    assert effect.modifies_protein_sequence is False
    if annotator == "transcript_model":
        outcome = effect.candidates[0].outcomes[0]
        assert outcome.baseline.start_codon_present is True
        assert outcome.mutant.start_codon_present is True
        assert outcome.baseline.protein_sequence == str(transcript.protein_sequence)
        assert outcome.mutant.protein_sequence == str(transcript.protein_sequence)


@pytest.mark.parametrize("transcript_id,alt_codon", [
    ("ENST00000003084", "GTG"),
    ("ENST00000357654", "GTG"),
    ("ENST00000361624", "CTG"),
    ("ENST00000361624", "TTG"),
])
def test_start_recognition_uses_the_transcripts_codon_table(
        transcript_id, alt_codon):
    genome = cached_release(81)
    transcript = genome.transcript_by_id(transcript_id)
    ref, alt = "ATG", alt_codon
    if transcript.on_backward_strand:
        ref, alt = reverse_complement(ref), reverse_complement(alt)
    variant = Variant(
        transcript.contig, min(transcript.start_codon_positions), ref, alt,
        genome=genome)

    effect = variant.effect_on_transcript(transcript, annotator="transcript_model")

    assert type(effect) is StartLoss
    assert effect.modifies_protein_sequence is True
    assert effect.candidates[0].outcomes[0].mutant.start_codon_present is False


@pytest.fixture(params=["exon_path", "splice_plan"])
def realize(request):
    if request.param == "exon_path":
        return realize_exon_path

    def splice_plan(transcript, layout):
        plan = SplicePlan(
            choices=(),
            kept_runs=tuple(run.key for run in build_exon_runs(transcript, layout)))
        product, unresolved = realize_splice_plan(transcript, layout, plan, ())
        assert unresolved == ()
        return product

    return splice_plan


@pytest.mark.parametrize("reverse", [False, True])
def test_both_layout_paths_translate_alternate_start_as_methionine(realize, reverse):
    transcript = Transcript(reverse)
    layout = GenomicLayout.from_transcript(
        transcript, flank=0, sequence_provider=provider(reverse))
    variant = point(18, "T", "G") if reverse else point(1, "A", "C")
    baseline = realize(transcript, layout)
    mutant = realize(transcript, layout.apply_point_variant(variant))

    effect = classify_products(variant, transcript, baseline, mutant)

    assert type(effect) is AlternateStartCodon
    assert (effect.ref_codon, effect.alt_codon) == ("ATG", "CTG")
    assert mutant.cdna_sequence == "CTGCAAGCTTAG"
    assert mutant.protein_sequence == "MQA"
    assert mutant.start_codon_present is True


@pytest.mark.parametrize("reverse", [False, True])
def test_deleted_start_is_still_start_loss(realize, reverse):
    transcript = Transcript(reverse)
    layout = GenomicLayout.from_transcript(
        transcript, flank=0, sequence_provider=provider(reverse))
    variant = structural("DEL", *(16, 18) if reverse else (1, 3))
    baseline = realize(transcript, layout)
    mutant = realize(transcript, layout.apply_structural_variant(variant))

    effect = classify_products(variant, transcript, baseline, mutant)

    assert type(effect) is StartLoss
    assert mutant.start_codon_present is False
    assert mutant.protein_sequence == ""


def test_alternate_start_is_compared_to_patient_baseline(realize):
    transcript = Transcript()
    layout = GenomicLayout.from_transcript(
        transcript, flank=0, sequence_provider=provider())
    layout = layout.apply_point_variant(point(1, "A", "C"))
    baseline = realize(transcript, layout)
    variant = point(1, "C", "T")
    mutant = realize(transcript, layout.apply_point_variant(variant))

    effect = classify_products(variant, transcript, baseline, mutant)

    assert type(effect) is AlternateStartCodon
    assert (effect.ref_codon, effect.alt_codon) == ("CTG", "TTG")
    assert baseline.protein_sequence == mutant.protein_sequence == "MQA"


def test_alternate_start_does_not_hide_downstream_protein_change(realize):
    transcript = Transcript()
    layout = GenomicLayout.from_transcript(
        transcript, flank=0, sequence_provider=provider())
    baseline = realize(transcript, layout)
    variant = point(1, "A", "C")
    mutant_layout = layout.apply_point_variant(variant).apply_point_variant(
        point(6, "C", "A"))
    mutant = realize(transcript, mutant_layout)

    effect = classify_products(variant, transcript, baseline, mutant)

    assert type(effect) is Substitution
    assert (effect.aa_ref, effect.aa_alt) == ("Q", "K")
    assert mutant.protein_sequence == "MKA"
