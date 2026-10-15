"""Source-declared initiation survives ORF variant and RNA-product annotation (#597)."""

from dataclasses import replace
from types import SimpleNamespace

import pytest

from varcode import VariantCollection, annotate_orf_rna_product
from varcode.effect_hypotheses import RealizedTranscriptProduct, SplicePlan
from varcode.effects import Failure, StartLoss, Substitution
from varcode.germline import GermlineContext
from varcode.genomic_layout import GenomicLayout
from varcode.orf import _ORFTranscript, _reference_provider
from varcode.transcript_layout import build_exon_runs, realize_exon_path, realize_splice_plan

from .test_extra_orfs import synthetic

orfdata = pytest.importorskip("orfdata", minversion="0.27.1")


def declared_orf(start, reverse=False, split=False):
    transcript = synthetic(reverse, split)
    transcript.sequence = transcript.sequence[:3] + start + transcript.sequence[6:]
    return replace(transcript.record(), start_codon=start).bind(transcript)


@pytest.mark.parametrize("start", ["ATG", "CTG", "GTG", "ACG"])
@pytest.mark.parametrize("reverse,split", [(False, False), (True, False), (False, True), (True, True)])
def test_downstream_missense_retains_source_initiation(start, reverse, split):
    reference = declared_orf(start, reverse, split)
    variant = reference.transcript.variant(6, "A", "G")
    annotation = variant.effects(extra_orfs=[reference]).orf_annotations[0]
    assert isinstance(annotation.effect, Substitution)
    assert annotation.effect.short_description == "p.K2E"
    assert len(annotation.outcomes) == 1
    outcome = annotation.outcomes[0]
    assert outcome.baseline.protein_sequence == "MKPG"
    assert outcome.mutant.protein_sequence == "MEPG"
    assert outcome.baseline.start_codon_present and outcome.mutant.start_codon_present
    assert outcome.mutant.initiation_context.start_codon == start
    assert outcome.mutant.reference_evidence["mutant_translation_evidence"] == "not_assessed"
    assert annotation.to_dict()["model"]["start_codon"] == start


@pytest.mark.parametrize("start", ["GTG", "ACG"])
@pytest.mark.parametrize("offset,ref,alt,effect,protein", [
    (6, "AAA", "GTG", "Substitution", "MVPG"),
    (6, "AAA", "ACG", "Substitution", "MTPG"),
    (8, "A", "G", "Silent", "MKPG"),
    (6, "A", "T", "PrematureStop", "M"),
    (9, "CCC", "", "Deletion", "MKG"),
    (9, "C", "", "FrameShift", "MKPGK"),
    (15, "T", "C", "StopLoss", "MKPGQMK"),
])
def test_declared_start_does_not_change_internal_translation(start, offset, ref, alt, effect, protein):
    reference = declared_orf(start)
    variant = reference.transcript.variant(offset, ref, alt)
    annotation = variant.effects(extra_orfs=[reference]).orf_annotations[0]
    assert type(annotation.effect).__name__ == effect
    assert annotation.outcomes[0].baseline.protein_sequence == "MKPG"
    assert annotation.outcomes[0].mutant.protein_sequence == protein


@pytest.mark.parametrize("start", ["GTG", "ACG"])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("alternate", ["ATG", "GCG", ""])
def test_changed_or_deleted_declared_start_keeps_reference_product(start, reverse, alternate):
    reference = declared_orf(start, reverse)
    variant = reference.transcript.variant(3, start, alternate)
    annotation = variant.effects(extra_orfs=[reference]).orf_annotations[0]
    outcome = annotation.outcomes[0]
    assert outcome.baseline.protein_sequence == "MKPG"
    assert outcome.baseline.start_codon_present
    if alternate == "ATG":
        assert type(annotation.effect).__name__ == "AlternateStartCodon"
        assert annotation.effect.ref_codon == start and annotation.effect.alt_codon == alternate
        assert outcome.mutant.protein_sequence == "MKPG"
    else:
        assert isinstance(annotation.effect, StartLoss)
        assert not outcome.mutant.start_codon_present
        assert outcome.mutant.protein_sequence == ""
    assert outcome.mutant.reference_evidence["mutant_translation_evidence"] == "not_assessed"


@pytest.mark.parametrize("start,alternate", [("ATG", "GTG"), ("GTG", "ACG"), ("ACG", "GTG")])
def test_undeclared_non_aug_start_is_not_enabled(start, alternate):
    reference = declared_orf(start)
    variant = reference.transcript.variant(3, start, alternate)
    annotation = variant.effects(extra_orfs=[reference]).orf_annotations[0]
    assert isinstance(annotation.effect, StartLoss)
    assert annotation.outcomes[0].baseline.protein_sequence == "MKPG"
    assert not annotation.outcomes[0].mutant.start_codon_present


@pytest.mark.parametrize("start", ["GTG", "ACG"])
def test_patient_baseline_and_joint_haplotype_keep_declared_start(start):
    reference = declared_orf(start)
    transcript = reference.transcript
    first, second = transcript.variant(6, "A", "G"), transcript.variant(7, "A", "C")
    annotation = first.effects(extra_orfs=[reference], germline=GermlineContext.from_variants(
        [second])).orf_annotations[0]
    assert {o.baseline.protein_sequence for o in annotation.outcomes} == {"MKPG", "MTPG"}
    assert {o.mutant.protein_sequence for o in annotation.outcomes} == {"MEPG", "MAPG"}
    resolver = SimpleNamespace(in_cis=lambda *args, **kwargs: True)
    result = VariantCollection([first, second]).effects(extra_orfs=[reference], phase_resolver=resolver)
    joint = [a for a in result.orf_annotations if len(a.variants) == 2]
    assert len(joint) == 1
    assert joint[0].outcomes[0].mutant.protein_sequence == "MAPG"


@pytest.mark.parametrize("start", ["GTG", "ACG"])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("splice_plan", [False, True])
def test_both_splice_realization_paths_keep_split_declared_start(start, reverse, splice_plan):
    reference = declared_orf(start, reverse, split=True)
    view = _ORFTranscript(reference)
    layout = GenomicLayout.from_transcript(view, flank=0, sequence_provider=_reference_provider(view))
    variant = reference.transcript.variant(6, "A", "G")
    mutant_layout = layout.apply_point_variant(variant)
    if splice_plan:
        plan = SplicePlan(choices=(), kept_runs=tuple(run.key for run in build_exon_runs(view, mutant_layout)))
        mutant, unresolved = realize_splice_plan(view, mutant_layout, plan, ())
        assert unresolved == ()
    else:
        mutant = realize_exon_path(view, mutant_layout)
    assert mutant.start_codon_present
    assert mutant.protein_sequence == "MEPG"


@pytest.mark.parametrize("start", ["GTG", "ACG"])
@pytest.mark.parametrize("reverse,split", [(False, False), (True, False), (False, True), (True, True)])
def test_mapped_rna_product_uses_declared_start(start, reverse, split):
    reference = declared_orf(start, reverse, split)
    transcript = reference.transcript
    sequence = transcript.sequence[:6] + "G" + transcript.sequence[7:]
    product = RealizedTranscriptProduct(transcript, sequence, None, evidence={
        "base_origins": [(transcript.contig, p, transcript.strand) for p in transcript.positions],
        "rna_support": 10})
    annotation = annotate_orf_rna_product(reference, product)
    assert annotation.status == "conditional_translation"
    assert isinstance(annotation.effect, Substitution)
    assert annotation.effect.short_description == "p.K2E"
    assert annotation.product.protein_sequence == "MEPG"
    assert annotation.product.start_codon_present
    assert annotation.product.evidence["rna_support"] == 10
    assert annotation.product.reference_evidence["mutant_translation_evidence"] == "not_assessed"


@pytest.mark.parametrize("unsupported", ["initiation", "genetic_code"])
def test_unsupported_translation_instruction_retains_explicit_failure(unsupported):
    if unsupported == "initiation":
        reference = declared_orf("TAA")
        reason = "Unsupported declared initiation codon"
    else:
        reference = declared_orf("GTG")
        reference = replace(reference.record, translation_table=3).bind(reference.transcript)
        reason = "Unsupported explicit NCBI translation table: 3"
    variant = reference.transcript.variant(6, "A", "G")
    with pytest.raises(ValueError, match=reason):
        variant.effects(extra_orfs=[reference])
    effects = variant.effects(extra_orfs=[reference], raise_on_error=False)
    annotation = effects.orf_annotations[0]
    assert any(isinstance(e, Failure) and e.orf_annotation is annotation for e in effects)
    assert isinstance(annotation.reference, orfdata.ReferenceORF)
    assert reason in annotation.unresolved_reason
    assert annotation.outcomes == ()


@pytest.mark.parametrize("start", ["GTG", "ACG"])
def test_mapped_rna_with_changed_start_retains_conditional_start_loss(start):
    reference = declared_orf(start)
    transcript = reference.transcript
    product = RealizedTranscriptProduct(transcript, transcript.sequence[:3] + "GCG" + transcript.sequence[6:],
        None, evidence={"base_origins": [(transcript.contig, p, transcript.strand) for p in transcript.positions]})
    annotation = annotate_orf_rna_product(reference, product)
    assert annotation.reference.protein_sequence == "MKPG"
    assert annotation.status == "conditional_translation"
    assert isinstance(annotation.effect, StartLoss)
    assert annotation.product.protein_sequence == ""
    assert not annotation.product.start_codon_present
    assert annotation.product.reference_evidence["mutant_translation_evidence"] == "not_assessed"
