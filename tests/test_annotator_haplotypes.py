"""The selected annotator owns joint predictions through the collection API."""

import pytest
from pyensembl import cached_release

from varcode import (
    EffectCollection, GermlineContext, MutantTranscript, StructuralVariant, Variant,
    VariantCollection, get_annotator, use_annotator,
)
from varcode.effects.effect_classes import (
    ComplexSubstitution, Deletion, Failure, HaplotypeEffect, PrematureStop,
    Silent, Unresolved,
)


@pytest.fixture
def transcript():
    return cached_release(81).transcript_by_id("ENST00000003084")


@pytest.fixture
def variants():
    genome = cached_release(81)
    # Each allele alone is missense; together GTA becomes TAA.
    return (
        Variant("7", 117530914, "G", "T", genome),
        Variant("7", 117530915, "T", "A", genome),
    )


class Phase:
    phase_source = "test_phase"

    def __init__(self, answer=True, observed=None):
        self.answer = answer
        self.observed = observed or {}

    def in_cis(self, left, right, transcript=None):
        return self.answer

    def mutant_transcript(self, variant, transcript):
        return self.observed.get((variant, transcript.id))


def joint_effects(effects, transcript):
    return [e for e in effects if getattr(e, "transcript", None) == transcript
            and len(getattr(e, "variants", ())) > 1]


@pytest.mark.parametrize("name", ["fast", "protein_diff", "transcript_model"])
def test_selected_model_composes_cis_codons_and_keeps_individuals(name, variants, transcript):
    with use_annotator(name):
        effects = VariantCollection(variants).effects(phase_resolver=Phase())
    joint, = joint_effects(effects, transcript)
    expected_stop = (transcript.spliced_offset(variants[0].start)
                     - min(transcript.start_codon_spliced_offsets)) // 3
    assert joint.mutant_protein_sequence == transcript.protein_sequence[:expected_stop]
    assert type(joint) is (PrematureStop if name == "transcript_model" else HaplotypeEffect)
    assert joint.variants == variants
    assert joint.phase_source == "test_phase"
    assert joint.annotator == effects.annotator == name
    assert joint.annotator_version == effects.annotator_version
    for variant in variants:
        individual, = [e for e in effects if e.transcript == transcript
                       and e.variant == variant
                       and len(getattr(e, "variants", ())) < 2]
        assert len(individual.mutant_protein_sequence) == len(transcript.protein_sequence)


@pytest.mark.parametrize("name", ["fast", "transcript_model"])
@pytest.mark.parametrize("phase", [False, None])
def test_trans_and_unknown_phase_do_not_form_a_joint_prediction(name, phase, variants, transcript):
    effects = VariantCollection(variants).effects(
        annotator=name, phase_resolver=Phase(phase))
    assert not joint_effects(effects, transcript)


class Declining:
    name = "declining"
    version = "test-1"

    def annotate_on_transcript(self, variant, transcript):
        return NotImplemented


@pytest.mark.parametrize("explicit_hook", [False, True])
def test_declining_plugin_keeps_unresolved_group_without_default_prediction(
        explicit_hook, variants, transcript):
    annotator = Declining()
    if explicit_hook:
        annotator.annotate_haplotype = lambda *args, **kwargs: NotImplemented
    effects = VariantCollection(variants).effects(
        annotator=annotator, phase_resolver=Phase())
    joint, = joint_effects(effects, transcript)
    assert isinstance(joint, Unresolved)
    assert joint.mechanism == "unsupported_haplotype"
    assert "declining" in joint.reason
    assert joint.variants == variants
    assert joint.phase_source == "test_phase"
    assert joint.annotator == effects.annotator == "declining"
    assert joint.annotator_version == "test-1"
    assert not any(isinstance(e, HaplotypeEffect) for e in effects)


def test_plugin_receives_whole_group_and_original_context(variants, transcript):
    context = GermlineContext.from_variants([variants[0]])
    resolver = Phase()
    calls = []

    class Plugin(Declining):
        def annotate_haplotype(self, members, target, germline_ctx=None, phase_resolver=None):
            calls.append((members, target, germline_ctx, phase_resolver))
            return Silent(members[0], target, aa_pos=0, aa_ref="M")

    effects = VariantCollection(variants).effects(
        annotator=Plugin(), phase_resolver=resolver, germline=context)
    joint, = joint_effects(effects, transcript)
    assert type(joint) is Silent
    assert (variants, transcript, context, resolver) in calls


@pytest.mark.parametrize("name", ["fast", "protein_diff"])
def test_default_does_not_discard_joint_patient_baseline(name, variants, transcript):
    context = GermlineContext.from_variants([
        Variant("7", 117531101, "T", "C", cached_release(81))])
    effects = VariantCollection(variants).effects(
        annotator=name, germline=context, phase_resolver=Phase())
    joint, = joint_effects(effects, transcript)
    assert isinstance(joint, Unresolved)
    assert "germline context" in joint.reason


def test_model_uses_germline_near_every_member(transcript):
    genome = cached_release(81)
    members = (Variant("7", 117531100, "T", "A", genome),
               Variant("7", 117534280, "T", "G", genome))
    position = 117534281
    offset = transcript.spliced_offset(position)
    ref = transcript.sequence[offset]
    alt = "A" if ref != "A" else "C"
    germline = Variant("7", position, ref, alt, genome)
    effects = VariantCollection(members).effects(
        annotator="transcript_model", phase_resolver=Phase(),
        germline=GermlineContext.from_variants([germline]))
    joint, = joint_effects(effects, transcript)
    for candidate in joint.candidates:
        for outcome in candidate.outcomes:
            assert outcome.baseline.cdna_sequence[offset] == alt
            assert (germline, "cis") in outcome.hypothesis.phase
            for variant in members:
                index = transcript.spliced_offset(variant.start)
                assert outcome.baseline.cdna_sequence[index] == variant.ref
                assert outcome.mutant.cdna_sequence[index] == variant.alt


@pytest.mark.parametrize("name", ["fast", "protein_diff", "transcript_model"])
def test_rna_observation_is_preserved_with_its_provenance(name, variants, transcript):
    observed = MutantTranscript(
        reference_transcript=transcript, cdna_sequence="ATGGCTTAA",
        mutant_protein_sequence="MA", annotator_name="rna_source")
    resolver = Phase(observed={(variants[1], transcript.id): observed})
    effects = VariantCollection(variants).effects(annotator=name, phase_resolver=resolver)
    joint, = joint_effects(effects, transcript)
    if name == "transcript_model":
        assert joint.observed_mutant_transcripts == ((variants[1], observed),)
        assert joint.mutant_protein_sequence != observed.mutant_protein_sequence
    else:
        assert joint.mutant_transcript is observed
        assert joint.mutant_protein_sequence == "MA"
    assert observed.annotator_name == "rna_source"


@pytest.mark.parametrize("name", ["fast", "protein_diff", "transcript_model"])
@pytest.mark.parametrize("kind", ["splice", "sv"])
def test_splice_and_sv_groups_are_delegated_and_never_silently_dropped(name, kind, transcript):
    genome = cached_release(81)
    coding = Variant("7", 117531100, "T", "A", genome)
    if kind == "splice":
        second = Variant("7", 117531115, "G", "A", genome)
    else:
        exon = transcript.exons[4]
        second = StructuralVariant(
            "7", exon.start - 300, "DEL", end=exon.end + 300, genome=genome)
    effects = VariantCollection([coding, second]).effects(
        annotator=name, phase_resolver=Phase())
    joint, = joint_effects(effects, transcript)
    assert joint.variants == (coding, second)
    if name == "transcript_model":
        assert type(joint) is (Deletion if kind == "splice" else ComplexSubstitution)
        assert joint.candidates
        assert len(joint.mutant_protein_sequence) < len(transcript.protein_sequence)
    else:
        assert isinstance(joint, Unresolved)
        assert joint.mechanism == "unsupported_haplotype"


@pytest.mark.parametrize("name", ["fast", "protein_diff"])
def test_conflicting_point_edits_remain_an_unresolved_group(name, variants, transcript):
    conflict = Variant("7", variants[0].start, "G", "A", cached_release(81))
    members = (variants[0], conflict)
    effects = VariantCollection(members).effects(annotator=name, phase_resolver=Phase())
    joint, = joint_effects(effects, transcript)
    assert isinstance(joint, Unresolved)
    assert joint.variants == members


def test_joint_error_policy_and_invalid_return(variants, transcript):
    class Broken(Declining):
        def annotate_haplotype(self, *args, **kwargs):
            raise ValueError("joint prediction failed")

    collection = VariantCollection(variants)
    with pytest.raises(ValueError, match="joint prediction failed"):
        collection.effects(annotator=Broken(), phase_resolver=Phase(), raise_on_error=True)
    effects = collection.effects(
        annotator=Broken(), phase_resolver=Phase(), raise_on_error=False)
    joint, = joint_effects(effects, transcript)
    assert isinstance(joint, Failure)
    assert joint.variants == variants
    invalid = Declining()
    invalid.annotate_haplotype = lambda *args, **kwargs: None
    with pytest.raises(TypeError, match="MutationEffect or NotImplemented"):
        collection.effects(annotator=invalid, phase_resolver=Phase(), raise_on_error=False)


def test_model_conflicting_edits_follow_error_policy(variants, transcript):
    conflict = Variant("7", variants[0].start, "G", "A", cached_release(81))
    collection = VariantCollection([variants[0], conflict])
    with pytest.raises(ValueError, match="does not match realized layout"):
        collection.effects(annotator="transcript_model", phase_resolver=Phase())
    effects = collection.effects(
        annotator="transcript_model", phase_resolver=Phase(), raise_on_error=False)
    joint, = joint_effects(effects, transcript)
    assert isinstance(joint, Failure)
    assert set(joint.variants) == {variants[0], conflict}


@pytest.mark.parametrize("name", ["fast", "protein_diff", "transcript_model"])
def test_reverse_strand_group_keeps_both_edits(name):
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000357654")
    variants = (Variant("17", 43082563, "T", "A", genome),
                Variant("17", 43082570, "C", "A", genome))
    effects = VariantCollection(variants).effects(
        annotator=name, phase_resolver=Phase(), germline=GermlineContext.empty())
    joint, = joint_effects(effects, transcript)
    cdna = (joint.candidates[0].outcomes[0].mutant.cdna_sequence
            if name == "transcript_model" else joint.mutant_transcript.cdna_sequence)
    assert [i for i, (ref, alt) in enumerate(zip(transcript.sequence, cdna))
            if ref != alt] == [4309, 4316]
    assert cdna[4309] == cdna[4316] == "T"


@pytest.mark.parametrize("name", ["fast", "transcript_model", "declining"])
def test_joint_membership_and_provenance_survive_json(name, variants, transcript):
    annotator = Declining() if name == "declining" else name
    observed = MutantTranscript(
        reference_transcript=transcript, cdna_sequence="ATGGCTTAA",
        mutant_protein_sequence="MA", annotator_name="rna_source")
    effects = VariantCollection(variants).effects(
        annotator=annotator,
        phase_resolver=Phase(observed={(variants[1], transcript.id): observed}))
    original, = joint_effects(effects, transcript)
    restored = EffectCollection.from_json(effects.to_json())
    joint, = joint_effects(restored, transcript)
    assert joint.variants == original.variants
    assert joint.phase_source == original.phase_source
    assert joint.annotator == restored.annotator == name
    assert joint.annotator_version == original.annotator_version
    if name == "declining":
        assert joint.reason == original.reason
    elif name == "transcript_model":
        (variant, model), = joint.observed_mutant_transcripts
        assert variant == variants[1]
        assert model.mutant_protein_sequence == "MA"
        assert model.annotator_name == "rna_source"
    else:
        assert joint.mutant_protein_sequence == "MA"
