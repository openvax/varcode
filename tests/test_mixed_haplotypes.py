"""Mixed inherited/somatic groups use one patient baseline (#500)."""

import pytest
from pyensembl import cached_release

from varcode import (
    EffectCollection, GermlineAlleleOverlap, GermlineContext, HypothesisLimit,
    TranscriptModelEffectAnnotator, Variant, VariantCollection,
    predict_transcript_model_effect,
)


class Phase:
    phase_source = "mixed_haplotype_test"

    def __init__(self, *links):
        self.links = {frozenset((left, right)): answer
                      for left, right, answer in links}

    def in_cis(self, left, right, transcript=None):
        return self.links.get(frozenset((left, right)))


@pytest.fixture
def transcript():
    return cached_release(81).transcript_by_id("ENST00000003084")


@pytest.fixture(params=["C", "TGCC", ""], ids=["snv", "insertion", "deletion"])
def alleles(request, transcript):
    novel = Variant("7", 117531100, "T", "A", genome=81)
    inherited = Variant("7", 117531101, "T", request.param, genome=81)
    offset = transcript.spliced_offset(117531101)
    baseline = (transcript.sequence[:offset] + request.param
                + transcript.sequence[offset + 1:])
    # The novel edit precedes the inherited edit, so this offset is unchanged.
    offset = transcript.spliced_offset(novel.start)
    mutant = baseline[:offset] + "A" + baseline[offset + 1:]
    return novel, inherited, baseline, mutant


def outcomes(effect):
    return [outcome for candidate in effect.candidates
            for outcome in candidate.outcomes]


@pytest.mark.parametrize("inherited_first", [False, True])
@pytest.mark.parametrize("resolved", [False, True])
def test_mixed_edits_are_applied_once(alleles, transcript, inherited_first, resolved):
    novel, inherited, baseline, mutant = alleles
    members = ((inherited, novel, inherited) if inherited_first
               else (novel, inherited, inherited))
    resolver = Phase((novel, inherited, True)) if resolved else None
    effect = predict_transcript_model_effect(
        members, transcript, germline_variants=(inherited, inherited),
        phase_resolver=resolver, max_phase_hypotheses=1)

    assert effect.variant == novel
    assert effect.variants == members
    outcome, = outcomes(effect)
    assert outcome.baseline.cdna_sequence == baseline
    assert outcome.mutant.cdna_sequence == mutant
    assert outcome.hypothesis.phase == ((inherited, "cis"),)
    assert outcome.hypothesis.evidence["somatic_variants"] == (novel,)
    assert outcome.hypothesis.evidence["phase_state"] == "phased"
    assert outcome.hypothesis.evidence["phase_source"] == (
        resolver.phase_source if resolved else None)


@pytest.mark.parametrize("relation", [True, False, None])
def test_other_germline_phase_is_preserved(alleles, transcript, relation):
    novel, inherited, baseline, mutant = alleles
    position = 117531107
    offset = transcript.spliced_offset(position)
    ref = transcript.sequence[offset]
    other = Variant("7", position, ref, "A" if ref != "A" else "C", genome=81)
    # Evidence through the inherited member must propagate to the novel edit.
    resolver = Phase((inherited, other, relation))
    effect = predict_transcript_model_effect(
        (novel, inherited), transcript, germline_variants=(inherited, other),
        phase_resolver=resolver)

    realized = outcomes(effect)
    assert len(realized) == (2 if relation is None else 1)
    expected_sides = {"cis", "trans"} if relation is None else {
        "cis" if relation else "trans"}
    assert {dict(o.hypothesis.phase)[other] for o in realized} == expected_sides
    # Locate the other site after the inherited indel has shifted it.
    shifted = offset + len(baseline) - len(transcript.sequence)
    for outcome in realized:
        phase = dict(outcome.hypothesis.phase)
        assert phase[inherited] == "cis"
        base, mut = baseline, mutant
        if phase[other] == "cis":
            base = base[:shifted] + other.alt + base[shifted + 1:]
            mut = mut[:shifted] + other.alt + mut[shifted + 1:]
        assert outcome.baseline.cdna_sequence == base
        assert outcome.mutant.cdna_sequence == mut
        assert outcome.hypothesis.phase_probability == (0.5 if relation is None else 1)
        assert outcome.hypothesis.evidence["phase_state"] == (
            "unknown" if relation is None else "phased")


def test_haplotype_hook_keeps_repeated_inherited_alleles(alleles, transcript):
    novel, inherited, baseline, mutant = alleles
    members = (inherited, inherited, novel)
    effect = TranscriptModelEffectAnnotator().annotate_haplotype(
        members, transcript, germline_ctx=GermlineContext.from_variants([inherited]))
    assert effect.variant == novel
    assert effect.variants == members
    outcome, = outcomes(effect)
    assert outcome.baseline.cdna_sequence == baseline
    assert outcome.mutant.cdna_sequence == mutant


@pytest.mark.parametrize("alt", ["TGAA", ""], ids=["insertion", "deletion"])
def test_novel_indels_use_the_same_patient_baseline(alleles, transcript, alt):
    _, inherited, baseline, _ = alleles
    novel = Variant("7", 117531100, "T", alt, genome=81)
    effect = predict_transcript_model_effect(
        (inherited, novel, inherited), transcript, germline_variants=(inherited,))
    outcome, = outcomes(effect)
    offset = transcript.spliced_offset(117531100)
    assert outcome.baseline.cdna_sequence == baseline
    assert outcome.mutant.cdna_sequence == baseline[:offset] + alt + baseline[offset + 1:]


@pytest.mark.parametrize("relation", [True, False, None])
def test_collection_only_composes_supported_cis_group(alleles, transcript, relation):
    novel, inherited, baseline, mutant = alleles
    resolver = Phase((novel, inherited, relation))
    effects = VariantCollection([inherited, novel, inherited]).effects(
        annotator="transcript_model", phase_resolver=resolver,
        germline=GermlineContext.from_variants([inherited, inherited]))
    local = [e for e in effects if e.transcript == transcript]
    joints = [e for e in local if len(getattr(e, "variants", ())) > 1]
    assert any(isinstance(e, GermlineAlleleOverlap) for e in local)
    if relation is not True:
        assert not joints
        individual, = [e for e in local if e.variant == novel]
        assert {dict(o.hypothesis.phase)[inherited] for o in outcomes(individual)} == (
            {"trans"} if relation is False else {"cis", "trans"})
        return
    joint, = joints
    assert joint.variant == novel
    assert joint.phase_source == resolver.phase_source
    outcome, = outcomes(joint)
    assert outcome.baseline.cdna_sequence == baseline
    assert outcome.mutant.cdna_sequence == mutant
    restored = EffectCollection.from_json(EffectCollection([joint]).to_json())[0]
    assert restored.variants == joint.variants
    assert restored.variant == novel
    assert restored.phase_source == resolver.phase_source


def test_all_inherited_members_establish_no_new_allele(alleles, transcript):
    _, inherited, _, _ = alleles
    other = Variant("7", 117531100, "T", "A", genome=81)
    members = (inherited, other, inherited)
    effect = predict_transcript_model_effect(
        members, transcript, germline_variants=(inherited, other, inherited))
    assert isinstance(effect, GermlineAlleleOverlap)
    assert effect.variants == members
    assert effect.is_loh is None
    assert effect.modifies_coding_sequence is False
    assert effect.modifies_protein_sequence is False


def test_normalized_equivalent_alleles_share_baseline(transcript):
    novel = Variant("7", 117531100, "T", "A", genome=81)
    inherited = Variant("7", 117531101, "T", "TGCC", genome=81)
    unpadded = Variant("7", 117531101, "", "GCC", genome=81)
    effect = predict_transcript_model_effect(
        (novel, unpadded, inherited), transcript,
        germline_variants=(inherited, unpadded))
    outcome, = outcomes(effect)
    assert outcome.hypothesis.phase == ((inherited, "cis"),)
    assert outcome.hypothesis.evidence["somatic_variants"] == (novel,)
    assert len(outcome.baseline.cdna_sequence) == len(transcript.sequence) + 3
    assert len(outcome.mutant.cdna_sequence) == len(transcript.sequence) + 3


def test_mixed_haplotype_retains_phase_limit(transcript):
    novel = Variant("7", 117531100, "T", "A", genome=81)
    inherited = Variant("7", 117531101, "T", "C", genome=81)
    other = Variant("7", 117531102, "G", "A", genome=81)
    effect = predict_transcript_model_effect(
        (inherited, novel, inherited), transcript,
        germline_variants=(inherited, other), max_phase_hypotheses=1)
    assert isinstance(effect, HypothesisLimit)
    assert effect.variant == novel
    assert effect.phase.cis == (inherited,)
    assert effect.phase.unphased == (other,)
    assert effect.reference_effect is not None


def test_reverse_strand_mixed_haplotype():
    transcript = cached_release(81).transcript_by_id("ENST00000357654")
    inherited = Variant("17", 43082570, "C", "A", genome=81)
    novel = Variant("17", 43082563, "T", "A", genome=81)
    effect = predict_transcript_model_effect(
        (inherited, novel, inherited), transcript, germline_variants=(inherited,))
    outcome, = outcomes(effect)
    baseline = list(transcript.sequence)
    baseline[transcript.spliced_offset(inherited.start)] = "T"
    assert outcome.baseline.cdna_sequence == "".join(baseline)
    baseline[transcript.spliced_offset(novel.start)] = "T"
    assert outcome.mutant.cdna_sequence == "".join(baseline)
