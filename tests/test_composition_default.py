"""The six original #421 failures compose through the ordinary public API."""

from types import SimpleNamespace

import pytest
from pyensembl import cached_release

from varcode import EffectCollection, GermlineContext, StructuralVariant, Variant, VariantCollection


@pytest.fixture
def transcript():
    return cached_release(81).transcript_by_id("ENST00000003084")


def outcomes(effect):
    return [outcome for candidate in effect.candidates for outcome in candidate.outcomes]


def canonical_baseline(outcome):
    return all(option.mechanism == "normal"
               for _, option in outcome.hypothesis.baseline_splice_plan.choices)


def delete_exons(transcript, sequence, first, last):
    start = transcript.spliced_offset(transcript.exons[first - 1].start)
    end = transcript.spliced_offset(transcript.exons[last - 1].end) + 1
    return sequence[:start] + sequence[end:]


def span_deletion(transcript):
    return StructuralVariant(
        "7", transcript.exons[4].start - 100,
        end=transcript.exons[5].end + 100, sv_type="DEL", genome=transcript.genome)


@pytest.mark.parametrize("kind", ["donor", "deletion"])
def test_default_splice_and_deletion_keep_patient_baselines(transcript, kind):
    somatic = (Variant("7", 117531115, "G", "A", transcript.genome)
               if kind == "donor" else span_deletion(transcript))
    position = 117531113 if kind == "donor" else transcript.exons[4].end - 1
    reference_base = transcript.sequence[transcript.spliced_offset(position)]
    alternate = "C" if reference_base != "C" else "A"
    germline = Variant("7", position, reference_base, alternate, transcript.genome)
    context = GermlineContext.from_variants([germline])
    result = somatic.effect_on_transcript(transcript, germline=context)
    offset = transcript.spliced_offset(germline.start)
    inherited = transcript.sequence[:offset] + alternate + transcript.sequence[offset + 1:]
    relevant = [o for o in outcomes(result) if canonical_baseline(o)]
    assert {o.baseline.cdna_sequence for o in relevant} == {transcript.sequence, inherited}
    expected_mutant = delete_exons(transcript, transcript.sequence, 4 if kind == "donor" else 5, 4 if kind == "donor" else 6)
    skipped = [o for o in relevant if o.mutant.cdna_sequence == expected_mutant]
    assert {o.baseline.cdna_sequence for o in skipped} == {transcript.sequence, inherited}
    assert len({o.effect.original_protein_sequence for o in skipped}) == 2
    restored = EffectCollection.from_json(EffectCollection([result]).to_json())[0]
    assert [o.hypothesis.baseline_splice_plan for o in outcomes(restored)] == [
        o.hypothesis.baseline_splice_plan for o in outcomes(result)]


@pytest.mark.parametrize("kind", ["donor", "deletion"])
def test_default_cis_splice_and_deletion_groups_keep_coding_edits(transcript, kind):
    # Coding edit outside the potentially deleted exons remains in every
    # resolved product retaining its exon.
    position = transcript.exons[6].start + 10
    offset = transcript.spliced_offset(position)
    ref = transcript.sequence[offset]
    alt = "A" if ref != "A" else "C"
    coding = Variant("7", position, ref, alt, transcript.genome)
    structural = (Variant("7", 117531115, "G", "A", transcript.genome)
                  if kind == "donor" else span_deletion(transcript))
    members = (structural, coding)
    effects = VariantCollection(members).effects(
        phase_resolver=SimpleNamespace(in_cis=lambda *args, **kwargs: True),
        transcript_filter=lambda t: t.id == transcript.id)
    joint, = [e for e in effects if len(e.variants) == 2]
    assert joint.variants == members
    edited = transcript.sequence[:offset] + alt + transcript.sequence[offset + 1:]
    expected = delete_exons(transcript, edited, 4 if kind == "donor" else 5, 4 if kind == "donor" else 6)
    assert expected in {o.mutant.cdna_sequence for o in outcomes(joint)}
    assert all(o.baseline.cdna_sequence == transcript.sequence for o in outcomes(joint))
    assert len(effects) == 3  # both individuals and their joint comparison


@pytest.mark.parametrize("start,end,strength", [
    (117531117, 117531125, "weak_intronic"),
    (117531114, 117531116, "lost"),
])
def test_default_deletion_crossing_splice_window_keeps_mechanisms(transcript, start, end, strength):
    variant = StructuralVariant("7", start, end=end, sv_type="DEL", genome=transcript.genome)
    result = variant.effect_on_transcript(transcript)
    choices = [option for o in outcomes(result) for key, option in o.hypothesis.splice_plan.choices
               if key.side == "donor" and key.exon_number == 4]
    assert {option.mechanism for option in choices} == {
        "normal", "exon_skip", "intron_retention", "cryptic"}
    assert {option.evidence["strength"] for option in choices} == {strength}
    assert delete_exons(transcript, transcript.sequence, 4, 4) in {
        o.mutant.cdna_sequence for o in outcomes(result)}


def test_default_filtering_and_priority_preserve_alternatives_to_normal_splicing(transcript):
    from varcode import effect_priority
    from varcode.effects import Deletion, Intronic

    variant = Variant("7", 117531117, "A", "C", transcript.genome)
    effect = variant.effect_on_transcript(transcript)
    assert isinstance(effect, Intronic)  # preferred normal-splicing comparison
    assert effect.modifies_protein_sequence is False
    changing = [c for c in effect.candidates if isinstance(c.effect, Deletion)]
    assert changing
    collection = EffectCollection([effect])
    assert list(collection.drop_silent_and_noncoding()) == [effect]
    assert list(collection.drop_silent_and_noncoding(keep_unresolved=False)) == [effect]
    assert list(collection.filter_by_effect_priority(Deletion)) == [effect]
    assert effect.highest_priority_candidate in changing
    missense = Variant("7", 117531100, "T", "A", transcript.genome).effect_on_transcript(transcript)
    assert effect_priority(effect) > effect_priority(missense)
    assert EffectCollection([missense, effect]).top_priority_effect() is effect


def test_original_donor_case_keeps_germline_in_last_codon(transcript):
    donor = Variant("7", 117531115, "G", "A", transcript.genome)
    position = min(transcript.stop_codon_positions) - 3
    offset = transcript.spliced_offset(position)
    ref = transcript.sequence[offset]
    alt = "A" if ref != "A" else "C"
    inherited = Variant("7", position, ref, alt, transcript.genome)
    result = donor.effect_on_transcript(
        transcript, germline=GermlineContext.from_variants([inherited]))
    patient = transcript.sequence[:offset] + alt + transcript.sequence[offset + 1:]
    expected = {
        (baseline, delete_exons(transcript, baseline, 4, 4))
        for baseline in (transcript.sequence, patient)}
    actual = {(outcome.baseline.cdna_sequence, outcome.mutant.cdna_sequence)
              for outcome in outcomes(result)}
    assert expected <= actual
    skipped = [outcome for outcome in outcomes(result)
               if (outcome.baseline.cdna_sequence, outcome.mutant.cdna_sequence) in expected]
    assert len({outcome.mutant.protein_sequence for outcome in skipped}) == 2
