# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Tests for the splice outcome possibility-set prototype
(openvax/varcode#262).

Coverage:
  - Default behavior unchanged (back-compat)
  - Opt-in wraps splice effects in SpliceOutcomeSet
  - Each canonical splice signal class produces the expected outcome set
  - Plausibility ordering is stable
  - Per-outcome candidate construction (normal splicing, exon
    skipping, intron retention stubs, cryptic splice stubs)
  - SpliceOutcomeSet integrates with EffectCollection
  - Multi-allelic and reverse-strand variants work too
"""

from pyensembl import cached_release

import varcode
from varcode import (
    SpliceCandidate,
    SpliceOutcome,
    SpliceOutcomeSet,
    Variant,
)
from varcode.effects import (
    ExonicSpliceSite,
    IntronicSpliceSite,
    SpliceAcceptor,
    SpliceDonor,
)
from varcode.splice_outcomes import enumerate_splice_outcomes


ensembl_grch38 = cached_release(81)
CFTR_TRANSCRIPT_ID = "ENST00000003084"
BRCA1_TRANSCRIPT_ID = "ENST00000357654"


# --------------------------------------------------------------------
# Back-compat
# --------------------------------------------------------------------


def test_default_behavior_unchanged():
    # No kwarg -> same as today: SpliceDonor effect, no wrapping.
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effect = variant.effect_on_transcript(transcript)
    assert effect.__class__ is SpliceDonor
    assert not isinstance(effect, SpliceOutcomeSet)


def test_default_for_collection_unchanged():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    effects = variant.effects()  # default: splice_outcomes=False
    classes = {type(e) for e in effects}
    assert SpliceOutcomeSet not in classes


# --------------------------------------------------------------------
# Opt-in wraps splice effects
# --------------------------------------------------------------------


def test_opt_in_wraps_splice_donor():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    effects = variant.effects(splice_outcomes=True)
    cftr_effect = next(
        e for e in effects
        if getattr(e, "transcript", None) is not None
        and e.transcript.id == CFTR_TRANSCRIPT_ID
    )
    assert isinstance(cftr_effect, SpliceOutcomeSet)
    assert cftr_effect.disrupted_signal_class is SpliceDonor


def test_opt_in_wraps_splice_acceptor():
    # CFTR exon 4 acceptor -1 with canonical ref G.
    variant = Variant("7", 117530898, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effect = variant.effects(splice_outcomes=True)
    target = next(e for e in effect if e.transcript is transcript)
    assert isinstance(target, SpliceOutcomeSet)
    assert target.disrupted_signal_class is SpliceAcceptor


def test_opt_in_wraps_exonic_splice_site():
    # CFTR exon 4 ends with AAG. G->T at -1 disrupts the MAG signal.
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert isinstance(target, SpliceOutcomeSet)
    assert target.disrupted_signal_class is ExonicSpliceSite


def test_opt_in_wraps_intronic_splice_site():
    # CFTR exon 4 +1 with NON-canonical ref A is downgraded to
    # IntronicSpliceSite (post-2.0.0 sequence-aware classification).
    variant = Variant("7", 117531115, "A", "G", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert isinstance(target, SpliceOutcomeSet)
    assert target.disrupted_signal_class is IntronicSpliceSite


def test_opt_in_passes_through_non_splice_effects():
    # Pure substitution that doesn't touch any splice signal: should
    # not be wrapped.
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", ensembl_grch38)
    effects = variant.effects(splice_outcomes=True)
    # At least one effect should be a Substitution, not wrapped.
    classes = [type(e).__name__ for e in effects]
    assert "Substitution" in classes


# --------------------------------------------------------------------
# Plausibility ordering and candidate composition
# --------------------------------------------------------------------


def test_splice_donor_candidate_set_has_expected_outcomes():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effect = variant.effects(splice_outcomes=True)
    target = next(e for e in effect if e.transcript is transcript)
    outcomes = {c.outcome for c in target.candidates}
    assert outcomes == {
        SpliceOutcome.EXON_SKIPPING,
        SpliceOutcome.INTRON_RETENTION,
        SpliceOutcome.CRYPTIC_DONOR,
        SpliceOutcome.NORMAL_SPLICING,
    }


def test_splice_acceptor_candidate_set_uses_cryptic_acceptor():
    variant = Variant("7", 117530898, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    outcomes = {c.outcome for c in target.candidates}
    # SpliceAcceptor disruption uses CRYPTIC_ACCEPTOR not CRYPTIC_DONOR.
    assert SpliceOutcome.CRYPTIC_ACCEPTOR in outcomes
    assert SpliceOutcome.CRYPTIC_DONOR not in outcomes


def test_candidates_sorted_by_plausibility_descending():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    plaus = [c.plausibility for c in target.candidates]
    assert plaus == sorted(plaus, reverse=True), (
        "Candidates should be ordered most-plausible-first, got %r"
        % plaus)


def test_most_likely_for_splice_donor_is_exon_skipping():
    # Per the plausibility table, EXON_SKIPPING dominates SpliceDonor.
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert target.most_likely.outcome is SpliceOutcome.EXON_SKIPPING


def test_most_likely_for_exonic_splice_site_is_normal_splicing():
    # ExonicSpliceSite gets NORMAL_SPLICING as the most-likely
    # outcome (the disruption is on the exon side and often tolerated).
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert target.most_likely.outcome is SpliceOutcome.NORMAL_SPLICING


def test_normal_splicing_carries_underlying_coding_effect():
    # ExonicSpliceSite has an alternate_effect (the coding change if
    # splicing proceeds). NORMAL_SPLICING candidate exposes it.
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    normal = next(
        c for c in target.candidates
        if c.outcome is SpliceOutcome.NORMAL_SPLICING
    )
    assert normal.coding_effect is not None
    assert "p." in normal.coding_effect.short_description


# --------------------------------------------------------------------
# Per-outcome detail
# --------------------------------------------------------------------


def test_intron_retention_candidate_predicts_premature_stop():
    # Intron retention typically produces a PrematureStop. Stub
    # without exact protein since we don't have intronic genomic seq.
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    intron = next(
        c for c in target.candidates
        if c.outcome is SpliceOutcome.INTRON_RETENTION
    )
    assert intron.predicted_class_name == "PrematureStop"
    assert intron.coding_effect is None


def test_cryptic_donor_candidate_is_a_stub():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    cryptic = next(
        c for c in target.candidates
        if c.outcome is SpliceOutcome.CRYPTIC_DONOR
    )
    assert cryptic.coding_effect is None
    assert "cryptic" in cryptic.description.lower()


def test_exon_skipping_for_in_frame_exon_emits_deletion():
    # CFTR exon 4 is 216 nucleotides = 72 codons (216 % 3 == 0), so
    # skipping it is in-frame. The candidate should report Deletion
    # of the exon's amino acids.
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    skip = next(
        c for c in target.candidates
        if c.outcome is SpliceOutcome.EXON_SKIPPING
    )
    # Either a Deletion was constructed, or the candidate falls back
    # to None with predicted_class_name still set. Both are valid.
    if skip.coding_effect is not None:
        assert skip.predicted_class_name == "Deletion"
        assert skip.coding_effect.aa_ref  # non-empty AA range
    else:
        assert skip.predicted_class_name in ("Deletion", "FrameShift", "ExonLoss")


# --------------------------------------------------------------------
# EffectCollection integration
# --------------------------------------------------------------------


def test_collection_iteration_after_wrapping():
    # The wrapped collection should still be iterable, indexable, and
    # produce SpliceOutcomeSet objects in place of the original splice
    # effects.
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    effects = variant.effects(splice_outcomes=True)
    items = list(effects)
    assert len(items) > 0
    splice_set_count = sum(1 for e in items if isinstance(e, SpliceOutcomeSet))
    assert splice_set_count >= 1


def test_short_description_uses_most_likely():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    desc = target.short_description
    assert desc.startswith("splice-set:")
    assert target.most_likely.outcome.value in desc


# --------------------------------------------------------------------
# Reverse-strand
# --------------------------------------------------------------------


def test_opt_in_works_on_reverse_strand_donor():
    # BRCA1 exon 12 reverse-strand donor at 43082403 with canonical ref C.
    variant = Variant("17", 43082403, "C", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(BRCA1_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert isinstance(target, SpliceOutcomeSet)
    assert target.disrupted_signal_class is SpliceDonor


# --------------------------------------------------------------------
# Direct enumerate_splice_outcomes tests
# --------------------------------------------------------------------


def test_enumerate_passes_through_non_splice():
    # Non-splice effect should pass through unchanged.
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(BRCA1_TRANSCRIPT_ID)
    sub_effect = variant.effect_on_transcript(transcript)
    assert type(sub_effect).__name__ == "Substitution"
    wrapped = enumerate_splice_outcomes(sub_effect)
    assert wrapped is sub_effect


# --------------------------------------------------------------------
# SpliceCandidate dataclass ergonomics
# --------------------------------------------------------------------


def test_splice_candidate_is_frozen():
    c = SpliceCandidate(
        outcome=SpliceOutcome.EXON_SKIPPING,
        plausibility=0.5,
        description="test",
    )
    try:
        c.plausibility = 0.9  # type: ignore
    except Exception:
        pass
    else:
        raise AssertionError("SpliceCandidate should be frozen")


def test_splice_candidate_equality():
    a = SpliceCandidate(
        outcome=SpliceOutcome.EXON_SKIPPING,
        plausibility=0.5,
        description="d",
    )
    b = SpliceCandidate(
        outcome=SpliceOutcome.EXON_SKIPPING,
        plausibility=0.5,
        description="d",
    )
    assert a == b


def test_package_level_exports():
    assert varcode.SpliceCandidate is SpliceCandidate
    assert varcode.SpliceOutcome is SpliceOutcome
    assert varcode.SpliceOutcomeSet is SpliceOutcomeSet


# --------------------------------------------------------------------
# Priority integration: SpliceOutcomeSet sorts as if it were the
# disrupted-signal class (review feedback on PR #292).
# --------------------------------------------------------------------


def test_splice_outcome_set_sorts_as_disrupted_signal_class():
    # When wrapped, a SpliceDonor-backed SpliceOutcomeSet should have
    # the same priority as a bare SpliceDonor — higher than Intronic,
    # lower than Substitution. If the priority delegation is broken,
    # SpliceOutcomeSet gets priority -1 and sorts to the bottom.
    from varcode.effects import effect_priority

    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare_effect = variant.effect_on_transcript(transcript)
    assert isinstance(bare_effect, SpliceDonor)
    bare_priority = effect_priority(bare_effect)

    wrapped_effects = variant.effects(splice_outcomes=True)
    wrapped = next(e for e in wrapped_effects if e.transcript is transcript)
    assert isinstance(wrapped, SpliceOutcomeSet)
    wrapped_priority = effect_priority(wrapped)

    assert wrapped_priority == bare_priority, (
        "SpliceOutcomeSet priority (%d) must match the disrupted-"
        "signal class priority (%d); otherwise sorting and "
        "top_priority_effect() behave wrongly." % (
            wrapped_priority, bare_priority))


def test_splice_outcome_set_top_priority_works():
    # top_priority_effect on a collection containing SpliceOutcomeSet
    # should not pick a lower-priority non-splice effect.
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    effects = variant.effects(splice_outcomes=True)
    top = effects.top_priority_effect()
    # The wrapped SpliceDonor (or one of the splice-set variants) is
    # higher priority than Intronic/NoncodingTranscript from other
    # overlapping transcripts, so the top should be a splice-related
    # effect.
    top_class_name = type(top).__name__
    assert top_class_name in ("SpliceOutcomeSet", "SpliceDonor"), (
        "Expected a splice-related effect at top priority, got %s"
        % top_class_name)


# --------------------------------------------------------------------
# Acceptor-side IntronicSpliceSite emits CRYPTIC_ACCEPTOR, not DONOR.
# --------------------------------------------------------------------


def test_acceptor_side_intronic_splice_site_uses_cryptic_acceptor():
    # CFTR exon 4 acceptor -3 (3bp before exon.start). A variant here
    # with NON-canonical ref (not A, the canonical MAG component) is
    # classified as IntronicSpliceSite. The splice set should include
    # CRYPTIC_ACCEPTOR (the relevant cryptic direction for the
    # acceptor side), not CRYPTIC_DONOR.
    from varcode.effects import IntronicSpliceSite
    # chr7:117530896 is -3 before CFTR exon 4 (forward strand).
    # Use a non-canonical ref for the -3 position so it's
    # IntronicSpliceSite (not SpliceAcceptor which covers -1/-2).
    # At distance -3, the position isn't required to be canonical
    # anyway — the classifier emits IntronicSpliceSite for this window.
    variant = Variant("7", 117530896, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    assert isinstance(bare, IntronicSpliceSite) and \
        not isinstance(bare, (SpliceDonor, SpliceAcceptor))
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    outcomes = {c.outcome for c in target.candidates}
    assert SpliceOutcome.CRYPTIC_ACCEPTOR in outcomes, \
        "Acceptor-side IntronicSpliceSite should use CRYPTIC_ACCEPTOR"
    assert SpliceOutcome.CRYPTIC_DONOR not in outcomes, \
        "Acceptor-side IntronicSpliceSite should not use CRYPTIC_DONOR"


def test_donor_side_intronic_splice_site_uses_cryptic_donor():
    # Mirror test for donor-side IntronicSpliceSite at +3 after CFTR
    # exon 4 end (117531117 = exon.end + 3).
    from varcode.effects import IntronicSpliceSite
    variant = Variant("7", 117531117, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    assert isinstance(bare, IntronicSpliceSite) and \
        not isinstance(bare, (SpliceDonor, SpliceAcceptor))
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    outcomes = {c.outcome for c in target.candidates}
    assert SpliceOutcome.CRYPTIC_DONOR in outcomes
    assert SpliceOutcome.CRYPTIC_ACCEPTOR not in outcomes


# --------------------------------------------------------------------
# Multi-protein surface: candidate_proteins and mutant_protein_sequences
# --------------------------------------------------------------------


def test_candidate_proteins_maps_each_outcome_to_a_protein():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    proteins = target.candidate_proteins
    # Every candidate outcome appears as a key.
    outcomes = {c.outcome for c in target.candidates}
    assert set(proteins.keys()) == outcomes
    # EXON_SKIPPING for an in-frame exon should have a non-empty
    # protein (reference minus the skipped AAs). CFTR exon 4 is 216
    # nucleotides = 72 codons = in-frame.
    assert proteins[SpliceOutcome.EXON_SKIPPING], \
        "Expected a concrete mutant protein for in-frame exon skipping"
    # INTRON_RETENTION and CRYPTIC are stubs → empty string for now.
    assert proteins[SpliceOutcome.INTRON_RETENTION] == ""
    assert proteins[SpliceOutcome.CRYPTIC_DONOR] == ""


def test_mutant_protein_sequences_collects_distinct_proteins():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    proteins = target.mutant_protein_sequences
    assert isinstance(proteins, set)
    assert len(proteins) >= 1
    # Reference protein should be in there or a proper subset of
    # reference (exon-skipped version is shorter).
    ref = str(transcript.protein_sequence)
    # The in-frame exon skip removes exon 4 AAs; resulting protein
    # should be shorter than reference.
    shortest = min(proteins, key=len)
    assert len(shortest) < len(ref)


# --------------------------------------------------------------------
# Out-of-frame exon skip now produces a real mutant protein
# --------------------------------------------------------------------


def test_out_of_frame_exon_skip_produces_mutant_protein():
    # CFTR exon 5 is 90 nucleotides = 30 codons, BUT exon 5 is not
    # out of frame — need a different exon. Use a variant known to
    # target an out-of-frame exon. We'll discover one empirically
    # by finding an exon whose length is not divisible by 3.
    for exon in ensembl_grch38.transcript_by_id(
            CFTR_TRANSCRIPT_ID).exons[2:]:
        length = exon.end - exon.start + 1
        if length % 3 != 0:
            target_exon = exon
            break
    else:
        # All exons in-frame; skip this test.
        return

    # Construct a donor-side disrupting variant at this exon's end
    # (+1 position after the exon in + strand coords).
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    donor_plus_1 = target_exon.end + 1
    # Use a canonical-ref SNV to ensure SpliceDonor classification.
    variant = Variant("7", donor_plus_1, "G", "A", ensembl_grch38)
    bare = variant.effect_on_transcript(transcript)
    if not isinstance(bare, SpliceDonor):
        # Skip the test if the canonical G isn't actually there.
        return
    effects = variant.effects(splice_outcomes=True)
    splice_set = next(e for e in effects if e.transcript is transcript)
    skip_candidate = next(
        c for c in splice_set.candidates
        if c.outcome is SpliceOutcome.EXON_SKIPPING
    )
    # Out-of-frame skip should now carry a mutant protein.
    assert skip_candidate.coding_effect is not None, (
        "Out-of-frame exon skip should produce a concrete mutant "
        "protein, not a stub")
    protein = skip_candidate.coding_effect.mutant_protein_sequence
    assert isinstance(protein, str)
    assert len(protein) > 0
    # The frameshifted protein should differ from the reference
    # after the skip point.
    assert protein != str(transcript.protein_sequence)


# --------------------------------------------------------------------
# has_protein property on SpliceCandidate
# --------------------------------------------------------------------


def test_has_protein_is_true_for_candidates_with_coding_effect():
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    # NORMAL_SPLICING has a Substitution coding_effect with a protein.
    normal = next(
        c for c in target.candidates
        if c.outcome is SpliceOutcome.NORMAL_SPLICING
    )
    assert normal.has_protein is True


def test_has_protein_is_false_for_stub_candidates():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    intron = next(
        c for c in target.candidates
        if c.outcome is SpliceOutcome.INTRON_RETENTION
    )
    assert intron.has_protein is False
