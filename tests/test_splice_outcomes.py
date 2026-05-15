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
(openvax/varcode#262), updated for the unified MultiOutcomeEffect
shape (openvax/varcode#382): ``SpliceOutcomeSet.candidates`` is a
``tuple[EffectCandidate, ...]`` like every other ``MultiOutcomeEffect``
subclass, with ``candidate.evidence["splice_outcome"]`` carrying the
biological outcome and ``candidate.effect`` always a real
``MutationEffect`` (concrete coding effect or placeholder).
"""

import pytest
from pyensembl import cached_release

import varcode
from varcode import (
    EffectCandidate,
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


def _outcome_of(candidate):
    return candidate.evidence.get("splice_outcome")


def _candidate_for(splice_set, outcome):
    return next(
        c for c in splice_set.candidates if _outcome_of(c) is outcome)


# --------------------------------------------------------------------
# Back-compat
# --------------------------------------------------------------------


def test_default_behavior_unchanged():
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
    variant = Variant("7", 117530898, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert isinstance(target, SpliceOutcomeSet)
    assert target.disrupted_signal_class is SpliceAcceptor


def test_opt_in_wraps_exonic_splice_site():
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert isinstance(target, SpliceOutcomeSet)
    assert target.disrupted_signal_class is ExonicSpliceSite


def test_opt_in_wraps_intronic_splice_site():
    variant = Variant("7", 117531115, "A", "G", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert isinstance(target, SpliceOutcomeSet)
    assert target.disrupted_signal_class is IntronicSpliceSite


def test_opt_in_passes_through_non_splice_effects():
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", ensembl_grch38)
    effects = variant.effects(splice_outcomes=True)
    classes = [type(e).__name__ for e in effects]
    assert "Substitution" in classes


# --------------------------------------------------------------------
# Probability ordering and candidate composition
# --------------------------------------------------------------------


def test_splice_donor_candidate_set_has_expected_outcomes():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effect = variant.effects(splice_outcomes=True)
    target = next(e for e in effect if e.transcript is transcript)
    outcomes = {_outcome_of(c) for c in target.candidates}
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
    outcomes = {_outcome_of(c) for c in target.candidates}
    assert SpliceOutcome.CRYPTIC_ACCEPTOR in outcomes
    assert SpliceOutcome.CRYPTIC_DONOR not in outcomes


def test_candidates_sorted_by_probability_descending():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    probs = [c.probability for c in target.candidates]
    assert probs == sorted(probs, reverse=True), (
        "Candidates should be ordered most-plausible-first, got %r" % probs)


def test_most_likely_for_splice_donor_is_exon_skipping():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert _outcome_of(target.most_likely_candidate) is SpliceOutcome.EXON_SKIPPING


def test_most_likely_for_exonic_splice_site_is_normal_splicing():
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert _outcome_of(target.most_likely_candidate) is SpliceOutcome.NORMAL_SPLICING


def test_normal_splicing_carries_underlying_coding_effect():
    # ExonicSpliceSite has an alternate_effect (the coding change if
    # splicing proceeds). The NORMAL_SPLICING candidate's .effect IS
    # that coding change.
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    normal = _candidate_for(target, SpliceOutcome.NORMAL_SPLICING)
    assert "p." in normal.effect.short_description


# --------------------------------------------------------------------
# Per-outcome detail
# --------------------------------------------------------------------


def test_intron_retention_candidate_is_placeholder_without_provider():
    from varcode.effects.effect_classes import PredictedIntronRetention
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    intron = _candidate_for(target, SpliceOutcome.INTRON_RETENTION)
    assert isinstance(intron.effect, PredictedIntronRetention)
    assert intron.effect.mutant_protein_sequence is None
    # Placeholder builds are marked so consumers (and
    # ``alternate_effect``) can tell them apart from real effects of
    # the same class.
    assert intron.evidence.get("placeholder") is True


def test_resolved_candidate_does_not_carry_placeholder_marker():
    """When a real coding effect is built (e.g. EXON_SKIPPING with
    a computable protein on an in-frame exon), the
    ``placeholder`` evidence flag must NOT be set."""
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    skip = _candidate_for(target, SpliceOutcome.EXON_SKIPPING)
    # Exon 3 is in-frame so we get a concrete Deletion, not a
    # placeholder.
    assert type(skip.effect).__name__ == "Deletion"
    assert "placeholder" not in skip.evidence


def test_cryptic_donor_candidate_is_placeholder_without_provider():
    from varcode.effects.effect_classes import PredictedCrypticSpliceSite
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    cryptic = _candidate_for(target, SpliceOutcome.CRYPTIC_DONOR)
    assert isinstance(cryptic.effect, PredictedCrypticSpliceSite)
    assert "cryptic" in cryptic.evidence["description"].lower()


def test_exon_skipping_for_in_frame_exon_emits_deletion():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    skip = _candidate_for(target, SpliceOutcome.EXON_SKIPPING)
    # The inner effect is either a concrete coding effect (Deletion/
    # Substitution/ComplexSubstitution) or the ExonLoss placeholder.
    assert type(skip.effect).__name__ in (
        "Deletion", "FrameShift", "Substitution",
        "ComplexSubstitution", "ExonLoss")


# --------------------------------------------------------------------
# EffectCollection integration
# --------------------------------------------------------------------


def test_collection_iteration_after_wrapping():
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
    assert _outcome_of(target.most_likely_candidate).value in desc


# --------------------------------------------------------------------
# Reverse-strand
# --------------------------------------------------------------------


def test_opt_in_works_on_reverse_strand_donor():
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
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(BRCA1_TRANSCRIPT_ID)
    sub_effect = variant.effect_on_transcript(transcript)
    assert type(sub_effect).__name__ == "Substitution"
    wrapped = enumerate_splice_outcomes(sub_effect)
    assert wrapped is sub_effect


# --------------------------------------------------------------------
# Package-level exports
# --------------------------------------------------------------------


def test_package_level_exports():
    assert varcode.SpliceOutcome is SpliceOutcome
    assert varcode.SpliceOutcomeSet is SpliceOutcomeSet


# --------------------------------------------------------------------
# MultiOutcomeEffect protocol
# --------------------------------------------------------------------


def test_splice_outcome_set_is_a_multi_outcome_effect():
    from varcode import MultiOutcomeEffect
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert isinstance(target, MultiOutcomeEffect)
    assert hasattr(target, "candidates") and len(target.candidates) > 0
    assert all(isinstance(c, EffectCandidate) for c in target.candidates)
    assert target.most_likely_candidate is target.candidates[0]
    assert target.priority_class is target.disrupted_signal_class


def test_multi_outcome_effect_exported_at_package_level():
    from varcode import MultiOutcomeEffect
    from varcode import effects
    assert MultiOutcomeEffect is effects.MultiOutcomeEffect
    assert issubclass(SpliceOutcomeSet, MultiOutcomeEffect)


# --------------------------------------------------------------------
# Explicit accessor names (#382): most_likely_candidate vs.
# most_likely_effect vs. highest_impact_candidate vs.
# highest_impact_effect.
# --------------------------------------------------------------------


def test_most_likely_candidate_returns_effect_candidate():
    """most_likely_candidate is the wrapped form (provenance + effect)."""
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert isinstance(target.most_likely_candidate, EffectCandidate)


def test_most_likely_effect_returns_inner_mutation_effect():
    """most_likely_effect peels off the wrapper for callers that
    only need the inner Effect (e.g. ``most_likely_effect.short_description``).
    """
    from varcode.effects import MutationEffect
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert isinstance(target.most_likely_effect, MutationEffect)
    assert not isinstance(target.most_likely_effect, EffectCandidate)
    # And it's the same effect as inside most_likely_candidate.
    assert target.most_likely_effect is target.most_likely_candidate.effect


def test_highest_impact_candidate_returns_most_disruptive():
    """highest_impact_candidate picks by effect_priority, not by
    probability. For a SpliceDonor-backed SpliceOutcomeSet, the
    EXON_SKIPPING / INTRON_RETENTION candidates resolve to concrete
    coding effects (Deletion / PrematureStop placeholder etc.) while
    NORMAL_SPLICING is a placeholder Intronic — the highest-impact
    one should never be the Intronic placeholder.
    """
    from varcode.effects import effect_priority, Intronic
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    top = target.highest_impact_candidate
    assert isinstance(top, EffectCandidate)
    # It's the candidate with the highest priority among all.
    expected_priority = max(
        effect_priority(c.effect) for c in target.candidates)
    assert effect_priority(top.effect) == expected_priority
    # Intronic (the NORMAL_SPLICING placeholder) has the lowest
    # priority — should never win.
    assert type(top.effect) is not Intronic


def test_highest_impact_effect_unwraps_to_mutation_effect():
    from varcode.effects import MutationEffect
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert isinstance(target.highest_impact_effect, MutationEffect)
    assert target.highest_impact_effect is target.highest_impact_candidate.effect


def test_effects_unwraps_candidate_tuple():
    """The .effects helper returns just the inner MutationEffects."""
    from varcode.effects import MutationEffect
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    inner = target.effects
    assert isinstance(inner, tuple)
    assert len(inner) == len(target.candidates)
    assert all(isinstance(e, MutationEffect) for e in inner)
    assert inner == tuple(c.effect for c in target.candidates)


# --------------------------------------------------------------------
# Serialization invariants (#382): no duplicate `candidates` /
# `_candidates` keys; to_dict doesn't mutate self.
# --------------------------------------------------------------------


def test_to_dict_has_single_candidates_key_no_dup():
    """The serialized form must not carry both ``candidates`` (the
    init param / property name) AND ``_candidates`` (the private
    storage slot) — that would duplicate the payload."""
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    d = target.to_dict()
    assert "candidates" in d
    assert "_candidates" not in d


def test_to_dict_does_not_mutate_self():
    """``to_dict`` must build the serialized payload without
    mutating ``self._candidates`` mid-call. A reader observing
    ``self.candidates`` during serialization should always see the
    runtime form (enum, not string)."""
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    before = tuple(target.candidates)
    target.to_dict()
    after = tuple(target.candidates)
    assert before == after
    # Spot-check that the in-memory form still has enum objects, not
    # the stringified wire form.
    for c in target.candidates:
        outcome = c.evidence.get("splice_outcome")
        assert isinstance(outcome, SpliceOutcome)


# --------------------------------------------------------------------
# Unified EffectCandidate.effect contract (#339 / #382).
# --------------------------------------------------------------------


def test_candidate_effect_is_always_a_mutation_effect():
    from varcode.effects.effect_classes import MutationEffect
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    for candidate in target.candidates:
        assert isinstance(candidate.effect, MutationEffect)
        assert isinstance(candidate.effect.short_description, str)


def test_candidates_carry_splice_outcome_in_evidence():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    tags = {c.evidence["splice_outcome"] for c in target.candidates}
    assert SpliceOutcome.EXON_SKIPPING in tags
    assert SpliceOutcome.INTRON_RETENTION in tags


def test_candidates_carry_probability():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    for candidate in target.candidates:
        assert candidate.probability is not None
        assert 0.0 <= candidate.probability <= 1.0


def test_intron_retention_candidate_effect_is_placeholder_class():
    from varcode.effects.effect_classes import PredictedIntronRetention
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    ir = _candidate_for(target, SpliceOutcome.INTRON_RETENTION)
    assert isinstance(ir.effect, PredictedIntronRetention)
    assert ir.effect.mutant_protein_sequence is None
    assert "intron-retention" in ir.effect.short_description


def test_cryptic_donor_candidate_effect_is_placeholder_class():
    from varcode.effects.effect_classes import PredictedCrypticSpliceSite
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    try:
        cryptic = _candidate_for(target, SpliceOutcome.CRYPTIC_DONOR)
    except StopIteration:
        pytest.skip("Plausibility table doesn't include CRYPTIC_DONOR here.")
    assert isinstance(cryptic.effect, PredictedCrypticSpliceSite)
    assert cryptic.effect.direction == "donor"
    assert cryptic.effect.mutant_protein_sequence is None


# --------------------------------------------------------------------
# _make_splice_candidate aliasing guard (#382).
# --------------------------------------------------------------------


def test_make_splice_candidate_rejects_overwriting_shared_mutant_transcript():
    """The internal helper mutates the inner effect's
    ``mutant_transcript`` field. If a future caller passes a shared
    ``coding_effect`` that already has a ``mutant_transcript`` set
    AND also passes a different ``mutant_transcript``, the
    assertion fires rather than silently clobbering the shared
    state."""
    from unittest.mock import MagicMock
    from varcode.splice_outcomes import _make_splice_candidate
    from varcode import MutantTranscript

    fake_splice_effect = MagicMock()
    fake_splice_effect.transcript = None
    fake_splice_effect.variant = MagicMock()

    pre_owned = MagicMock(spec=MutantTranscript)
    coding = MagicMock()
    coding.mutant_transcript = pre_owned  # pretend it's already assigned

    other = MagicMock(spec=MutantTranscript)
    with pytest.raises(AssertionError):
        _make_splice_candidate(
            fake_splice_effect,
            SpliceOutcome.EXON_SKIPPING,
            0.5,
            coding_effect=coding,
            mutant_transcript=other,
        )


# --------------------------------------------------------------------
# Boundary-codon reconstruction (#298).
# --------------------------------------------------------------------


def test_mid_codon_in_frame_exon_skip_reshapes_boundary():
    variant = Variant("7", 117536674, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    if not isinstance(bare, SpliceDonor):
        pytest.skip(
            "Canonical donor G not present at 117536674; "
            "classifier emitted %s." % type(bare).__name__)
    effects = variant.effects(splice_outcomes=True)
    splice_set = next(e for e in effects if e.transcript is transcript)
    skip = _candidate_for(splice_set, SpliceOutcome.EXON_SKIPPING)
    effect_class = type(skip.effect).__name__
    assert effect_class in ("Deletion", "Substitution", "ComplexSubstitution"), (
        "Expected Deletion/Substitution/ComplexSubstitution, got %s"
        % effect_class)


def test_codon_aligned_exon_skip_still_pure_deletion():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    splice_set = next(e for e in effects if e.transcript is transcript)
    skip = _candidate_for(splice_set, SpliceOutcome.EXON_SKIPPING)
    assert type(skip.effect).__name__ == "Deletion"


# --------------------------------------------------------------------
# Serialization round-trip
# --------------------------------------------------------------------


def test_splice_outcome_set_json_round_trip():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)

    payload = target.to_json()
    restored = SpliceOutcomeSet.from_json(payload)

    assert type(restored) is SpliceOutcomeSet
    assert len(restored.candidates) == len(target.candidates)
    assert restored.disrupted_signal_class is target.disrupted_signal_class
    assert _outcome_of(restored.most_likely_candidate) is _outcome_of(target.most_likely_candidate)
    for rt_c, og_c in zip(restored.candidates, target.candidates):
        assert _outcome_of(rt_c) is _outcome_of(og_c)
        assert rt_c.probability == og_c.probability
        assert type(rt_c.effect) is type(og_c.effect)


def test_splice_outcome_set_dict_round_trip_is_idempotent():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)

    d = target.to_dict()
    restored = SpliceOutcomeSet.from_dict(d)
    assert restored.to_dict() == d


def test_splice_outcome_set_round_trip_preserves_priority_class():
    from varcode.effects import effect_priority
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    restored = SpliceOutcomeSet.from_json(target.to_json())
    assert effect_priority(restored) == effect_priority(target)


def test_non_splice_effects_are_not_multi_outcome():
    from varcode import MultiOutcomeEffect
    from varcode.effects import Substitution, Silent, Intronic, MutationEffect
    for cls in (Substitution, Silent, Intronic, MutationEffect):
        assert not issubclass(cls, MultiOutcomeEffect), (
            "%s should not be a MultiOutcomeEffect" % cls.__name__)


# --------------------------------------------------------------------
# Priority integration.
# --------------------------------------------------------------------


def test_splice_outcome_set_sorts_as_disrupted_signal_class():
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

    assert wrapped_priority == bare_priority


def test_splice_outcome_set_top_priority_works():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    effects = variant.effects(splice_outcomes=True)
    top = effects.top_priority_effect()
    top_class_name = type(top).__name__
    assert top_class_name in ("SpliceOutcomeSet", "SpliceDonor")


# --------------------------------------------------------------------
# Acceptor-side IntronicSpliceSite emits CRYPTIC_ACCEPTOR.
# --------------------------------------------------------------------


def test_acceptor_side_intronic_splice_site_uses_cryptic_acceptor():
    from varcode.effects import IntronicSpliceSite
    variant = Variant("7", 117530896, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    assert isinstance(bare, IntronicSpliceSite) and \
        not isinstance(bare, (SpliceDonor, SpliceAcceptor))
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    outcomes = {_outcome_of(c) for c in target.candidates}
    assert SpliceOutcome.CRYPTIC_ACCEPTOR in outcomes
    assert SpliceOutcome.CRYPTIC_DONOR not in outcomes


def test_donor_side_intronic_splice_site_uses_cryptic_donor():
    from varcode.effects import IntronicSpliceSite
    variant = Variant("7", 117531117, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    assert isinstance(bare, IntronicSpliceSite) and \
        not isinstance(bare, (SpliceDonor, SpliceAcceptor))
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    outcomes = {_outcome_of(c) for c in target.candidates}
    assert SpliceOutcome.CRYPTIC_DONOR in outcomes
    assert SpliceOutcome.CRYPTIC_ACCEPTOR not in outcomes


# --------------------------------------------------------------------
# Multi-protein surface
# --------------------------------------------------------------------


def test_candidate_proteins_maps_each_outcome_to_a_protein():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    proteins = target.candidate_proteins
    outcomes = {_outcome_of(c) for c in target.candidates}
    assert set(proteins.keys()) == outcomes
    assert proteins[SpliceOutcome.EXON_SKIPPING], \
        "Expected a concrete mutant protein for in-frame exon skipping"
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
    ref = str(transcript.protein_sequence)
    shortest = min(proteins, key=len)
    assert len(shortest) < len(ref)


# --------------------------------------------------------------------
# Out-of-frame exon skip produces a real mutant protein.
# --------------------------------------------------------------------


def test_out_of_frame_exon_skip_produces_mutant_protein():
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    target_exon = None
    for exon in transcript.exons[2:]:
        length = exon.end - exon.start + 1
        if length % 3 != 0:
            target_exon = exon
            break
    if target_exon is None:
        pytest.skip("No out-of-frame exon found in CFTR beyond exon 2")
    donor_plus_1 = target_exon.end + 1
    variant = Variant("7", donor_plus_1, "G", "A", ensembl_grch38)
    bare = variant.effect_on_transcript(transcript)
    if not isinstance(bare, SpliceDonor):
        pytest.skip(
            "Canonical donor G not present at %d; classifier emitted %s "
            "rather than SpliceDonor." % (donor_plus_1, type(bare).__name__))
    effects = variant.effects(splice_outcomes=True)
    splice_set = next(e for e in effects if e.transcript is transcript)
    skip = _candidate_for(splice_set, SpliceOutcome.EXON_SKIPPING)
    protein = skip.effect.mutant_protein_sequence
    assert isinstance(protein, str)
    assert len(protein) > 0
    assert protein != str(transcript.protein_sequence)


# --------------------------------------------------------------------
# mutant_transcript exposure on the inner effect (#305).
# --------------------------------------------------------------------


def test_exon_skipping_candidate_exposes_mutant_transcript():
    from varcode import MutantTranscript
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    exon_skip = _candidate_for(target, SpliceOutcome.EXON_SKIPPING)
    mt = exon_skip.effect.mutant_transcript
    assert isinstance(mt, MutantTranscript)
    assert mt.cdna_sequence is not None
    assert len(mt.edits) == 1
    edit = mt.edits[0]
    assert edit.alt_bases == ""
    assert edit.cdna_end > edit.cdna_start
    assert mt.annotator_name == "splice_outcomes"


def test_stub_candidates_have_no_mutant_transcript():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    intron = _candidate_for(target, SpliceOutcome.INTRON_RETENTION)
    assert getattr(intron.effect, "mutant_transcript", None) is None


# --------------------------------------------------------------------
# Intron-retention with a genomic_sequence provider (#296).
# --------------------------------------------------------------------


def _synthetic_intron_provider(stop_offset=30, pad_with="A"):
    def provider(contig, start, end):
        length = end - start + 1
        core = (pad_with * stop_offset) + "TAA" + (pad_with * length)
        return core[:length]
    return provider


def test_intron_retention_with_provider_populates_mutant_transcript():
    from varcode import MutantTranscript
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    assert isinstance(bare, SpliceDonor)
    splice_set = enumerate_splice_outcomes(
        bare, genomic_sequence=_synthetic_intron_provider())
    intron = _candidate_for(splice_set, SpliceOutcome.INTRON_RETENTION)
    mt = intron.effect.mutant_transcript
    assert isinstance(mt, MutantTranscript)
    assert len(mt.cdna_sequence) > len(str(transcript.sequence))
    assert mt.mutant_protein_sequence is not None
    assert len(mt.mutant_protein_sequence) < len(
        str(transcript.protein_sequence))
    # Effect is now a real classified coding effect, not the
    # PredictedIntronRetention placeholder.
    from varcode.effects.effect_classes import PredictedIntronRetention
    assert not isinstance(intron.effect, PredictedIntronRetention)


def test_intron_retention_without_provider_still_stub():
    from varcode.effects.effect_classes import PredictedIntronRetention
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    splice_set = enumerate_splice_outcomes(bare)  # no provider
    intron = _candidate_for(splice_set, SpliceOutcome.INTRON_RETENTION)
    assert getattr(intron.effect, "mutant_transcript", None) is None
    assert isinstance(intron.effect, PredictedIntronRetention)


def _synthetic_cryptic_donor_provider(motif_offset=65, motif="CAGGTAAGT"):
    def provider(contig, start, end):
        length = end - start + 1
        seq = bytearray(b"A" * length)
        for i, c in enumerate(motif):
            if motif_offset + i < length:
                seq[motif_offset + i] = ord(c)
        return seq.decode()
    return provider


def test_cryptic_donor_with_provider_picks_up_implanted_motif():
    from varcode import MutantTranscript
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    splice_set = enumerate_splice_outcomes(
        bare, genomic_sequence=_synthetic_cryptic_donor_provider(
            motif_offset=65))
    cryptic = _candidate_for(splice_set, SpliceOutcome.CRYPTIC_DONOR)
    mt = cryptic.effect.mutant_transcript
    assert isinstance(mt, MutantTranscript)
    delta = len(mt.cdna_sequence) - len(str(transcript.sequence))
    assert delta > 0
    desc = cryptic.evidence["description"]
    assert "motif score" in desc
    assert "Cryptic donor" in desc


def test_cryptic_acceptor_with_provider_picks_up_implanted_motif():
    variant = Variant("7", 117530898, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    if not isinstance(bare, SpliceAcceptor):
        pytest.skip("Variant didn't classify as SpliceAcceptor")

    def provider(contig, start, end):
        length = end - start + 1
        seq = bytearray(b"A" * length)
        motif = b"CAGG"
        pos = 40
        for i, c in enumerate(motif):
            if pos + i < length:
                seq[pos + i] = c
        return seq.decode()

    splice_set = enumerate_splice_outcomes(bare, genomic_sequence=provider)
    cryptic = _candidate_for(splice_set, SpliceOutcome.CRYPTIC_ACCEPTOR)
    assert cryptic.effect.mutant_transcript is not None
    assert "Cryptic acceptor" in cryptic.evidence["description"]


def test_cryptic_without_provider_still_stub():
    from varcode.effects.effect_classes import PredictedCrypticSpliceSite
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    splice_set = enumerate_splice_outcomes(bare)
    cryptic = _candidate_for(splice_set, SpliceOutcome.CRYPTIC_DONOR)
    assert isinstance(cryptic.effect, PredictedCrypticSpliceSite)
    assert getattr(cryptic.effect, "mutant_transcript", None) is None


def test_cryptic_with_no_motif_in_range_is_stub():
    from varcode.effects.effect_classes import PredictedCrypticSpliceSite
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)

    def provider(contig, start, end):
        return "N" * (end - start + 1)

    splice_set = enumerate_splice_outcomes(bare, genomic_sequence=provider)
    cryptic = _candidate_for(splice_set, SpliceOutcome.CRYPTIC_DONOR)
    assert isinstance(cryptic.effect, PredictedCrypticSpliceSite)
    assert getattr(cryptic.effect, "mutant_transcript", None) is None


def test_intron_retention_provider_can_raise_without_crashing():
    from varcode.effects.effect_classes import PredictedIntronRetention
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)

    def failing_provider(contig, start, end):
        raise IOError("FASTA not available")
    splice_set = enumerate_splice_outcomes(
        bare, genomic_sequence=failing_provider)
    intron = _candidate_for(splice_set, SpliceOutcome.INTRON_RETENTION)
    assert isinstance(intron.effect, PredictedIntronRetention)
