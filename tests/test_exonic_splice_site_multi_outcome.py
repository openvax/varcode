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

"""Tests for the ExonicSpliceSite → MultiOutcomeEffect retrofit and
the SpliceOutcomeSet accessor surface
(openvax/varcode#299, #391, #392).

The retrofit makes ExonicSpliceSite a MultiOutcomeEffect — it
exposes the ``candidates`` / ``most_likely`` / ``priority_class``
protocol uniformly with SpliceOutcomeSet. The SpliceOutcomeSet
side exposes ``effect_if_splicing_unchanged``, which resolves to
the NormalSplicing candidate's coding_effect (the opt-in shape's
analogue of the legacy ``ExonicSpliceSite.alternate_effect``).
"""

from pyensembl import cached_release

from varcode import MultiOutcomeEffect, NormalSplicing, Variant
from varcode.effects import ExonicSpliceSite


ensembl_grch38 = cached_release(81)
CFTR_TRANSCRIPT_ID = "ENST00000003084"


# ====================================================================
# ExonicSpliceSite is now a MultiOutcomeEffect
# ====================================================================


def test_exonic_splice_site_is_multi_outcome_effect():
    # CFTR exon 4 ends in ...AAG. A G->T at the last exonic base is a
    # canonical ExonicSpliceSite case.
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effect = variant.effect_on_transcript(transcript)
    assert isinstance(effect, ExonicSpliceSite)
    assert isinstance(effect, MultiOutcomeEffect), (
        "ExonicSpliceSite should pass isinstance(e, MultiOutcomeEffect) "
        "after #299 Part 1 retrofit")


def test_exonic_splice_site_candidates_include_self_and_alternate():
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effect = variant.effect_on_transcript(transcript)
    assert effect.alternate_effect is not None
    # 2-candidate form: splice-disruption (self) + coding consequence.
    # Each entry is an EffectCandidate wrapping the inner effect.
    inners = tuple(c.effect for c in effect.candidates)
    assert inners == (effect, effect.alternate_effect)


def test_exonic_splice_site_most_likely_is_self():
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effect = variant.effect_on_transcript(transcript)
    # most_likely_effect is the inner MutationEffect, which IS the
    # splice-disruption outcome (this effect itself).
    assert effect.most_likely_effect is effect
    # most_likely_candidate wraps it in an EffectCandidate.
    assert effect.most_likely_candidate.effect is effect


def test_exonic_splice_site_priority_class_is_self_class():
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effect = variant.effect_on_transcript(transcript)
    assert effect.priority_class is ExonicSpliceSite


def test_exonic_splice_site_alternate_effect_still_first_class():
    # Back-compat: the attribute access that downstream code uses
    # today must keep working. The retrofit adds MultiOutcomeEffect
    # machinery alongside it, not in place of it.
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effect = variant.effect_on_transcript(transcript)
    alt = effect.alternate_effect
    assert alt is not None
    # alternate_effect is a MutationEffect (not a SpliceCandidate).
    from varcode.effects import MutationEffect
    assert isinstance(alt, MutationEffect)


# ====================================================================
# SpliceOutcomeSet exposes effect_if_splicing_unchanged
# ====================================================================


def test_effect_if_splicing_unchanged_matches_alternate_effect():
    # ExonicSpliceSite wrapped under splice_outcomes=True becomes a
    # SpliceOutcomeSet. Its .effect_if_splicing_unchanged should
    # return the same coding_effect that ExonicSpliceSite.alternate_effect
    # would.
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)

    bare = variant.effect_on_transcript(transcript)
    assert isinstance(bare, ExonicSpliceSite)
    bare_alt = bare.alternate_effect

    wrapped_effects = variant.effects(splice_outcomes=True)
    wrapped = next(e for e in wrapped_effects if e.transcript is transcript)
    # Same underlying coding-effect class; don't require object
    # identity since the wrapped path may rebuild it.
    assert type(wrapped.effect_if_splicing_unchanged) is type(bare_alt)
    assert (
        wrapped.effect_if_splicing_unchanged.short_description
        == bare_alt.short_description)


def test_effect_if_splicing_unchanged_none_for_intronic_donor():
    # SpliceDonor-backed SpliceOutcomeSet: the NormalSplicing candidate
    # exists but its coding_effect is None (intronic variant, no
    # underlying coding change). effect_if_splicing_unchanged should be None.
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    wrapped_effects = variant.effects(splice_outcomes=True)
    wrapped = next(e for e in wrapped_effects if e.transcript is transcript)
    normal = next(
        c for c in wrapped.candidates
        if isinstance(c.effect, NormalSplicing))
    assert normal.effect.coding_effect is None
    assert wrapped.effect_if_splicing_unchanged is None


# ====================================================================
# Unified downstream pattern: isinstance(e, MultiOutcomeEffect) works
# across both shapes.
# ====================================================================


def test_multi_outcome_effect_check_works_on_both_shapes():
    # Same variant expressed two ways; both match the MultiOutcomeEffect
    # isinstance check so downstream code doesn't need to special-case.
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)

    bare = variant.effect_on_transcript(transcript)
    assert isinstance(bare, MultiOutcomeEffect)

    wrapped_effects = variant.effects(splice_outcomes=True)
    wrapped = next(e for e in wrapped_effects if e.transcript is transcript)
    assert isinstance(wrapped, MultiOutcomeEffect)

    # The "if splicing unchanged" coding consequence is reachable on
    # both shapes (under their respective names — ExonicSpliceSite
    # still carries the legacy ``alternate_effect``; SpliceOutcomeSet
    # uses ``effect_if_splicing_unchanged``).
    assert bare.alternate_effect is not None
    assert wrapped.effect_if_splicing_unchanged is not None
    assert (
        type(bare.alternate_effect)
        is type(wrapped.effect_if_splicing_unchanged))
