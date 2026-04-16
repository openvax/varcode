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

"""Tests for the ExonicSpliceSite → MultiOutcomeEffect retrofit
(openvax/varcode#299 Part 1).

The retrofit is a class-hierarchy change only:

* ExonicSpliceSite inherits from MultiOutcomeEffect and exposes
  the ``candidates`` / ``most_likely`` / ``priority_class``
  protocol.
* SpliceOutcomeSet gains a back-compat ``alternate_effect``
  property that resolves to the NORMAL_SPLICING candidate's
  coding_effect.

No behaviour change — the existing test suite (tests/test_splice_*)
continues to lock in byte-for-byte output.
"""

from pyensembl import cached_release

from varcode import MultiOutcomeEffect, SpliceOutcome, Variant
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
    assert effect.candidates == (effect, effect.alternate_effect)


def test_exonic_splice_site_most_likely_is_self():
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effect = variant.effect_on_transcript(transcript)
    assert effect.most_likely is effect


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
# SpliceOutcomeSet gets alternate_effect back-compat property
# ====================================================================


def test_splice_outcome_set_alternate_effect_resolves_to_normal_splicing():
    # ExonicSpliceSite wrapped under splice_outcomes=True becomes a
    # SpliceOutcomeSet. Its .alternate_effect should return the same
    # coding_effect that ExonicSpliceSite.alternate_effect would.
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)

    bare = variant.effect_on_transcript(transcript)
    assert isinstance(bare, ExonicSpliceSite)
    bare_alt = bare.alternate_effect

    wrapped_effects = variant.effects(splice_outcomes=True)
    wrapped = next(e for e in wrapped_effects if e.transcript is transcript)
    # Same underlying coding-effect class; don't require object
    # identity since the wrapped path may rebuild it.
    assert type(wrapped.alternate_effect) is type(bare_alt)
    assert wrapped.alternate_effect.short_description == bare_alt.short_description


def test_splice_outcome_set_alternate_effect_none_when_no_normal_splicing():
    # SpliceDonor-backed SpliceOutcomeSet: the NORMAL_SPLICING candidate
    # exists but its coding_effect is None (intronic variant, no
    # underlying coding change). alternate_effect should be None.
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    wrapped_effects = variant.effects(splice_outcomes=True)
    wrapped = next(e for e in wrapped_effects if e.transcript is transcript)
    # Sanity: this is a SpliceDonor-backed set with a NORMAL_SPLICING
    # candidate that has no coding_effect.
    normal = next(
        c for c in wrapped.candidates
        if c.outcome is SpliceOutcome.NORMAL_SPLICING
    )
    assert normal.coding_effect is None
    assert wrapped.alternate_effect is None


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

    # alternate_effect is callable on both shapes post-retrofit.
    assert bare.alternate_effect is not None
    assert wrapped.alternate_effect is not None
    assert type(bare.alternate_effect) is type(wrapped.alternate_effect)
