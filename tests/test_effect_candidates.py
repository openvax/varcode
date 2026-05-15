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

"""Tests for the unified :class:`EffectCandidate` type (openvax/varcode#299).

The type is deliberately minimal (dataclass + one helper), so tests
focus on (a) the contract — probability bounds, defaults, frozen
semantics — and (b) the ``MultiOutcomeEffect.candidates`` harmonized
accessor lifting existing ``candidates`` into the new shape.
"""

import pytest
from pyensembl import cached_release

from varcode import EffectCandidate, Variant, candidates_from_effects
from varcode.effects.effect_classes import (
    ExonicSpliceSite,
    MultiOutcomeEffect,
)

ensembl_grch38 = cached_release(81)


def test_effect_candidate_defaults():
    o = EffectCandidate(effect=object())
    assert o.probability is None
    assert o.source == "varcode"
    assert o.evidence == {}


def test_effect_candidate_rejects_removed_description_kwarg():
    """The pre-rename ``Outcome`` carried a ``description`` field that
    was a pure passthrough to ``effect.short_description``. The rename
    dropped it; lock that in by pinning the TypeError for any caller
    still using the old kwarg."""
    with pytest.raises(TypeError):
        EffectCandidate(effect=object(), description="old kwarg")


def test_effect_candidate_probability_bounds():
    EffectCandidate(effect=object(), probability=0.0)
    EffectCandidate(effect=object(), probability=1.0)
    EffectCandidate(effect=object(), probability=None)
    with pytest.raises(ValueError):
        EffectCandidate(effect=object(), probability=-0.1)
    with pytest.raises(ValueError):
        EffectCandidate(effect=object(), probability=1.5)


def test_effect_candidate_is_frozen():
    o = EffectCandidate(effect=object(), source="test")
    with pytest.raises(Exception):
        # dataclasses.FrozenInstanceError is a subclass of AttributeError
        # in some Python versions — catch the broad shape.
        o.source = "mutated"


def test_effect_candidate_short_description_passthrough():
    class _FakeEffect:
        short_description = "p.L101del"
    o = EffectCandidate(effect=_FakeEffect())
    assert o.short_description == "p.L101del"


def test_candidates_from_effects_tags_source():
    class _C:
        short_description = "c1"
    outcomes = candidates_from_effects((_C(), _C()), source="test_source")
    assert len(outcomes) == 2
    assert all(o.source == "test_source" for o in outcomes)
    assert all(o.probability is None for o in outcomes)


def test_exonic_splice_site_exposes_outcomes():
    """Real integration: an SNV at the last base of an exon yields
    ``ExonicSpliceSite``, which is a ``MultiOutcomeEffect``. Its
    ``.candidates`` should return two :class:`EffectCandidate` entries (the
    splice-disruption outcome and the coding-change alternate)."""
    # CFTR exon 3 ends at 117531114 (last exon base).
    variant = Variant("7", 117531114, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id("ENST00000003084")
    effect = variant.effect_on_transcript(transcript)
    assert isinstance(effect, MultiOutcomeEffect)
    assert isinstance(effect, ExonicSpliceSite)

    outcomes = effect.candidates
    assert len(outcomes) == 2
    assert all(isinstance(o, EffectCandidate) for o in outcomes)
    # First outcome: the splice-disruption classification (the
    # ExonicSpliceSite itself).
    assert outcomes[0].effect is effect
    # Second outcome: the coding-change alternate.
    assert outcomes[1].effect is effect.alternate_effect
    # Both default to varcode-source and unscored probability.
    assert all(o.source == "varcode" for o in outcomes)
    assert all(o.probability is None for o in outcomes)


def test_sv_outcomes_carry_sv_type_in_evidence():
    """:class:`StructuralVariantEffect.candidates` attaches the
    ``sv_type`` to each outcome's evidence so consumers iterating
    outcomes can filter by SV kind uniformly (#339)."""
    from varcode import StructuralVariant, StructuralVariantAnnotator
    transcript = ensembl_grch38.transcript_by_id("ENST00000003084")
    sv = StructuralVariant(
        contig="7",
        start=transcript.start + 100,
        end=transcript.start + 5_000,
        sv_type="DEL",
        genome=ensembl_grch38)
    effect = StructuralVariantAnnotator().annotate_on_transcript(sv, transcript)
    outcomes = effect.candidates
    assert len(outcomes) >= 1
    for o in outcomes:
        assert o.evidence["sv_type"] == "DEL"


def test_uniform_iteration_sv_and_splice_outcomes():
    """The point of #339: downstream code iterating
    ``outcome.effect.short_description`` works across SV and splice
    outcomes without any ``isinstance`` branching."""
    from varcode import StructuralVariant, StructuralVariantAnnotator
    from varcode.effects.effect_classes import MutationEffect

    transcript = ensembl_grch38.transcript_by_id("ENST00000003084")

    sv = StructuralVariant(
        contig="7",
        start=transcript.start + 100,
        end=transcript.start + 5_000,
        sv_type="DEL",
        genome=ensembl_grch38)
    sv_effect = StructuralVariantAnnotator().annotate_on_transcript(
        sv, transcript)

    splice_variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    splice_effects = splice_variant.effects(splice_outcomes=True)
    splice_effect = next(
        e for e in splice_effects if e.transcript is transcript)

    all_outcomes = list(sv_effect.candidates) + list(splice_effect.candidates)
    # One-hop read across both outcome kinds:
    for o in all_outcomes:
        assert isinstance(o.effect, MutationEffect)
        assert isinstance(o.effect.short_description, str)
        # mutant_protein_sequence is the canonical access point; may
        # be None for placeholder effects but the attribute exists.
        _ = o.effect.mutant_protein_sequence


def test_effect_candidate_round_trips_via_json():
    """``EffectCandidate`` now inherits :class:`DataclassSerializable`, so
    ``to_json`` / ``from_json`` round-trip the full outcome — including
    a polymorphic :class:`MutationEffect` ``effect`` field — without
    any custom serialization code (#343)."""
    variant = Variant("7", 117531114, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id("ENST00000003084")
    real_effect = variant.effect_on_transcript(transcript)
    o = EffectCandidate(
        effect=real_effect,
        probability=0.75,
        source="spliceai",
        evidence={"ds_ag": 0.12})
    rt = EffectCandidate.from_json(o.to_json())
    assert rt.probability == 0.75
    assert rt.source == "spliceai"
    assert rt.evidence == {"ds_ag": 0.12}
    # effect round-trips polymorphically — same class as the original.
    assert type(rt.effect) is type(real_effect)


def test_effect_candidate_accepts_external_scorer_shape():
    """An external predictor (SpliceAI-style) can construct an
    ``EffectCandidate`` with its own probability and evidence dict. This
    pins the interchange contract — varcode doesn't ship a scorer,
    but the type stays usable by one."""
    class _FakeEffect:
        short_description = "splice-donor"
    scored = EffectCandidate(
        effect=_FakeEffect(),
        probability=0.87,
        source="spliceai",
        evidence={"ds_ag": 0.02, "ds_al": 0.01, "ds_dg": 0.87, "ds_dl": 0.03})
    assert scored.probability == 0.87
    assert scored.source == "spliceai"
    assert scored.evidence["ds_dg"] == 0.87
