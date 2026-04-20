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

"""Tests for the unified :class:`Outcome` type (openvax/varcode#299).

The type is deliberately minimal (dataclass + one helper), so tests
focus on (a) the contract — probability bounds, defaults, frozen
semantics — and (b) the ``MultiOutcomeEffect.outcomes`` harmonized
accessor lifting existing ``candidates`` into the new shape.
"""

import pytest
from pyensembl import cached_release

from varcode import Outcome, Variant, outcomes_from_candidates
from varcode.effects.effect_classes import (
    ExonicSpliceSite,
    MultiOutcomeEffect,
)

ensembl_grch38 = cached_release(81)


def test_outcome_defaults():
    o = Outcome(effect=object())
    assert o.probability is None
    assert o.source == "varcode"
    assert o.evidence == {}


def test_outcome_probability_bounds():
    Outcome(effect=object(), probability=0.0)
    Outcome(effect=object(), probability=1.0)
    Outcome(effect=object(), probability=None)
    with pytest.raises(ValueError):
        Outcome(effect=object(), probability=-0.1)
    with pytest.raises(ValueError):
        Outcome(effect=object(), probability=1.5)


def test_outcome_is_frozen():
    o = Outcome(effect=object(), source="test")
    with pytest.raises(Exception):
        # dataclasses.FrozenInstanceError is a subclass of AttributeError
        # in some Python versions — catch the broad shape.
        o.source = "mutated"


def test_outcome_short_description_passthrough():
    class _FakeEffect:
        short_description = "p.L101del"
    o = Outcome(effect=_FakeEffect())
    assert o.short_description == "p.L101del"


def test_outcomes_from_candidates_tags_source():
    class _C:
        short_description = "c1"
    outcomes = outcomes_from_candidates((_C(), _C()), source="test_source")
    assert len(outcomes) == 2
    assert all(o.source == "test_source" for o in outcomes)
    assert all(o.probability is None for o in outcomes)


def test_exonic_splice_site_exposes_outcomes():
    """Real integration: an SNV at the last base of an exon yields
    ``ExonicSpliceSite``, which is a ``MultiOutcomeEffect``. Its
    ``.outcomes`` should return two :class:`Outcome` entries (the
    splice-disruption outcome and the coding-change alternate)."""
    # CFTR exon 3 ends at 117531114 (last exon base).
    variant = Variant("7", 117531114, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id("ENST00000003084")
    effect = variant.effect_on_transcript(transcript)
    assert isinstance(effect, MultiOutcomeEffect)
    assert isinstance(effect, ExonicSpliceSite)

    outcomes = effect.outcomes
    assert len(outcomes) == 2
    assert all(isinstance(o, Outcome) for o in outcomes)
    # First outcome: the splice-disruption classification (the
    # ExonicSpliceSite itself).
    assert outcomes[0].effect is effect
    # Second outcome: the coding-change alternate.
    assert outcomes[1].effect is effect.alternate_effect
    # Both default to varcode-source and unscored probability.
    assert all(o.source == "varcode" for o in outcomes)
    assert all(o.probability is None for o in outcomes)


def test_outcome_accepts_external_scorer_shape():
    """An external predictor (SpliceAI-style) can construct an
    ``Outcome`` with its own probability and evidence dict. This
    pins the interchange contract — varcode doesn't ship a scorer,
    but the type stays usable by one."""
    class _FakeEffect:
        short_description = "splice-donor"
    scored = Outcome(
        effect=_FakeEffect(),
        probability=0.87,
        source="spliceai",
        evidence={"ds_ag": 0.02, "ds_al": 0.01, "ds_dg": 0.87, "ds_dl": 0.03})
    assert scored.probability == 0.87
    assert scored.source == "spliceai"
    assert scored.evidence["ds_dg"] == 0.87
