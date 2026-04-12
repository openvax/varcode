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
Regression tests for https://github.com/openvax/varcode/issues/227

EffectCollection previously had `sort_key=None` by default, so the
returned effects were in whatever order the gene/transcript loop
produced them. Users expect effects to be ordered by priority (highest
first) so that iterating over a collection surfaces the most severe
effects up front.
"""

from varcode import Variant
from varcode.effects.effect_ordering import effect_priority


def test_effect_collection_is_sorted_by_priority_descending():
    # Use a substitution that produces a mix of Substitution (high),
    # ExonicSpliceSite (lower), Intronic (lower), NoncodingTranscript
    # (lowest) across many transcripts. Previously these came back in
    # insertion order; now they should be ordered by priority, highest
    # first.
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", "GRCh38")
    effects = variant.effects()
    assert len(effects) > 1, "Test variant should produce multiple effects"

    priorities = [effect_priority(e) for e in effects]
    assert priorities == sorted(priorities, reverse=True), (
        "Expected effects in descending priority order, got %r" % priorities
    )


def test_effect_collection_top_priority_effect_matches_first_element():
    # When the collection is sorted by priority, the first element should
    # match top_priority_effect() (modulo tie-breaking, which uses the same
    # priority ordering).
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", "GRCh38")
    effects = variant.effects()
    top = effects.top_priority_effect()
    assert effect_priority(effects[0]) == effect_priority(top), \
        "First effect should have the same priority as top_priority_effect()"


def test_effect_collection_explicit_sort_key_still_works():
    # Users who pass their own sort_key should still get that ordering.
    from varcode.effects.effect_collection import EffectCollection
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", "GRCh38")
    effects_list = list(variant.effects())
    # Sort by class name alphabetically
    custom = EffectCollection(
        effects=effects_list,
        sort_key=lambda e: e.__class__.__name__,
    )
    names = [e.__class__.__name__ for e in custom]
    assert names == sorted(names), \
        "Custom sort_key should be respected, got %r" % names
