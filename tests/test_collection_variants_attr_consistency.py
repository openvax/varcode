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
Regression tests for https://github.com/openvax/varcode/issues/220

VariantCollection.variants was assigned the raw input list before
Collection.__init__ sorted / de-duplicated elements, so iterating a
VariantCollection produced one order while `.variants` produced
another. The same pattern existed in EffectCollection.effects.
"""

import varcode
from varcode import Variant, VariantCollection
from varcode.effects.effect_collection import EffectCollection


def _unsorted_variants():
    return [
        Variant(contig="1", start=4, ref="A", alt="T"),
        Variant(contig="1", start=1, ref="A", alt="T"),
        Variant(contig="1", start=3, ref="A", alt="T"),
        Variant(contig="1", start=2, ref="A", alt="T"),
    ]


def test_variant_collection_variants_matches_iteration_with_default_sort():
    vc = VariantCollection(variants=_unsorted_variants())
    iterated = [v.start for v in vc]
    via_attr = [v.start for v in vc.variants]
    assert iterated == via_attr, (
        "vc.variants should match iteration order, got %r vs %r"
        % (via_attr, iterated)
    )


def test_variant_collection_variants_matches_iteration_with_sort_key_none_distinct_false():
    # Explicitly opt out of sorting and deduplication. Iteration and
    # .variants should both preserve the original input order.
    variants = _unsorted_variants()
    vc = VariantCollection(
        variants=variants, sort_key=None, distinct=False
    )
    iterated = [v.start for v in vc]
    via_attr = [v.start for v in vc.variants]
    assert iterated == via_attr
    # And the order should match the input.
    assert iterated == [4, 1, 3, 2]


def test_variant_collection_variants_matches_iteration_with_distinct_true():
    # When distinct=True and sort_key=None, the base Collection uses
    # set() which loses order. Regardless of what order it ends up in,
    # iteration and .variants should agree.
    variants = _unsorted_variants()
    vc = VariantCollection(
        variants=variants, sort_key=None, distinct=True
    )
    iterated = [v.start for v in vc]
    via_attr = [v.start for v in vc.variants]
    assert iterated == via_attr


def test_variant_collection_len_matches_variants_attr():
    vc = VariantCollection(variants=_unsorted_variants())
    assert len(vc) == len(vc.variants)


def test_effect_collection_effects_matches_iteration():
    # Same bug existed in EffectCollection.effects — produce an effect
    # collection and verify the two access paths agree.
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", "GRCh38")
    effects = variant.effects()
    iterated = [e.__class__.__name__ for e in effects]
    via_attr = [e.__class__.__name__ for e in effects.effects]
    assert iterated == via_attr, (
        "effects.effects should match iteration order, got %r vs %r"
        % (via_attr, iterated)
    )
