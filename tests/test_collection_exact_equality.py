"""Exact equality includes source records and structural INFO (#670)."""

import pytest

from varcode import StructuralVariant, Variant, VariantCollection


def collection(variant, value=1, source="caller"):
    return VariantCollection([variant], source_to_metadata_dict={
        source: {variant: {"info": {"DP": value}}}})


@pytest.mark.parametrize("structural", [False, True])
def test_exact_equality_uses_current_record_storage(structural):
    variant = (StructuralVariant("1", 10, sv_type="DEL", end=20, genome=81, info={"SVLEN": -10})
               if structural else Variant("1", 10, "A", "C", genome=81))
    first = collection(variant)
    assert first.exactly_equal(first)
    assert first.exactly_equal(collection(variant))
    assert first == collection(variant, value=2)
    assert not first.exactly_equal(collection(variant, value=2))
    assert not first.exactly_equal(collection(variant, source="another caller"))
    assert not first.exactly_equal(None)
    assert not first.exactly_equal([variant])
    assert first.exactly_equal(first.filter(lambda _: True))
    assert first.exactly_equal(first.groupby(lambda _: "key")["key"])
    if structural:
        changed = StructuralVariant("1", 10, sv_type="DEL", end=20, genome=81, info={"SVLEN": -11})
        assert not first.exactly_equal(collection(changed))


def test_exact_equality_preserves_order_and_empty_collection_behavior():
    first = Variant("1", 10, "A", "C", genome=81)
    second = Variant("1", 11, "A", "C", genome=81)
    left = VariantCollection([first, second], sort_key=None, distinct=False)
    right = VariantCollection([second, first], sort_key=None, distinct=False)
    assert not left.exactly_equal(right)
    assert VariantCollection([]).exactly_equal(VariantCollection([]))
