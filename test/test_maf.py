from varcode import load_variants, VariantCollection, Variant
from nose.tools import eq_

def test_maf():
    variant_collection_from_maf = load_variants("data/tcga_ov.head.maf")
    eq_(variant_collection_from_maf.reference_name, "GRCh37")
    expected_variants = [
        Variant(1, 1650797, "A", "G"),
        Variant(1, 231401797, "A", "C"),
        Variant(1, 23836447, "C", "A"),
        Variant(11,124617502, "C", "G"),
    ]
    eq_(len(variant_collection_from_maf), len(expected_variants))
    for v1, v2 in zip(expected_variants, variant_collection_from_maf):
        eq_(v1, v2)
