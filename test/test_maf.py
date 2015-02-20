from nose.tools import eq_
from pyensembl import EnsemblRelease
from varcode import load_maf, VariantCollection, Variant

def test_maf():
    ensembl = EnsemblRelease(75)
    variant_collection_from_maf = load_maf("data/tcga_ov.head.maf")
    expected_variants = [
        Variant(1, 1650797, "A", "G", ensembl),
        Variant(1, 231401797, "A", "C", ensembl),
        Variant(1, 23836447, "C", "A", ensembl),
        Variant(11,124617502, "C", "G", ensembl),
    ]
    eq_(len(variant_collection_from_maf), len(expected_variants))
    for v_expect, v_maf in zip(expected_variants, variant_collection_from_maf):
        eq_(v_expect, v_maf)
        gene_name = v_maf.info['Hugo_Symbol']
        assert any(gene.name == gene_name for gene in v_maf.genes()), \
            "Expected gene name %s but got %s" % (gene_name, v_maf.genes())
