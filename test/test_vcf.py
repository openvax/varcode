from varcode import VariantCollection, VariantAnnotator

VCF_FILENAME = "somatic_hg19_14muts.vcf"

def test_vcf_reference_name():
    variants = VariantCollection(VCF_FILENAME)
    # the raw reference name can be a file path to the hg19 FASTA file
    assert "hg19" in variants.raw_reference_name, \
        "Expected hg19 reference, got %s" % (variants.raw_reference_name,)

    # after normalization, hg19 should be remapped to GRCh37
    assert variants.reference_name == "GRCh37"

def test_vcf_number_entries():
    # there are 14 mutations listed in the VCF, make sure they are all parsed
    variants = VariantCollection(VCF_FILENAME)
    assert len(variants) == 14, \
        "Expected 14 mutations, got %d" % (len(variants),)

def _check_effect_gene_name(variant, effect):
    expected_gene_name = variant.info['GE']
    gene_names = [gene.name for gene in effect.genes]
    assert expected_gene_name in gene_names, \
        "Expected gene name %s for variant %s, got %s" % (
            expected_gene_name, variant, gene_names)

def test_vcf_gene_names():
    variants = VariantCollection(VCF_FILENAME)
    for variant, effect in variants.variant_effects():
        yield (_check_effect_gene_name, variant, effect)
