"""
Any variants which are encountered in the wild and either cause Varcode
to crash or return an incorrect annotation should be added to this
test module.
"""

from pyensembl import EnsemblRelease
from varcode import Variant

ensembl75 = EnsemblRelease(75)

# variants which have previously resulted in raised exceptions
# during effect annotation
should_not_crash_variants = [
    # error message:
    # "Couldn't find position 92979124 on any exon of ENST00000540033"
    Variant(
        contig=1,
        pos=92979092,
        ref="ATATATATATATATATATATATATATATATATG",
        alt="A",
        ensembl=ensembl75),
    # error message:
    # "Expect non-silent stop-loss variant to cause longer variant protein"
    # "" but got len(original) = 653, len(variant) = 653"
    Variant(
        contig=1,
        pos=167385324,
        ref="TAA",
        alt="T",
        ensembl=ensembl75),
    # error message:
    # "Variant which span 5' UTR and CDS not supported"
    Variant(
        contig=19,
        pos=44351166,
        ref="GGGAGAT",
        alt="G",
        ensembl=ensembl75)
]

def try_effect_annotation(variant):
    effect = variant.top_effect()
    assert effect is not None

def test_crashing_variants():
    for variant in should_not_crash_variants:
        yield (try_effect_annotation, variant)