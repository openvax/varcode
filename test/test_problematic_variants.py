"""
Any variants which are encountered in the wild and either cause Varcode
to crash or return an incorrect annotation should be added to this
test module.
"""

from pyensembl import EnsemblRelease
from varcode import Variant

ensembl75 = EnsemblRelease(75)
ensembl78 = EnsemblRelease(78)

# variants which have previously resulted in raised exceptions
# during effect annotation
should_not_crash_variants = [
    # error message:
    # "Couldn't find position 92979124 on any exon of ENST00000540033"
    Variant(
        contig=1,
        start=92979092,
        ref="ATATATATATATATATATATATATATATATATG",
        alt="A",
        ensembl=ensembl75),
    # error message:
    # "Expect non-silent stop-loss variant to cause longer variant protein"
    # "" but got len(original) = 653, len(variant) = 653"
    Variant(
        contig=1,
        start=167385324,
        ref="TAA",
        alt="T",
        ensembl=ensembl75),
    # error message:
    # "Variant which span 5' UTR and CDS not supported"
    Variant(
        contig=19,
        start=44351166,
        ref="GGGAGAT",
        alt="G",
        ensembl=ensembl75),
    # error message:
    # "Can't have ref = '' and alt = 'E' at aa_pos = 445, cds_pos = 1335"
    Variant(
        contig=1,
        start=1684347,
        ref="",
        alt="CCT",
        ensembl=ensembl75),
    Variant(
        contig=11,
        start=47640416,
        ref="",
        alt="TCTTT",
        ensembl=ensembl75),
    Variant(
        contig=12,
        start=98880902,
        ref="A",
        alt="",
        ensembl=ensembl75),
    Variant(
        contig=19,
        start=52803670,
        ref="TG",
        alt="",
        ensembl=ensembl75),
    Variant(
        contig=1,
        start=109792735,
        ref="",
        alt="CGC",
        ensembl=ensembl75),
    # error message:
    # "expected ref 'GATGTCGG' at offset 1412 of ENST00000297524...CDS has 'G'"
    Variant(
        contig=8,
        start=87226635,
        ref="CCGACATC",
        alt="",
        ensembl=ensembl75),
    # error message:
    # "Can't have empty aa_ref and aa_alt"
    Variant(
        contig=8,
        start=141488566,
        ref="T",
        alt="C",
        ensembl=ensembl78),
    # error message:
    # "len(aa_alt) = 0"
    Variant(
        contig=11,
        start=57741870,
        ref="G",
        alt="C",
        ensembl=ensembl78),
]

def try_effect_annotation(variant):
    effect = variant.top_effect()
    assert effect is not None

def test_crashing_variants():
    for variant in should_not_crash_variants:
        yield (try_effect_annotation, variant)
