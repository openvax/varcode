"""
Any variants which are encountered in the wild and either cause Varcode
to crash or return an incorrect annotation should be added to this
test module.
"""

import pytest
from varcode import Variant

from .common import check_effect_properties

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
        genome="GRCh37"),
    # error message:
    # "Expect non-silent stop-loss variant to cause longer variant protein"
    # "" but got len(original) = 653, len(variant) = 653"
    Variant(
        contig=1,
        start=167385324,
        ref="TAA",
        alt="T",
        genome="GRCh37"),
    # error message:
    # "Variant which span 5' UTR and CDS not supported"
    Variant(
        contig=19,
        start=44351166,
        ref="GGGAGAT",
        alt="G",
        genome="GRCh37"),
    # error message:
    # "Can't have ref = '' and alt = 'E' at aa_pos = 445, cds_pos = 1335"
    Variant(
        contig=1,
        start=1684347,
        ref="",
        alt="CCT",
        genome="GRCh37"),
    Variant(
        contig=11,
        start=47640416,
        ref="",
        alt="TCTTT",
        genome="GRCh37"),
    Variant(
        contig=12,
        start=98880902,
        ref="A",
        alt="",
        genome="GRCh37"),
    Variant(
        contig=19,
        start=52803670,
        ref="TG",
        alt="",
        genome="GRCh37"),
    Variant(
        contig=1,
        start=109792735,
        ref="",
        alt="CGC",
        genome="GRCh37"),
    # error message:
    # "expected ref 'GATGTCGG' at offset 1412 of ENST00000297524...CDS has 'G'"
    Variant(
        contig=8,
        start=87226635,
        ref="CCGACATC",
        alt="",
        genome="GRCh37"),
    # error message: "Can't have empty aa_ref and aa_alt"
    Variant(
        contig=8,
        start=141488566,
        ref="T",
        alt="C",
        genome="GRCh38"),
    # error message: "len(aa_alt) = 0"
    Variant(
        contig=11,
        start=57741870,
        ref="G",
        alt="C",
        genome="GRCh38"),
    # error message: "IndexError: string index out of range"
    Variant(
        contig=11,
        start=63676705,
        ref="T", alt="",
        genome="GRCh37"),
    # AssertionError: aa_ref and aa_alt can't both be empty string
    Variant(
        contig=1,
        start=56962223,
        ref='C',
        alt='T',
        genome="GRCh37"),
    # AssertionError: aa_ref and aa_alt can't both be empty string
    Variant(
        contig=1,
        start=56962223,
        ref="C",
        alt="T",
        genome="GRCh37"),
    # AssertionError: aa_ref and aa_alt can't both be empty string
    Variant(
        contig=1,
        start=151314663,
        ref="C",
        alt="T",
        genome="GRCh37"),
    # AssertionError: aa_ref and aa_alt can't both be empty string
    Variant(
        contig=1,
        start=153409535,
        ref="C",
        alt="T",
        genome="GRCh37"),
    # AssertionError: aa_ref and aa_alt can't both be empty string
    Variant(
        contig=10,
        start=105791994,
        ref="C",
        alt="T",
        genome="GRCh37"),
    # Expected frameshift_insertion to be before stop codon
    # for Variant(contig=1, start=109925189, ref=., alt=A, genome=GRCh38)
    # on transcript_id=ENST00000329608
    # len(protein) = 554, aa_pos = 554
    Variant(
        contig=1,
        start=109925189,
        ref="",
        alt="A",
        genome="GRCh38"),
    Variant(
        contig=7,
        start=117120188,
        ref="A",
        alt="AAGT",
        genome="GRCh37"),
    # had problems with end coordinate loading this one from a MAF but also
    # want to make sure it doesn't cause other trouble
    Variant(
        contig=1,
        start=109461324,
        ref="GG",
        alt="TT",
        genome="GRCh37")
]


@pytest.mark.parametrize(['variant'], [(v,) for v in should_not_crash_variants])
def test_crashing_variants(variant):
    effect = variant.effects().top_priority_effect()
    check_effect_properties(effect)