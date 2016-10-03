# Copyright (c) 2015. Mount Sinai School of Medicine
#
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
The different effects get tested in an ad-hoc manner throughout the
unit test but the goal of this test module is to make sure that there is
at least one test for each effect class
"""

from varcode import Variant
from varcode.effects import (
    IncompleteTranscript,
    NoncodingTranscript,
    FivePrimeUTR,
    ThreePrimeUTR,
    Intronic,
    Silent,
    Insertion,
    Deletion,
    Substitution,
    StopLoss,
    StartLoss,
    AlternateStartCodon,
    PrematureStop,
    FrameShift,
    ExonLoss,
    ExonicSpliceSite,
    FrameShiftTruncation,
    # TODO: SpliceDonor, SpliceReceptor
)
from pyensembl import ensembl_grch37, cached_release

from .common import expect_effect

# tried using more recent releases but found that many of them
# are very specific to Ensembl data between releases 77-81
ensembl_grch38 = cached_release(81)

def test_incomplete():
    # transcript EGFR-009 (ENST00000450046 in Ensembl 78)
    # has an incomplete 3' end
    # chrom. 7 starting at 55,109,723
    # first exon begins: ATCATTCCTTTGGGCCTAGGA

    # change the first nucleotide of the 5' UTR A>T
    variant = Variant("7", 55109723, "A", "T", ensembl=ensembl_grch38)
    expect_effect(
        variant,
        transcript_id="ENST00000450046",
        effect_class=IncompleteTranscript,
        modifies_coding_sequence=False,
        modifies_protein_sequence=False)


def test_noncoding_polymorphic_pseudogene():
    # variant in MROH5-001, which is a polymorphic pseudogene
    variant = Variant("8", 142458077, "C", "T", ensembl_grch37)
    expect_effect(
        variant,
        transcript_id="ENST00000430863",
        effect_class=NoncodingTranscript,
        modifies_coding_sequence=False,
        modifies_protein_sequence=False)

def test_start_loss():
    # transcript EGFR-005 (ENST00000420316 in Ensembl 77)
    # location: chrom 7 @ 55,019,034-55,156,951 forward strand

    # CDS starts at position 244 of the first exon,
    # which is 55,019,034 + 244 of chr7 = 55019278
    # change the first two nucleotides of the 5' UTR AT>GG
    # making what used to be a start codon into GGG (Glycine)
    variant = Variant("7", 55019278, "AT", "GG", ensembl=ensembl_grch38)
    expect_effect(
        variant,
        transcript_id="ENST00000420316",
        effect_class=StartLoss,
        modifies_coding_sequence=True,
        modifies_protein_sequence=True)

def test_alternate_start_codon():
    # transcript EGFR-005 (ENST00000420316 in Ensembl 77)
    # location: chrom 7 @ 55,019,034-55,156,951 forward strand

    # CDS starts at position 244 of the first exon,
    # which is 55,019,034 + 244 of chr7 = 55019278
    # change the first nucleotide of the 5' UTR A>T
    # making what used to be the standard start codon ATG into TTG,
    # which normally codes for Leucine but can be used as an alternate
    # start codon
    variant = Variant("7", 55019278, "A", "T", ensembl=ensembl_grch38)
    expect_effect(
        variant,
        transcript_id="ENST00000420316",
        effect_class=AlternateStartCodon,
        modifies_coding_sequence=True,
        modifies_protein_sequence=False)

def test_stop_loss():
    # transcript MAST2-001 (ENST00000361297 in Ensembl 75)
    # location: chrom 1 @ 46,501,738 forward strand

    # change G>C in stop codon, should result in stop-loss mutation
    # causing an elongated protein
    variant = Variant("1", 46501738, "G", "C", ensembl=ensembl_grch37)
    expect_effect(
        variant,
        transcript_id="ENST00000361297",
        effect_class=StopLoss,
        modifies_coding_sequence=True,
        modifies_protein_sequence=True)

def test_stop_loss_from_larger_deletion_before_stop_codon():
    # transcript MAST2-001 (ENST00000361297 in Ensembl 75)
    # location: chrom 1 @ 46,501,733 forward strand

    # delete stop codon and the codon before it,
    # should result in stop-loss mutation
    # causing an elongated protein
    variant = Variant("1", 46501733, "ACATAG", "", ensembl=ensembl_grch37)
    expect_effect(
        variant,
        transcript_id="ENST00000361297",
        effect_class=StopLoss,
        modifies_coding_sequence=True,
        modifies_protein_sequence=True)

def test_stop_loss_from_larger_deletion_after_stop_codon():
    # transcript MAST2-001 (ENST00000361297 in Ensembl 75)
    # location: chrom 1 @ 46,501,736 forward strand

    # delete stop codon and the codon after it,
    # should result in stop-loss mutation
    # causing an elongated protein
    variant = Variant("1", 46501736, "TAGCAG", "", ensembl=ensembl_grch37)
    expect_effect(
        variant,
        transcript_id="ENST00000361297",
        effect_class=StopLoss,
        modifies_coding_sequence=True,
        modifies_protein_sequence=True)

def test_stop_loss_from_out_of_frame_deletion_in_stop_codon():
    # transcript MAST2-001 (ENST00000361297 in Ensembl 75)
    # location: chrom 1 @ 46,501,736 forward strand
    #
    # delete first two nucleotides of stop codon
    # TAG CAG... -> GCA G...
    #
    # should result in stop-loss mutation causing an elongated protein
    variant = Variant("1", 46501736, "TA", "", ensembl=ensembl_grch37)
    expect_effect(
        variant,
        transcript_id="ENST00000361297",
        effect_class=StopLoss,
        modifies_coding_sequence=True,
        modifies_protein_sequence=True)


def test_silent_from_in_frame_deletion_in_stop_codon():
    # transcript MAST2-001 (ENST00000361297 in Ensembl 75)
    # location: chrom 1 @ 46,501,737 forward strand
    #
    # delete first two nucleotides of stop codon
    # T_AG C_AG... -> TAG...
    #
    # should result in stop-loss mutation causing an elongated protein
    variant = Variant("1", 46501737, "AGC", "", ensembl=ensembl_grch37)
    expect_effect(
        variant,
        transcript_id="ENST00000361297",
        effect_class=Silent,
        modifies_coding_sequence=True,
        modifies_protein_sequence=False)

def test_silent_from_out_of_frame_deletion_in_stop_codon():
    # transcript MAST2-001 (ENST00000361297 in Ensembl 75)
    # location: chrom 1 @ 46,501,738 forward strand
    #
    # delete first two nucleotides of stop codon
    # TA_G C_AG... -> TAA G...
    #
    # should result in stop-loss mutation causing an elongated protein
    variant = Variant("1", 46501738, "GC", "", ensembl=ensembl_grch37)
    expect_effect(
        variant,
        transcript_id="ENST00000361297",
        effect_class=Silent,
        modifies_coding_sequence=True,
        modifies_protein_sequence=False)

def test_silent_from_out_of_frame_insertion_in_stop_codon():
    # transcript MAST2-001 (ENST00000361297 in Ensembl 75)
    # location: chrom 1 @ 46,501,737 forward strand
    #
    # insert "A" after first two nucleotides of stop codon
    # TAG CAG... -> TAA GCA G...
    #
    # should result in stop-loss mutation causing an elongated protein
    variant = Variant("1", 46501737, "AG", "AAG", ensembl=ensembl_grch37)
    expect_effect(
        variant,
        transcript_id="ENST00000361297",
        effect_class=Silent,
        modifies_coding_sequence=True,
        modifies_protein_sequence=False)

def test_stop_gain():
    # transcript BBRCA1-001 ENST00000357654 (looked up Ensembl 79)
    # Chromosome 17: 43,044,295-43,125,370 reverse strand.
    #
    # Let's look at exon #12
    # ENSE00003527960 43,082,575  43,082,404  start_phase = 0
    #
    # We can insert a reverse complement stop codon near the beginning since
    # the phase is 0.
    variant = Variant(
        "17",
        43082575 - 6,
        ref="",
        alt="CTA",
        ensembl=ensembl_grch38)
    expect_effect(
        variant,
        transcript_id="ENST00000357654",
        effect_class=PrematureStop,
        modifies_coding_sequence=True,
        modifies_protein_sequence=True)

def test_stop_gain_with_extra_amino_acids():
    # transcript BBRCA1-001 ENST00000357654 (looked up Ensembl 79)
    # Chromosome 17: 43,044,295-43,125,370 reverse strand.
    #
    # Let's look at exon #12
    # ENSE00003527960 43,082,575  43,082,404  start_phase = 0
    #
    # We can insert a reverse complement AAA followed by a stop codon
    # at near beginning since the phase is 0.
    variant = Variant(
        "17",
        43082575 - 6,
        ref="",
        alt="CTAAAA",
        ensembl=ensembl_grch38)
    expect_effect(
        variant,
        transcript_id="ENST00000357654",
        effect_class=PrematureStop,
        modifies_coding_sequence=True,
        modifies_protein_sequence=True)

def test_exon_loss():
    # transcript BBRCA1-001 ENST00000357654 (looked up Ensembl 79)
    # Chromosome 17: 43,044,295-43,125,370 reverse strand.
    #
    # Deleting exon #12
    # ENSE00003527960 43,082,575  43,082,404  start_phase = 0
    #
    variant = Variant(
        "17",
        43082404,
        ref="".join([
            "CTTTTTCTGATGTGCTTTGTTCTGGATTTCGCAGGTCCTCAAGGGCAGAAGAGTCACTTATGATG",
            "GAAGGGTAGCTGTTAGAAGGCTGGCTCCCATGCTGTTCTAACACAGCTTCAGTAATTAGATTAGT",
            "TAAAGTGATGTGGTGTTTTCTGGCAAACTTGTACACGAGCAT"
        ]),
        alt="",
        ensembl=ensembl_grch38)
    expect_effect(
        variant,
        transcript_id="ENST00000357654",
        effect_class=ExonLoss,
        modifies_coding_sequence=True,
        modifies_protein_sequence=True)


def test_exonic_splice_site():
    # transcript BBRCA1-001 ENST00000357654 (looked up Ensembl 79)
    # Chromosome 17: 43,044,295-43,125,370 reverse strand.
    #
    # Deleting last nucleotide of exon #12
    # ENSE00003527960 43,082,575  43,082,404  start_phase = 0
    #
    variant = Variant(
        "17",
        43082404,
        ref="C",
        alt="",
        ensembl=ensembl_grch38)
    expect_effect(
        variant,
        transcript_id="ENST00000357654",
        effect_class=ExonicSpliceSite,
        modifies_coding_sequence=True,
        modifies_protein_sequence=True)

def test_deletion():
    # transcript BBRCA1-001 ENST00000357654 (looked up Ensembl 79)
    # Chromosome 17: 43,044,295-43,125,370 reverse strand.
    #
    # Deleting second to last codon of exon #12
    # ENSE00003527960 43,082,575  43,082,404  start_phase = 0
    #
    variant = Variant(
        "17",
        43082404 + 4,
        ref="TTC",
        alt="",
        ensembl=ensembl_grch38)
    expect_effect(
        variant,
        transcript_id="ENST00000357654",
        effect_class=Deletion,
        modifies_coding_sequence=True,
        modifies_protein_sequence=True)


def test_deletion_of_last_amino_acid():
    # transcript CFTR-001
    # last two codons of coding sequence are CTT|TAG
    variant = Variant("7", 117667103, "CTTT", "T", "GRCh38")
    expect_effect(
        variant,
        transcript_id="ENST00000003084",
        effect_class=Deletion,
        modifies_coding_sequence=True,
        modifies_protein_sequence=True,
        aa_ref="L")

def test_insertion():
    # transcript BBRCA1-001 ENST00000357654 (looked up Ensembl 79)
    # Chromosome 17: 43,044,295-43,125,370 reverse strand.
    #
    # Insert codon after first two codons of exon #12
    # ENSE00003527960 43,082,575  43,082,404  start_phase = 0
    #
    variant = Variant(
        "17",
        43082575 - 6,
        ref="",
        alt="AAA",
        ensembl=ensembl_grch38)
    expect_effect(
        variant,
        transcript_id="ENST00000357654",
        effect_class=Insertion,
        modifies_coding_sequence=True,
        modifies_protein_sequence=True)

def test_frameshift_near_start_of_BRCA1_001():
    #
    # Insertion of genomic "A" after second codon of coding sequence.
    #
    # Transcript: BRCA1-001 (ENST00000357654)
    # Manually annotated using Ensembl release 85
    #
    # Original mRNA coding sequnce:
    #   ATG GAT TTA TCT GCT CTT CGC GTT GAA GAA GTA CAA
    #   -M- -D- -L- -S- -A- -L- -A- -V- -E- -E- -V- -Q-
    #
    # After variant:
    #   ATG GAT TTT ATC TGC TCT TCG CGT TGA
    #   -M- -D- -F- -I- -C- -S- -S- -R-  *
    variant = Variant("17", 43124096 - 6, ref="", alt="A", ensembl=ensembl_grch38)
    expect_effect(
        variant,
        transcript_id="ENST00000357654",
        effect_class=FrameShift,
        modifies_coding_sequence=True,
        modifies_protein_sequence=True,
        aa_alt="FICSSR")

def test_frameshift_truncation():
    # transcript BBRCA1-001 ENST00000357654 (looked up Ensembl 84)
    # Chromosome 17: 43,044,295-43,125,370 reverse strand.
    #
    # Out of frame insertion after first codon of exon #11
    # creates immediate "TAG" stop codon
    # Inspired by rs786202998 from dbSNP, turns GAA -> TGA|A
    variant = Variant(
        "17",
        43091031,
        ref="",
        alt="A",
        ensembl=ensembl_grch38)
    expect_effect(
        variant,
        transcript_id="ENST00000357654",
        effect_class=FrameShiftTruncation,
        modifies_coding_sequence=True,
        modifies_protein_sequence=True)


def test_frameshift_truncation_in_exon_12_of_BRCA1_001():
    # transcript BBRCA1-001 ENST00000357654 (looked up Ensembl 79)
    # Chromosome 17: 43,044,295-43,125,370 reverse strand.
    #
    # Out of frame insertion after first two codons of exon #12
    # ENSE00003527960 43,082,575  43,082,404  start_phase = 0
    #
    # Original mRNA sequence for exon #12:
    #   CAG AGG GAT ACC ATG
    #   -Q- -R- -D- -T- -M-
    # After variant:
    #   CAG AGG TGA
    #   -Q- -R-  *
    variant = Variant(
        "17",
        43082575 - 6,
        ref="",
        alt="A",
        ensembl=ensembl_grch38)
    expect_effect(
        variant,
        transcript_id="ENST00000357654",
        effect_class=FrameShiftTruncation,
        modifies_coding_sequence=True,
        modifies_protein_sequence=True,
        aa_alt="")

def test_substitution():
    # transcript BBRCA1-001 ENST00000357654 (looked up Ensembl 79)
    # Chromosome 17: 43,044,295-43,125,370 reverse strand.
    #
    # Substitute second codon of exon #12 AGG > CCC (amino acid R>P)
    # ENSE00003527960 43,082,575  43,082,404  start_phase = 0
    #
    variant = Variant(
        "17",
        43082575 - 5,
        ref="CCT",
        alt="GGG",
        ensembl=ensembl_grch38)
    expect_effect(
        variant,
        transcript_id="ENST00000357654",
        effect_class=Substitution,
        modifies_coding_sequence=True,
        modifies_protein_sequence=True)

def test_silent():
    # transcript BBRCA1-001 ENST00000357654 (looked up Ensembl 79)
    # Chromosome 17: 43,044,295-43,125,370 reverse strand.
    #
    # Substitute second codon of exon #12 AGG > AGA (amino acid R>R silent)
    # ENSE00003527960 43,082,575  43,082,404  start_phase = 0
    #
    variant = Variant(
        "17",
        43082575 - 5,
        ref="CCT",
        alt="TCT",
        ensembl=ensembl_grch38)
    expect_effect(
        variant,
        transcript_id="ENST00000357654",
        effect_class=Silent,
        modifies_coding_sequence=True,
        modifies_protein_sequence=False)

def test_silent_stop_codons():
    silent_stop_codon_variants = {
        "ENST00000290524": Variant(
            contig=1,
            start=151314663,
            ref="C",
            alt="T",
            ensembl=ensembl_grch37),
        "ENST00000368725": Variant(
            contig=1,
            start=153409535,
            ref="C",
            alt="T",
            ensembl=ensembl_grch37),
        "ENST00000353479": Variant(
            contig=10,
            start=105791994,
            ref="C",
            alt="T",
            ensembl=ensembl_grch37),
    }
    for transcript_id, variant in silent_stop_codon_variants.items():
        yield (
            expect_effect,
            variant,
            transcript_id,
            Silent)

def test_five_prime_utr():
    # transcript BBRCA1-001 ENST00000357654 (looked up Ensembl 79)
    # Chromosome 17: 43,044,295-43,125,370 reverse strand.
    #
    # Exon #1 is the beginning of the 5' UTR
    # 1 ENSE00001871077 43,125,370  43,125,271  -   -   length=100
    # Sequence:
    # GAGCTCGCTGAGACTTCCTGGACGGGGGACAGGCTGTGGGGTTTCTCAGATAACTGGGCC
    # CCTGCGCTCAGGAGGCCTTCACCCTCTGCTCTGGGTAAAG
    variant = Variant(
        "17",
        43125370 - 2,
        ref="CTC",
        alt="",
        ensembl=ensembl_grch38)
    expect_effect(
        variant,
        transcript_id="ENST00000357654",
        effect_class=FivePrimeUTR,
        modifies_coding_sequence=False,
        modifies_protein_sequence=False)

def test_three_prime_utr():
    # transcript BBRCA1-001 ENST00000357654 (looked up Ensembl 79)
    # Chromosome 17: 43,044,295-43,125,370 reverse strand.
    #
    # Exon #23 contains the 3' UTR
    # 23  ENSE00001814242 43,045,802  43,044,295  1   -   length=1,508
    # Sequence end with:
    # ...CACTTCCA
    variant = Variant(
        "17",
        43044295,
        ref="TGG",
        alt="",
        ensembl=ensembl_grch38)
    expect_effect(
        variant,
        transcript_id="ENST00000357654",
        effect_class=ThreePrimeUTR,
        modifies_coding_sequence=False,
        modifies_protein_sequence=False)

def test_intronic():
    # transcript BBRCA1-001 ENST00000357654 (looked up Ensembl 79)
    # Chromosome 17: 43,044,295-43,125,370 reverse strand.
    #
    # Insertion 20bp before exon #12 should be consider "Intronic"
    # #12 ENSE00003527960 43,082,575  43,082,404  start_phase = 0
    variant = Variant(
        "17",
        43082575 + 20,
        ref="",
        alt="CCC",
        ensembl=ensembl_grch38)
    expect_effect(
        variant,
        transcript_id="ENST00000357654",
        effect_class=Intronic,
        modifies_coding_sequence=False,
        modifies_protein_sequence=False)

def test_disruptive_insertion_stop_gain():
    variant = Variant(
        "7",
        117182143,
        ref="",
        alt="GTA",
        ensembl=ensembl_grch37)
    expect_effect(
        variant,
        transcript_id="ENST00000003084",
        effect_class=PrematureStop,
        modifies_coding_sequence=True,
        modifies_protein_sequence=True)

def test_substitution_nonfinal_stop_gain():
    variant = Variant(
        "7",
        117182145,
        ref="ACAG",
        alt="TAAC",
        ensembl=ensembl_grch37)
    expect_effect(
        variant,
        transcript_id="ENST00000003084",
        effect_class=PrematureStop,
        modifies_coding_sequence=True,
        modifies_protein_sequence=True)

def test_insertion_nonfinal_stop_gain():
    variant = Variant(
        "7",
        117182144,
        ref="",
        alt="TAGGTT",
        ensembl=ensembl_grch37)
    expect_effect(
        variant,
        transcript_id="ENST00000003084",
        effect_class=PrematureStop,
        modifies_coding_sequence=True,
        modifies_protein_sequence=True)
