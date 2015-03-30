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

from varcode import (
    Variant,
    #
    # transcript effects
    #
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
    # TODO: SpliceDonor, SpliceReceptor
)
from pyensembl import EnsemblRelease

ensembl_grch37 = EnsemblRelease(75)
ensembl_grch38 = EnsemblRelease(79)

def expect_effect(variant, transcript_id, effect_class):
    transcript = variant.ensembl.transcript_by_id(transcript_id)
    effect = variant.effect_on_transcript(transcript)
    assert isinstance(effect, effect_class), \
        "Expected %s on %s to be %s, got %s" % (
            variant, transcript, effect_class.__name__, effect)

def test_incomplete():
    # transcript EGFR-009 (ENST00000450046 in Ensembl 78)
    # has an incomplete 3' end
    # chrom. 7 starting at 55,109,723
    # first exon begins: ATCATTCCTTTGGGCCTAGGA

    # change the first nucleotide of the 5' UTR A>T
    variant = Variant("7", 55109723, "A", "T", ensembl=ensembl_grch38)
    expect_effect(variant, "ENST00000450046", IncompleteTranscript)


def test_start_loss():
    # transcript EGFR-005 (ENST00000420316 in Ensembl 77)
    # location: chrom 7 @ 55,019,034-55,156,951 forward strand

    # CDS starts at position 244 of the first exon,
    # which is 55,019,034 + 244 of chr7 = 55019278
    # change the first two nucleotides of the 5' UTR AT>GG
    # making what used to be a start codon into GGG (Glycine)
    variant = Variant("7", 55019278, "AT", "GG", ensembl=ensembl_grch38)
    expect_effect(variant, "ENST00000420316", StartLoss)

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
    expect_effect(variant, "ENST00000420316", AlternateStartCodon)


def test_stop_loss():
    # transcript MAST2-001 (ENST00000361297 in Ensembl 75)
    # location: chrom 1 @ 46,501,738 forward strand

    # change G>C in stop codon, should result in stop-loss mutation
    # causing an elongated protein
    variant = Variant("1", 46501738, "G", "C", ensembl=ensembl_grch37)
    expect_effect(variant, "ENST00000361297", StopLoss)

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
    expect_effect(variant, "ENST00000357654", PrematureStop)

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
    expect_effect(variant, "ENST00000357654", PrematureStop)
