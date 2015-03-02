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
    NoncodingTranscript,
    IncompleteTranscript,
    FivePrimeUTR,
    ThreePrimeUTR,
    Intronic,
    Silent,
    Insertion,
    Deletion,
    Substitution,
    StopLoss,
    StartLoss,
    PrematureStop,
    FrameShift,
    # TODO: SpliceDonor, SpliceReceptor
)
from pyensembl import EnsemblRelease

ensembl = EnsemblRelease(release=78)

def test_incomplete():
    # transcript EGFR-009 (ENST00000450046 in Ensembl 77)
    # has an incomplete 3' end
    # chrom. 7 starting at 55,109,723
    # first exon begins: ATCATTCCTTTGGGCCTAGGA

    # change the first nucleotide of the 5' UTR A>T
    variant = Variant("7", 55109723, "A", "T", ensembl=ensembl)

    transcript = ensembl.transcript_by_id("ENST00000450046")
    effect = variant.transcript_effect(transcript)
    assert isinstance(effect, IncompleteTranscript), \
        "Expected %s on %s to be IncompleteTranscript, got %s" % (
            variant, transcript, effect)

def test_start_loss():
    # transcript EGFR-005 (ENST00000420316 in Ensembl 77)
    # location: chrom 7 @ 55,019,034-55,156,951 forward strand

    # CDS starts at position 244 of the first exon,
    # which is 55,019,034 + 244 of chr7 = 55019278
    # change the first nucleotide of the 5' UTR A>T
    # making what used to be a start codon into TTG (Leucine)
    variant = Variant("7", 55019278, "A", "T", ensembl=ensembl)
    transcript = ensembl.transcript_by_id("ENST00000420316")
    effect = variant.transcript_effect(transcript)
    assert isinstance(effect, StartLoss), \
        "Expected StartLoss, got %s for %s on %s" % (
            effect, variant, transcript, )
