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
Regression tests for misclassification bugs involving variants in or near
the stop codon:

* https://github.com/openvax/varcode/issues/250 - SNV in stop codon reported
  as Insertion instead of StopLoss, when the 3' UTR starts with a stop codon.
* https://github.com/openvax/varcode/issues/205 - Same root cause as #250
  on a different transcript.
* https://github.com/openvax/varcode/issues/201 - Insertion before the stop
  codon that produces an identical protein sequence is reported as an
  Insertion rather than Silent.
"""

from varcode import Variant
from varcode.effects import Silent, StopLoss


# -----------------------------------------------------------------------
# Issue #250: SNV in stop codon reported as Insertion
#
# The stop codon TGA is changed to a sense codon, but the 3' UTR begins
# with another stop codon (TGA...), so translation terminates immediately.
# Previously this was classified as Insertion; the correct answer is
# StopLoss with a single additional amino acid.
# -----------------------------------------------------------------------


def test_250_snv_in_stop_codon_with_immediate_utr_stop_is_stoploss():
    # chr16:17549386 A>G, Lrrc74b-204 / ENSMUST00000232637 on GRCm38.
    # Transcript is on the reverse strand, so the + strand A>G is a T>C
    # on the transcript, changing the stop codon TGA to CGA (Arg).
    # The 3' UTR begins with another TGA, so the new protein gains
    # exactly one amino acid (R).
    variant = Variant(
        contig=16,
        start=17549386,
        ref="A",
        alt="G",
        genome="GRCm38",
    )
    transcript = variant.ensembl.transcript_by_id("ENSMUST00000232637")
    effect = variant.effect_on_transcript(transcript)
    assert effect.__class__ is StopLoss, \
        "Expected StopLoss, got %s" % effect.__class__.__name__
    assert effect.aa_alt == "R", \
        "Expected aa_alt='R', got %r" % effect.aa_alt


def test_250b_second_snv_in_stop_codon_with_immediate_utr_stop_is_stoploss():
    # chr7:126830599 A>T from the #250 follow-up comment.
    # Same bug on 4930451I11Rik-201 / ENSMUST00000061695.
    variant = Variant(
        contig=7,
        start=126830599,
        ref="A",
        alt="T",
        genome="GRCm38",
    )
    transcript = variant.ensembl.transcript_by_id("ENSMUST00000061695")
    effect = variant.effect_on_transcript(transcript)
    assert effect.__class__ is StopLoss, \
        "Expected StopLoss, got %s" % effect.__class__.__name__


# -----------------------------------------------------------------------
# Issue #205: Same underlying bug
# -----------------------------------------------------------------------


def test_205_stop_codon_snv_is_stoploss_when_utr_has_consecutive_stops():
    # chr11:21319856 T>C, Vps54-201 / ENSMUST00000006221 on GRCm38.
    # CDS ends with TGA; 3' UTR starts with TGA. SNV changes first TGA to
    # CGA (Arg). Reporter says the correct annotation is StopLoss.
    variant = Variant(
        contig=11,
        start=21319856,
        ref="T",
        alt="C",
        genome="GRCm38",
    )
    transcript = variant.ensembl.transcript_by_id("ENSMUST00000006221")
    effect = variant.effect_on_transcript(transcript)
    assert effect.__class__ is StopLoss, \
        "Expected StopLoss, got %s" % effect.__class__.__name__
    assert effect.aa_alt == "R", \
        "Expected aa_alt='R', got %r" % effect.aa_alt


# -----------------------------------------------------------------------
# Issue #201: Insertion before stop codon that yields identical protein
#
# Insertion of ATATAA after the last C of "...TTC | ATC TGA" gives
# "...TTC ATA TAA ATC TGA". The new TAA terminates translation after F I,
# matching the original protein ending F I. The proteins are identical;
# the annotation should be Silent, not Insertion.
# -----------------------------------------------------------------------


def test_201_synonymous_insertion_before_stop_is_silent():
    variant = Variant(
        contig=1,
        start=100484695,
        ref="C",
        alt="CATATAA",
        genome="GRCm38",
    )
    transcript = variant.ensembl.transcript_by_id("ENSMUST00000086738")
    effect = variant.effect_on_transcript(transcript)
    assert effect.__class__ is Silent, \
        "Expected Silent, got %s (%s)" % (
            effect.__class__.__name__,
            getattr(effect, "short_description", effect))
