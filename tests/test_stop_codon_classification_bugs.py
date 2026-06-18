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
* https://github.com/openvax/varcode/issues/394 - In-frame deletion that
  removes the stop codon of a transcript with no 3' UTR sequence crashes
  while constructing a StopLoss with an empty aa_alt; it should be a
  Deletion. (The earlier #246 only fixed the non-empty-UTR variant.)
"""

from varcode import Variant
from varcode.effects import Silent, StopLoss, Deletion
from varcode.effects.effect_prediction_coding_in_frame import (
    predict_in_frame_coding_effect,
)


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


# -----------------------------------------------------------------------
# Issue #394: in-frame deletion of the stop codon on a transcript with an
# empty 3' UTR.
#
# The reported variant is a large in-frame deletion in MAPK3 on GRCh37
# (Ensembl 75) that removes the stop codon. It overlaps several MAPK3
# transcripts, which gives us both sides of the fix from a single variant:
#
#   * ENST00000395199 (MAPK3-006) has NO 3' UTR sequence, so there's
#     nothing to translate into once the stop codon is gone — aa_alt ends
#     up empty and StopLoss can't be constructed. Correct answer: a
#     C-terminal Deletion. This is the case that used to raise
#     "If no amino acids added by StopLoss then it should be Silent".
#   * ENST00000403394 has an 804nt 3' UTR, so the same deletion reads
#     through and adds residues — that one must stay a StopLoss, proving
#     the fix didn't over-broaden the Deletion fallback.
#
# The reporter's start coordinate was off by one against the reference
# genome; the equivalent variant that matches GRCh37 ref starts at
# 30128151 (normalized to a pure 18nt deletion at 30128152).
# -----------------------------------------------------------------------


def _mapk3_stop_deletion_variant():
    return Variant(
        contig="16",
        start=30128151,
        ref="GGGATGCCTACGTGCCCCC",
        alt="G",
        genome="GRCh37",
    )


def _coding_effect(effect):
    """Unwrap a SpliceOutcomeSet to the coding effect if splicing is
    unchanged; otherwise return the effect itself."""
    return getattr(effect, "effect_if_splicing_unchanged", effect)


def test_394_stop_deletion_with_empty_3p_utr_is_deletion():
    variant = _mapk3_stop_deletion_variant()
    transcript = variant.ensembl.transcript_by_id("ENST00000395199")
    # precondition for the bug: this transcript has no 3' UTR sequence
    assert transcript.three_prime_utr_sequence == ""

    # used to raise ValueError before the fix
    effect = _coding_effect(variant.effect_on_transcript(transcript))
    assert effect.__class__ is Deletion, \
        "Expected Deletion, got %s (%s)" % (
            effect.__class__.__name__,
            getattr(effect, "short_description", effect))
    assert effect.aa_ref == "GGT", \
        "Expected aa_ref='GGT', got %r" % effect.aa_ref
    assert effect.aa_alt == "", \
        "Expected empty aa_alt, got %r" % effect.aa_alt
    assert effect.short_description == "p.GGT354del", \
        "Expected p.GGT354del, got %r" % effect.short_description


def test_394_full_effects_call_does_not_raise():
    # The crux of the bug: predicting effects across *all* MAPK3
    # transcripts (including ENST00000395199) must not raise.
    variant = _mapk3_stop_deletion_variant()
    effects = variant.effects()
    assert len(effects) > 0


def test_394_same_deletion_with_nonempty_3p_utr_stays_stoploss():
    # Same variant, different transcript: a non-empty 3' UTR means the
    # stop-loss reads through and adds residues, so this must remain a
    # StopLoss (guards against over-broadening the Deletion fallback).
    variant = _mapk3_stop_deletion_variant()
    transcript = variant.ensembl.transcript_by_id("ENST00000403394")
    assert len(transcript.three_prime_utr_sequence) > 0

    effect = _coding_effect(variant.effect_on_transcript(transcript))
    assert effect.__class__ is StopLoss, \
        "Expected StopLoss, got %s (%s)" % (
            effect.__class__.__name__,
            getattr(effect, "short_description", effect))
    assert len(effect.aa_alt) > 0, \
        "StopLoss should have added residues, got aa_alt=%r" % effect.aa_alt


def test_394_predict_in_frame_empty_3p_utr_is_deletion_unit():
    # Splice-free unit test of the exact branch in
    # predict_in_frame_coding_effect: deleting the final residue codon plus
    # the stop codon (in-frame) on a transcript with no 3' UTR.
    variant = _mapk3_stop_deletion_variant()
    transcript = variant.ensembl.transcript_by_id("ENST00000395199")
    assert transcript.three_prime_utr_sequence == ""

    start = transcript.first_start_codon_spliced_offset
    sequence_from_start_codon = str(transcript.sequence[start:])
    # length of the CDS including the stop codon
    cds_plus_stop_len = (len(transcript.protein_sequence) + 1) * 3

    # delete the last residue codon + stop codon (6nt, stays in-frame)
    cds_offset = cds_plus_stop_len - 6
    trimmed_cdna_ref = sequence_from_start_codon[cds_offset:cds_offset + 6]

    effect = predict_in_frame_coding_effect(
        variant=variant,
        transcript=transcript,
        trimmed_cdna_ref=trimmed_cdna_ref,
        trimmed_cdna_alt="",
        sequence_from_start_codon=sequence_from_start_codon,
        cds_offset=cds_offset)
    assert effect.__class__ is Deletion, \
        "Expected Deletion, got %s (%s)" % (
            effect.__class__.__name__,
            getattr(effect, "short_description", effect))
    assert effect.aa_ref == "T"
    assert effect.aa_alt == ""
