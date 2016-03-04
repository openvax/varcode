# Copyright (c) 2016. Mount Sinai School of Medicine
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

from __future__ import print_function, division, absolute_import

from .frameshift_coding_effect import frameshift_coding_effect
from .in_frame_coding_effect import in_frame_coding_effect


def coding_effect(
        ref,
        alt,
        transcript_offset,
        variant,
        transcript):
    """
    Given a minimal ref/alt nucleotide string pair and an offset into a given
    transcript, determine the coding effect of this nucleotide substitution
    onto the translated protein.

    Parameters
    ----------
    ref : str
        Reference nucleotides we expect to find in the transcript's CDS

    alt : str
        Alternate nucleotides we're replacing the reference with

    transcript_offset : int
        Offset into the full transcript sequence of the ref->alt substitution

    transcript : Transcript

    variant : Variant
    """
    if not transcript.complete:
        raise ValueError(
            ("Can't annotate coding effect for %s"
             " on incomplete transcript %s" % (variant, transcript)))

    sequence = transcript.sequence

    # reference nucleotides found on the transcript, if these don't match
    # what we were told to expect from the variant then raise an exception
    ref_nucleotides_from_transcript = str(
        sequence[transcript_offset:transcript_offset + len(ref)])

    # Make sure that the reference sequence agrees with what we expected
    # from the VCF
    assert ref_nucleotides_from_transcript == ref, \
        "%s: expected ref '%s' at offset %d of %s, transcript has '%s'" % (
            variant,
            ref,
            transcript_offset,
            transcript,
            ref_nucleotides_from_transcript)

    cds_start_offset = transcript.first_start_codon_spliced_offset
    cds_stop_offset = transcript.last_stop_codon_spliced_offset

    cds_len = cds_stop_offset + cds_start_offset + 1

    if cds_len < 3:
        raise ValueError(
            "Coding sequence for %s is too short: '%s'" % (
                transcript,
                transcript.sequence[cds_start_offset:cds_stop_offset + 1]))

    if len(ref) == 0 and transcript.strand == "-":
        # By convention, genomic insertions happen *after* their base 1 position on
        # a chromosome. On the reverse strand, however, an insertion has to go
        # before the nucleotide at some transcript offset.
        # Example:
        #    chromosome sequence:
        #        TTT|GATCTCGTA|CCC
        #    transcript on reverse strand:
        #        CCC|ATGCTCTAG|TTT
        #    where the CDS is emphasized:
        #            ATGCTCTAG
        # If we have a genomic insertion g.6insATT
        # the genomic sequence becomes:
        #       TTT|GAT_ATT_CTCGTA|CCC
        # (insert the "ATT" after the "T" at position 6)
        # On the reverse strand this becomes:
        #       CCC|ATGCTC_TTA_TAG|TTT
        # (insert the "ATT" *before* the "T" at position 10)
        #
        # To preserve the interpretation of the start offset as the base
        # before the insertion, need to subtract one
        cds_offset = transcript_offset - cds_start_offset - 1
    else:
        cds_offset = transcript_offset - cds_start_offset

    assert cds_offset < cds_len, \
        "Expected CDS offset (%d) < |CDS| (%d) for %s on %s" % (
            cds_offset, cds_len, variant, transcript)

    sequence_from_start_codon = str(sequence[cds_start_offset:])

    # is this an in-frame mutations?
    if (len(ref) - len(alt)) % 3 == 0:
        return in_frame_coding_effect(
            ref=ref,
            alt=alt,
            cds_offset=cds_offset,
            sequence_from_start_codon=sequence_from_start_codon,
            variant=variant,
            transcript=transcript)
    else:
        return frameshift_coding_effect(
            ref=ref,
            alt=alt,
            cds_offset=cds_offset,
            sequence_from_start_codon=sequence_from_start_codon,
            variant=variant,
            transcript=transcript)
