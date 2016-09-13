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

"""
Helper functions for determine effect annotation for a variant
"""

from __future__ import print_function, division, absolute_import

from ..nucleotides import PURINE_NUCLEOTIDES, AMINO_NUCLEOTIDES

def variant_overlaps_interval(
        variant_start,
        n_ref_bases,
        interval_start,
        interval_end):
    """
    Does a variant overlap a given interval on the same chromosome?

    Parameters
    ----------
    variant_start : int
        Inclusive base-1 position of variant's starting location
        (or location before an insertion)

    n_ref_bases : int
        Number of reference bases affect by variant (used to compute
        end coordinate or determine whether variant is an insertion)

    interval_start : int
        Interval's inclusive base-1 start position

    interval_end : int
        Interval's inclusive base-1 end position
    """

    if n_ref_bases == 0:
        # insertions only overlap intervals which start before and
        # end after the insertion point, they must be fully contained
        # by the other interval
        return interval_start <= variant_start and interval_end >= variant_start
    variant_end = variant_start + n_ref_bases
    """
    if self._changes_exonic_splice_site(
            strand_ref,
            strand_alt,)
    """
    # overlap means other interval starts before this variant ends
    # and the interval ends after this variant starts
    return interval_start <= variant_end and interval_end >= variant_start


def matches_exon_end_pattern(seq):
    """Does the end of the nucleotide string `seq` match the canonical splice
    signal for the 3' end of an exon: "MAG", where M is either amino base.
    """
    if len(seq) < 3:
        return False
    return seq[-3] in AMINO_NUCLEOTIDES and seq[-2] == "A" and seq[-1] == "G"

def changes_exonic_splice_site(
        transcript_offset,
        transcript,
        transcript_ref,
        transcript_alt,
        exon_start_offset,
        exon_end_offset,
        exon_number):
    """Does the given exonic mutation of a particular transcript change a
    splice site?

    Parameters
    ----------
    transcript_offset : int
        Offset from start of transcript of first reference nucleotide
        (or the last nucleotide before an insertion)

    transcript : pyensembl.Transcript

    transcript_ref : str
        Reference nucleotides

    transcript_alt : alt
        Alternate nucleotides

    exon_start_offset : int
        Start offset of exon relative to beginning of transcript

    exon_end_offset : int
        End offset of exon relative to beginning of transcript

    exon_number : int
        Which exon in the order they form the transcript
    """
    # first we're going to make sure the variant doesn't disrupt the
    # splicing sequences we got from Divina et. al's
    #   Ab initio prediction of mutation-induced cryptic
    #   splice-site activation and exon skipping
    # (http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2947103/)
    #
    # 5' splice site: MAG|GURAGU consensus
    #   M is A or C; R is purine; | is the exon-intron boundary
    #
    # 3' splice site: YAG|R
    #
    if exon_number > 1 and transcript_offset == exon_start_offset:
        # if this is any exon past the first, check to see if it lost
        # the purine on its left side
        #
        # the 3' splice site sequence has just a single purine on
        # the exon side
        if len(transcript_ref) > 0 and transcript_ref[0] in PURINE_NUCLEOTIDES:
            if len(transcript_alt) > 0:
                if transcript_alt[0] not in PURINE_NUCLEOTIDES:
                    return True
            else:
                # if the mutation is a deletion, are there ref nucleotides
                # afterward?
                offset_after_deletion = transcript_offset + len(transcript_ref)
                if len(transcript.sequence) > offset_after_deletion:
                    next_base = transcript.sequence[offset_after_deletion]
                    if next_base not in PURINE_NUCLEOTIDES:
                        return True

    if exon_number < len(transcript.exons):
        # if the mutation affects an exon whose right end gets spliced
        # to a next exon, check if the variant alters the exon side of
        # 5' consensus splicing sequence
        #
        # splicing sequence:
        #   MAG|GURAGU
        # M is A or C; R is purine; | is the exon-intron boundary
        #
        # TODO: check for overlap of two intervals instead of just
        # seeing if the mutation starts inside the exonic splice site
        if variant_overlaps_interval(
                variant_start=transcript_offset,
                n_ref_bases=len(transcript_ref),
                interval_start=exon_end_offset - 2,
                interval_end=exon_end_offset):
            end_of_reference_exon = transcript.sequence[
                exon_end_offset - 2:exon_end_offset + 1]

            if matches_exon_end_pattern(end_of_reference_exon):
                # if the last three nucleotides conform to the consensus
                # sequence then treat any deviation as an ExonicSpliceSite
                # mutation
                end_of_variant_exon = end_of_reference_exon
                if matches_exon_end_pattern(end_of_variant_exon):
                    # end of exon matches splicing signal, check if it still
                    # does after the mutation
                    return True
