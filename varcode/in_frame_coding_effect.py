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
Effect annotation for variants which modify the coding sequence without
changing the reading frame.
"""

import logging

from .effects import (
    IncompleteTranscript,
    Silent,
    Insertion,
    Deletion,
    Substitution,
    ComplexSubstitution,
    PrematureStop,
    AlternateStartCodon,
    StartLoss,
    StopLoss,
    FrameShift,
    FrameShiftTruncation,
    ThreePrimeUTR,
)
from .string_helpers import trim_shared_flanking_strings
from .translate import transcript_protein_sequence, START_CODONS, translate

def _choose_annotation(
        aa_pos,
        aa_ref,
        aa_alt,
        transcript,
        variant):
    """Assumption: mutations which modify the start are not handled here."""

    assert len(aa_ref) > 0 or len(aa_alt) > 0

    aa_ref, aa_alt, shared_prefix, shared_suffix = \
        trim_shared_flanking_strings(
            aa_ref,
            aa_alt)

    if len(aa_ref) == len(aa_alt) == 0:
        shared_amino_acids = shared_prefix + shared_suffix
        return Silent(
            variant=variant,
            transcript=transcript,
            aa_pos=aa_pos,
            aa_ref=shared_amino_acids)

    # index of first amino acid which is different from the reference
    aa_pos += len(shared_prefix)

    if aa_pos == len(transcript.protein_sequence):
        # if non-silent mutation is at the end of the protein then
        # should be a stop-loss
        assert aa_ref == "", \
            "Expected end of coding sequence for %s, got '%s'" % (
                aa_ref)
        return StopLoss(
            variant,
            transcript,
            aa_pos=aa_pos,
            aa_alt=aa_alt)
    elif len(aa_alt) == 0:
        return Deletion(
            variant,
            transcript,
            aa_pos=aa_pos,
            aa_ref=aa_ref)
    elif len(aa_ref) == 0:
        return Insertion(
            variant,
            transcript,
            aa_pos=aa_pos,
            aa_alt=aa_alt)
    elif len(aa_alt) == len(aa_ref) == 1:
        # simple substitution e.g. p.V600E
        return Substitution(
            variant,
            transcript,
            aa_pos=aa_pos,
            aa_ref=aa_ref,
            aa_alt=aa_alt)
    else:
        return ComplexSubstitution(
            variant,
            transcript,
            aa_pos=aa_pos,
            aa_ref=aa_ref,
            aa_alt=aa_alt)

def in_frame_coding_effect(
        ref,
        alt,
        cds_offset,
        sequence_from_start_codon,
        transcript,
        variant):
    """Coding effect of an in-frame nucleotide change

    Parameters
    ----------
    ref : str
        Reference nucleotides

    alt : str
        Nucleotides to insert in place of the reference nucleotides

    cds_offset : int
        Index of first ref nucleotide, starting from 0 = beginning of coding
        sequence. If variant is a pure insertion (no ref nucleotides) then this
        argument indicates the offset *after* which to insert the `alt`
        nucleotides.

    sequence_from_start_codon : Bio.Seq
        Transcript sequence from the CDS start codon (including the 3' UTR).
        This sequence includes the 3' UTR since a mutation may delete the stop
        codon and we'll have to translate past the normal end of the CDS to
        determine the new protein sequence.

    transcript : Transcript

    variant : Variant
    """

    n_ref_nucleotides = len(ref)

    if n_ref_nucleotides == 0:
        # since insertions interact with the "base counting, one start"
        # coordinates used by VCF & Ensembl to create lots of special cases
        # let's handle all logic for in-frame insertions separately
        return in_frame_insertion_effect(
            inserted_nucleotides=alt,
            cds_offset_before_insertion=cds_offset,
            transcript=transcript,
            variant=variant)

    original_protein_sequence = transcript_protein_sequence(transcript)

    first_ref_codon_index = int(cds_offset / 3)

    assert first_ref_codon_index <= len(original_protein_sequence), \
        ("Unexpected mutation at offset %d (5' UTR starts at %d"
         " while annotating %s on %s") % (
         first_ref_codon_index,
         len(transcript.protein_sequence))

    last_ref_codon_index = int((cds_offset + n_ref_nucleotides - 1) / 3)

    assert last_ref_codon_index >= first_ref_codon_index, \
        ("Expected first_ref_amino_acid_index (%d) <="
         "last_ref_amino_acid_index (%d) while annotating %s on %s") % (
         first_ref_codon_index,
         last_ref_codon_index,
         variant,
         transcript)

    # codons in the reference sequence
    ref_codons = sequence_from_start_codon[
        first_ref_codon_index * 3:last_ref_codon_index * 3 + 3]

    # We construct the new codons by taking the unmodified prefix
    # of the first ref codon, the unmodified suffix of the last ref codon
    # and sticking the alt nucleotides in between.
    # Since this is supposed to be an in-frame mutation, the concatenated
    # nucleotide string is expected to have a length that is a multiple of
    # three.

    # which nucleotide of the codon got changed?
    offset_in_first_ref_codon = cds_offset % 3
    prefix = ref_codons[:offset_in_first_ref_codon]

    offset_in_last_ref_codon = (cds_offset + len(ref) - 1) % 3
    if offset_in_last_ref_codon == 0:
        suffix = ref_codons[-2:]
    elif offset_in_last_ref_codon == 1:
        suffix = ref_codons[-1:]
    else:
        suffix = ""

    mutant_codons = prefix + alt + suffix

    assert len(mutant_codons) % 3 == 0, \
        "Expected in-frame mutation but got %s (length = %d)" % (
            mutant_codons, len(mutant_codons))

    original_protein_subsequence = transcript.protein_sequence[
        first_ref_codon_index:last_ref_codon_index + 1]

    mutant_protein_subsequence = translate(
        mutant_codons,
        first_codon_is_start=(first_ref_codon_index == 0))

    if first_ref_codon_index == 0:
        if len(mutant_codons) == 3 and mutant_codons in START_CODONS:
            assert len(original_protein_subsequence) == 1
            assert len(mutant_protein_subsequence) == 1
            # if we changed the starting codon treat then
            # this is technically a Silent mutation but may cause
            # alternate starts or other effects
            return AlternateStartCodon(
                variant,
                transcript,
                aa_ref=original_protein_subsequence,
                aa_alt=mutant_protein_subsequence)
        elif mutant_codons[:3] not in START_CODONS:
            # if we changed a start codon to something else then
            # we no longer know where the protein begins (or even in
            # what frame).
            # TODO: use the Kozak consensus sequence or a predictive model
            # to identify the most likely start site
            return StartLoss(
                variant=variant,
                transcript=transcript,
                aa_alt=mutant_protein_subsequence)

    return _choose_annotation(
        aa_pos=first_ref_codon_index,
        aa_ref=original_protein_subsequence,
        aa_alt=mutant_protein_subsequence,
        variant=variant,
        transcript=transcript)


def in_frame_insertion_effect(
        inserted_nucleotides,
        cds_offset_before_insertion,
        transcript,
        variant):

    assert len(inserted_nucleotides) > 0, \
        "Expected len(inserted_nucleotides) > 0 for %s on %s" % (
            transcript, variant)
    assert False
    """
    assert cds_offset_before_insertion < cds_len, \
        "Expected CDS offset (%d) < |CDS| (%d) for %s on %s" % (
            cds_offset, cds_len, variant, transcript)

    if start_codon:
        # IF insertion into start codon, StartLoss
        pass

    if stop_codon:
        # IF insertion into stop codon, StopLoss (scan forward to next stop?)
        pass

    pass
    """
