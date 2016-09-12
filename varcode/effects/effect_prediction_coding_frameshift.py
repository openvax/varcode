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
Effect annotation for variants which modify the coding sequence and change
reading frame.
"""

from __future__ import print_function, division, absolute_import

from ..string_helpers import trim_shared_prefix

from .effect_classes import (
    FrameShift,
    FrameShiftTruncation,
    StartLoss,
    StopLoss,
    Silent
)
from .mutate import substitute
from .translate import translate

def create_frameshift_effect(
        mutated_codon_index,
        sequence_from_mutated_codon,
        variant,
        transcript):
    """
    Determine frameshift effect within a coding sequence (possibly affecting
    either the start or stop codons, or anythign in between)

    Parameters
    ----------
    mutated_codon_index : int
        Codon offset (starting from 0 = start codon) of first non-reference
        amino acid in the variant protein

    sequence_from_mutated_codon: Bio.Seq
        Sequence of mutated cDNA, starting from first mutated codon, until
        the end of the transcript

    variant : Variant

    transcript : transcript
    """

    assert transcript.protein_sequence is not None, \
        "Expect transcript %s to have protein sequence" % transcript

    original_protein_sequence = transcript.protein_sequence
    original_protein_length = len(original_protein_sequence)

    mutant_protein_suffix = translate(
        nucleotide_sequence=sequence_from_mutated_codon,
        first_codon_is_start=False,
        to_stop=True,
        truncate=True)

    if mutated_codon_index == 0:
        # TODO: scan through sequence_from_mutated_codon for
        # Kozak sequence + start codon to choose the new start
        return StartLoss(variant=variant, transcript=transcript)

    # the frameshifted sequence may contain some amino acids which are
    # the same as the original protein!
    _, mutant_protein_suffix, unchanged_amino_acids = trim_shared_prefix(
        ref=original_protein_sequence[mutated_codon_index:],
        alt=mutant_protein_suffix)
    n_unchanged_amino_acids = len(unchanged_amino_acids)
    offset_to_first_different_amino_acid = mutated_codon_index + n_unchanged_amino_acids
    if offset_to_first_different_amino_acid >= original_protein_length:
        # frameshift is either extending the protein or leaving it unchanged
        if len(mutant_protein_suffix) == 0:
            # miraculously, this frameshift left the protein unchanged,
            # most likely by turning one stop codon into another stop codon
            aa_ref = original_protein_sequence[-n_unchanged_amino_acids:]
            return Silent(
                variant=variant,
                transcript=transcript,
                aa_pos=mutated_codon_index,
                aa_ref=aa_ref)
        else:
            # When all the amino acids are the same as the original, we either
            # have the original protein or we've extended it.
            # If we've extended it, it means we must have lost our stop codon.
            return StopLoss(
                variant=variant,
                transcript=transcript,
                extended_protein_sequence=mutant_protein_suffix)
    # original amino acid at the mutated codon before the frameshift occurred
    aa_ref = original_protein_sequence[offset_to_first_different_amino_acid]

    # TODO: what if all the shifted amino acids were the same and the protein
    # ended up the same length? Add a Silent case?
    if len(mutant_protein_suffix) == 0:
        # if a frameshift doesn't create any new amino acids, then
        # it must immediately have hit a stop codon
        return FrameShiftTruncation(
            variant=variant,
            transcript=transcript,
            stop_codon_offset=offset_to_first_different_amino_acid)
    return FrameShift(
        variant=variant,
        transcript=transcript,
        aa_mutation_start_offset=offset_to_first_different_amino_acid,
        shifted_sequence=str(mutant_protein_suffix))

def cdna_codon_sequence_after_insertion_frameshift(
        sequence_from_start_codon,
        cds_offset_before_insertion,
        inserted_nucleotides):
    """
    Returns index of mutated codon and nucleotide sequence starting at the first
    mutated codon.
    """
    # special logic for insertions
    coding_sequence_after_insertion = \
        sequence_from_start_codon[cds_offset_before_insertion + 1:]

    if cds_offset_before_insertion % 3 == 2:
        # insertion happens after last nucleotide in a codon,
        # doesn't disrupt the existing codon from cds_offset-2 to cds_offset
        mutated_codon_index = cds_offset_before_insertion // 3 + 1
        nucleotides_before = ""
    elif cds_offset_before_insertion % 3 == 1:
        # insertion happens after 2nd nucleotide of a codon
        # codon positions:
        #   1) cds_offset - 1
        #   2) cds_offset
        #    <----- Insertsion
        #   3) cds_offset + 1
        mutated_codon_index = cds_offset_before_insertion // 3
        # the first codon in the returned sequence will contain two reference
        # nucleotides before the insertion
        nucleotides_before = sequence_from_start_codon[
            cds_offset_before_insertion - 1:cds_offset_before_insertion + 1]
    elif cds_offset_before_insertion % 3 == 0:
        # insertion happens after 1st nucleotide of a codon
        # codon positions:
        #   1) cds_offset
        #    <----- Insertsion
        #   2) cds_offset + 1
        #   3) cds_offset + 2
        mutated_codon_index = cds_offset_before_insertion // 3
        # the first codon in the returned sequence will contain one reference
        # nucleotide before the insertion
        nucleotides_before = sequence_from_start_codon[cds_offset_before_insertion]
    sequence_from_mutated_codon = (
        nucleotides_before +
        inserted_nucleotides +
        coding_sequence_after_insertion)
    return mutated_codon_index, sequence_from_mutated_codon


def cdna_codon_sequence_after_deletion_or_substitution_frameshift(
        sequence_from_start_codon,
        cds_offset,
        trimmed_cdna_ref,
        trimmed_cdna_alt):
    """
    Logic for any frameshift which isn't an insertion.

    We have insertions as a special case since our base-inclusive
    indexing means something different for insertions:
       cds_offset = base before insertion
    Whereas in this case:
      cds_offset = first reference base affected by a variant

    Returns index of first modified codon and sequence from that codon
    onward.
    """
    mutated_codon_index = cds_offset // 3
    # get the sequence starting from the first modified codon until the end
    # of the transcript.
    sequence_after_mutated_codon = \
        sequence_from_start_codon[mutated_codon_index * 3:]

    # the variant's ref nucleotides should start either 0, 1, or 2 nucleotides
    # into `sequence_after_mutated_codon`
    offset_into_mutated_codon = cds_offset % 3

    sequence_from_mutated_codon = substitute(
        sequence=sequence_after_mutated_codon,
        offset=offset_into_mutated_codon,
        ref=trimmed_cdna_ref,
        alt=trimmed_cdna_alt)
    return mutated_codon_index, sequence_from_mutated_codon

def predict_frameshift_coding_effect(
        variant,
        transcript,
        trimmed_cdna_ref,
        trimmed_cdna_alt,
        cds_offset,
        sequence_from_start_codon):
    """
    Coding effect of a frameshift mutation.

    Parameters
    ----------
    variant : Variant

    transcript : Transcript

    trimmed_cdna_ref : nucleotide sequence
        Reference nucleotides in the coding sequence of the given transcript.

    trimmed_cdna_alt : nucleotide sequence
        Alternate nucleotides introduced by mutation

    cds_offset : int
        Offset into the CDS of first ref nucleotide. For insertions, this
        is the offset of the last ref nucleotide before the insertion.

    sequence_from_start_codon : nucleotide sequence
        Nucleotides of the coding sequence and 3' UTR

    """
    if len(trimmed_cdna_ref) != 0:
        mutated_codon_index, sequence_from_mutated_codon = \
            cdna_codon_sequence_after_deletion_or_substitution_frameshift(
                sequence_from_start_codon=sequence_from_start_codon,
                cds_offset=cds_offset,
                trimmed_cdna_ref=trimmed_cdna_ref,
                trimmed_cdna_alt=trimmed_cdna_alt)
    else:
        mutated_codon_index, sequence_from_mutated_codon = \
            cdna_codon_sequence_after_insertion_frameshift(
                sequence_from_start_codon=sequence_from_start_codon,
                cds_offset_before_insertion=cds_offset,
                inserted_nucleotides=trimmed_cdna_alt)
    return create_frameshift_effect(
        mutated_codon_index=mutated_codon_index,
        sequence_from_mutated_codon=sequence_from_mutated_codon,
        variant=variant,
        transcript=transcript)
