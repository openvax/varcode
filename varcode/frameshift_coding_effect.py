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
Effect annotation for variants which modify the coding sequence and change
reading frame.
"""

from .effects import (
    FrameShift,
    FrameShiftTruncation,
    StartLoss,
    StopLoss,
    Silent
)
from .mutate import substitute
from .string_helpers import trim_shared_prefix
from .translate import translate

def _frameshift(
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
    elif mutated_codon_index == original_protein_length:
        return StopLoss(
            variant=variant,
            transcript=transcript,
            extended_protein_sequence=mutant_protein_suffix)

    # the frameshifted sequence may contain some amino acids which are
    # the same as the original protein!
    _, mutant_protein_suffix, unchanged_amino_acids = trim_shared_prefix(
        ref=original_protein_sequence[mutated_codon_index:],
        alt=mutant_protein_suffix)
    n_unchanged_amino_acids = len(unchanged_amino_acids)
    mutation_start_position = mutated_codon_index + n_unchanged_amino_acids
    if mutation_start_position >= original_protein_length:
        # frameshift is either extending the protein or leaving it unchanged
        if len(mutant_protein_suffix) == 0:
            # miraculously, this frameshift left the protein unchanged,
            # most likely by turning one stop codon into another stop codon
            aa_ref = original_protein_sequence[mutated_codon_index]
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
    aa_ref = original_protein_sequence[mutation_start_position]

    # TODO: what if all the shifted amino acids were the same and the protein
    # ended up the same length? Add a Silent case?
    if len(mutant_protein_suffix) == 0:
        # if a frameshift doesn't create any new amino acids, then
        # it must immediately have hit a stop codon
        return FrameShiftTruncation(
            variant=variant,
            transcript=transcript,
            stop_codon_offset=mutation_start_position,
            aa_ref=aa_ref)
    return FrameShift(
        variant=variant,
        transcript=transcript,
        aa_mutation_start_offset=mutation_start_position,
        aa_ref=aa_ref,
        shifted_sequence=mutant_protein_suffix)


def frameshift_coding_effect(
        ref,
        alt,
        cds_offset,
        sequence_from_start_codon,
        variant,
        transcript):
    """
    Coding effect of a frameshift mutation.

    Parameters
    ----------
    ref : nucleotide sequence
        Reference nucleotides

    alt : nucleotide sequence
        Alternate nucleotides introduced by mutation

    cds_offset : int
        Offset into the CDS of first ref nucleotide. For insertions, this
        is the offset of the last ref nucleotide before the insertion.

    sequence_from_start_codon : nucleotide sequence
        Nucleotides of the coding sequence and 3' UTR

    variant : Variant

    transcript : Transcript
    """

    if len(ref) != 0:
        # Logic for any frameshift which isn't an insertion.
        # We have reat insertions as a special case since our base-inclusive
        # indexing means something different for insertions:
        #   cds_offset = base before insertion
        # Whereas in this case:
        #   cds_offset = first reference base affected by a variant
        mutated_codon_index = int(cds_offset / 3)
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
            ref=ref,
            alt=alt)
    else:
        # special logic for insertions
        coding_sequence_after_insertion = \
            sequence_from_start_codon[cds_offset + 1:]
        if cds_offset % 3 == 2:
            # insertion happens after last nucleotide in a codon,
            # doesn't disrupt the existing codon from cds_offset-2 to cds_offset
            mutated_codon_index = int(cds_offset / 3) + 1
            sequence_from_mutated_codon = (
                alt + coding_sequence_after_insertion)
        elif cds_offset % 3 == 1:
            # insertion happens after 2nd nucleotide of a codon
            # codon positions:
            #   1) cds_offset - 1
            #   2) cds_offset
            #    <----- Insertsion
            #   3) cds_offset + 1
            mutated_codon_index = int(cds_offset / 3)
            nucleotides_before = sequence_from_start_codon[
                cds_offset - 1:cds_offset + 1]
            nucleotide_after = sequence_from_start_codon[cds_offset + 1]
            sequence_from_mutated_codon = (
                nucleotides_before +
                alt +
                nucleotide_after +
                coding_sequence_after_insertion)
        elif cds_offset % 3 == 0:
            # insertion happens after 1st nucleotide of a codon
            # codon positions:
            #   1) cds_offset
            #    <----- Insertsion
            #   2) cds_offset + 1
            #   3) cds_offset + 2
            mutated_codon_index = int(cds_offset / 3)
            nucleotide_before = sequence_from_start_codon[cds_offset]
            nucleotides_after = sequence_from_start_codon[
                cds_offset + 1:cds_offset + 3]
            sequence_from_mutated_codon = (
                nucleotide_before +
                alt +
                nucleotides_after +
                coding_sequence_after_insertion)
    return _frameshift(
        mutated_codon_index=mutated_codon_index,
        sequence_from_mutated_codon=sequence_from_mutated_codon,
        variant=variant,
        transcript=transcript)
