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
Effect annotation for variants which modify the coding sequence without
changing the reading frame.
"""
from __future__ import division, absolute_import, print_function

from six.moves import range

from ..string_helpers import trim_shared_flanking_strings
from .effect_classes import (
    Silent,
    Insertion,
    Deletion,
    Substitution,
    ComplexSubstitution,
    PrematureStop,
    AlternateStartCodon,
    StartLoss,
    StopLoss,
)
from .translate import START_CODONS, STOP_CODONS, translate

def translate_mutant_codons(
        mutant_codons,
        first_ref_codon_index,
        last_ref_codon_index,
        reference_protein_length,
        transcript):

    mutant_protein_subsequence = translate(
        mutant_codons,
        first_codon_is_start=(first_ref_codon_index == 0))

    mutant_codons_contain_stop = contains_stop_codon(mutant_codons)

    if not mutant_codons_contain_stop and (
            last_ref_codon_index >= reference_protein_length):
        # if the mutant codons didn't contain a stop but did mutate the last
        # reference codon then the translated sequence might involve the 3' UTR
        three_prime_utr = transcript.three_prime_utr_sequence
        n_utr_codons = len(three_prime_utr) // 3
        # trim the 3' UTR sequence to have a length that is a multiple of 3
        truncated_utr_sequence = three_prime_utr[:n_utr_codons * 3]
        translated_utr = translate(
            truncated_utr_sequence, first_codon_is_start=False)
        # combine the in-frame mutant codons with the truncated sequence of
        # the 3' UTR
        mutant_codons += truncated_utr_sequence
        mutant_protein_subsequence += translated_utr
    return mutant_protein_subsequence

def contains_stop_codon(mutant_codons):
    """
    Given a sequence of codons (expected to have length multiple of three),
    are any of them a stop codon?
    """
    n_mutant_codons = len(mutant_codons) // 3
    return STOP_CODONS.intersection(
        mutant_codons[3 * i:3 * i + 3]
        for i in range(n_mutant_codons))

def choose_in_frame_effect_annotation(
        variant,
        transcript,
        first_ref_codon_index,
        last_ref_codon_index,
        mutant_codons):
    """Choose a coding effect annotation for in-frame mutations which do
    not affect the start codon and do not introduce a premature stop codon.
    This function encompasses all the logic which does not need to look at the
    specific nucleotides which created each amino acid (can deal only with
    amino acid sequences).

    Parameters
    ----------
    variant : Variant

    transcript : Transcript

    first_ref_codon_index : int
        Inclusive (starting from 0) amino acid position of the first ref
        amino acid which is changed by the mutation. Might include
        synonymous substitutions.

    last_ref_codon_index : int

    mutant_codons : str
        cDNA nucleotide sequence of mutated codons
    """

    modifies_start_codon = (first_ref_codon_index == 0)

    if modifies_start_codon and (mutant_codons[:3] not in START_CODONS):
        # if we changed a start codon to something else then
        # we no longer know where the protein begins (or even in
        # what frame).
        # TODO: use the Kozak consensus sequence or a predictive model
        # to identify the most likely start site
        return StartLoss(
            variant=variant,
            transcript=transcript)

    original_protein_sequence = transcript.protein_sequence

    aa_ref = original_protein_sequence[
        first_ref_codon_index:last_ref_codon_index + 1]

    reference_protein_length = len(original_protein_sequence)

    aa_alt = translate_mutant_codons(
        mutant_codons=mutant_codons,
        first_ref_codon_index=first_ref_codon_index,
        last_ref_codon_index=last_ref_codon_index,
        reference_protein_length=reference_protein_length,
        transcript=transcript)

    if modifies_start_codon and (aa_ref == aa_alt):
        # Substitution between start codons gets special treatment since,
        # though superficially synonymous, this could still potentially
        # cause a start loss / change in reading frame and might be worth
        # closer scrutiny
        return AlternateStartCodon(
            variant=variant,
            transcript=transcript,
            aa_ref=aa_ref,
            ref_codon=transcript.sequence[:3],
            alt_codon=mutant_codons)

    aa_ref, aa_alt, shared_prefix, shared_suffix = \
        trim_shared_flanking_strings(
            aa_ref,
            aa_alt)

    # index of first amino acid which is different from the reference
    aa_mutation_start_offset = first_ref_codon_index + len(shared_prefix)
    last_ref_amino_acid_index = aa_mutation_start_offset + len(aa_ref)

    mutant_codons_contain_stop = contains_stop_codon(mutant_codons)

    if mutant_codons_contain_stop:
        # if the new coding sequence contains a stop codon, then this is a
        # PrematureStop mutation
        n_remaining_amino_acids_in_ref = (
            reference_protein_length - aa_mutation_start_offset)
        if len(aa_alt) < n_remaining_amino_acids_in_ref:
            # only call this mutation a premature stop if it decreases
            # the length of the protein
            return PrematureStop(
                variant=variant,
                transcript=transcript,
                aa_mutation_start_offset=aa_mutation_start_offset,
                aa_ref=aa_ref,
                aa_alt=aa_alt)

    if len(aa_ref) == len(aa_alt) == 0:
        return Silent(
            variant=variant,
            transcript=transcript,
            aa_pos=aa_mutation_start_offset,
            aa_ref=shared_prefix + shared_suffix)
    elif last_ref_amino_acid_index == reference_protein_length and (
            not mutant_codons_contain_stop):
        # if non-silent mutation is at the end of the protein then
        # should be a stop-loss
        return StopLoss(
            variant,
            transcript,
            extended_protein_sequence=aa_alt)
    elif len(aa_alt) == 0:
        return Deletion(
            variant,
            transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref)
    elif len(aa_ref) == 0:
        return Insertion(
            variant,
            transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_alt=aa_alt)
    elif len(aa_alt) == len(aa_ref) == 1:
        # simple substitution e.g. p.V600E
        return Substitution(
            variant,
            transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt=aa_alt)
    else:
        return ComplexSubstitution(
            variant,
            transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt=aa_alt)


def get_codons(
        ref,
        alt,
        cds_offset,
        sequence_from_start_codon,
        transcript_id,
        protein_sequence,
        variant):
    """
    Returns reference and mutated codons
    """
    # index (starting from 0) of first affected reference codon
    first_ref_codon_index = cds_offset // 3

    # which nucleotide of the first codon got changed?
    offset_in_first_ref_codon = cds_offset % 3

    if len(ref) == 0:
        # inserting inside a reference codon
        # include an extra codon at the end of the reference so that if we
        # insert a stop before a stop, we can return Silent
        insertion_ref_codon = sequence_from_start_codon[
            first_ref_codon_index * 3:first_ref_codon_index * 3 + 6]
        if offset_in_first_ref_codon == 2:
            # if insertion happens between codons then we don't actually
            # need to select any reference codons for modification
            last_ref_codon_index = first_ref_codon_index
        else:
            last_ref_codon_index = first_ref_codon_index + 1
        # split the reference codon into nucleotides before/after insertion
        prefix = insertion_ref_codon[:offset_in_first_ref_codon + 1]
        suffix = insertion_ref_codon[offset_in_first_ref_codon + 1:]
        mutant_codons = prefix + alt + suffix
    else:
        assert first_ref_codon_index <= len(protein_sequence), \
            ("Unexpected mutation at offset %d (5' UTR starts at %d)"
             " while annotating %s on %s") % (
                 first_ref_codon_index,
                 len(protein_sequence),
                 variant,
                 transcript_id)
        n_ref_nucleotides = len(ref)
        last_ref_codon_index = int((cds_offset + n_ref_nucleotides - 1) / 3)

        assert last_ref_codon_index >= first_ref_codon_index, \
            ("Expected first_ref_codon_index (%d) <= "
             "last_ref_codon_index (%d) while annotating %s on %s") % (
                first_ref_codon_index,
                last_ref_codon_index,
                variant,
                transcript_id)
        # codons in the reference sequence
        ref_codons = sequence_from_start_codon[
            first_ref_codon_index * 3:last_ref_codon_index * 3 + 3]

        # We construct the new codons by taking the unmodified prefix
        # of the first ref codon, the unmodified suffix of the last ref codon
        # and sticking the alt nucleotides in between.
        # Since this is supposed to be an in-frame mutation, the concatenated
        # nucleotide string is expected to have a length that is a multiple of
        # three.
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
    return first_ref_codon_index, last_ref_codon_index, mutant_codons

def predict_in_frame_coding_effect(
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
    first_ref_codon_index, last_ref_codon_index, mutant_codons = get_codons(
        ref=ref,
        alt=alt,
        cds_offset=cds_offset,
        sequence_from_start_codon=sequence_from_start_codon,
        transcript_id=transcript.id,
        protein_sequence=transcript.protein_sequence,
        variant=variant)

    return choose_in_frame_effect_annotation(
        variant=variant,
        transcript=transcript,
        first_ref_codon_index=first_ref_codon_index,
        last_ref_codon_index=last_ref_codon_index,
        mutant_codons=mutant_codons)
