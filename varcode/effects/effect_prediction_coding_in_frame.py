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
from .translate import translate_in_frame_mutation, START_CODONS

def get_codons(
        variant,
        trimmed_cdna_ref,
        trimmed_cdna_alt,
        sequence_from_start_codon,
        cds_offset):
    """
    Returns indices of first and last reference codons affected by the variant,
    as well as the actual sequence of the mutated codons which replace those
    reference codons.

    Parameters
    ----------
    variant : Variant

    trimmed_cdna_ref : str
        Trimmed reference cDNA nucleotides affected by the variant

    trimmed_cdna_alt : str
        Trimmed alternate cDNA nucleotides which replace the reference

    sequence_from_start_codon : str
        cDNA nucleotide coding sequence

    cds_offset : int
        Integer offset into the coding sequence where ref is replace with alt
    """
    # index (starting from 0) of first affected reference codon
    ref_codon_start_offset = cds_offset // 3
    # which nucleotide of the first codon got changed?
    nucleotide_offset_into_first_ref_codon = cds_offset % 3

    n_ref_nucleotides = len(trimmed_cdna_ref)
    if n_ref_nucleotides == 0:
        if nucleotide_offset_into_first_ref_codon == 2:
            # if we're inserting between codons
            ref_codon_end_offset = ref_codon_start_offset
        else:
            # inserting inside a reference codon
            ref_codon_end_offset = ref_codon_start_offset + 1
        ref_codons = sequence_from_start_codon[
            ref_codon_start_offset * 3:ref_codon_end_offset * 3]
        # split the reference codon into nucleotides before/after insertion
        prefix = ref_codons[:nucleotide_offset_into_first_ref_codon + 1]
        suffix = ref_codons[nucleotide_offset_into_first_ref_codon + 1:]
    else:
        ref_codon_end_offset = (cds_offset + n_ref_nucleotides - 1) // 3 + 1
        # codons in the reference sequence
        ref_codons = sequence_from_start_codon[
            ref_codon_start_offset * 3:ref_codon_end_offset * 3]

        # We construct the new codons by taking the unmodified prefix
        # of the first ref codon, the unmodified suffix of the last ref codon
        # and sticking the alt nucleotides in between.
        # Since this is supposed to be an in-frame mutation, the concatenated
        # nucleotide string is expected to have a length that is a multiple of
        # three.
        prefix = ref_codons[:nucleotide_offset_into_first_ref_codon]

        offset_in_last_ref_codon = (cds_offset + n_ref_nucleotides - 1) % 3

        if offset_in_last_ref_codon == 0:
            suffix = ref_codons[-2:]
        elif offset_in_last_ref_codon == 1:
            suffix = ref_codons[-1:]
        else:
            suffix = ""
    mutant_codons = prefix + trimmed_cdna_alt + suffix
    assert len(mutant_codons) % 3 == 0, \
        "Expected in-frame mutation but got %s (length = %d)" % (
            mutant_codons, len(mutant_codons))
    return ref_codon_start_offset, ref_codon_end_offset, mutant_codons

def predict_in_frame_coding_effect(
        variant,
        transcript,
        trimmed_cdna_ref,
        trimmed_cdna_alt,
        sequence_from_start_codon,
        cds_offset):
    """Coding effect of an in-frame nucleotide change

    Parameters
    ----------
    variant : Variant

    transcript : Transcript

    trimmed_cdna_ref : str
        Reference nucleotides from the coding sequence of the transcript

    trimmed_cdna_alt : str
        Nucleotides to insert in place of the reference nucleotides

    sequence_from_start_codon : Bio.Seq or str
        Transcript sequence from the CDS start codon (including the 3' UTR).
        This sequence includes the 3' UTR since a mutation may delete the stop
        codon and we'll have to translate past the normal end of the CDS to
        determine the new protein sequence.

    cds_offset : int
        Index of first ref nucleotide, starting from 0 = beginning of coding
        sequence. If variant is a pure insertion (no ref nucleotides) then this
        argument indicates the offset *after* which to insert the alt
        nucleotides.
    """
    ref_codon_start_offset, ref_codon_end_offset, mutant_codons = get_codons(
        variant=variant,
        trimmed_cdna_ref=trimmed_cdna_ref,
        trimmed_cdna_alt=trimmed_cdna_alt,
        sequence_from_start_codon=sequence_from_start_codon,
        cds_offset=cds_offset)

    mutation_affects_start_codon = (ref_codon_start_offset == 0)

    if mutation_affects_start_codon and mutant_codons[:3] not in START_CODONS:
        # if we changed a start codon to something else then
        # we no longer know where the protein begins (or even in
        # what frame).
        # TODO: use the Kozak consensus sequence or a predictive model
        # to identify the most likely start site
        return StartLoss(
            variant=variant,
            transcript=transcript)

    # rely on Ensembl's annotation of the protein sequence since we can't
    # easily predict whether the starting nucleotide is a methionine
    # (most common) or leucine
    aa_ref = transcript.protein_sequence[ref_codon_start_offset:ref_codon_end_offset]

    reference_protein_length = len(transcript.protein_sequence)

    aa_alt, mutant_stop_codon_index, using_three_prime_utr = \
        translate_in_frame_mutation(
            transcript=transcript,
            ref_codon_start_offset=ref_codon_start_offset,
            ref_codon_end_offset=ref_codon_end_offset,
            mutant_codons=mutant_codons)

    mutant_codons_contain_stop = mutant_stop_codon_index != -1

    # trim shared subsequences at the start and end of reference
    # and mutated amino acid sequences
    aa_ref, aa_alt, shared_prefix, shared_suffix = \
        trim_shared_flanking_strings(
            aa_ref,
            aa_alt)

    n_aa_ref = len(aa_ref)
    n_aa_alt = len(aa_alt)
    n_aa_shared = len(shared_prefix)

    is_insertion = (ref_codon_start_offset == ref_codon_end_offset)

    # index of first amino acid which is different from the reference
    aa_mutation_start_offset = (
        ref_codon_start_offset + n_aa_shared + is_insertion)

    if mutant_codons_contain_stop:
        mutant_stop_codon_index += n_aa_shared

    if mutation_affects_start_codon and (aa_ref == aa_alt):
            # Substitution between start codons gets special treatment since,
            # though superficially synonymous, this could still potentially
            # cause a start loss / change in reading frame and might be worth
            # closer scrutiny
            return AlternateStartCodon(
                variant=variant,
                transcript=transcript,
                ref_codon=transcript.sequence[:3],
                alt_codon=mutant_codons[:3])

    n_ref_amino_acids_after_mutated_site = (
        reference_protein_length - aa_mutation_start_offset - 1)

    if mutant_codons_contain_stop and (
            n_aa_alt <= n_ref_amino_acids_after_mutated_site):
        # if the new coding sequence contains a stop codon, then this is a
        # PrematureStop mutation if it decreases the length of the protein
        return PrematureStop(
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt=aa_alt)

    if (aa_mutation_start_offset > reference_protein_length) or (
            n_aa_ref == n_aa_alt == 0):
        # if inserted nucleotides go after original stop codon or if nothing
        # is changed in the amino acid sequence then this is a Silent variant
        return Silent(
            variant=variant,
            transcript=transcript,
            aa_pos=aa_mutation_start_offset,
            aa_ref=shared_prefix + shared_suffix)
    elif using_three_prime_utr:
        # if non-silent mutation is at the end of the protein then
        # should be a stop-loss
        return StopLoss(
            variant,
            transcript,
            extended_protein_sequence=aa_alt)
    elif n_aa_alt == 0:
        return Deletion(
            variant,
            transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref)
    elif n_aa_ref == 0:
        return Insertion(
            variant,
            transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_alt=aa_alt)
    elif n_aa_ref == n_aa_alt == 1:
        # simple substitution e.g. p.V600E
        return Substitution(
            variant,
            transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt=aa_alt)
    else:
        # multiple amino acids were substituted e.g. p.VQQ39FF
        return ComplexSubstitution(
            variant,
            transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt=aa_alt)
