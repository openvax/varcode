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

"""Helpers for cDNA -> protein translation.

TODO: generalize this to work with the mitochondrial codon table.
"""

from __future__ import division, absolute_import, print_function

from six.moves import range
from Bio.Data import CodonTable
from Bio.Seq import Seq

DNA_CODON_TABLE = CodonTable.standard_dna_table.forward_table
START_CODONS = set(CodonTable.standard_dna_table.start_codons)
STOP_CODONS = set(CodonTable.standard_dna_table.stop_codons)

def translate_codon(codon, aa_pos):
    """Translate a single codon into a single amino acid or stop '*'

    Parameters
    ----------
    codon : str
        Expected to be of length 3
    aa_pos : int
        Codon/amino acid offset into the protein (starting from 0)
    """
    # not handling rare Leucine or Valine starts!
    if aa_pos == 0 and codon in START_CODONS:
        return "M"
    elif codon in STOP_CODONS:
        return "*"
    else:
        return DNA_CODON_TABLE[codon]

def translate(
        nucleotide_sequence,
        first_codon_is_start=True,
        to_stop=True,
        truncate=False):
    """Translates cDNA coding sequence into amino acid protein sequence.

    Should typically start with a start codon but allowing non-methionine
    first residues since the CDS we're translating might have been affected
    by a start loss mutation.

    The sequence may include the 3' UTR but will stop translation at the first
    encountered stop codon.

    Parameters
    ----------
    nucleotide_sequence : BioPython Seq
        cDNA sequence

    first_codon_is_start : bool
        Treat the beginning of nucleotide_sequence (translates methionin)

    truncate : bool
        Truncate sequence if it's not a multiple of 3 (default = False)
    Returns BioPython Seq of amino acids
    """
    if not isinstance(nucleotide_sequence, Seq):
        nucleotide_sequence = Seq(nucleotide_sequence)

    if truncate:
        # if sequence isn't a multiple of 3, truncate it so BioPython
        # doesn't complain
        n_nucleotides = int(len(nucleotide_sequence) / 3) * 3
        nucleotide_sequence = nucleotide_sequence[:n_nucleotides]
    else:
        n_nucleotides = len(nucleotide_sequence)

    assert n_nucleotides % 3 == 0, \
        ("Expected nucleotide sequence to be multiple of 3"
         " but got %s of length %d") % (
            nucleotide_sequence,
            n_nucleotides)

    # passing cds=False to translate since we may want to deal with premature
    # stop codons
    protein_sequence = nucleotide_sequence.translate(to_stop=to_stop, cds=False)

    if first_codon_is_start and (
            len(protein_sequence) == 0 or protein_sequence[0] != "M"):
        if nucleotide_sequence[:3] in START_CODONS:
            # TODO: figure out when these should be made into methionines
            # and when left as whatever amino acid they normally code for
            # e.g. Leucine start codons
            # See: DOI: 10.1371/journal.pbio.0020397
            return "M" + protein_sequence[1:]
        else:
            raise ValueError(
                ("Expected first codon of %s to be start codon"
                 " (one of %s) but got %s") % (
                    protein_sequence[:10],
                    START_CODONS,
                    nucleotide_sequence))

    return protein_sequence


def find_first_stop_codon(nucleotide_sequence):
    """
    Given a sequence of codons (expected to have length multiple of three),
    return index of first stop codon, or -1 if none is in the sequence.
    """
    n_mutant_codons = len(nucleotide_sequence) // 3
    for i in range(n_mutant_codons):
        codon = nucleotide_sequence[3 * i:3 * i + 3]
        if codon in STOP_CODONS:
            return i
    return -1

def translate_in_frame_mutation(
        transcript,
        ref_codon_start_offset,
        ref_codon_end_offset,
        mutant_codons):
    """
    Returns:
        - mutant amino acid sequence
        - offset of first stop codon in the mutant sequence (or -1 if there was none)
        - boolean flag indicating whether any codons from the 3' UTR were used

    Parameters
    ----------
    transcript : pyensembl.Transcript
        Reference transcript to which a cDNA mutation should be applied.

    ref_codon_start_offset : int
        Starting (base 0) integer offset into codons (character triplets) of the
        transcript's reference coding sequence.

    ref_codon_end_offset : int
        Final (base 0) integer offset into codons of the transcript's
        reference coding sequence.

    mutant_codons : str
        Nucleotide sequence to replace the reference codons with
        (expected to have length that is a multiple of three)
    """
    mutant_stop_codon_index = find_first_stop_codon(mutant_codons)

    using_three_prime_utr = False

    if mutant_stop_codon_index != -1:
        mutant_codons = mutant_codons[:3 * mutant_stop_codon_index]
    elif ref_codon_end_offset > len(transcript.protein_sequence):
        # if the mutant codons didn't contain a stop but did mutate the
        # true reference stop codon then the translated sequence might involve
        # the 3' UTR
        three_prime_utr = transcript.three_prime_utr_sequence
        n_utr_codons = len(three_prime_utr) // 3
        # trim the 3' UTR sequence to have a length that is a multiple of 3
        truncated_utr_sequence = three_prime_utr[:n_utr_codons * 3]

        # note the offset of the first stop codon in the combined
        # nucleotide sequence of both the end of the CDS and the 3' UTR
        first_utr_stop_codon_index = find_first_stop_codon(truncated_utr_sequence)

        if first_utr_stop_codon_index > 0:
            # if there is a stop codon in the 3' UTR sequence and it's not the
            # very first codon
            using_three_prime_utr = True
            n_mutant_codons_before_utr = len(mutant_codons) // 3
            mutant_stop_codon_index = n_mutant_codons_before_utr + first_utr_stop_codon_index
            # combine the in-frame mutant codons with the truncated sequence of
            # the 3' UTR
            mutant_codons += truncated_utr_sequence[:first_utr_stop_codon_index * 3]
        elif first_utr_stop_codon_index == -1:
            # if there is no stop codon in the 3' UTR sequence
            using_three_prime_utr = True
            mutant_codons += truncated_utr_sequence

    amino_acids = translate(
        mutant_codons,
        first_codon_is_start=(ref_codon_start_offset == 0))

    return amino_acids, mutant_stop_codon_index, using_three_prime_utr
