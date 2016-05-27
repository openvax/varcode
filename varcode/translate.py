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
