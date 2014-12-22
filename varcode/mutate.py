# Copyright (c) 2014. Mount Sinai School of Medicine
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

import logging
import math
from collections import namedtuple

from Bio.Seq import Seq


def mutate_split(sequence, position, ref, alt):
    """
    Return the prefix, mutated region, and suffix of
    applying a mutation to a sequence.

    Parameters
    ----------
    sequence : sequence
        String of amino acids or DNA bases

    position : int
        Position in the sequence, starts from 0

    ref : sequence or str
        What do we expect to find at the position?

    alt : sequence or str
        Alternate allele to insert
    """

    transcript_ref_base = sequence[position:position+len(ref)]
    if ref != '.':
        assert str(transcript_ref_base) == ref, \
            "Transcript ref base %s at position %d != given reference %s" % \
            (transcript_ref_base, position, ref)
    prefix = sequence[:position]
    suffix = sequence[position+len(ref):]
    return prefix, alt, suffix

def mutate(sequence, position, ref, alt):
    """
    Mutate a sequence by substituting given `alt` at instead of `ref` at the
    given `position`.

    Parameters
    ----------
    sequence : sequence
        String of amino acids or DNA bases

    position : int
        Position in the sequence, starts from 0

    ref : sequence or str
        What do we expect to find at the position?

    alt : sequence or str
        Alternate allele to insert
    """
    prefix, alt, suffix = mutate_split(sequence, position, ref, alt)
    return prefix + alt + suffix


def gene_mutation_description(pos, ref, alt):
    if ref == alt:
        return "g.%d %s=%s" % (pos, ref, alt)
    elif len(ref) == 0 or alt.startswith(ref):
        return "g.%d ins%s" % (pos + len(ref),  alt[len(ref):])
    elif len(alt) == 0 or ref.startswith(alt):
        return "g.%d_%d del%s" % (pos + len(alt), pos + len(ref), ref[len(alt):])
    else:
        return "g.%d %s>%s" % (pos, ref, alt)

def protein_mutation_description(
        original_protein,
        mutated_protein,
        aa_position,
        n_deleted,
        n_inserted,
        frameshift,
        early_stop):

    ref_start = aa_position

    assert len(original_protein) > ref_start, \
        "Protein sequence too short (%d residues) for index %d" % \
        (len(original_protein), ref_start)

    ref_stop = max(aa_position + n_deleted, 1)
    aa_ref = original_protein[ref_start : ref_stop]
    mut_start = aa_position
    mut_stop = max(aa_position + n_inserted, 1)
    aa_mut = mutated_protein[mut_start : mut_stop]

    if frameshift:
        return "%s%dfs" % (aa_ref[0], aa_position+1)
    elif early_stop:
        return "%s%d%s*" % \
            (aa_ref, aa_position+1, aa_mut)
    elif n_inserted == 0:
        return "%s%ddel" % (aa_ref, aa_position+1)
    elif n_deleted == 0:
        return "%dins%s" % (aa_position + 1, aa_mut)
    else:
        if aa_ref == aa_mut:
            logging.warn("Not a mutation: ref and alt = %s at position %d",  aa_ref, aa_position)
        return "%s%d%s" % \
            (aa_ref, aa_position+1, aa_mut)

Mutation = \
    namedtuple(
        "Mutation",
        (
            "seq",    # mutated region sequence
            "start",  # start position of region in the protein
            "stop",   # stop position of region in the protein
            "mutation_start", # where in the region is the first mutated AA?
            "n_removed",  # how many wildtype residues removed?
            "n_inserted",  # how many new residues in the seq?
            "annot", # mutation annotation i.e. "V600E"
        ))


def mutate_protein_from_transcript(
        transcript_seq,
        position,
        dna_ref,
        dna_alt,
        padding = None):
    """
    Mutate a sequence by inserting the allele into the genomic transcript
    and translate to protein sequence

    Parameters
    ----------
    transcript_seq :  sequence
        Protein sequence we're going to mutate

    position : int
        Position in `transcript_seq`, starts from 0

    dna_ref : sequence or str
        What's already supposed to be at the given position

    dna_alt : sequence or str
        Alternate substring to insert

    padding : int, optional
        Number of wildtype amino acids to keep left and right of the
        mutation. Default is to return whole mutated string.
    """

    # turn any character sequence into a BioPython sequence
    transcript_seq = Seq(str(transcript_seq))

    transcript_ref_base = transcript_seq[position:position+len(dna_ref)]

    if dna_ref != '.':
        assert str(transcript_ref_base) == dna_ref, \
            "Transcript reference base %s at position %d != reference %s" % \
            (transcript_ref_base, position, dna_ref)

    original_protein = transcript_seq.translate()
    n_original_protein = len(original_protein)

    mutated_dna = mutate(transcript_seq, position, dna_ref, dna_alt)
    mutated_protein = mutated_dna.translate()
    n_mutated_protein = len(mutated_protein)

    if str(original_protein) == str(mutated_protein):
        # if protein product is unmodified then
        # this is a silent mutation
        return Mutation(
            seq = "",
            start = 0,
            stop = 0,
            mutation_start = 0,
            n_removed = 0,
            n_inserted = 0,
            annot = "SILENT"
        )

    aa_position = int(position / 3)  # genomic position to codon position

    n_dna_ref = len(dna_ref)
    n_dna_alt = len(dna_alt)

    # is this a frameshift mutation?
    if abs(n_dna_ref - n_dna_alt) % 3 != 0:
        # frameshifts 'delete' the rest of the origin protein string
        # and 'insert' the frameshifted residues.
        # These numbers, since they encompass the entire protein,
        # will typically be larger than n_wildtype_deleted/n_mutant_inserted
        # on the region struct below, since that restricts the residues counts
        # to those within a particular region
        n_aa_deleted = n_original_protein - aa_position
        # careful, this doesn't include stop codons, will
        # have to update it later
        n_aa_inserted = n_mutated_protein - aa_position
        frameshift = True
    else:
        n_aa_deleted = int(math.ceil(n_dna_ref / 3.0))

        n_aa_inserted = int(math.ceil(n_dna_alt / 3.0))
        frameshift = False

    if padding is None:
        start_pos = 0
        end_pos = n_mutated_protein
    else:
        start_pos = max(0, aa_position - padding)
        end_pos = min(n_mutated_protein, aa_position + padding + 1)


    # move padding to not include stop codons
    prefix_stop_codon = str(mutated_protein[:aa_position]).rfind("*")
    #  stop codon found to the left
    if prefix_stop_codon != -1:
        start_pos = max(start_pos, prefix_stop_codon + 1)

    suffix_stop_codon = str(mutated_protein[aa_position:]).find("*")
    # stop codon found to the right
    early_stop = False
    if suffix_stop_codon != -1:
        if suffix_stop_codon < n_aa_inserted:
            n_aa_inserted = suffix_stop_codon
            early_stop = True
        end_pos = min(end_pos, suffix_stop_codon + aa_position)


    seq_region = mutated_protein[start_pos : end_pos]
    mutation_start_pos_in_region = aa_position - start_pos

    annot = \
        protein_mutation_description(
            original_protein,
            mutated_protein,
            aa_position,
            n_aa_deleted,
            n_aa_inserted,
            frameshift,
            early_stop)

    return Mutation(
        seq = str(seq_region),
        start = start_pos,
        stop = end_pos,
        mutation_start = mutation_start_pos_in_region,
        n_removed = n_aa_deleted,
        n_inserted = n_aa_inserted,
        annot = annot)
