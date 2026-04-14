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

Codon tables live in :mod:`varcode.effects.codon_tables` as plain
Python dicts/frozensets — no BioPython dependency for codon lookup.
Mitochondrial transcripts are routed to NCBI translation table 2 via
:func:`codon_table_for_transcript`; everything else uses table 1.
"""

from .codon_tables import (
    STANDARD,
    codon_table_for_transcript,
    translate_sequence,
)

# Back-compat module-level constants — callers that imported these
# names directly still get the standard (nuclear) table. Code paths
# that need mitochondrial handling look up the per-transcript table
# through codon_table_for_transcript().
DNA_CODON_TABLE = STANDARD.forward_table
START_CODONS = STANDARD.start_codons
STOP_CODONS = STANDARD.stop_codons


def translate_codon(codon, aa_pos, codon_table=STANDARD):
    """Translate a single codon into a single amino acid or stop '*'.

    Parameters
    ----------
    codon : str
        Expected to be of length 3.
    aa_pos : int
        Codon/amino acid offset into the protein (starting from 0).
    codon_table : CodonTable
        Defaults to the standard nuclear table. Pass
        :data:`VERTEBRATE_MITOCHONDRIAL` for mt transcripts.
    """
    # not handling rare Leucine or Valine starts!
    if aa_pos == 0 and codon in codon_table.start_codons:
        return "M"
    elif codon in codon_table.stop_codons:
        return "*"
    else:
        return codon_table.forward_table[codon]


def translate(
        nucleotide_sequence,
        first_codon_is_start=True,
        to_stop=True,
        truncate=False,
        codon_table=STANDARD):
    """Translate a cDNA coding sequence into an amino acid string.

    Parameters
    ----------
    nucleotide_sequence : str or convertible-to-str
        cDNA sequence.
    first_codon_is_start : bool
        Treat the first codon as a start codon (translates to 'M' even
        when the codon is a non-ATG alternate start like TTG/GTG).
    to_stop : bool
        Stop translation at the first stop codon.
    truncate : bool
        Truncate the sequence to a multiple of 3 before translating.
    codon_table : CodonTable
        Defaults to standard. For mt transcripts pass the vertebrate
        mitochondrial table.

    Returns
    -------
    str
        Amino acid sequence (previously returned a ``Bio.Seq.Seq``;
        now a plain string — all in-repo callers already use string
        semantics).
    """
    seq = str(nucleotide_sequence)
    if truncate:
        n_nucleotides = (len(seq) // 3) * 3
        seq = seq[:n_nucleotides]
    else:
        n_nucleotides = len(seq)

    assert n_nucleotides % 3 == 0, \
        ("Expected nucleotide sequence to be multiple of 3"
         " but got %s of length %d") % (seq[:30], n_nucleotides)

    protein_sequence = translate_sequence(
        seq, codon_table=codon_table, to_stop=to_stop)

    if first_codon_is_start and (
            len(protein_sequence) == 0 or protein_sequence[0] != "M"):
        if seq[:3].upper() in codon_table.start_codons:
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
                    sorted(codon_table.start_codons),
                    seq[:30]))

    return protein_sequence


def find_first_stop_codon(nucleotide_sequence, codon_table=STANDARD):
    """
    Given a sequence of codons (expected to have length multiple of three),
    return the codon offset of the first stop codon, or -1 if none is in
    the sequence. The meaning of "stop" depends on ``codon_table`` — in
    the vertebrate mitochondrial table, AGA and AGG are stops.
    """
    seq = str(nucleotide_sequence).upper()
    n_mutant_codons = len(seq) // 3
    stops = codon_table.stop_codons
    for i in range(n_mutant_codons):
        codon = seq[3 * i:3 * i + 3]
        if codon in stops:
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

    The codon table is selected automatically from the transcript's
    contig — mitochondrial transcripts use NCBI table 2 (AGA/AGG as
    stops, TGA as Trp, ATA as Met); everything else uses the standard
    table.

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
    codon_table = codon_table_for_transcript(transcript)

    mutant_stop_codon_index = find_first_stop_codon(
        mutant_codons, codon_table=codon_table)

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
        first_utr_stop_codon_index = find_first_stop_codon(
            truncated_utr_sequence, codon_table=codon_table)

        if first_utr_stop_codon_index >= 0:
            # if there is a stop codon in the 3' UTR sequence (including the
            # case where it is the very first codon immediately after the
            # disrupted original stop)
            using_three_prime_utr = True
            n_mutant_codons_before_utr = len(mutant_codons) // 3
            mutant_stop_codon_index = n_mutant_codons_before_utr + first_utr_stop_codon_index
            # combine the in-frame mutant codons with the truncated sequence of
            # the 3' UTR up to the first UTR stop codon (empty when the UTR
            # itself starts with a stop)
            mutant_codons += truncated_utr_sequence[:first_utr_stop_codon_index * 3]
        else:
            # if there is no stop codon in the 3' UTR sequence
            using_three_prime_utr = True
            mutant_codons += truncated_utr_sequence

    amino_acids = translate(
        mutant_codons,
        first_codon_is_start=(ref_codon_start_offset == 0),
        codon_table=codon_table)

    return amino_acids, mutant_stop_codon_index, using_three_prime_utr
