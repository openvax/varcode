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

from __future__ import print_function, division, absolute_import

import logging

from .effects import (
    Silent,
    Insertion,
    Deletion,
    Substitution,
    ComplexSubstitution,
    PrematureStop,
    StartLoss,
    FrameShift,
    FrameShiftTruncation,
    IncompleteTranscript,
)
from .string_helpers import trim_shared_flanking_strings

from Bio.Seq import Seq, CodonTable

def mutate(sequence, position, variant_ref, variant_alt):
    """
    Mutate a sequence by substituting given `alt` at instead of `ref` at the
    given `position`.

    Parameters
    ----------
    sequence : sequence
        String of amino acids or DNA bases

    position : int
        Position in the sequence, starts from 0

    variant_ref : sequence or str
        What do we expect to find at the position?

    variant_alt : sequence or str
        Alternate sequence to insert
    """
    n_variant_ref = len(variant_ref)
    sequence_ref = sequence[position:position + n_variant_ref]
    assert str(sequence_ref) == str(variant_ref), \
        "Reference %s at position %d != expected reference %s" % \
        (sequence_ref, position, variant_ref)
    prefix = sequence[:position]
    suffix = sequence[position + n_variant_ref:]
    return prefix + variant_alt + suffix

def infer_coding_effect(
        ref,
        alt,
        cds_offset,
        transcript,
        variant):
    """
    Given a minimal ref/alt nucleotide string pair and an offset into a given
    transcript, determine the coding effect of this nucleotide substitution
    onto the translated protein.

    Parameters
    ----------
    ref : str

    alt : str

    cds_offset : int

    transcript : Transcript

    variant : Variant
    """
    # Don't need a pyfaidx.Sequence object here, just convert it to the an str
    cds_seq = str(transcript.coding_sequence)

    # past this point we know that we're somewhere in the coding sequence
    cds_ref = cds_seq[cds_offset:cds_offset + len(ref)]

    # Make sure that the reference sequence agrees with what we expected
    # from the VCF
    # TODO: check that the ref allele is correct for UTR by looking
    # at transcript.sequence instead of transcript.coding_sequence
    assert cds_ref == ref, \
        "Expected ref '%s', got '%s' in %s (offset %d)" % (
            ref,
            cds_ref,
            variant,
            cds_offset)

    # turn cDNA sequence into a BioPython sequence, translate
    # to amino acids. For the original CDS make sure that it starts with
    # a start codon and ends with a stop codon. Don't include the stop codon
    # in the translated sequence.
    try:
        original_protein = str(Seq(cds_seq).translate(cds=True, to_stop=True))
    except CodonTable.TranslationError as error:
        # coding sequence may contain premature stop codon or may have
        # an incorrectly annotated frame
        logging.warning(
            "Translation error in coding sequence for %s" % transcript)
        logging.warning(error)
        return IncompleteTranscript(variant, transcript)

    assert len(original_protein) > 0, \
        "Translated protein sequence of %s is empty" % (transcript,)
    assert original_protein[0] == "M", \
        "Expected coding sequence of %s to start with M, instead got '%s'" % (
            transcript, original_protein[0])

    variant_cds_seq = mutate(cds_seq, cds_offset, ref, alt)

    # In case sequence isn't a multiple of 3, then truncate it
    # TODO: if we get a frameshift by the end of a CDS (e.g. in the stop codon)
    # then we should use some of the 3' UTR to finish translation.
    truncated_variant_cds_seq = variant_cds_seq[:int(len(variant_cds_seq) / 3) * 3]

    # Can't be sure that the variant is a complete CDS, so passing cds=False
    # in the case of a frameshift, the transcript might not actually contain
    # a stop codon. So, we're for it manually (by passing to_stop=False).
    variant_protein = str(Seq(truncated_variant_cds_seq).translate(
        cds=False, to_stop=False))

    assert len(variant_protein) > 0, \
        "Protein sequence empty for variant %s on transcript %s" % (
            variant, transcript)

    # genomic position to codon position
    aa_pos = int(cds_offset / 3)

    if variant_protein[0] != "M":
        assert aa_pos == 0, \
            "Unexpected change to start codon (M>%s) when aa_pos=%s" % (
                variant_protein[0], aa_pos)
        return StartLoss(
            variant,
            transcript,
            aa_pos=aa_pos,
            aa_alt=variant_protein[0])

    # variant_protein sequence includes stop codon, whereas original_protein
    # doesn't
    if variant_protein[-1] == "*" and original_protein == variant_protein[:-1]:
        return Silent(
            variant,
            transcript,
            aa_pos=aa_pos,
            aa_ref=variant_protein[aa_pos])

    stop_codon_index = variant_protein.find("*")

    # if contained stop codon, truncate sequence before it
    if stop_codon_index > -1:
        variant_protein = variant_protein[:stop_codon_index]

    n_cdna_ref = len(ref)
    n_cdna_alt = len(alt)

    if n_cdna_ref == 0:
        last_aa_ref_pos = aa_pos
        aa_ref = ""
    else:
        last_aa_ref_pos = int((cds_offset + n_cdna_ref - 1) / 3)
        aa_ref = original_protein[aa_pos:last_aa_ref_pos + 1]
        assert len(aa_ref) > 0, \
            "len(aa_ref) = 0 for variant %s on transcript %s (aa_pos=%d:%d)" % (
                variant, transcript, aa_pos, last_aa_ref_pos)

    # is this a premature stop codon?
    if stop_codon_index == aa_pos:
        return PrematureStop(
            variant,
            transcript,
            cds_offset,
            aa_ref)

    if abs(n_cdna_ref - n_cdna_alt) % 3 != 0:
        shifted_sequence = variant_protein[aa_pos:]

        # frameshift may still preserve some of the same codons
        # so trim any shared amino acids
        for i, aa_mutant in enumerate(shifted_sequence):
            if original_protein[aa_pos + i] != aa_mutant:
                break
            aa_pos += 1
        shifted_sequence = shifted_sequence[i:]

        if len(shifted_sequence) == 0:
            return FrameShiftTruncation(
                variant=variant,
                transcript=transcript,
                aa_pos=aa_pos,
                aa_ref=aa_ref)
        else:
            return FrameShift(
                variant=variant,
                transcript=transcript,
                aa_pos=aa_pos,
                aa_ref=aa_ref,
                shifted_sequence=shifted_sequence)

    if n_cdna_alt == 0:
        last_aa_alt_pos = aa_pos
        aa_alt = ""
    else:
        last_aa_alt_pos = int((cds_offset + n_cdna_alt - 1) / 3)
        aa_alt = variant_protein[aa_pos:last_aa_alt_pos + 1]
        assert len(aa_alt) > 0, \
            "len(aa_alt) = 0 for variant %s on transcript %s (aa_pos=%d:%d)" % (
                variant, transcript, aa_pos, last_aa_ref_pos)

    assert aa_alt != aa_ref, \
        "Unexpected silent mutation for variant %s on transcript %s (aa=%s)" % (
            variant, transcript, aa_ref)

    # in case of simple insertion like FY>FYGL or deletions FYGL > FY,
    # get rid of the shared prefixes/suffixes
    aa_ref, aa_alt, prefix, _ = trim_shared_flanking_strings(aa_ref, aa_alt)
    aa_pos += len(prefix)

    # Deletion e.g. p.389delQQ
    if len(aa_alt) == 0:
        assert len(aa_ref) > 0, \
            "Can't have empty ref and alt for variant %s on transcript %s" % (
                variant, transcript)
        n_kept = len(aa_alt)
        # if aa_ref = FYPQ and aa_alt = F, then deleted = YPQ
        deleted = aa_ref[n_kept:]

        return Deletion(
            variant,
            transcript,
            aa_pos=aa_pos + n_kept,
            aa_ref=deleted)

    # Insertion, e.g. p.37insA
    elif aa_alt.startswith(aa_ref):
        assert len(aa_alt) > 0, \
            "Can't have empty ref and alt for variant %s on transcript %s" % (
                variant, transcript)
        inserted = aa_alt[len(aa_ref):]

        return Insertion(
            variant, transcript,
            aa_pos=aa_pos,
            aa_alt=inserted)

    # simple substitution e.g. p.V600E
    elif len(aa_ref) == 1 and len(aa_alt) == 1:
        return Substitution(
            variant,
            transcript,
            aa_pos=aa_pos,
            aa_ref=aa_ref,
            aa_alt=aa_alt)

    # substitution which involes multiple amino acids
    # Example: p.V600EEQ, p.IL49AQY
    else:
        return ComplexSubstitution(
            variant,
            transcript,
            aa_pos=aa_pos,
            aa_ref=aa_ref,
            aa_alt=aa_alt)
