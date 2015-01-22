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

from common import reverse_complement
from variant import Variant
from transcript_mutation_effects import (
    NoncodingTranscript,
    IncompleteTranscript,
    FivePrimeUTR,
    ThreePrimeUTR,
    Intronic,
    Silent,
    Insertion,
    Deletion,
    Substitution,
    PrematureStop,
    StartLoss,
    FrameShift
)

from Bio.Seq import Seq
from pyensembl.transcript import Transcript
from pyensembl.biotypes import is_coding_biotype

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
    sequence_ref = sequence[position:position+n_variant_ref]
    assert str(sequence_ref) == str(variant_ref), \
        "Reference %s at position %d != expected reference %s" % \
        (sequence_ref, position, variant_ref)
    prefix = sequence[:position]
    suffix = sequence[position+n_variant_ref:]
    return prefix + variant_alt + suffix

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

def overlaps_any_exon(variant, transcript):
    return any(
        exon.overlaps(
            contig=variant.contig,
            start=variant.pos,
            end=variant.end_pos)
        for exon in transcript.exons)

def group_by(records, field_name):
    groups = {}
    for record in records:
        value = getattr(record, field_name)
        if value in groups:
            groups[value].append(record)
        else:
            groups[value] = [record]
    return groups

def apply_variant_to_transcript(variant, transcript):
    """
    Generate a transcript effect (such as FrameShift) by applying a genomic
    variant to a particular transcript.

    Parameters
    ----------
    variant : Variant
        Genomic variant

    transcript :  Transcript
        Transcript we're going to mutate
    """

    if not isinstance(variant, Variant):
        raise TypeError(
            "Expected %s : %s to have type Variant" % (
                variant, type(variant)))

    if not isinstance(transcript, Transcript):
        raise TypeError(
            "Expected %s : %s to have type Transcript" % (
                transcript, type(transcript)))

    if not transcript.complete:
        return IncompleteTranscript(variant, transcript)

    if not is_coding_biotype(transcript.biotype):
        return NoncodingTranscript(variant, transcript)

    if not overlaps_any_exon(variant, transcript):
        return Intronic(variant, transcript)

    if transcript.on_backward_strand:
        ref = reverse_complement(variant.ref)
        alt = reverse_complement(variant.alt)
    else:
        ref = variant.ref
        alt = variant.alt

    # offsets into the spliced transcript
    offsets = [
        transcript.spliced_offset(variant.pos),
        transcript.spliced_offset(variant.end_pos)
    ]


    start_offset_with_utr5 = min(offsets)
    end_offset_with_utr5 = max(offsets)

    assert start_offset_with_utr5 >= 0, \
        "Position %d is before start of transcript %s" % (
            start_offset_with_utr5, transcript)
    assert end_offset_with_utr5 >= 0, \
        "Position %d is before start of transcript %s" % (
            end_offset_with_utr5, transcript)

    utr5_length = transcript.first_start_codon_spliced_offset
    if utr5_length > start_offset_with_utr5:
        # TODO: what do we do if the variant spans the beginning of
        # the coding sequence?
        assert utr5_length > end_offset_with_utr5, \
            "Variant which span the 5' UTR and CDS not yet supported: %s" % (
                variant)
        return FivePrimeUTR(variant, transcript)

    # get offsets into coding sequence by subtracting off
    # untranslated region lengths
    # TODO: move subtraction of 5' UTR length into
    # pyensembl.Transcript, call the method "coding_offset"
    cds_start_offset = start_offset_with_utr5 - utr5_length
    cds_end_offset = end_offset_with_utr5 - utr5_length

    # Don't need a pyfaidx.Sequence object here, just convert it to the an str
    cds_seq = str(transcript.coding_sequence)

    if cds_start_offset >= len(cds_seq) and cds_end_offset >= len(cds_seq):
        return ThreePrimeUTR(variant, transcript)

    # past this point we know that we're somewhere in the coding sequence
    cds_ref = cds_seq[cds_start_offset:cds_end_offset+1]

    # Make sure that the reference sequence agrees with what we expected
    # from the VCF
    # TODO: check that the ref allele is correct for UTR by looking
    # at transcript.sequence instead of transcript.coding_sequence
    assert cds_ref == ref, \
        "Expected ref '%s', got '%s' in %s (offset %d:%d)" % (
            ref,
            cds_ref,
            variant,
            cds_start_offset,
            cds_end_offset)

    # turn cDNA sequence into a BioPython sequence, translate
    # to amino acids. For the original CDS make sure that it starts with
    # a start codon and ends with a stop codon. Don't include the stop codon
    # in the translated sequence.
    original_protein = str(Seq(cds_seq).translate(cds=True, to_stop=True))

    assert len(original_protein) > 0, \
        "Translated protein sequence of %s is empty" % (transcript,)
    assert original_protein[0] == "M", \
        "Expected coding sequence of %s to start with M, instead got '%s'" % (
            transcript, original_protein[0])

    variant_cds_seq = mutate(cds_seq, cds_start_offset, ref, alt)

    # In case sequence isn't a multiple of 3, then truncate it
    # TODO: if we get a frameshift by the end of a CDS (e.g. in the stop codon)
    # then we should use some of the 3' UTR to finish translation.
    truncated_variant_cds_seq = variant_cds_seq[:len(variant_cds_seq) / 3 * 3]

    # Can't be sure that the variant is a complete CDS, so passing cds=False
    # in the case of a frameshift, the transcript might not actually contain
    # a stop codon. So, we're for it manually (by passing to_stop=False).
    variant_protein = str(Seq(truncated_variant_cds_seq).translate(
        cds=False, to_stop=False))

    assert len(variant_protein) > 0, \
        "Protein sequence empty for variant %s on transcript %s" % (
            variant, transcript)

    aa_position = int(cds_start_offset / 3) # genomic position to codon position

    if variant_protein[0] != "M":
        assert aa_position == 0, \
            "Unexpected change to start codon (M>%s) when aa_pos=%s" % (
                variant_protein[0], aa_position)
        return StartLoss(
            variant,
            transcript,
            cds_start_offset,
            aa_alt=variant_protein[0])

    # variant_protein sequence includes stop codon, whereas original_protein
    # doesn't
    if variant_protein[-1] == "*" and original_protein == variant_protein[:-1]:
        return Silent(
            variant,
            transcript,
            cds_start_offset,
            aa_ref=variant_protein[aa_position])

    stop_codon_index = variant_protein.find("*")

    # if contained stop codon, truncate sequence before it
    if stop_codon_index > -1:
        variant_protein = variant_protein[:stop_codon_index]

    n_cdna_ref = len(ref)
    n_cdna_alt = len(alt)

    last_aa_ref_position = (cds_start_offset + n_cdna_ref) / 3
    aa_ref = original_protein[aa_position:last_aa_ref_position+1]
    assert len(aa_ref) > 0, \
        "len(aa_ref) = 0 for variant %s on transcript %s (aa_pos=%d:%d)" % (
            variant, transcript, aa_position, last_aa_ref_position)

    # is this a premature stop codon?
    if stop_codon_index == aa_position:
        return PrematureStop(
            variant,
            transcript,
            cds_start_offset,
            aa_ref)

    if abs(n_cdna_ref - n_cdna_alt) % 3 != 0:
        shifted_sequence = variant_protein[aa_position:]
        return FrameShift(
            variant, transcript, cds_start_offset,
            aa_ref, shifted_sequence)

    last_aa_alt_position = (cds_start_offset + n_cdna_alt) / 3
    aa_alt = variant_protein[aa_position:last_aa_alt_position+1]

    assert len(aa_alt) > 0, \
        "len(aa_alt) = 0 for variant %s on transcript %s (aa_pos=%d:%d)" % (
            variant, transcript, aa_position, last_aa_ref_position)

    assert aa_alt != aa_ref, \
        "Unexpected silent mutation for variant %s on transcript %s (aa=%s)" % (
            variant, transcript, aa_ref)

    # Deletion e.g. FYPQ > F
    if aa_ref.startswith(aa_alt):
        assert len(aa_ref) > 0, \
            "Can't have empty ref and alt for variant %s on transcript %s" % (
                variant, transcript)
        # if aa_ref = FYPQ and aa_alt = F, then deleted = YPQ
        n_kept = len(aa_alt)
        deleted = aa_ref[n_kept:]

        return Deletion(
            variant,
            transcript,
            cds_start_offset,
            n_kept,
            deleted)

    # Insertion, e.g. Q > QLLQ
    elif aa_alt.startswith(aa_ref):
        assert len(aa_alt) > 0, \
            "Can't have empty ref and alt for variant %s on transcript %s" % (
                variant, transcript)
        inserted = aa_alt[len(aa_ref):]
        return Insertion(
            variant, transcript,
            cds_start_offset,
            aa_ref,
            inserted)
    else:
        return Substitution(
            variant,
            transcript,
            cds_start_offset,
            aa_ref,
            aa_alt)

def apply_variant(ensembl, variant):
    """
    Determine the effects of a variant on any transcripts it overlaps,
    return the list of overlapping Gene objects, and a dictionary
    mapping gene IDs to a list of transcript mutation effects (e.g. FrameShift).

    Parameters
    ----------

    ensembl : Ensembl
        Ensembl release which lets us query which genes/transcripts are at a
        particular locus.

    variant : Variant
    """

    contig = variant.contig
    pos = variant.pos
    end_pos = variant.end_pos
    ref = variant.ref
    alt = variant.alt

    overlapping_genes = ensembl.genes_at_locus(contig, pos, end_pos)

    if len(overlapping_genes) == 0:
        return [], {}
    else:
        overlapping_transcripts = ensembl.transcripts_at_locus(
                contig, pos, end_pos)

        assert len(overlapping_transcripts) > 0, \
            "No transcripts found for mutation %s:%d %s>%s" % (
                contig, pos, ref, alt)

    # group transcripts by their gene ID
    overlapping_transcript_groups = group_by(
        overlapping_transcripts, field_name='gene_id')

    variant_effect_groups = {}
    for gene_id, transcripts in overlapping_transcript_groups.iteritems():
        effects = []
        for transcript in transcripts:
            effect = apply_variant_to_transcript(variant, transcript)
            effects.append(effect)
        variant_effect_groups[gene_id] = effects
    return overlapping_genes, variant_effect_groups
