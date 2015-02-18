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

from __future__ import print_function, division, absolute_import

import logging
import math
from collections import namedtuple

from .common import reverse_complement, trim_shared_flanking_strings
from .variant import Variant
from .transcript_mutation_effects import (
    NoncodingTranscript,
    IncompleteTranscript,
    FivePrimeUTR,
    ThreePrimeUTR,
    Intronic,
    Silent,
    Insertion,
    Deletion,
    Substitution,
    ComplexSubstitution,
    PrematureStop,
    StartLoss,
    FrameShift,
    FrameShiftTruncation
)

from Bio.Seq import Seq, CodonTable
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



def overlaps_any_exon(variant, transcript):
    return any(
        exon.overlaps(
            contig=variant.contig,
            start=variant.pos,
            end=variant.end_pos)
        for exon in transcript.exons)

def first_transcript_offset(start_pos, end_pos, transcript):
    """
    Given a variant and a particular transcript,
    return the start/end offsets of the variant relative to the
    chromosomal positions of the transcript.
    """
    # ensure that start_pos:end_pos overlap with transcript positions
    assert start_pos <= end_pos, \
        "start_pos %d shouldn't be greater than end_pos %d" % (
            start_pos, end_pos)
    assert start_pos <= transcript.end, \
        "Range %d:%d starts after transcript %s (%d:%d)" % (
            start_pos,
            end_pos,
            transcript,
            transcript.start,
            transcript.end)
    assert end_pos >= transcript.start, \
        "Range %d:%d ends before transcript %s (%d:%d)" % (
            start_pos,
            end_pos,
            transcript,
            transcript.start,
            transcript.end)
    # trim the start position to the beginning of the transcript
    if start_pos < transcript.start:
        start_pos = transcript.start
    # trim the end position to the end of the transcript
    if end_pos > transcript.end:
        end_pos = transcript.end

    # offsets into the spliced transcript
    offsets = [
        transcript.spliced_offset(start_pos),
        transcript.spliced_offset(end_pos)
    ]

    start_offset_with_utr5 = min(offsets)

    assert start_offset_with_utr5 >= 0, \
        "Position %d is before start of transcript %s" % (
            start_offset_with_utr5, transcript)

    return start_offset_with_utr5

def normalize_variant_fields(variant, transcript):
    """
    Given a Variant (with fields `ref`, `alt`, `pos` and `end_pos`),
    compute the offset the variant position relative
    to the position and direction of the given transcript.
    Also reverse the nucleotide strings if the transcript is on the
    backwards strand and trim any shared substrings from the beginning
    or end of `ref` and `alt`.

    Example:
    If Variant(ref="ATGGC", alt="TCGC", pos=3000, end_pos=3005) is normalized
    relative to Transcript(start=2000, end=4000) then the result of this
    function will be:
        ("G", "C", 1002)
    where the elements of this tuple are (ref, alt, offset).
    """
    if transcript.on_backward_strand:
        ref = reverse_complement(variant.ref)
        alt = reverse_complement(variant.alt)
    else:
        ref = variant.ref
        alt = variant.alt

    start_offset_with_utr5 = first_transcript_offset(
        variant.pos, variant.end_pos, transcript)

    # in case nucleotide strings share prefix (e.g. ref="C", alt="CC")
    # bump the offsets and make the strings disjoint (ref="", alt="C")
    ref, alt, prefix, _ = trim_shared_flanking_strings(ref, alt)
    n_same = len(prefix)
    start_offset_with_utr5 += n_same
    return ref, alt, start_offset_with_utr5

def infer_coding_effect(
        ref,
        alt,
        cds_offset,
        variant,
        transcript):
    # Don't need a pyfaidx.Sequence object here, just convert it to the an str
    cds_seq = str(transcript.coding_sequence)

    # past this point we know that we're somewhere in the coding sequence
    cds_ref = cds_seq[cds_offset:cds_offset+len(ref)]

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

    aa_pos = int(cds_offset / 3) # genomic position to codon position

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
        aa_ref = original_protein[aa_pos:last_aa_ref_pos+1]
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
            if original_protein[aa_pos+i] != aa_mutant:
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
        aa_alt = variant_protein[aa_pos:last_aa_alt_pos+1]
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

def infer_exonic_effect(variant, transcript):
    ref, alt, start_offset_with_utr5 = normalize_variant_fields(
        variant, transcript)

    utr5_length = min(transcript.start_codon_spliced_offsets)

    if utr5_length > start_offset_with_utr5:
        # TODO: what do we do if the variant spans the beginning of
        # the coding sequence?
        assert utr5_length > start_offset_with_utr5 + len(ref), \
            "Variant which span the 5' UTR and CDS not yet supported: %s" % (
                variant)
        return FivePrimeUTR(variant, transcript)

    utr3_offset = max(transcript.stop_codon_spliced_offsets) + 1

    if start_offset_with_utr5 >= utr3_offset:
        return ThreePrimeUTR(variant, transcript)

    cds_offset = start_offset_with_utr5 - utr5_length
    return infer_coding_effect(
        ref,
        alt,
        cds_offset,
        variant,
        transcript)

def find_nearest_exon(start, end, exons):
    """
    Finds nearest exon to an intronic variant, returns the
    distance to that exon, along with the Exon object itself.
    """
    best_distance = float("inf")
    best_exon = None
    for exon in exons:
        if exon.start > end:
            # distance from end of variant to start of exon
            # (when variant starts before the exon)
            distance = exon.start - end
        elif exon.end < start:
            distance = start - exon.end
        else:
            assert exon.start <= end and exon.end >= start, \
                "Expected interval %d:%d to overlap with exon %s" % (
                    start, end, exon)
            return 0, exon

        if best_distance > distance:
            best_distance = distance
            best_exon = exon

    return best_exon, best_exon


def infer_intronic_effect(variant, transcript, nearest_exon, distance_to_exon):
    """
    Infer effect of variant which does not overlap any exon of
    the given transcript.
    """
    assert distance_to_exon > 0, \
        "Expected intronic effect to have distance_to_exon > 0, got %d" % (
            distance_to_exon,)
    before_forward_exon = (
        transcript.strand == "+" and
        variant.pos < nearest_exon.start)
    before_backward_exon = (
        transcript.strand == "-" and
        variant.end_pos > nearest_exon.end)
    before_exon = before_forward_exon or before_backward_exon

    if distance_to_exon <= 2:
        if before_exon:
            # 2 last nucleotides of intron before exon are the splice acceptor
            # site
            effect_class = SpliceAcceptor
        else:
            # 2 first nucleotides of intron after exon are the splice donor
            # site
            effect_class = SpliceDonor
    elif not before_exon and distance_to_exon <= 6:
        # variants in nucleotides 3-6 at start of intron aren't as certain
        # to cause problems as nucleotides 1-2 but still implicated in
        # alternative splicing
        effect_class = IntronicSpliceSite
    elif before_exon and distance_to_exon <= 4:
        # nucleotides -4 and -3 before exon are part of the 3' splicing motif
        # but allow for more degeneracy than the -2, -1 nucleotides
        effect_class = IntronicSpliceSite
    else:
        assert distance_to_exon > 6, \
            "Looks like we didn't cover all possible splice site mutations"
        # intronic mutation unrelated to splicing
        effect_class = Intronic
    return effect_class(
            variant=variant,
            transcript=transcript,
            nearest_exon=nearest_exon,
            distance_to_exon=distance_to_exon)

def infer_transcript_effect(variant, transcript):
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

    # check for non-coding transcripts first, since
    # every non-coding transcript is "incomplete".
    if not is_coding_biotype(transcript.biotype):
        return NoncodingTranscript(variant, transcript)

    if not transcript.complete:
        return IncompleteTranscript(variant, transcript)

    distance_to_exon, nearest_exon = find_nearest_exon(
        start=variant.pos,
        end=variant.end_pos,
        exons=transcript.exons)

    if distance_to_exon > 0:
        return infer_intronic_effect(
            variant=variant,
            transcript=transcript,
            nearest_exon=nearest_exon,
            distance_to_exon=distance_to_exon)

    # TODO: exonic splice site mutations
    return infer_exonic_effect(variant, transcript)
