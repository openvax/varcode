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

# from Bio.Seq import Seq
from Bio.Data import CodonTable

from .effects import (
    Silent,
    Insertion,
    Deletion,
    Substitution,
    ComplexSubstitution,
    PrematureStop,
    AlternateStartCodon,
    StartLoss,
    StopLoss,
    FrameShift,
    FrameShiftTruncation,
    ThreePrimeUTR,
)
from .mutate import substitute, insert_after
from .string_helpers import trim_shared_flanking_strings

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
        to_stop=True):
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

    Returns BioPython Seq of amino acids
    """
    if len(nucleotide_sequence) < 3:
        raise ValueError("Sequence '%s' is too short to translate" % (
            nucleotide_sequence))

    # if sequence isn't a multiple of 3, truncate it so BioPython
    # doesn't complain
    truncated_cds_len = int(len(nucleotide_sequence) / 3) * 3
    truncated_cds_seq = nucleotide_sequence[:truncated_cds_len]

    # passing cds=False to translate since we may want to deal with premature
    # stop codons
    protein_sequence = truncated_cds_seq.translate(to_stop=to_stop, cds=False)

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
                 truncated_cds_seq))

    return protein_sequence

def transcript_protein_sequence(transcript):
    cds_start_offset = transcript.first_start_codon_spliced_offset
    cds_stop_offset = transcript.last_stop_codon_spliced_offset
    cds_len = cds_stop_offset - cds_start_offset + 1

    original_protein = transcript.protein_sequence

    if not original_protein:
        # Ensembl should have given us a protein sequence for every
        # transcript but it's possible that we're trying to annotate a
        # transcript whose biotype isn't included in the protein sequence FASTA
        logging.warn("No protein sequence for %s in Ensembl", transcript)

        original_protein = translate(
            transcript.sequence[cds_start_offset:cds_stop_offset + 1])

    # subtract 3 for the stop codon and divide by 3 since
    # 3 nucleotides = 1 codon = 1 amino acid
    expected_protein_length = int((cds_len - 3) / 3)
    if len(original_protein) != expected_protein_length:
        raise ValueError(
            "Expected protein sequence of %s to be %d amino acids"
            "but got %d : %s" % (
                transcript,
                expected_protein_length,
                len(original_protein),
                original_protein))
    return original_protein

def frameshift_insertion_effect(
        cds_offset_before_insertion,
        inserted_nucleotides,
        sequence_from_start_codon,
        variant,
        transcript):
    """
    Assumption:
        The insertion is happening after the start codon and before the stop
        codon of this coding sequence.
    """
    # TODO:
    # get the protein sequence until the first modified codon
    # translate the insert + cds_after_insertion
    # concatenate

    if cds_offset_before_insertion % 3 == 2:
        # if insertion happens after last nucleotide in a codons
        codon_index_before_insertion = int(cds_offset_before_insertion / 3)
    else:
        # if insertion happens after the 1st or 2nd nucleotide in a codon,
        # then it disrupts that codon
        codon_index_before_insertion = int(cds_offset_before_insertion / 3) - 1

    assert codon_index_before_insertion >= 0, \
        "Expected frameshift_insertion to be after start codon for %s on %s" % (
            variant, transcript)
    assert codon_index_before_insertion < len(original_protein_sequence) - 1, \
        "Expected frameshift_insertion to be before stop codon for %s on %s" % (
            variant, transcript)

    original_protein_sequence = transcript.protein_sequence

    if codon_index_before_insertion == len(original_protein_sequence) - 1
        # if insertion is into the stop codon then this is a stop-loss
    protein_before_insertion = \
        transcript.protein_sequence[:codon_index_before_insertion + 1]


    cds_offset_after_insertion = (codon_index_before_insertion + 1) * 3
    original_coding_sequence_after_insertion = \
        sequence_from_start_codon[cds_offset_after_insertion:]
    coding_sequence_after_insertion = \
        inserted_nucleotides + original_coding_sequence_after_insertion
    protein_after_insertion = translate(
        nucleotide_sequence=coding_sequence_after_insertion,
        first_codon_is_start=(codon_index_before_insertion == -1),
        to_stop=True)

    if len(protein_after_insertion) == 0:
        return FrameShiftTruncation()
    else:
        return FrameShift()

def in_frame_insertion_effect(
        cds_offset_before_insertion,
        inserted_nucleotides,
        sequence_from_start_codon,
        variant,
        transcript):
    pass


def insertion_effect(
        inserted_nucleotides,
        transcript_offset,
        transcript,
        variant):
    assert len(inserted_nucleotides) > 0, \
        "Expected len(inserted_nucleotides) > 0 for %s on %s" % (
            transcript, variant)

    assert cds_offset_before_insertion < cds_len, \
        "Expected CDS offset (%d) < |CDS| (%d) for %s on %s" % (
            cds_offset, cds_len, variant, transcript)

    if
    # IF insertion into start codon, StartLoss

    # IF insertion into stop codon, StopLoss (scan forward to next stop?)

    if len(inserted_nucleotides) % 3 != 0:
        return frameshift_insertion_effect()
    else:
        return in_frame_insertion_effect()

def in_frame_coding_effect(
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

    if len(ref) == 0:
        return in_frame_insertion_effect()

    n_ref_nucleotides = len(ref)
    n_alt_nucleotides = len(alt)

    original_protein_sequence = transcript.protein_sequence
    first_ref_amino_acid_index = int(cds_offset / 3)
    assert first_ref_amino_acid_index <= len(original_protein_sequence), \
        ("Unexpected mutation at offset %d (5' UTR starts at %d"
         " while annotating %s on %s") % (
         first_ref_amino_acid_index,
         len(transcript.protein_sequence))

    last_ref_amino_acid_index = int((cds_offset + n_ref_nucleotides - 1) / 3)

    assert last_ref_amino_acid_index >= first_ref_amino_acid_index, \
        ("Expected first_ref_amino_acid_index (%d) <="
         "last_ref_amino_acid_index (%d) while annotating %s on %s") % (
         first_ref_amino_acid_index,
         last_ref_amino_acid_index,
         variant,
         transcript)


    # codon in the reference sequence
    ref_codons = str(cds_seq_with_utr3[
        first_ref_amino_acid_index * 3:last_ref_amino_acid_index * 3 + 3])
    # which nucleotide of the codon got changed?
    codon_offset = cds_offset % 3
    mutant_codon = (
        ref_codon[:codon_offset] + alt + ref_codon[codon_offset + 1:])
    assert len(mutant_codon) == 3, \
        "Expected codon to have length 3, got %s (length = %d)" % (
            mutant_codon, len(mutant_codon))
    if aa_pos == 0:
        if mutant_codon in START_CODONS:
            # if we changed the starting codon treat then
            # this is technically a Silent mutation but may cause
            # alternate starts or other effects
            return AlternateStartCodon(
                variant,
                transcript,
                ref_codon,
                mutant_codon)
        else:
            # if we changed a start codon to something else then
            # we no longer know where the protein begins (or even in
            # what frame).
            # TODO: use the Kozak consensus sequence or a predictive model
            # to identify the most likely start site
            return StartLoss(
                variant=variant,
                transcript=transcript,
                aa_alt=translate_codon(mutant_codon, 0))

    original_amino_acid = translate_codon(ref_codon, aa_pos)
    mutant_amino_acid = translate_codon(mutant_codon, aa_pos)

    if original_amino_acid == mutant_amino_acid:
        return Silent(
            variant,
            transcript,
            aa_pos=aa_pos,
            aa_ref=original_amino_acid)
    elif aa_pos == len(transcript.protein_sequence):
        # if non-silent mutation is at the end of the protein then
        # should be a stop-loss
        assert original_amino_acid == "*"
        # if mutatin is at the end of the protein and both
        assert mutant_amino_acid != "*"
        return StopLoss(
            variant,
            transcript,
            aa_pos=aa_pos,
            aa_alt=mutant_amino_acid)
    else:
        # simple substitution e.g. p.V600E
        return Substitution(
            variant,
            transcript,
            aa_pos=aa_pos,
            aa_ref=original_amino_acid,
            aa_alt=mutant_amino_acid)



def coding_effect(
        ref,
        alt,
        transcript_offset,
        transcript,
        variant):
    """
    Given a minimal ref/alt nucleotide string pair and an offset into a given
    transcript, determine the coding effect of this nucleotide substitution
    onto the translated protein.

    Parameters
    ----------
    ref : str
        Reference nucleotides we expect to find in the transcript's CDS

    alt : str
        Alternate nucleotides we're replacing the reference with

    transcript_offset : int
        Offset into the full transcript sequence of the ref->alt substitution

    transcript : Transcript

    variant : Variant
    """
    if not transcript.complete:
        raise ValueError(
            ("Can't annotate coding effect for %s"
             " on incomplete transcript %s" % (variant, transcript)))

    sequence = transcript.sequence

    # reference nucleotides found on the transcript, if these don't match
    # what we were told to expect from the variant then raise an exception
    ref_nucleotides_from_transcript = \
        sequence[transcript_offset:transcript_offset + len(ref)]

    # Make sure that the reference sequence agrees with what we expected
    # from the VCF
    assert ref_nucleotides_from_transcript == ref, \
        "%s: expected ref '%s' at offset %d of %s, transcript has '%s'" % (
            variant,
            ref,
            transcript_offset,
            transcript,
            ref_nucleotides_from_transcript)

    cds_start_offset = transcript.first_start_codon_spliced_offset
    cds_stop_offset = transcript.last_stop_codon_spliced_offset

    cds_len = cds_start_offset - cds_stop_offset + 1

    if cds_len < 3:
        raise ValueError(
            "Coding sequence for %s is too short: '%s'" % (
                transcript,
                transcript.sequence[cds_start_offset:cds_stop_offset + 1]))

    if len(ref) == 0 and transcript.strand == "-":
        # By convention, genomic insertions happen *after* their base 1 position on
        # a chromosome. On the reverse strand, however, an insertion has to go
        # before the nucleotide at some transcript offset.
        # Example:
        #    chromosome sequence:
        #        TTT|GATCTCGTA|CCC
        #    transcript on reverse strand:
        #        CCC|ATGCTCTAG|TTT
        #    where the CDS is emphasized:
        #            ATGCTCTAG
        # If we have a genomic insertion g.6insATT
        # the genomic sequence becomes:
        #       TTT|GAT_ATT_CTCGTA|CCC
        # (insert the "ATT" after the "T" at position 6)
        # On the reverse strand this becomes:
        #       CCC|ATGCTC_TTA_TAG|TTT
        # (insert the "ATT" *before* the "T" at position 10)
        #
        # To preserve the interpretation of the start offset as the base
        # before the insertion, need to subtract one
        cds_offset = transcript_offset - cds_start_offset - 1
    else:
        cds_offset = transcript_offset - cds_start_offset

    assert cds_offset < cds_len, \
        "Expected CDS offset (%d) < |CDS| (%d) for %s on %s" % (
            cds_offset, cds_len, variant, transcript)

    # did the mutation disrupt the start codon?

    # did the mutation disrupt the stop codon?

    # since insertions create lots of special cases for "base counting"
    # inclusive indexing, let's handle all insertion logic in its own
    # function
    if len(ref) == 0:
        return insertion_effect(
            inserted_nucleotides=alt,
            cds_offset_before_insertion=cds_offset,
            transcript=transcript,
            variant=variant)

    # past this point, we know the mutation is not a start-loss, a stop-loss
    # or a pure insertion

    sequence_after_start_codon = sequence[cds_start_offset:]

    # is this an in-frame mutations?
    if (len(ref) - len(alt)) % 3 == 0:
        return in_frame_coding_effect(
            ref,
            alt,
            cds_offset,
            sequence_after_start_codon,
            transcript,
            variant)
    else:



    #
    # Further up, we set the CDS offset for insertions on the reverse strand to
    # have an offset one less than they otherwise would, which makes the
    # insertion go to the correct location.
    if len(ref) == 0:
        variant_sequence = insert_after(
            sequence_after_start_codon,
            cds_offset,
            alt)
    else:
        variant_sequence = substitute(
            sequence_after_start_codon,
            cds_offset,
            ref,
            alt)

    variant_protein = translate(variant_sequence)

    if len(variant_protein) == 0:
        raise ValueError(
            "Translated mutant protein sequence of %s is empty" % (transcript,))

    # genomic position to codon position
    aa_pos = int(cds_offset / 3)

    if original_protein == variant_protein:
        original_start_codon = sequence_after_start_codon[:3]
        variant_start_codon = variant_sequence[:3]
        if original_start_codon != variant_start_codon:
            # mutation is silent on the amino acid sequence but
            # uses a different start codon, which may cause the transcript
            # to not be translated or translated in unexpected ways
            return AlternateStartCodon(
                variant,
                transcript,
                original_start_codon,
                variant_start_codon)
        elif aa_pos < len(original_protein):
            aa_ref = original_protein[aa_pos]
        elif aa_pos == len(original_protein):
            aa_ref = "*"
        elif aa_pos > len(original_protein):
            # We got into this function because the mutation was expected
            # to start in the coding sequence
            # If the first affected amino acid is after the end of the original
            # protein then it's most likely that the stop codon used to
            # terminate translation was actually a selenocysteine.
            # TODO: look up selenocysteine annotations and pass them
            # to translate.
            if cds_seq[:len(original_protein) * 3 + 3].endswith("TGA"):
                logging.info(
                    ("Possible selenocysteine codon (TGA)"
                     " at position %d of %s") % (
                        aa_pos * 3,
                        transcript))
                return ThreePrimeUTR(variant, transcript)
            logging.warn(
                ("Unexpected aa_pos = %d  for len(protein) = %d"
                 " in 3' UTR of %s for %s"),
                aa_pos,
                len(original_protein),
                transcript,
                variant)
            aa_ref = "?"
        return Silent(
            variant,
            transcript,
            aa_pos=aa_pos,
            aa_ref=aa_ref)

    if aa_pos == 0 and (
            variant_protein[0] != original_protein[0] or
            len(original_protein) > len(variant_protein)):
        # if change is in first codon of the protein and it either
        # changes the amino acid or truncates the protein, consider that
        # a start loss
        return StartLoss(
            variant=variant,
            transcript=transcript,
            aa_alt=variant_protein[0])
    elif aa_pos == len(variant_protein):
        # is this a premature stop codon?
        last_codon = variant_cds_seq[aa_pos * 3:aa_pos * 3 + 3]
        if last_codon not in STOP_CODONS:
            # if protein ends at the mutation point but there wasn't a stop
            # codon there?
            logging.warn(
                ("Truncated protein doesn't end with stop codon for %s"
                " on %s, original len = %d, mutant len = %d") % (
                    variant,
                    transcript,
                    len(original_protein),
                    len(variant_protein)))
        return PrematureStop(
            variant,
            transcript,
            cds_offset,
            aa_ref=original_protein[aa_pos])
    elif aa_pos == len(original_protein):
        # if mutation begins at the stop codon of this protein and isn't silent
        if len(variant_protein) == len(original_protein):
            logging.info(
                "Expected non-silent stop-loss variant to cause longer "
                "protein but got len(original) = len(variant) = %d for "
                "%s, transcript probably lacks 3' UTR" % (
                    len(variant_protein),
                    transcript))
        aa_alt = variant_protein[aa_pos:]
        return StopLoss(
            variant,
            transcript,
            aa_pos=aa_pos,
            aa_alt=aa_alt)

    elif aa_pos >= len(original_protein):
        # we hit an early stop codon which, in some individuals,
        # is mutated into an amino acid codon
        if transcript.biotype == "polymorphic_pseudogene":
            return ThreePrimeUTR(variant, transcript)
        # Selenocysteine hijack the TGA stop codon
        # See: http://en.wikipedia.org/wiki/Selenocysteine
        elif cds_seq[:len(original_protein) * 3 + 3].endswith("TGA"):
            logging.info(
                "Possible selenocysteine codon (TGA) at position %d of %s" % (
                    aa_pos * 3,
                    transcript))
            return ThreePrimeUTR(variant, transcript)
        else:
            raise ValueError(
                ("Expected aa_pos (%d) < |protein| (%d)"
                 " for %s on %s (CDS offset = %d/%d)" % (
                    aa_pos,
                    len(original_protein),
                    variant,
                    transcript,
                    cds_offset,
                    len(cds_seq))))

    frameshift = False

    # does the mutation shift the open reading frame?
    if abs(len(ref) - len(alt)) % 3 != 0:
        frameshift = True
        aa_alt = variant_protein[aa_pos:]

    # the position of deleted amino acids on the variant protein
    # will be from aa_pos:aa_pos, where aa_pos is the last position before
    # the deleted residues
    elif len(alt) == 0:
        aa_alt = ""
    # insertions happen after cds_offset, so we need slightly different logic
    # for them than a substitution, whose variant nucleotides start
    #  *at* cds_offset
    elif len(ref) == 0:
        last_aa_alt_pos = int((cds_offset + len(alt)) / 3)
        aa_alt = variant_protein[aa_pos:last_aa_alt_pos + 1]
    # if not a frameshift, insertion, deletion, or premature stop,
    # then pull out the new or modified amino acids into `aa_alt`
    # and determine the type of mutation later
    else:
        last_aa_alt_pos = int((cds_offset + len(alt) - 1) / 3)
        aa_alt = variant_protein[aa_pos:last_aa_alt_pos + 1]

    assert len(alt) == 0 or len(aa_alt) > 0, \
            "len(aa_alt) = 0 for variant %s on transcript %s (aa_pos=%d)" % (
                variant, transcript, aa_pos)
    last_aa_ref_pos = int((cds_offset + max(0, len(ref) - 1)) / 3)
    aa_ref = original_protein[aa_pos:last_aa_ref_pos + 1]
    assert len(aa_ref) > 0, \
        "len(aa_ref) = 0 for variant %s on transcript %s (aa_pos=%d:%d)" % (
            variant, transcript, aa_pos, last_aa_ref_pos)

    # in case of simple insertion like FY>FYGL or deletions FYGL > FY,
    # get rid of the shared prefixes/suffixes
    aa_ref, aa_alt, prefix, suffix = trim_shared_flanking_strings(
        aa_ref, aa_alt)
    aa_pos += len(prefix)

    if frameshift:
        aa_ref = original_protein[aa_pos]
        # if a frameshift doesn't create any new amino acids, then
        # it must immediately have hit a stop codon
        if len(aa_alt) == 0:
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
                shifted_sequence=aa_alt)

    # Deletion e.g. p.389delQQ
    if len(aa_alt) == 0:
        assert len(aa_ref) > 0, \
            ("Can't have empty aa_ref and aa_alt for variant %s on"
             " transcript %s, shared prefix = '%s', shared suffix = '%s'") % (
             variant, transcript, prefix, suffix)
        return Deletion(
            variant,
            transcript,
            aa_pos=aa_pos,
            aa_ref=aa_ref)

    # Insertion, e.g. p.37insA
    elif len(aa_ref) == 0:
        assert len(aa_alt) > 0, \
            ("Can't have ref = '' and alt = '%s' at aa_pos = %d, cds_pos = %d"
             " for variant %s on transcript %s with shared prefix ='%s',"
             " shared suffix = '%s'") % (
                aa_alt,
                aa_pos,
                cds_offset,
                variant,
                transcript,
                prefix,
                suffix)
        return Insertion(
            variant, transcript,
            aa_pos=aa_pos,
            aa_alt=aa_alt)
    elif len(aa_ref) == len(aa_alt) == 1:
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
