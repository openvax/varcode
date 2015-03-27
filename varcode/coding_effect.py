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

def translate(cds_seq):
    """Translates cDNA coding sequence into amino acid protein sequence.

    Should typically start with a start codon but allowing non-methionine
    first residues since the CDS we're translating might have been affected
    by a start loss mutation.

    The sequence may include the 3' UTR but will stop translation at the first
    encountered stop codon.

    Parameters
    ----------
    cds_seq : BioPython Seq
        cDNA coding sequence

    Returns BioPython Seq of amino acids
    """
    # In case sequence isn't a multiple of 3, then truncate it
    truncated_cds_len = int(len(cds_seq) / 3) * 3
    truncated_cds_seq = cds_seq[:truncated_cds_len]

    # turn cDNA sequence into a BioPython sequence, translate
    # to amino acids.
    # passing cds=False since we may want to deal with premature
    # stop codons
    protein = truncated_cds_seq.translate(to_stop=True, cds=False)

    if protein[0] != "M" and cds_seq[:3] in START_CODONS:
        # TODO: figure out when these should be made into methionines
        # and when left as whatever amino acid they normally code for
        # e.g. Leucine start codons
        # See: DOI: 10.1371/journal.pbio.0020397
        return "M" + protein[1:]
    return protein

def snv_coding_effect(
        ref,
        alt,
        cds_seq,
        cds_offset,
        transcript,
        variant):
    """Coding effect of a single nucleotide substitution"""
    aa_pos = int(cds_offset / 3)
    assert aa_pos <= len(transcript.protein_sequence)

    # codon in the reference sequence
    ref_codon = str(cds_seq[aa_pos * 3:aa_pos * 3 + 3])
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
    transcript_ref = sequence[transcript_offset:transcript_offset + len(ref)]

    # Make sure that the reference sequence agrees with what we expected
    # from the VCF
    assert transcript_ref == ref, \
        "%s: expected ref '%s' at offset %d of %s, transcript has '%s'" % (
            variant,
            ref,
            transcript_offset,
            transcript,
            transcript_ref)

    cds_start_offset = min(transcript.start_codon_spliced_offsets)
    cds_stop_offset = max(transcript.stop_codon_spliced_offsets)

    cds_seq = sequence[cds_start_offset:cds_stop_offset + 1]

    if len(ref) == 0 and transcript.strand == "-":
        # for insertions the CDS offset is supposed to point to the
        # nucleotide immediately before the insertion, but on the reverse
        # strand this is actually the nucleotide immediately after the
        # insertion
        # Need to adjust this by moving the CDS offset back one
        cds_offset = transcript_offset - cds_start_offset - 1
    else:
        cds_offset = transcript_offset - cds_start_offset

    assert cds_offset < len(cds_seq), \
        "Expected CDS offset (%d) < |CDS| (%d) for %s on %s" % (
            cds_offset, len(cds_seq), variant, transcript)

    if len(cds_seq) < 3:
        raise ValueError("Coding sequence for %s is too short: '%s'" % (
            transcript, cds_seq))

    original_protein = transcript.protein_sequence

    if not original_protein:
        # Ensembl should have given us a protein sequence for every
        # transcript but it's possible that we're trying to annotate a
        # transcript whose biotype isn't included in the protein sequence FASTA
        logging.warn("No protein sequence for %s in Ensembl", transcript)
        original_protein = translate(cds_seq)

    # subtract 3 for the stop codon and divide by 3 since
    # 3 nucleotides = 1 codon = 1 amino acid
    expected_protein_length = int((len(cds_seq) - 3) / 3)
    if len(original_protein) != expected_protein_length:
        raise ValueError(
            "Expected protein sequence of %s to be %d amino acids"
            "but got %d : %s" % (
                transcript,
                expected_protein_length,
                len(original_protein),
                original_protein))

    # genomic position to codon position
    aa_pos = int(cds_offset / 3)

    # special case simplified logic for an SNV
    # TODO: generalize for all in-frame substitutions
    if len(ref) == len(alt) == 1:
        return snv_coding_effect(
            ref,
            alt,
            cds_seq,
            cds_offset,
            transcript,
            variant)

    transcript_after_start_codon = sequence[cds_start_offset:]

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
    # Further up, we set the CDS offset for insertions on the reverse strand to
    # have an offset one less than they otherwise would, which makes the
    # insertion go to the correct location.
    if len(ref) == 0:
        variant_cds_seq = insert_after(
            transcript_after_start_codon, cds_offset, alt)
    else:
        variant_cds_seq = substitute(
            transcript_after_start_codon,
            cds_offset,
            ref,
            alt)

    variant_protein = translate(variant_cds_seq)

    if len(variant_protein) == 0:
        raise ValueError(
            "Translated mutant protein sequence of %s is empty" % (transcript,))

    if original_protein == variant_protein:
        original_start_codon = cds_seq[:3]
        variant_start_codon = variant_cds_seq[:3]
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
