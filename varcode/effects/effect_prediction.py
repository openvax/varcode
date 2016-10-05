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

from __future__ import print_function, division, absolute_import
import logging

from Bio.Seq import reverse_complement
from pyensembl import Transcript

from ..common import groupby_field

from .transcript_helpers import interval_offset_on_transcript
from .effect_helpers import changes_exonic_splice_site
from .effect_collection import EffectCollection
from .effect_prediction_coding import predict_variant_coding_effect_on_transcript
from .effect_classes import (
    Failure,
    Intergenic,
    Intragenic,
    NoncodingTranscript,
    IncompleteTranscript,
    FivePrimeUTR,
    ThreePrimeUTR,
    Intronic,
    IntronicSpliceSite,
    SpliceAcceptor,
    SpliceDonor,
    StartLoss,
    ExonLoss,
    ExonicSpliceSite,
)


logger = logging.getLogger(__name__)


def predict_variant_effects(variant, raise_on_error=False):
    """Determine the effects of a variant on any transcripts it overlaps.
    Returns an EffectCollection object.

    Parameters
    ----------
    variant : Variant

    raise_on_error : bool
        Raise an exception if we encounter an error while trying to
        determine the effect of this variant on a transcript, or simply
        log the error and continue.
    """
    # if this variant isn't overlapping any genes, return a
    # Intergenic effect
    # TODO: look for nearby genes and mark those as Upstream and Downstream
    # effects
    if len(variant.gene_ids) == 0:
        effects = [Intergenic(variant)]
    else:
        # list of all MutationEffects for all genes & transcripts
        effects = []

        # group transcripts by their gene ID
        transcripts_grouped_by_gene = groupby_field(variant.transcripts, 'gene_id')

        # want effects in the list grouped by the gene they come from
        for gene_id in sorted(variant.gene_ids):
            if gene_id not in transcripts_grouped_by_gene:
                # intragenic variant overlaps a gene but not any transcripts
                gene = variant.ensembl.gene_by_id(gene_id)
                effects.append(Intragenic(variant, gene))
            else:
                # gene ID  has transcripts overlapped by this variant
                for transcript in transcripts_grouped_by_gene[gene_id]:
                    if raise_on_error:
                        effect = predict_variant_effect_on_transcript(
                            variant=variant,
                            transcript=transcript)
                    else:
                        effect = predict_variant_effect_on_transcript_or_failure(
                            variant=variant,
                            transcript=transcript)
                    effects.append(effect)
    return EffectCollection(effects)


def predict_variant_effect_on_transcript_or_failure(variant, transcript):
    """
    Try predicting the effect of a variant on a particular transcript but
    suppress raised exceptions by converting them into `Failure` effect
    values.
    """
    try:
        return predict_variant_effect_on_transcript(
            variant=variant,
            transcript=transcript)
    except (AssertionError, ValueError) as error:
        logger.warn(
            "Encountered error annotating %s for %s: %s",
            variant,
            transcript,
            error)
        return Failure(variant, transcript)

def predict_variant_effect_on_transcript(variant, transcript):
        """Return the transcript effect (such as FrameShift) that results from
        applying this genomic variant to a particular transcript.

        Parameters
        ----------
        transcript :  Transcript
            Transcript we're going to apply mutation to.
        """

        if transcript.__class__ is not Transcript:
            raise TypeError(
                "Expected %s : %s to have type Transcript" % (
                    transcript, type(transcript)))

        # check for non-coding transcripts first, since
        # every non-coding transcript is "incomplete".
        if not transcript.is_protein_coding:
            return NoncodingTranscript(variant, transcript)

        if not transcript.complete:
            return IncompleteTranscript(variant, transcript)

        # since we're using inclusive base-1 coordinates,
        # checking for overlap requires special logic for insertions
        is_insertion = variant.is_insertion

        # determine if any exons are deleted, and if not,
        # what is the closest exon and how far is this variant
        # from that exon (overlapping the exon = 0 distance)
        completely_lost_exons = []

        # list of which (exon #, Exon) pairs this mutation overlaps
        overlapping_exon_numbers_and_exons = []

        distance_to_nearest_exon = float("inf")

        start_in_exon = False
        end_in_exon = False

        nearest_exon = None

        variant_start = variant.trimmed_base1_start
        variant_end = variant.trimmed_base1_end

        for i, exon in enumerate(transcript.exons):
            if variant_start <= exon.start and variant_end >= exon.end:
                completely_lost_exons.append(exon)

            if is_insertion and exon.strand == "+" and variant_end == exon.end:
                # insertions after an exon don't overlap the exon
                distance = 1
            elif is_insertion and exon.strand == "-" and variant_start == exon.start:
                distance = 1
            else:
                distance = exon.distance_to_interval(variant_start, variant_end)

            if distance == 0:
                overlapping_exon_numbers_and_exons.append((i + 1, exon))
                # start is contained in current exon
                if exon.start <= variant_start <= exon.end:
                    start_in_exon = True
                # end is contained in current exon
                if exon.end >= variant_end >= exon.start:
                    end_in_exon = True
            elif distance < distance_to_nearest_exon:
                    distance_to_nearest_exon = distance
                    nearest_exon = exon

        if len(overlapping_exon_numbers_and_exons) == 0:
            intronic_effect_class = choose_intronic_effect_class(
                variant=variant,
                nearest_exon=nearest_exon,
                distance_to_exon=distance_to_nearest_exon)
            return intronic_effect_class(
                variant=variant,
                transcript=transcript,
                nearest_exon=nearest_exon,
                distance_to_exon=distance_to_nearest_exon)
        elif len(completely_lost_exons) > 0 or (
                len(overlapping_exon_numbers_and_exons) > 1):
            # if spanning multiple exons, or completely deleted an exon
            # then consider that an ExonLoss mutation
            exons = [exon for (_, exon) in overlapping_exon_numbers_and_exons]
            return ExonLoss(variant, transcript, exons)

        assert len(overlapping_exon_numbers_and_exons) == 1

        exon_number, exon = overlapping_exon_numbers_and_exons[0]

        exonic_effect_annotation = exonic_transcript_effect(
            variant, exon, exon_number, transcript)

        # simple case: both start and end are in the same
        if start_in_exon and end_in_exon:
            return exonic_effect_annotation
        elif isinstance(exonic_effect_annotation, ExonicSpliceSite):
            # if mutation bleeds over into intro but even just
            # the exonic portion got annotated as an exonic splice site
            # then return it
            return exonic_effect_annotation

        return ExonicSpliceSite(
            variant=variant,
            transcript=transcript,
            exon=exon,
            alternate_effect=exonic_effect_annotation)

def choose_intronic_effect_class(
        variant,
        nearest_exon,
        distance_to_exon):
    """
    Infer effect of variant which does not overlap any exon of
    the given transcript.
    """
    assert distance_to_exon > 0, \
        "Expected intronic effect to have distance_to_exon > 0, got %d" % (
            distance_to_exon,)

    if nearest_exon.strand == "+":
        # if exon on positive strand
        start_before = variant.trimmed_base1_start < nearest_exon.start
        start_same = variant.trimmed_base1_start == nearest_exon.start
        before_exon = start_before or (variant.is_insertion and start_same)
    else:
        # if exon on negative strand
        end_after = variant.trimmed_base1_end > nearest_exon.end
        end_same = variant.trimmed_base1_end == nearest_exon.end
        before_exon = end_after or (variant.is_insertion and end_same)

    # distance cutoffs based on consensus splice sequences from
    # http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2947103/
    # 5' splice site: MAG|GURAGU consensus
    #   M is A or C; R is purine; | is the exon-intron boundary
    # 3' splice site: YAG|R
    if distance_to_exon <= 2:
        if before_exon:
            # 2 last nucleotides of intron before exon are the splice
            # acceptor site, typically "AG"
            return SpliceAcceptor
        else:
            # 2 first nucleotides of intron after exon are the splice donor
            # site, typically "GT"
            return SpliceDonor
    elif not before_exon and distance_to_exon <= 6:
        # variants in nucleotides 3-6 at start of intron aren't as certain
        # to cause problems as nucleotides 1-2 but still implicated in
        # alternative splicing
        return IntronicSpliceSite
    elif before_exon and distance_to_exon <= 3:
        # nucleotide -3 before exon is part of the 3' splicing
        # motif but allows for more degeneracy than the -2, -1 nucleotides
        return IntronicSpliceSite
    else:
        # intronic mutation unrelated to splicing
        return Intronic

def exonic_transcript_effect(variant, exon, exon_number, transcript):
    """Effect of this variant on a Transcript, assuming we already know
    that this variant overlaps some exon of the transcript.

    Parameters
    ----------
    variant : Variant

    exon : pyensembl.Exon
        Exon which this variant overlaps

    exon_number : int
        Index (starting from 1) of the given exon in the transcript's
        sequence of exons.

    transcript : pyensembl.Transcript
    """

    genome_ref = variant.trimmed_ref
    genome_alt = variant.trimmed_alt
    variant_start = variant.trimmed_base1_start
    variant_end = variant.trimmed_base1_end

    # clip mutation to only affect the current exon
    if variant_start < exon.start:
        # if mutation starts before current exon then only look
        # at nucleotides which overlap the exon
        assert len(genome_ref) > 0, "Unexpected insertion into intron"
        n_skip_start = exon.start - variant_start
        genome_ref = genome_ref[n_skip_start:]
        genome_alt = genome_alt[n_skip_start:]
        genome_start = exon.start
    else:
        genome_start = variant_start

    if variant_end > exon.end:
        # if mutation goes past exon end then only look at nucleotides
        # which overlap the exon
        n_skip_end = variant_end - exon.end
        genome_ref = genome_ref[:-n_skip_end]
        genome_alt = genome_alt[:len(genome_ref)]
        genome_end = exon.end
    else:
        genome_end = variant_end

    transcript_offset = interval_offset_on_transcript(
        genome_start, genome_end, transcript)

    if transcript.on_backward_strand:
        cdna_ref = reverse_complement(genome_ref)
        cdna_alt = reverse_complement(genome_alt)
    else:
        cdna_ref = genome_ref
        cdna_alt = genome_alt

    n_ref = len(cdna_ref)

    expected_ref = str(
        transcript.sequence[transcript_offset:transcript_offset + n_ref])

    if cdna_ref != expected_ref:
        raise ValueError(
            ("Found ref nucleotides '%s' in sequence"
             " of %s at offset %d (chromosome positions %d:%d)"
             " but variant %s has '%s'") % (
                 expected_ref,
                 transcript,
                 transcript_offset,
                 genome_start,
                 genome_end,
                 variant,
                 cdna_ref))

    utr5_length = min(transcript.start_codon_spliced_offsets)

    # does the variant start inside the 5' UTR?
    if utr5_length > transcript_offset:
        # does the variant end after the 5' UTR, within the coding region?
        if utr5_length < transcript_offset + n_ref:
            # TODO: we *might* lose the Kozak sequence or the start codon
            # but without looking at the modified sequence how can we tell
            # for sure that this is a start-loss variant?
            return StartLoss(variant, transcript)
        else:
            # if variant contained within 5' UTR
            return FivePrimeUTR(variant, transcript)

    utr3_offset = max(transcript.stop_codon_spliced_offsets) + 1

    if transcript_offset >= utr3_offset:
        return ThreePrimeUTR(variant, transcript)

    exon_start_offset = interval_offset_on_transcript(
        exon.start, exon.end, transcript)
    exon_end_offset = exon_start_offset + len(exon) - 1

    # Further below we're going to try to predict exonic splice site
    # modifications, which will take this effect_annotation as their
    # alternative hypothesis for what happens if splicing doesn't change.
    # If the mutation doesn't affect an exonic splice site, then
    # we'll just return this effect.
    coding_effect_annotation = predict_variant_coding_effect_on_transcript(
        variant=variant,
        transcript=transcript,
        trimmed_cdna_ref=cdna_ref,
        trimmed_cdna_alt=cdna_alt,
        transcript_offset=transcript_offset)

    if changes_exonic_splice_site(
            transcript=transcript,
            transcript_ref=cdna_ref,
            transcript_alt=cdna_alt,
            transcript_offset=transcript_offset,
            exon_start_offset=exon_start_offset,
            exon_end_offset=exon_end_offset,
            exon_number=exon_number):
        return ExonicSpliceSite(
            variant=variant,
            transcript=transcript,
            exon=exon,
            alternate_effect=coding_effect_annotation)
    return coding_effect_annotation
