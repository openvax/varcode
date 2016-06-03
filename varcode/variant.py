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
import json

from collections import namedtuple

from Bio.Seq import reverse_complement
from memoized_property import memoized_property
from pyensembl import (
    Transcript,
    cached_release,
    genome_for_reference_name,
    Genome,
    ensembl_grch38,
)
from pyensembl.locus import normalize_chromosome
from pyensembl.biotypes import is_coding_biotype
from typechecks import require_instance

from .coding_effect import coding_effect
from .common import groupby_field
from .effect_helpers import changes_exonic_splice_site
from .effect_collection import EffectCollection
from .effects import (
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
from .nucleotides import (
    normalize_nucleotide_string,
    STANDARD_NUCLEOTIDES,
    is_purine
)
from .string_helpers import trim_shared_flanking_strings
from .transcript_helpers import interval_offset_on_transcript


class Variant(object):
    __slots__ = (
        "contig",
        "start",
        "end",
        "ref",
        "alt",
        "ensembl",
        "original_ref",
        "original_alt",
        "original_start",
        "_transcripts",
        "_genes",)

    def __init__(
            self,
            contig,
            start,
            ref,
            alt,
            ensembl=ensembl_grch38,
            allow_extended_nucleotides=False,
            normalize_contig_name=True):
        """
        Construct a Variant object.

        Parameters
        ----------
        contig : str
            Chromosome that this variant is on

        start : int
            1-based position on the chromosome of first reference nucleotide

        ref : str
            Reference nucleotide(s)

        alt : str
            Alternate nucleotide(s)

        ensembl : Genome or EnsemblRelease
            Object used for determining gene/transcript annotations

        allow_extended_nucleotides : bool
            Extended nucleotides include 'Y' for pyrimidies or 'N' for any base

        normalize_contig_name : bool
            By default the contig name will be normalized by trimming a 'chr'
            prefix and converting all letters to upper-case. If we don't want
            this behavior then pass normalize_contig_name=False.
        """
        # user might supply Ensembl release as an integer, reference name,
        # or pyensembl.Genome object
        if isinstance(ensembl, Genome):
            self.ensembl = ensembl
        elif isinstance(ensembl, int):
            self.ensembl = cached_release(ensembl)
        elif isinstance(ensembl, str):
            self.ensembl = genome_for_reference_name(ensembl)
        else:
            raise TypeError(
                ("Expected ensembl to be an int, string, or pyensembl.Genome "
                 "instance, got %s : %s") % (type(ensembl), str(ensembl)))
        if normalize_contig_name:
            self.contig = normalize_chromosome(contig)
        else:
            self.contig = contig

        if ref != alt and (ref in STANDARD_NUCLEOTIDES) and (
                alt in STANDARD_NUCLEOTIDES):
            # Optimization for common case.
            self.original_ref = self.ref = ref
            self.original_alt = self.alt = alt
            self.original_start = self.start = self.end = int(start)
            return

        # we want to preserve the ref/alt/pos both as they appeared in the
        # original VCF or MAF file but also normalize variants to get rid
        # of shared prefixes/suffixes between the ref and alt nucleotide
        # strings e.g. g.10 CTT>T can be normalized into g.10delCT
        #
        # The normalized variant properties go into fields
        #    Variant.{original_ref, original_alt, original_pos}
        # whereas the trimmed fields are:
        #    Variant.{ref, alt, start, end}

        # the original entries must preserve the number of nucleotides in
        # ref and alt but we still want to normalize e.g. '-' and '.' into ''
        self.original_ref = normalize_nucleotide_string(
            ref,
            allow_extended_nucleotides=allow_extended_nucleotides)
        self.original_alt = normalize_nucleotide_string(
            alt,
            allow_extended_nucleotides=allow_extended_nucleotides)
        self.original_start = int(start)

        # normalize the variant by trimming any shared prefix or suffix
        # between ref and alt nucleotide sequences and then
        # offset the variant position in a strand-dependent manner
        (trimmed_ref, trimmed_alt, prefix, suffix) = (
            trim_shared_flanking_strings(self.original_ref, self.original_alt))

        self.ref = trimmed_ref
        self.alt = trimmed_alt

        if len(trimmed_ref) == 0:
            # insertions must be treated differently since the meaning of a
            # position for an insertion is:
            #   "insert the alt nucleotides after this position"
            #
            # Aside: what if both trimmed ref and alt strings are empty?
            # This means we had a "null" variant, probably from a VCF
            # generated by force-calling mutations which weren't actually
            # found in the sample.
            # Null variants are interepted as inserting zero nucleotides
            # after the whole reference sequence.
            #
            # Start and end both are base-1 nucleotide position before
            # insertion.
            self.start = self.original_start + max(0, len(prefix) - 1)
            self.end = self.start
        else:
            # for substitutions and deletions the [start:end] interval is
            # an inclusive selection of reference nucleotides
            self.start = self.original_start + len(prefix)
            self.end = self.start + len(trimmed_ref) - 1

    @property
    def reference_name(self):
        return self.ensembl.reference_name

    def __str__(self):
        return "Variant(contig='%s', start=%d, ref='%s', alt='%s', reference_name='%s')" % (
            self.contig,
            self.start,
            self.ref,
            self.alt,
            self.reference_name)

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return self.start

    def __lt__(self, other):
        """
        Variants are ordered by locus.
        """
        require_instance(other, Variant, name="variant")
        if self.contig == other.contig:
            return self.start < other.start
        return self.contig < other.contig

    # The identifying fields of a variant
    BasicFields = namedtuple(
        "BasicFields",
        "contig start end ref alt release")

    def fields(self):
        """
        All identifying fields of a variant (contig, pos, ref, alt, genome)
        in a single tuple. This makes for cleaner printing, hashing, and
        comparisons.
        """
        return Variant.BasicFields(
            self.contig,
            self.start,
            self.end,
            self.ref,
            self.alt,
            self.original_start,
            self.original_ref,
            self.original_alt,
            self.ensembl.release)

    def __eq__(self, other):
        if self is other:
            return True
        return (
            self.contig == other.contig and
            self.start == other.start and
            self.end == other.end and
            self.ref == other.ref and
            self.alt == other.alt and
            self.ensembl == other.ensembl)

    def __getstate__(self):
        # When pickle.load'ing, we want the original values vs. the normalized values.
        # Normalization will happen upon object creation.
        return [
            self.contig,
            self.original_start,
            self.original_ref,
            self.original_alt,
            self.ensembl.release
        ]

    def __setstate__(self, fields):
        # Special logic needed for JSON serialization, since EnsemblRelease
        # is not JSON serializable.
        # TODO: this does not work for Genome objects; only EnsemblRelease.
        ensembl = cached_release(fields.pop())
        fields.append(ensembl)

        self.__init__(*fields)

    def to_json(self):
        """
        Serialize to JSON.

        Returns
        ----------
        String giving the Variant instance serialized in JSON format.

        """
        return json.dumps(self.__getstate__())

    @classmethod
    def from_dict(cls, state):
        """
        Unserialize a Variant encoded as a Python dict.
        """
        instance = cls.__new__(cls)
        instance.__setstate__(state)
        return instance

    @classmethod
    def from_json(cls, serialized):
        """
        Unserialize a Variant encoded as a string of JSON.
        """
        return Variant.from_dict(json.loads(serialized))

    @property
    def short_description(self):
        """
        HGVS nomenclature for genomic variants
        More info: http://www.hgvs.org/mutnomen/
        """
        if self.ref == self.alt:
            # no change
            return "chr%s g.%d%s" % (
                self.contig, self.start, self.ref)
        elif len(self.ref) == 0:
            # insertions
            return "chr%s g.%d_%dins%s" % (
                self.contig,
                self.start,
                self.start + 1,
                self.alt)
        elif len(self.alt) == 0:
            # deletion
            return "chr%s g.%d_%ddel%s" % (
                self.contig, self.start, self.end, self.ref)
        else:
            # substitution
            return "chr%s g.%d%s>%s" % (
                self.contig, self.start, self.ref, self.alt)

    @memoized_property
    def transcripts(self):
        return self.ensembl.transcripts_at_locus(
            self.contig, self.start, self.end)

    @property
    def coding_transcripts(self):
        """
        Protein coding transcripts
        """
        return [
            transcript
            for transcript in self.transcripts
            if is_coding_biotype(transcript.biotype)
        ]

    @property
    def transcript_ids(self):
        return [transcript.id for transcript in self.transcripts]

    @property
    def transcript_names(self):
        return [transcript.name for transcript in self.transcripts]

    @memoized_property
    def genes(self):
        """
        Return Gene object for all genes which overlap this variant.
        """
        return self.ensembl.genes_at_locus(
            self.contig, self.start, self.end)

    @property
    def gene_ids(self):
        """
        Return IDs of all genes which overlap this variant. Calling
        this method is significantly cheaper than calling `Variant.genes()`,
        which has to issue many more queries to construct each Gene object.
        """
        return self.ensembl.gene_ids_at_locus(
            self.contig, self.start, self.end)

    @property
    def gene_names(self):
        """
        Return names of all genes which overlap this variant. Calling
        this method is significantly cheaper than calling `Variant.genes()`,
        which has to issue many more queries to construct each Gene object.
        """
        return self.ensembl.gene_names_at_locus(
            self.contig, self.start, self.end)

    @property
    def coding_genes(self):
        """
        Protein coding transcripts
        """
        return [
            gene for gene in self.genes
            if is_coding_biotype(gene.biotype)
        ]

    def effects(self, raise_on_error=True):
        """Determine the effects of a variant on any transcripts it overlaps.
        Returns an EffectCollection object.

        Parameters
        ----------
        raise_on_error : bool
            Raise an exception if we encounter an error while trying to
            determine the effect of this variant on a transcript, or simply
            log the error and continue.
        """
        gene_ids = self.gene_ids

        # if this variant isn't overlapping any genes, return a
        # Intergenic effect
        # TODO: look for nearby genes and mark those as Upstream and Downstream
        # effects
        if len(gene_ids) == 0:
            return EffectCollection([Intergenic(self)])

        # group transcripts by their gene ID
        transcripts_grouped_by_gene = groupby_field(
            self.transcripts,
            'gene_id')

        # list of all MutationEffects for all genes & transcripts
        effects = []

        # want effects in the list grouped by the gene they come from
        for gene_id in sorted(gene_ids):
            if gene_id not in transcripts_grouped_by_gene:
                # intragenic variant overlaps a gene but not any transcripts
                gene = self.ensembl.gene_by_id(gene_id)
                effects.append(Intragenic(self, gene))
            else:
                # gene ID  has transcripts overlapped by this variant
                for transcript in transcripts_grouped_by_gene[gene_id]:
                    try:
                        effect = self.effect_on_transcript(transcript)
                        effects.append(effect)
                    except (AssertionError, ValueError) as error:
                        if raise_on_error:
                            raise
                        else:
                            effects.append(Failure(self, transcript))
                            logging.warn(
                                "Encountered error annotating %s for %s: %s",
                                self,
                                transcript,
                                error)
        return EffectCollection(effects)

    def effect_on_transcript(self, transcript):
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
        if not is_coding_biotype(transcript.biotype):
            return NoncodingTranscript(self, transcript)

        if not transcript.complete:
            return IncompleteTranscript(self, transcript)

        # since we're using inclusive base-1 coordinates,
        # checking for overlap requires special logic for insertions
        insertion = len(self.ref) == 0

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
        for i, exon in enumerate(transcript.exons):
            if self.start <= exon.start and self.end >= exon.end:
                completely_lost_exons.append(exon)

            if insertion and exon.strand == "+" and self.end == exon.end:
                # insertions after an exon don't overlap the exon
                distance = 1
            elif insertion and exon.strand == "-" and self.start == exon.start:
                distance = 1
            else:
                distance = exon.distance_to_interval(self.start, self.end)

            if distance == 0:
                overlapping_exon_numbers_and_exons.append((i + 1, exon))
                # start is contained in current exon
                if exon.start <= self.start <= exon.end:
                    start_in_exon = True
                # end is contained in current exon
                if exon.end >= self.end >= exon.start:
                    end_in_exon = True
            elif distance < distance_to_nearest_exon:
                    distance_to_nearest_exon = distance
                    nearest_exon = exon

        if len(overlapping_exon_numbers_and_exons) == 0:
            intronic_effect_class = self._choose_intronic_effect_class(
                nearest_exon, distance_to_nearest_exon)
            return intronic_effect_class(
                variant=self,
                transcript=transcript,
                nearest_exon=nearest_exon,
                distance_to_exon=distance_to_nearest_exon)
        elif len(completely_lost_exons) > 0 or (
                len(overlapping_exon_numbers_and_exons) > 1):
            # if spanning multiple exons, or completely deleted an exon
            # then consider that an ExonLoss mutation
            exons = [exon for (_, exon) in overlapping_exon_numbers_and_exons]
            return ExonLoss(self, transcript, exons)

        assert len(overlapping_exon_numbers_and_exons) == 1

        exon_number, exon = overlapping_exon_numbers_and_exons[0]

        exonic_effect_annotation = self._exonic_transcript_effect(
            exon, exon_number, transcript)

        # simple case: both start and end are in the same
        if start_in_exon and end_in_exon:
            return exonic_effect_annotation
        elif isinstance(exonic_effect_annotation, ExonicSpliceSite):
            # if mutation bleeds over into intro but even just
            # the exonic portion got annotated as an exonic splice site
            # then return it
            return exonic_effect_annotation

        return ExonicSpliceSite(
            variant=self,
            transcript=transcript,
            exon=exon,
            alternate_effect=exonic_effect_annotation)

    def _choose_intronic_effect_class(
            self,
            nearest_exon,
            distance_to_exon):
        """
        Infer effect of variant which does not overlap any exon of
        the given transcript.
        """
        assert distance_to_exon > 0, \
            "Expected intronic effect to have distance_to_exon > 0, got %d" % (
                distance_to_exon,)

        before_forward_exon = (
            nearest_exon.strand == "+" and
            (self.start < nearest_exon.start or
                (self.ref == "" and self.start == nearest_exon.start)))

        before_backward_exon = (
            nearest_exon.strand == "-" and
            (self.end > nearest_exon.end or
                (self.ref == "" and self.end == nearest_exon.end)))

        before_exon = before_forward_exon or before_backward_exon

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

    def _exonic_transcript_effect(self, exon, exon_number, transcript):
        """Effect of this variant on a Transcript, assuming we already know
        that this variant overlaps some exon of the transcript.

        Parameters
        ----------
        exon : pyensembl.Exon
            Exon which this variant overlaps

        exon_number : int
            Index (starting from 1) of the given exon in the transcript's
            sequence of exons.

        transcript : pyensembl.Transcript
        """

        genome_ref = self.ref
        genome_alt = self.alt

        # clip mutation to only affect the current exon
        if self.start < exon.start:
            # if mutation starts before current exon then only look
            # at nucleotides which overlap the exon
            assert len(genome_ref) > 0, "Unexpected insertion into intron"
            n_skip_start = exon.start - self.start
            genome_ref = genome_ref[n_skip_start:]
            genome_alt = genome_alt[n_skip_start:]
            genome_start = exon.start
        else:
            genome_start = self.start

        if self.end > exon.end:
            # if mutation goes past exon end then only look at nucleotides
            # which overlap the exon
            n_skip_end = self.end - exon.end
            genome_ref = genome_ref[:-n_skip_end]
            genome_alt = genome_alt[:len(genome_ref)]
            genome_end = exon.end
        else:
            genome_end = self.end

        transcript_offset = interval_offset_on_transcript(
            genome_start, genome_end, transcript)

        if transcript.on_backward_strand:
            strand_ref = reverse_complement(genome_ref)
            strand_alt = reverse_complement(genome_alt)
        else:
            strand_ref = genome_ref
            strand_alt = genome_alt

        expected_ref = str(transcript.sequence[
            transcript_offset:transcript_offset + len(strand_ref)])

        if strand_ref != expected_ref:
            raise ValueError(
                ("Found ref nucleotides '%s' in sequence"
                 " of %s at offset %d (chromosome positions %d:%d)"
                 " but variant %s has '%s'") % (
                     expected_ref,
                     transcript,
                     transcript_offset,
                     genome_start,
                     genome_end,
                     self,
                     strand_ref))

        utr5_length = min(transcript.start_codon_spliced_offsets)

        # does the variant start inside the 5' UTR?
        if utr5_length > transcript_offset:
            # does the variant end after the 5' UTR, within the coding region?
            if utr5_length < transcript_offset + len(strand_ref):
                return StartLoss(self, transcript)
            else:
                # if variant contained within 5' UTR
                return FivePrimeUTR(self, transcript)

        utr3_offset = max(transcript.stop_codon_spliced_offsets) + 1

        if transcript_offset >= utr3_offset:
            return ThreePrimeUTR(self, transcript)

        exon_start_offset = interval_offset_on_transcript(
            exon.start, exon.end, transcript)
        exon_end_offset = exon_start_offset + len(exon) - 1

        # Further below we're going to try to predict exonic splice site
        # modifications, which will take this effect_annotation as their
        # alternative hypothesis for what happens if splicing doesn't change.
        # If the mutation doesn't affect an exonic splice site, then
        # we'll just return this effect.
        coding_effect_annotation = coding_effect(
            ref=strand_ref,
            alt=strand_alt,
            transcript_offset=transcript_offset,
            variant=self,
            transcript=transcript)

        if changes_exonic_splice_site(
                transcript=transcript,
                transcript_ref=strand_ref,
                transcript_alt=strand_alt,
                transcript_offset=transcript_offset,
                exon_start_offset=exon_start_offset,
                exon_end_offset=exon_end_offset,
                exon_number=exon_number):
            return ExonicSpliceSite(
                variant=self,
                transcript=transcript,
                exon=exon,
                alternate_effect=coding_effect_annotation)
        return coding_effect_annotation

    @property
    def is_insertion(self):
        """
        Does this variant represent the insertion of nucleotides into the
        reference genome?
        """
        # An insertion would appear in a VCF like C>CT, so that the
        # alternate allele starts with the reference nucleotides.
        # Since the nucleotide strings may be normalized in the constructor,
        # it's worth noting that the normalized form of this variant would be
        # ''>'T', so that 'T'.startswith('') still holds.
        return (len(self.ref) < len(self.alt)) and self.alt.startswith(self.ref)

    @property
    def is_deletion(self):
        """
        Does this variant represent the deletion of nucleotides from the
        reference genome?
        """
        # A deletion would appear in a VCF like CT>C, so that the
        # reference allele starts with the alternate nucleotides.
        # This is true even in the normalized case, where the alternate
        # nucleotides are an empty string.
        return (len(self.ref) > len(self.alt)) and self.ref.startswith(self.alt)

    @property
    def is_indel(self):
        """Is this variant either an insertion or deletion?"""
        return self.is_insertion or self.is_deletion

    @property
    def is_snv(self):
        """Is the variant a single nucleotide variant"""
        return (len(self.ref) == len(self.alt) == 1) and (self.ref != self.alt)

    @property
    def is_transition(self):
        """Is this variant and pyrimidine to pyrimidine change or purine to purine change"""
        return self.is_snv and is_purine(self.ref) == is_purine(self.alt)

    @property
    def is_transversion(self):
        """Is this variant a pyrimidine to purine change or vice versa"""
        return self.is_snv and is_purine(self.ref) != is_purine(self.alt)
