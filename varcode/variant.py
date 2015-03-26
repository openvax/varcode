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

from Bio.Seq import reverse_complement
from pyensembl import Transcript, EnsemblRelease
from pyensembl.locus import normalize_chromosome
from pyensembl.biotypes import is_coding_biotype
from typechecks import require_instance

from .coding_effect import coding_effect
from .common import groupby_field, memoize
from .effects import (
    TranscriptMutationEffect,
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
from .effect_ordering import top_priority_effect
from .nucleotides import normalize_nucleotide_string
from .string_helpers import trim_shared_flanking_strings
from .transcript_helpers import interval_offset_on_transcript

DEFAULT_ENSEMBL_RELEASE = EnsemblRelease()

class Variant(object):
    def __init__(self,
            contig,
            start,
            ref,
            alt,
            ensembl=DEFAULT_ENSEMBL_RELEASE,
            info=None,
            allow_extended_nucleotides=False):
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

        ensembl : EnsemblRelease
            Ensembl object used for determining gene/transcript annotations

        info : dict, optional
            Extra metadata about this variant
        """
        self.contig = normalize_chromosome(contig)

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
        self.original_ref = normalize_nucleotide_string(ref,
            allow_extended_nucleotides=allow_extended_nucleotides)
        self.original_alt = normalize_nucleotide_string(alt,
            allow_extended_nucleotides=allow_extended_nucleotides)
        self.original_start = int(start)

        # normalize the variant by trimming any shared prefix or suffix
        # between ref and alt nucleotide sequences and then
        # offset the variant position in a strand-dependent manner
        trimmed_ref, trimmed_alt, prefix, suffix = trim_shared_flanking_strings(
            self.original_ref, self.original_alt)

        self.ref = trimmed_ref
        self.alt = trimmed_alt

        # insertions must be treated differently since the meaning of a
        # position for an insertion is
        #   "insert the alt nucleotides after this position"
        if len(trimmed_ref) == 0:
            self.start = self.original_start + max(0, len(prefix) - 1)
            self.end = self.start
        else:
            # for substitutions and deletions the [start:end] interval is
            # an inclusive selection of reference nucleotides
            self.start = self.original_start + len(prefix)
            self.end = self.start + len(trimmed_ref) - 1

        require_instance(ensembl, EnsemblRelease, "ensembl")
        self.ensembl = ensembl

        self.info = {} if info is None else info

    @property
    def reference_name(self):
        return self.ensembl.reference_name

    def __str__(self):
        return "Variant(contig=%s, start=%d, ref=%s, alt=%s, genome=%s)" % (
            self.contig,
            self.start,
            self.ref if self.ref else ".",
            self.alt if self.alt else ".",
            self.reference_name,)

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash(self.fields())

    def __lt__(self, other):
        '''
        Variants are ordered by locus.
        '''
        require_instance(other, Variant, name="variant")
        if self.contig == other.contig:
            return self.start < other.start
        return self.contig < other.contig

    def fields(self):
        """
        All identifying fields of a variant (contig, pos, ref, alt, genome)
        in a single tuple. This makes for cleaner printing, hashing, and
        comparisons.
        """
        return (
            self.contig,
            self.start,
            self.ref,
            self.alt,
            self.reference_name)

    def __eq__(self, other):
        return (
            isinstance(other, Variant) and
            self.fields() == other.fields())

    @memoize
    def short_description(self):
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

    @memoize
    def transcripts(self):
        return self.ensembl.transcripts_at_locus(
            self.contig, self.start, self.end)

    @memoize
    def coding_transcripts(self):
        """
        Protein coding transcripts
        """
        return [
            transcript for transcript in self.transcripts()
            if is_coding_biotype(transcript.biotype)
        ]

    @memoize
    def genes(self):
        """
        Return Gene object for all genes which overlap this variant.
        """
        return self.ensembl.genes_at_locus(
            self.contig, self.start, self.end)

    @memoize
    def gene_ids(self):
        """
        Return IDs of all genes which overlap this variant. Calling
        this method is significantly cheaper than calling `Variant.genes()`,
        which has to issue many more queries to construct each Gene object.
        """
        return self.ensembl.gene_ids_at_locus(
            self.contig, self.start, self.end)

    @memoize
    def gene_names(self):
        """
        Return names of all genes which overlap this variant. Calling
        this method is significantly cheaper than calling `Variant.genes()`,
        which has to issue many more queries to construct each Gene object.
        """
        return self.ensembl.gene_names_at_locus(
            self.contig, self.start, self.end)

    @memoize
    def coding_genes(self):
        """
        Protein coding transcripts
        """
        return [
            gene for gene in self.genes()
            if is_coding_biotype(gene.biotype)
        ]

    @memoize
    def effects(self, raise_on_error=True):
        """Determine the effects of a variant on any transcripts it overlaps.
        Returns a VariantEffectCollection object.

        Parameters
        ----------
        raise_on_error : bool
            Raise an exception if we encounter an error while trying to
            determine the effect of this variant on a transcript, or simply
            log the error and continue.
        """
        gene_ids = self.gene_ids()

        # if this variant isn't overlapping any genes, return a
        # Intergenic effect
        # TODO: look for nearby genes and mark those as Upstream and Downstream
        # effects
        if len(gene_ids) == 0:
            return [Intergenic(self)]

        overlapping_transcripts = self.transcripts()
        # group transcripts by their gene ID
        transcripts_grouped_by_gene = groupby_field(
            overlapping_transcripts, 'gene_id')

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
        return effects

    @memoize
    def transcript_effect_dict(self, *args, **kwargs):
        """Dictionary mapping transcript IDs to their associated
        TranscriptMutationEffect objects.

        Arguments are passed on to Variant.effects(*args, **kwargs).
        """
        return {
            effect.transcript.id: effect
            for effect in self.effects(*args, **kwargs)
            if isinstance(effect, TranscriptMutationEffect)
        }

    @memoize
    def top_effect(self, *args, **kwargs):
        """Highest priority MutationEffect of all genes/transcripts overlapped
        by this variant. If this variant doesn't overlap anything, then this
        this method will return an Intergenic effect.
        """
        return top_priority_effect(self.effects(*args, **kwargs))

    @memoize
    def effect_on_transcript(self, transcript):
        """Return the transcript effect (such as FrameShift) that results from
        applying this genomic variant to a particular transcript.

        Parameters
        ----------
        transcript :  Transcript
            Transcript we're going to apply mutation to.
        """

        if not isinstance(transcript, Transcript):
            raise TypeError(
                "Expected %s : %s to have type Transcript" % (
                    transcript, type(transcript)))

        # check for non-coding transcripts first, since
        # every non-coding transcript is "incomplete".
        if not is_coding_biotype(transcript.biotype):
            return NoncodingTranscript(self, transcript)

        if not transcript.complete:
            return IncompleteTranscript(self, transcript)

        # determine if any exons are deleted, and if not,
        # what is the closest exon and how far is this variant
        # from that exon (overlapping the exon = 0 distance)
        lost_exons = []
        overlapping_exons = set([])
        distance_to_nearest_exon = float("inf")
        start_in_exon = False
        end_in_exon = False
        nearest_exon = None
        for exon in transcript.exons:
            if self.start <= exon.start and self.end >= exon.end:
                lost_exons.append(exon)

            distance = exon.distance_to_interval(self.start, self.end)
            if distance == 0:
                overlapping_exons.add(exon)
                # start is contained in current exon
                if exon.start <= self.start <= exon.end:
                    start_in_exon = True
                # end is contained in current exon
                if exon.end >= self.end >= exon.start:
                    end_in_exon = True
            elif distance < distance_to_nearest_exon:
                    distance_to_nearest_exon = distance
                    nearest_exon = exon

        if len(lost_exons) > 0:
            return ExonLoss(self, transcript, lost_exons)
        elif len(overlapping_exons) == 0:
            intronic_effect_class = self._choose_intronic_effect_class(
                nearest_exon, distance_to_nearest_exon)
            return intronic_effect_class(
                variant=self,
                transcript=transcript,
                nearest_exon=nearest_exon,
                distance_to_exon=distance_to_nearest_exon)

        # simple case: both start and end are in the same
        elif len(overlapping_exons) == 1 and start_in_exon and end_in_exon:
            return self._exonic_transcript_effect(transcript)
        # if spanning multiple exons, or only part of the variant is inside
        # an exon, then consider than an exonic splice site mutation
        else:
            return ExonicSpliceSite(self, transcript)

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
            self.start < nearest_exon.start)

        before_backward_exon = (
            nearest_exon.strand == "-" and
            self.end > nearest_exon.end)

        before_exon = before_forward_exon or before_backward_exon

        if distance_to_exon <= 2:
            if before_exon:
                # 2 last nucleotides of intron before exon are the splice acceptor
                # site, typically "AG"
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
        elif before_exon and distance_to_exon <= 4:
            # nucleotides -4 and -3 before exon are part of the 3' splicing motif
            # but allow for more degeneracy than the -2, -1 nucleotides
            return IntronicSpliceSite
        else:
            # intronic mutation unrelated to splicing
            return Intronic

    def _exonic_transcript_effect(self, transcript):
        """
        Effect of this variant on a Transcript, assuming we already know
        that this variant overlaps some exon of the transcript.
        """
        if transcript.on_backward_strand:
            strand_ref = reverse_complement(self.ref)
            strand_alt = reverse_complement(self.alt)
        else:
            strand_ref = self.ref
            strand_alt = self.alt

        offset_with_utr5 = interval_offset_on_transcript(
            self.start, self.end, transcript)

        utr5_length = min(transcript.start_codon_spliced_offsets)

        # does the variant start inside the 5' UTR?
        if utr5_length > offset_with_utr5:

            # does the variant end after the 5' UTR, within the coding region?
            if utr5_length < offset_with_utr5 + len(self.ref):
                return StartLoss(self, transcript)
            else:
                # if variant contained within 5' UTR
                return FivePrimeUTR(self, transcript)

        utr3_offset = max(transcript.stop_codon_spliced_offsets) + 1

        if offset_with_utr5 >= utr3_offset:
            return ThreePrimeUTR(self, transcript)

        return coding_effect(
            strand_ref,
            strand_alt,
            offset_with_utr5,
            variant=self,
            transcript=transcript)
