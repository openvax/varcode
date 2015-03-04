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
from pyensembl import Transcript, find_nearest_locus, EnsemblRelease
from pyensembl.locus import normalize_chromosome
from pyensembl.biotypes import is_coding_biotype
from typechecks import require_instance

from .coding_effect import infer_coding_effect
from .common import group_by
from .nucleotides import normalize_nucleotide_string
from .string_helpers import trim_shared_flanking_strings
from .transcript_helpers import interval_offset_on_transcript
from .effects import (
    NoncodingTranscript,
    IncompleteTranscript,
    FivePrimeUTR,
    ThreePrimeUTR,
    Intronic,
    IntronicSpliceSite,
    SpliceAcceptor,
    SpliceDonor,
)
from .variant_effect_collection import VariantEffectCollection


class Variant(object):
    def __init__(self,
            contig,
            pos,
            ref,
            alt,
            ensembl,
            info=None,
            allow_extended_nucleotides=False):
        """
        Construct a Variant object.

        Parameters
        ----------
        contig : str
            Chromosome that this variant is on

        pos : int
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
        self.ref = normalize_nucleotide_string(ref,
            allow_extended_nucleotides=allow_extended_nucleotides)
        self.alt = normalize_nucleotide_string(alt,
            allow_extended_nucleotides=allow_extended_nucleotides)
        self.pos = int(pos)
        self.end = self.pos + len(self.ref) - 1

        require_instance(ensembl, EnsemblRelease, "ensembl")
        self.ensembl = ensembl

        self.info = {} if info is None else info

    @property
    def reference_name(self):
        return self.ensembl.reference.reference_name

    def __str__(self):
        fields = self.fields()
        return "Variant(contig=%s, pos=%d, ref=%s, alt=%s, genome=%s)" % fields

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
            return self.pos < other.pos
        return self.contig < other.contig

    def fields(self):
        """
        All identifying fields of a variant (contig, pos, ref, alt, genome)
        in a single tuple. This makes for cleaner printing, hashing, and
        comparisons.
        """
        return (
            self.contig,
            self.pos,
            self.ref,
            self.alt,
            self.reference_name)

    def __eq__(self, other):
        return (
            isinstance(other, Variant) and
            self.fields() == other.fields())

    def short_description(self):
        chrom, pos, ref, alt = self.contig, self.pos, self.ref, self.alt
        if ref == alt:
            return "chr%s g.%d %s=%s" % (chrom, pos, ref, alt)
        elif len(ref) == 0 or alt.startswith(ref):
            return "chr%s g.%d ins%s" % (chrom, pos + len(ref), alt[len(ref):])
        elif len(alt) == 0 or ref.startswith(alt):
            return "chr%s g.%d_%d del%s" % (
                chrom, pos + len(alt), pos + len(ref), ref[len(alt):])
        else:
            return "chr%s g.%d %s>%s" % (chrom, pos, ref, alt)

    def transcripts(self):
        return self.ensembl.transcripts_at_locus(
            self.contig, self.pos, self.end)

    def coding_transcripts(self):
        """
        Protein coding transcripts
        """
        return [
            transcript for transcript in self.transcripts()
            if is_coding_biotype(transcript.biotype)
        ]

    def genes(self):
        """
        Return Gene object for all genes which overlap this variant.
        """
        return self.ensembl.genes_at_locus(
            self.contig, self.pos, self.end)

    def gene_names(self):
        """
        Return names of all genes which overlap this variant. Calling
        this method is significantly cheaper than calling `Variant.genes()`,
        which has to issue many more queries to construct each Gene object.
        """
        return self.ensembl.gene_names_at_locus(
            self.contig, self.pos, self.end)

    def coding_genes(self):
        """
        Protein coding transcripts
        """
        return [
            gene for gene in self.genes()
            if is_coding_biotype(gene.biotype)
        ]

    def effects(self, only_coding_transcripts=False, raise_on_error=True):
        """
        Determine the effects of a variant on any transcripts it overlaps.
        Returns a VariantEffectCollection object.
        """

        overlapping_transcripts = self.transcripts()

        if only_coding_transcripts:
            # remove non-coding transcripts
            overlapping_transcripts = [
                transcript
                for transcript
                in overlapping_transcripts
                if is_coding_biotype(transcript.biotype)
            ]

        if len(overlapping_transcripts) == 0:
            # intergenic variant
            return VariantEffectCollection(
                variant=self,
                gene_effect_groups={},
                errors={})

        # group transcripts by their gene ID
        overlapping_transcript_groups = group_by(
            overlapping_transcripts, field_name='gene_id')

        # dictionary from gene ID to list of transcript effects
        gene_transcript_effects_groups = {}

        # mapping from Transcript objects to errors encountered
        # while trying to annotated them
        errors = {}

        for (gene_id, transcripts) in overlapping_transcript_groups.items():
            effects = []
            for transcript in transcripts:
                try:
                    effects.append(self.transcript_effect(transcript))
                except (AssertionError, ValueError) as error:
                    if raise_on_error:
                        raise
                    else:
                        logging.warn(
                            "Encountered error annotating %s for %s: %s",
                            self,
                            transcript,
                            error)
                    errors[transcript] = error
            gene_transcript_effects_groups[gene_id] = effects

        return VariantEffectCollection(
            variant=self,
            gene_effect_groups=gene_transcript_effects_groups,
            errors=errors)

    def transcript_effect(self, transcript):
        """
        Return the transcript effect (such as FrameShift) that results from
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

        distance_to_exon, nearest_exon = find_nearest_locus(
            start=self.pos,
            end=self.end,
            loci=transcript.exons)

        if distance_to_exon > 0:
            return self._intronic_transcript_effect(
                transcript=transcript,
                nearest_exon=nearest_exon)

        # TODO: exonic splice site mutations
        return self._exonic_transcript_effect(transcript)

    def _intronic_transcript_effect(self, transcript, nearest_exon):
        """
        Infer effect of variant which does not overlap any exon of
        the given transcript.
        """
        distance_to_exon = nearest_exon.distance_to_interval(
            self.pos, self.end)

        assert distance_to_exon > 0, \
            "Expected intronic effect to have distance_to_exon > 0, got %d" % (
                distance_to_exon,)

        before_forward_exon = (
            transcript.strand == "+" and
            self.pos < nearest_exon.start)

        before_backward_exon = (
            transcript.strand == "-" and
            self.end > nearest_exon.end)

        before_exon = before_forward_exon or before_backward_exon

        if distance_to_exon <= 2:
            if before_exon:
                # 2 last nucleotides of intron before exon are the splice acceptor
                # site, typically "AG"
                effect_class = SpliceAcceptor
            else:
                # 2 first nucleotides of intron after exon are the splice donor
                # site, typically "GT"
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
                variant=self,
                transcript=transcript,
                nearest_exon=nearest_exon,
                distance_to_exon=distance_to_exon)

    def _offset_on_transcript(self, transcript):
        """
        Given a Variant (with fields ref, alt, pos, and end), compute the offset
        of the variant position relative to the position and direction of the
        transcript.

        Returns a tuple (offset, ref, alt), where ref and alt are this variant's
        ref and alt nucleotide strings with shared prefixes and suffixes removed
        and the bases reversed if the transcript is on the backward strand.

        For example,

        >> variant = Variant(ref="ATGGC", alt="TCGC", pos=3000, end=3005)
        >> variant.offset_on_transcript(Transcript(start=2000, end=4000))
        (1002, "G", "C")
        """
        if transcript.on_backward_strand:
            ref = reverse_complement(self.ref)
            alt = reverse_complement(self.alt)
        else:
            ref = self.ref
            alt = self.alt

        # in case nucleotide strings share prefix (e.g. ref="C", alt="CC")
        # bump the offsets and make the strings disjoint (ref="", alt="C")
        ref, alt, prefix, _ = trim_shared_flanking_strings(ref, alt)
        n_same = len(prefix)
        start_offset_with_utr5 = interval_offset_on_transcript(
            self.pos, self.end, transcript)
        offset = start_offset_with_utr5 + n_same
        return offset, ref, alt

    def _exonic_transcript_effect(self, transcript):
        """
        Effect of this variant on a Transcript, assuming we already know
        that this variant overlaps some exon of the transcript.
        """
        offset_with_utr5, ref, alt = self._offset_on_transcript(transcript)

        utr5_length = min(transcript.start_codon_spliced_offsets)

        if utr5_length > offset_with_utr5:
            # TODO: what do we do if the variant spans the beginning of
            # the coding sequence?
            if utr5_length < offset_with_utr5 + len(ref):
                raise ValueError(
                    "Variant which span 5' UTR and CDS not supported: %s" % (
                        self,))
            return FivePrimeUTR(self, transcript)

        utr3_offset = max(transcript.stop_codon_spliced_offsets) + 1

        if offset_with_utr5 >= utr3_offset:
            return ThreePrimeUTR(self, transcript)

        cds_offset = offset_with_utr5 - utr5_length

        return infer_coding_effect(
            ref,
            alt,
            cds_offset,
            variant=self,
            transcript=transcript)
