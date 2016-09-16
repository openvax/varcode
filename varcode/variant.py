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

from pyensembl import (
    cached_release,
    genome_for_reference_name,
    Genome,
    ensembl_grch38,
)
from pyensembl.locus import normalize_chromosome
from serializable import Serializable
from typechecks import require_instance

from .nucleotides import (
    normalize_nucleotide_string,
    STANDARD_NUCLEOTIDES,
    is_purine
)
from .string_helpers import trim_shared_flanking_strings
from .effects import (
    predict_variant_effects,
    predict_variant_effect_on_transcript
)

class Variant(Serializable):
    __slots__ = (
        "contig",
        "start",
        "end",
        "ref",
        "alt",
        "ensembl",
        "normalize_contig_name",
        "allow_extended_nucleotides",
        "original_contig",
        "original_ref",
        "original_alt",
        "original_start",
        "_transcripts",
        "_genes",
    )

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

        # first initialize the _genes and _transcripts fields we use to cache
        # lists of overlapping pyensembl Gene and Transcript objects
        self._genes = self._transcripts = None

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

        self.normalize_contig_name = normalize_contig_name
        self.allow_extended_nucleotides = allow_extended_nucleotides
        self.original_contig = contig
        self.contig = normalize_chromosome(contig) if normalize_contig_name else contig

        if ref != alt and ref in STANDARD_NUCLEOTIDES and alt in STANDARD_NUCLEOTIDES:
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

    def to_dict(self):
        """
        We want the original values (un-normalized) field values while
        serializing since normalization will happen in __init__.
        """
        return dict(
            contig=self.original_contig,
            start=self.original_start,
            ref=self.original_ref,
            alt=self.original_alt,
            ensembl=self.ensembl,
            allow_extended_nucleotides=self.allow_extended_nucleotides,
            normalize_contig_name=self.normalize_contig_name)

    @property
    def trimmed_ref(self):
        """
        Eventually the field Variant.ref will store the reference nucleotides
        as given in a VCF or MAF and trimming of any shared prefix/suffix
        between ref and alt will be done via the properties trimmed_ref
        and trimmed_alt.
        """
        return self.ref

    @property
    def trimmed_alt(self):
        """
        Eventually the field Variant.ref will store the reference nucleotides
        as given in a VCF or MAF and trimming of any shared prefix/suffix
        between ref and alt will be done via the properties trimmed_ref
        and trimmed_alt.
        """
        return self.alt

    @property
    def trimmed_base1_start(self):
        """
        Currently the field Variant.start carries the base-1 starting position
        adjusted by trimming any shared prefix between Variant.ref and
        Variant.alt. Eventually this trimming should be done more explicitly
        via trimmed_* properties.
        """
        return self.start

    @property
    def trimmed_base1_end(self):
        """
        Currently the field Variant.end carries the base-1 "last" position of
        this variant, adjusted by trimming any shared suffix between
        Variant.ref and Variant.alt. Eventually this trimming should be done
        more explicitly via trimmed_* properties.
        """
        return self.end

    @property
    def short_description(self):
        """
        HGVS nomenclature for genomic variants
        More info: http://www.hgvs.org/mutnomen/
        """
        if self.is_insertion:
            return "chr%s g.%d_%dins%s" % (
                self.contig,
                self.start,
                self.start + 1,
                self.alt)
        elif self.is_deletion:
            return "chr%s g.%d_%ddel%s" % (
                self.contig,
                self.start,
                self.end,
                self.ref)
        elif self.ref == self.alt:
            return "chr%s g.%d%s" % (self.contig, self.start, self.ref)
        else:
            # substitution
            return "chr%s g.%d%s>%s" % (
                self.contig,
                self.start,
                self.ref,
                self.alt)

    @property
    def transcripts(self):
        if self._transcripts is None:
            self._transcripts = self.ensembl.transcripts_at_locus(
                self.contig, self.start, self.end)
        return self._transcripts

    @property
    def coding_transcripts(self):
        """
        Protein coding transcripts
        """
        return [
            transcript
            for transcript in self.transcripts
            if transcript.is_protein_coding
        ]

    @property
    def transcript_ids(self):
        return [transcript.id for transcript in self.transcripts]

    @property
    def transcript_names(self):
        return [transcript.name for transcript in self.transcripts]

    @property
    def genes(self):
        """
        Return Gene object for all genes which overlap this variant.
        """
        if self._genes is None:
            self._genes = self.ensembl.genes_at_locus(
                self.contig, self.start, self.end)
        return self._genes

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
            if gene.is_protein_coding
        ]

    def effects(self, raise_on_error=True):
        return predict_variant_effects(
            variant=self, raise_on_error=raise_on_error)

    def effect_on_transcript(self, transcript):
        return predict_variant_effect_on_transcript(self, transcript)

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


def variant_ascending_position_sort_key(variant):
    """
    Sort key function used to sort variants in ascending order by
    chromosomal position.
    """
    return (variant.contig, variant.start)
