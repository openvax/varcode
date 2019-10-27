# Copyright (c) 2016-2019. Mount Sinai School of Medicine
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

from pyensembl.locus import normalize_chromosome
from serializable import Serializable
from typechecks import require_instance

from .nucleotides import (
    normalize_nucleotide_string,
    STANDARD_NUCLEOTIDES,
    is_purine
)
from .reference import  infer_genome, ensembl_to_ucsc_reference_names
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
        "original_genome",
        "original_reference_name",
        "genome",
        "normalize_contig_names",
        "original_genome_was_ucsc",
        "convert_ucsc_contig_names",
        "allow_extended_nucleotides",
        "original_contig",
        "original_ref",
        "original_alt",
        "original_start",
        "_transcripts",
        "_genes",
        "_gene_ids",
        "_gene_names",
    )

    # any keywords which are renamed or removed should be added to this
    # dictionary to preserve the ability to deserialize old representations
    # of Variant objects
    _SERIALIZABLE_KEYWORD_ALIASES = {
        "normalize_contig_name": "normalize_contig_names"
    }
    
    def __init__(
            self,
            contig,
            start,
            ref,
            alt,
            genome=None,
            ensembl=None,
            allow_extended_nucleotides=False,
            normalize_contig_names=True,
            convert_ucsc_contig_names=None):
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

        genome : Genome, EnsemblRelease, or str, or int
            Name of reference genome, Ensembl release number, or object
            derived from pyensembl.Genome. Default to latest available release
            of GRCh38

        ensembl : Genome, EnsemblRelease, or str, or int (DEPRECATED)
            Previous name used instead of 'genome', the two arguments should
            be mutually exclusive.

        allow_extended_nucleotides : bool
            Extended nucleotides include 'Y' for pyrimidies or 'N' for any base

        normalize_contig_names : bool
            By default the contig name will be normalized by converting integers
            to strings (e.g. 1 -> "1"), and converting any letters after "chr"
            to uppercase (e.g. "chrx" -> "chrX"). If you don't want
            this behavior then pass normalize_contig_name=False.

        convert_ucsc_contig_names : bool, optional
            Setting this argument to True causes UCSC chromosome names to be
            coverted, such as "chr1" to "1". If the default value (None) is used
            then it defaults to whether or not a UCSC genome was pass in for
            the 'genome' argument.
        """

        # first initialize the fields we use to cache lists of overlapping
        # pyensembl Gene and Transcript objects, or their properties such
        # as names/IDs
        self._genes = None
        self._transcripts = None
        self._gene_ids = None
        self._gene_names = None

        # store the options which affect how properties of this variant
        # may be changed/transformed
        self.normalize_contig_names = normalize_contig_names
        self.allow_extended_nucleotides = allow_extended_nucleotides

        # if genome not specified, try the old name 'ensembl'
        # if ensembl is also None, then default to "GRCh38"
        if genome is None and ensembl is None:
            genome = "GRCh38"
        elif genome is None:
            genome = ensembl


        # user might supply Ensembl release as an integer, reference name,
        # or pyensembl.Genome object
        self.original_genome = genome
        self.genome, self.original_genome_was_ucsc = infer_genome(genome)

        self.reference_name = self.genome.reference_name
        if self.original_genome_was_ucsc:
            self.original_reference_name = ensembl_to_ucsc_reference_names[
                self.reference_name]
        else:
            self.original_reference_name = self.reference_name

        self.original_contig = contig
        self.contig = normalize_chromosome(contig) if normalize_contig_names else contig

        if convert_ucsc_contig_names is None:
            self.convert_ucsc_contig_names = self.original_genome_was_ucsc
        else:
            self.convert_ucsc_contig_names = convert_ucsc_contig_names

        # trim off the starting "chr" from hg19 chromosome names to make them
        # match GRCh37, also convert "chrM" to "MT".
        if self.convert_ucsc_contig_names:
            self.contig = self._convert_ucsc_contig_name_to_ensembl(self.contig)

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
        (trimmed_ref, trimmed_alt, prefix, _) = \
            trim_shared_flanking_strings(self.original_ref, self.original_alt)

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
    def ensembl(self):
        """
        Deprecated alias for Variant.genome

        Returns
        -------
        pyensembl.Genome
        """
        return self.genome

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
            self.genome == other.genome)

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
            genome=self.original_genome,
            allow_extended_nucleotides=self.allow_extended_nucleotides,
            normalize_contig_names=self.normalize_contig_names,
            convert_ucsc_contig_names=self.convert_ucsc_contig_names)

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


    # dictionary from reference names to sets of valid contigs
    _reference_name_to_valid_contig_names = {}

    def _check_that_genome_has_contig(self):
        """
        Any annotation code which wants to look up genomic loci in a reference
        genome should first run this helper to make sure the normalized contig
        name is in the genome.
        """
        reference_name = self.reference_name

        if reference_name not in self._reference_name_to_valid_contig_names:
            self._reference_name_to_valid_contig_names[reference_name] = set(
                self.genome.contigs())

        valid_contigs =  \
            self._reference_name_to_valid_contig_names[reference_name]

        if self.contig not in valid_contigs:
            raise ValueError("Invalid contig name '%s' for reference '%s'" % (
                self.contig,
                self.reference_name))



    @property
    def transcripts(self):
        if self._transcripts is None:
            self._check_that_genome_has_contig()
            self._transcripts = self.genome.transcripts_at_locus(
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
            self._check_that_genome_has_contig()
            self._genes = self.genome.genes_at_locus(
                self.contig, self.start, self.end)
        return self._genes

    @property
    def gene_ids(self):
        """
        Return IDs of all genes which overlap this variant. Calling
        this method is significantly cheaper than calling `Variant.genes()`,
        which has to issue many more queries to construct each Gene object.
        """
        if self._gene_ids is None:
            self._check_that_genome_has_contig()
            self._gene_ids = self.genome.gene_ids_at_locus(
                self.contig, self.start, self.end)
        return self._gene_ids

    @property
    def gene_names(self):
        """
        Return names of all genes which overlap this variant. Calling
        this method is significantly cheaper than calling `Variant.genes()`,
        which has to issue many more queries to construct each Gene object.
        """
        if self._gene_names is None:
            self._check_that_genome_has_contig()
            self._gene_names = self.genome.gene_names_at_locus(
                self.contig, self.start, self.end)
        return self._gene_names

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

    def _convert_ucsc_contig_name_to_ensembl(self, contig):
        """
        Convert names such as "chr1" to "1" but don't convert the
        unplaced contigs such as "chrUn_KI270438v1" or alt contigs
        such as "chr6_apd_hap1".

        Parameters
        ----------
        contig : str

        Returns
        -------
        str
        """
        if contig.startswith("chr") and "_" not in contig:
            contig = self.contig[3:]
        if contig == "M":
            contig = "MT"
        return contig

    def clone_without_ucsc_data(self):
        """
        Clone this variant but discarding the original format of its genome
        and contig: useful if we want to mix hg19 and GRCh37 variants.

        Returns
        -------
        Variant
        """
        return Variant(
            contig=self.contig,
            start=self.original_start,
            ref=self.original_ref,
            alt=self.original_alt,
            genome=self.genome,
            allow_extended_nucleotides=self.allow_extended_nucleotides,
            normalize_contig_names=self.normalize_contig_names,
            convert_ucsc_contig_names=False)

def variant_ascending_position_sort_key(variant):
    """
    Sort key function used to sort variants in ascending order by
    chromosomal position.

    This is a function rather than method to make it easier to serialize.

    Returns
    -------
    tuple of (str, int)
    """
    return (variant.contig, variant.start)
