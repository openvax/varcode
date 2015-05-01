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

from __future__ import absolute_import

from collections import namedtuple, defaultdict, OrderedDict
import logging
import itertools
import re

import pandas
import pysam
import typechecks
import pyensembl
import typechecks

from .. import Locus
from . import Pileup, PileupElement, alignment_key, read_key

MatchingEvidence = namedtuple("MatchingEvidence", "ref alt other")

class PileupCollection(object):
    '''
    A collection of Pileup instances for some loci.

    Attributes
    ----------
    pileups (optional) : dict of Locus -> Pileup
        A map from 1-base Locus instances to the Pileup at that locus.

    parent (optional) : PileupCollection
        For PileupCollection instances that are derived (e.g. by filtering)
        from another instance, this attribute tracks the parent instance.
    '''
    def __init__(self, pileups=None, parent=None):
        '''
        Construct a new PileupCollection.
        '''
        if pileups is None:
            pileups = {}
        self.pileups = pileups
        self.parent = parent

    def pileup(self, locus):
        '''
        Given a 1-base locus, return the Pileup at that locus.

        Raises a KeyError if this PileupCollection does not have a Pileup at
        the specified locus.
        '''
        if len(locus.positions) != 1:
            raise ValueError("Not a single-base locus: %s" % locus)
        return self.pileups[locus]

    def at(self, *loci):
        '''
        Return a new PileupCollection instance including only pileups for 
        the specified loci.
        '''
        single_position_loci = []
        for locus in loci:
            for position in locus.positions:
                single_position_loci.append(
                    Locus.from_interbase_coordinates(locus.contig, position))
        pileups = dict(
            (locus, self.pileups[locus]) for locus in single_position_loci)
        return PileupCollection(pileups, self)

    def loci(self):
        '''
        Returns the loci included in this instance.
        '''
        return list(self.pileups)

    def reads(self):
        '''
        The reads in this PileupCollection. All reads will have an alignment
        that overlaps at least one of the included loci.

        Since SAM (and pysam) have no real notion of a "read", the returned
        instances are actually pysam.AlignedSegment instances, (i.e. 
        alignments). However, only one alignment will be returned by this
        method per read.
        
        Returns
        ----------
        List of pysam.AlignedSegment instances. If a particular read has more
        than one alignment in this PileupCollection (e.g. one primary and one
        secondary), then the alignment returned is the one with the highest
        mapping quality.
        '''
        # TODO: Optimize this.

        def alignment_precedence(pysam_alignment_record):
            return pysam_alignment_record.mapping_quality

        result = {}
        for pileup in self.pileups.values():
            for e in pileup.elements:
                key = read_key(e.alignment)
                if key not in result or (
                        alignment_precedence(e.alignment) >
                        alignment_precedence(result[key])):
                    result[key] = e.alignment
        return list(result.values())

    def num_reads(self):
        '''
        Returns the number of reads in this PileupCollection.
        '''
        return len(self.reads())

    def read_attribute(self, attribute):
        '''
        Query a read attribute across all reads in this PileupCollection.

        Parameters
        ----------
        attribute : string
            Attribute to query. See `PileupCollection.read_attributes` for
            possible attributes.

        Returns
        ----------
        pandas.Series instance giving values for the attribute in each read.
        The order of reads is fixed, so multiple calls to this function will
        return corresponding values. 
        '''
        return self.read_attributes([attribute])[attribute]

    def read_attributes(self, attributes=None):
        '''
        Collect read attributes across reads in this PileupCollection into a
        pandas.DataFrame.

        Valid attributes are the following properties of a pysam.AlignedSegment
        instance. See:
            http://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedSegment
        for the meaning of these attributes.

             * cigarstring
             * flag
             * inferred_length
             * is_duplicate
             * is_paired
             * is_proper_pair
             * is_qcfail
             * is_read1
             * is_read2
             * is_reverse
             * is_secondary
             * is_unmapped
             * mapping_quality
             * mate_is_reverse
             * mate_is_unmapped
             * next_reference_id
             * next_reference_start
             * query_alignment_end
             * query_alignment_length
             * query_alignment_qualities
             * query_alignment_sequence
             * query_alignment_start
             * query_length
             * query_name
             * reference_end
             * reference_id
             * reference_length
             * reference_start
             * template_length

        (Note: the above list is parsed into the _READ_ATTRIBUTE_NAMES class
        variable, so be careful when modifying it.)

        Additionally, for alignment "tags" (arbitrary key values associated
        with an alignment), a column of the form "TAG_{tag name}" is
        included.

        Finally, the column "pysam_alignment_record" gives the underlying
        `pysam.AlignedSegment` instances.

        Parameters
        ----------
        attributes (optional): list of strings
            List of columns to include. If unspecified, all columns are
            included in the result.

        Returns
        ----------
        pandas.DataFrame of read attributes.

        '''
        def include(attribute):
            return attributes is None or attribute in attributes

        reads = self.reads()

        possible_column_names = list(PileupCollection._READ_ATTRIBUTE_NAMES)
        result = OrderedDict(
            (name, [getattr(read, name) for read in reads])
            for name in PileupCollection._READ_ATTRIBUTE_NAMES
            if include(name))

        # Add tag columns.
        if reads:
            tag_dicts = [dict(x.get_tags()) for x in reads]
            tag_keys = set.union(
                *[set(item.keys()) for item in tag_dicts])
            for tag_key in sorted(tag_keys):
                column_name = "TAG_%s" % tag_key
                possible_column_names.append(column_name)
                if include(column_name):
                    result[column_name] = [d.get(tag_key) for d in tag_dicts]

        # Lastly, we include the underlying pysam alignment record.
        possible_column_names.append("pysam_alignment_record")
        if include("pysam_alignment_record"):
            result["pysam_alignment_record"] = reads

        # If particular attributes were requested, check that they're here.
        if attributes is not None:
            for attribute in attributes:
                if attribute not in result:
                    raise ValueError(
                        "No such attribute: %s. Valid attributes are: %s"
                        % (attribute, " ".join(possible_column_names)))
            assert set(attributes) == set(result)

        return pandas.DataFrame(result)

    # Parse the bulleted list of attributes in the previous method's docstring.
    _READ_ATTRIBUTE_NAMES = []
    for line in read_attributes.__doc__.split("\n"):
        match = re.match("\s+\*\s*(\w+)\s*", line)
        if match:
            _READ_ATTRIBUTE_NAMES.append(match.groups()[0])

    def group_by_allele(self, locus):
        '''
        Split the PileupCollection by the alleles suggested by the reads at the
        specified locus.

        If a read has an insertion immediately following the locus, then the
        insertion is included in the allele. For example, if locus is the
        1-base range [5,6), one allele might be "AGA", indicating that at
        locus 5 some read has an "A" followed by a 2-base insertion ("GA"). If
        a read has a deletion at the specified locus, the allele is the empty
        string.

        The given locus may include any number of bases. If the locus includes
        multiple bases, then the alleles consist of all bases aligning to that
        range in any read. Note that only sequences actually sequenced in a
        particular read are included. For example, if one read has "ATT" at a
        locus and another read has "GCC", then the alleles are "ATT" and
        "GCC", but not "GTT". That is, the bases in each allele are phased. For
        this reason, only reads that overlap the entire locus are included.

        If the locus is an empty interval (e.g. [5,5) ), then the alleles
        consist only of inserted bases. In this example, only bases inserted
        immediately after locus 5 would be included (but *not* the base
        actually at locus 5). In the previous insertion example, the allele
        would be "GA", indicating a 2-base insertion. Reads that have no
        insertion at that position (matches or deletions) would have the empty
        string as their allele.

        Parameters
        ----------
        locus : Locus
            The reference locus, encompassing 0 or more bases.

        Returns
        ----------
        A dict of string -> PileupCollection. The keys are nucleotide strings
        giving the bases sequenced at the locus, and the values are
        PileupCollection instances of the alignments that support that allele.

        '''
        read_to_allele = None
        loci = []
        if locus.positions:
            # Our locus includes at least one reference base.
            for position in locus.positions:
                base_position = Locus.from_interbase_coordinates(
                    locus.contig, position)
                loci.append(base_position)
                new_read_to_allele = {}
                for element in self.pileups[base_position]:
                    allele_prefix = ""
                    key = alignment_key(element.alignment)
                    if read_to_allele is not None:
                        try:
                            allele_prefix = read_to_allele[key]
                        except KeyError:
                            continue
                    allele = allele_prefix + element.bases
                    new_read_to_allele[key] = allele
                read_to_allele = new_read_to_allele
        else:
            # Our locus is between reference bases.
            position_before = Locus.from_interbase_coordinates(
                locus.contig, locus.start)
            loci.append(position_before)
            read_to_allele = {}
            for element in self.pileups[position_before]:
                allele = element.bases[1:]
                read_to_allele[alignment_key(element.alignment)] = allele

        split = defaultdict(lambda: PileupCollection(pileups={}, parent=self))
        for locus in loci:
            pileup = self.pileups[locus]
            for e in pileup.elements:
                key = read_to_allele.get(alignment_key(e.alignment))
                if key is not None:
                    if locus in split[key].pileups:
                        split[key].pileups[locus].append(e)
                    else:
                        split[key].pileups[locus] = Pileup(locus, [e])

        # Sort by number of reads (descending). Break ties with the
        # lexicographic ordering of the allele string.
        def sorter(pair):
            (allele, pileup_collection) = pair
            return (-1 * pileup_collection.num_reads(), allele)
        return OrderedDict(sorted(split.items(), key=sorter))

    def allele_summary(self, locus, score=lambda x: x.num_reads()):
        '''
        Convenience method to summarize the evidence for each of the alleles
        present at a locus. Applies a score function to the PileupCollection
        associated with each allele.

        See also `PileupCollection.group_by_allele`.

        Parameters
        ----------
        locus : Locus
            The reference locus, encompassing 0 or more bases.

        score (optional) : PileupCollection -> object
            Function to apply to summarize the evidence for each allele.
            Default: count number of reads.

        Returns
        ----------
        List of (allele, score) pairs.
        '''
        return [
            (allele, score(x))
            for (allele, x) in self.group_by_allele(locus).items()
        ]

    def group_by_match(self, variant):
        '''
        Given a variant, split the PileupCollection based on whether it the
        data supports the reference allele, the alternate allele, or neither.

        Parameters
        ----------
        variant : Variant
            The variant. Must have fields 'locus', 'ref', and 'alt'.

        Returns
        ----------
        A MatchingEvidence named tuple with fields (ref, alt, other),
        each of which is a string -> PileupCollection dict mapping alleles
        to the PileupCollection of evidence supporting them.
        '''
        if len(variant.ref) != len(variant.locus.positions):
            logging.warning(
                "Ref is length %d but locus has %d bases in variant: %s" %
                (len(variant.ref), len(variant.locus.positions), str(variant)))

        alleles_dict = self.group_by_allele(variant.locus)
        single_base_loci = [
            Locus.from_interbase_coordinates(variant.locus.contig, position)
            for position in variant.locus.positions
        ]
        empty_pileups = dict(
            (locus, Pileup(locus=locus, elements=[]))
            for locus in single_base_loci)
        empty_collection = PileupCollection(pileups=empty_pileups, parent=self)

        ref = {variant.ref: alleles_dict.pop(variant.ref, empty_collection)}
        alt = {variant.alt: alleles_dict.pop(variant.alt, empty_collection)}
        other = alleles_dict

        # TODO: consider end of read issues for insertions
        return MatchingEvidence(ref, alt, other)

    def match_summary(self, variant, score=lambda x: x.num_reads()):
        '''
        Convenience method to summarize the evidence for and against a variant
        using a user-specified score function.

        See also `PileupCollection.group_by_match`.

        Parameters
        ----------
        variant : Variant
            The variant. Must have fields 'locus', 'ref', and 'alt'.

        score (optional) : PileupCollection -> object
            Function to apply to summarize the evidence for each allele.
            Default: count number of reads.

        Returns
        ----------
        List of (allele, score) pairs. This list will always have at least two
        elements. The first pair in the list is the reference allele. The
        second pair is the alternate. The subsequent items give the "third"
        alleles (neither ref nor alt), if any.
        '''
        split = self.group_by_match(variant)

        def name(allele_to_pileup_collection):
            return ",".join(allele_to_pileup_collection)

        def aggregate_and_score(pileup_collections):
            merged = PileupCollection.merge(*pileup_collections)
            return score(merged)

        result = [
            (name(split.ref), aggregate_and_score(split.ref.values())),
            (name(split.alt), aggregate_and_score(split.alt.values())),
        ]
        result.extend(
            (allele, score(collection))
            for (allele, collection) in split.other.items())
        return result

    def filter(self,
            drop_duplicates=False,
            drop_improper_mate_pairs=False,
            min_mapping_quality=None,
            min_base_quality=None,
            filters=None):
        '''
        Return a new PileupCollection that includes only pileup elements
        satisfying the specified criteria.

        Parameters
        ----------
        drop_duplicates (optional, default False) : boolean
            Remove alignments with the is_duplicate flag.

        drop_improper_mate_pairs (optional, default False) : boolean
            Retain only alignments that have mapped mate pairs, where one
            alignment in the pair is on the forward strand and the other
            is on the reverse.

        min_mapping_quality (optional) : int
            If specified, retain only alignments with mapping quality >= the
            specified threshold.

        min_base_quality (optional) : int
            If specified, retain only pileup elements where the base quality
            for the bases aligning to the pileup's locus are >= the specified
            threshold.

        filters (optional) : list of PileupElement -> bool functions
            User-specified filter functions to apply. This will be called on 
            each PileupElement, and should return True if the element should
            be retained.

        Returns
        ----------
        A new PileupCollection that includes the subset of the pileup elements
        matching all the specified filters.
        '''

        if filters is None:
            filters = []
        if drop_duplicates:
            filters.append(lambda e: not e.alignment.is_duplicate)
        if drop_improper_mate_pairs:
            filters.append(lambda e: e.alignment.is_proper_pair)
        if min_mapping_quality is not None:
            filters.append(
                lambda e: e.alignment.mapping_quality >= min_mapping_quality)
        if min_base_quality is not None:
            filters.append(
                lambda e: e.min_base_quality >= min_base_quality)
        pileups = OrderedDict(
            (locus, pileup.filter(filters))
            for (locus, pileup)
            in self.pileups.items())
        return PileupCollection(pileups=pileups, parent=self)

    def merge(self, *others):
        '''
        Return a new PileupCollection that is the union of self and the other
        specified collections.
        '''
        new_pileups = {}
        for collection in (self,) + others:
            for (locus, pileup) in collection.pileups.items():
                if locus in new_pileups:
                    new_pileups[locus].update(pileup)
                else:
                    new_pileups[locus] = Pileup(locus, pileup.elements)
        return PileupCollection(new_pileups, parent=self)

    @staticmethod
    def from_bam(pysam_samfile, loci):
        '''
        Create a PileupCollection for a set of loci from a BAM file.

        Parameters
        ----------
        pysam_samfile : `pysam.csamfile.Samfile` instance, or filename string
            to a BAM file. The BAM file must be indexed.

        loci : list of Locus instances
            Loci to collect pileups for.

        Returns
        ----------
        PileupCollection instance containing pileups for the specified loci.
        All alignments in the BAM file are included (e.g. duplicate reads,
        secondary alignments, etc.). See `PileupCollection.filter` if these
        need to be removed. 
        '''

        close_on_completion = False
        if typechecks.is_string(pysam_samfile):
            pysam_samfile = pysam.Samfile(pysam_samfile)
            close_on_completion = True        

        try:
            # Map from pyensembl normalized chromosome names used in Variant to
            # the names used in the BAM file.
            chromosome_name_map = {}
            for name in pysam_samfile.references:
                normalized = pyensembl.locus.normalize_chromosome(name)
                chromosome_name_map[normalized] = name

            result = PileupCollection({})

            # Optimization: we sort variants so our BAM reads are localized.
            locus_iterator = itertools.chain.from_iterable(
                (Locus.from_interbase_coordinates(locus_interval.contig, pos)
                    for pos
                    in locus_interval.positions)
                for locus_interval in sorted(loci))
            for locus in locus_iterator:
                result.pileups[locus] = Pileup(locus, [])
                try:
                    chromosome = chromosome_name_map[locus.contig]
                except KeyError:
                    logging.warn("No such contig in bam: %s" % locus.contig)
                    continue
                columns = pysam_samfile.pileup(
                    chromosome,
                    locus.position,
                    locus.position + 1,  # exclusive, 0-indexed
                    truncate=True,
                    stepper="nofilter")
                try:
                    column = next(columns)
                except StopIteration:
                    # No reads align to this locus.
                    continue

                # Note that storing the pileups here is necessary, since the
                # subsequent assertion will invalidate our column.
                pileups = column.pileups
                assert list(columns) == []  # column is invalid after this.
                for pileup_read in pileups:
                    if not pileup_read.is_refskip:
                        element = PileupElement.from_pysam_alignment(
                            locus, pileup_read)
                        result.pileups[locus].append(element)
            return result
        finally:
            if close_on_completion:
                pysam_samfile.close()

