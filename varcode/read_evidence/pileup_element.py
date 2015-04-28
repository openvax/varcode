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

from . import alignment_key

class PileupElement(object):
    '''
    A PileupElement represents the segment of an alignment that aligns to a
    particular base in the reference.

    Attributes
    ----------
    locus : Varcode.Locus
        The reference locus. Must be length 1, i.e. a single base.

    offset_start : int
        0-based start offset into the alignment sequence, inclusive

    offset_end : int
        0-based end offset into the alignment sequence, exclusive

    alignment : pysam.AlignedSegment
        pysam alignment instance

    alignment_key : tuple
        value computed from the alignment instance that uniquely specifies its
        properties. Used for comparisons since pysam.AlignedSegment instances
        do not support a useful notion of equality (they compare using object
        identity). See `read_evidence.alignment_key` for the implementation of
        this key.
    '''
    def __init__(self, locus, offset_start, offset_end, alignment):
        '''
        Construct a PileupElement object.
        '''
        assert offset_end >= offset_start, \
            "offset_start=%d > offset_end=%d" % (offset_start, offset_end)
        self.locus = locus
        self.offset_start = offset_start
        self.offset_end = offset_end
        self.alignment = alignment
        self.alignment_key = alignment_key(self.alignment)

    def fields(self):
        '''
        Fields that should be considered for our notion of object equality.
        '''
        return (
            self.locus, self.offset_start, self.offset_end, self.alignment_key)

    def __eq__(self, other):
        return hasattr(other, "fields") and self.fields() == other.fields()

    def __hash__(self):
        return hash(self.fields())

    @property
    def bases(self):
        '''
        The sequenced bases in the alignment that align to this locus in the
        genome, as a string.

        Empty string in the case of a deletion. String of length > 1 if there
        is an insertion here.
        '''
        sequence = self.alignment.query_sequence
        assert self.offset_end <= len(sequence), \
            "End offset=%d > sequence length=%d. CIGAR=%s. SEQUENCE=%s" % (
                self.offset_end,
                len(sequence),
                self.alignment.cigarstring,
                sequence)
        return sequence[self.offset_start:self.offset_end]

    @property
    def base_qualities(self):
        '''
        The phred-scaled base quality scores corresponding to `self.bases`, as
        a list.
        '''
        return self.alignment.query_qualities[
            self.offset_start:self.offset_end]

    @property
    def min_base_quality(self):
        '''
        The minimum of the base qualities. In the case of a deletion, in which
        case there are no bases in this PileupElement, the minimum is taken
        over the sequenced bases immediately before and after the deletion.
        '''
        try:
            return min(self.base_qualities)
        except ValueError:
            # We are mid-deletion. We return the minimum of the adjacent bases.
            assert self.offset_start == self.offset_end
            adjacent_qualities = [
                self.alignment.query_qualities[offset]
                for offset in [self.offset_start - 1, self.offset_start]
                if 0 <= offset < len(self.alignment.query_qualities)
            ]
            return min(adjacent_qualities)

    @staticmethod
    def from_pysam_alignment(locus, pileup_read):
        '''
        Factory function to create a new PileupElement from a pysam
        `PileupRead`.

        Parameters
        ----------
        locus : varcode.Locus
            Reference locus for which to construct a PileupElement. Must
            include exactly one base.

        pileup_read : pysam.calignmentfile.PileupRead
            pysam PileupRead instance. Its alignment must overlap the locus.

        Returns
        ----------
        PileupElement

        '''
        assert not pileup_read.is_refskip, (
            "Can't create a PileupElement in a refskip (typically an intronic "
            "gap in an RNA alignment)")

        # Pysam has an `aligned_pairs` method that gives a list of
        # (offset, locus) pairs indicating the correspondence between bases in
        # the alignment and reference loci. Here we use that to compute
        # offset_start and offset_end.
        #
        # This is slightly tricky in the case of insertions and deletions.
        # Here are examples of the desired logic.
        #
        # Target locus = 1000
        #
        # (1) Simple case: matching bases.
        #
        # OFFSET           LOCUS
        # 0                999
        # 1                1000
        # 2                1001
        #
        # DESIRED RESULT: offset_start=1, offset_end=2.
        #
        #
        # (2) A 1 base insertion at offset 2.
        #
        # OFFSET           LOCUS
        # 0                999
        # 1                1000
        # 2                None
        # 3                1001
        #
        # DESIRED RESULT: offset_start = 1, offset_end=3.
        #
        #
        # (3) A 2 base deletion at loci 1000 and 1001.
        #
        # OFFSET           LOCUS
        # 0                999
        # None             1000
        # None             1001
        # 1                1002
        #
        # DESIRED RESULT: offset_start = 1, offset_end=1.
        #
        offset_start = None
        offset_end = len(pileup_read.alignment.query_sequence)
        # TODO: doing this with get_blocks() may be faster.
        for (offset, position) in pileup_read.alignment.aligned_pairs:
            if offset is not None and position is not None:
                if position == locus.position:
                    offset_start = offset
                elif position > locus.position:
                    offset_end = offset
                    break
        if offset_start is None:
            offset_start = offset_end
        
        assert pileup_read.is_del == (offset_end - offset_start == 0), \
            "Deletion=%s but | [%d,%d) |=%d for locus %d in: \n%s" % (
                pileup_read.is_del,
                offset_start,
                offset_end,
                offset_end - offset_start,
                locus.position,
                pileup_read.alignment.aligned_pairs)

        assert offset_end >= offset_start
        result = PileupElement(
            locus, offset_start, offset_end, pileup_read.alignment)
        return result   

