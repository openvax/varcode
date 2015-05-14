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

import typechecks
from collections import namedtuple

from pyensembl import EnsemblRelease, ensembl_grch38
from pyensembl.locus import normalize_chromosome

class Locus(object):
    '''
    A genomic interval in 0-indexed interbase coordinates.

    See this blog post for a discussion on coordinate systems:
        http://alternateallele.blogspot.com/2012/03/genome-coordinate-conventions.html
    '''
    def __init__(
        self,
        contig,
        start=None,
        end=None,
        ensembl=ensembl_grch38):

        typechecks.require_string(contig)
        typechecks.require_integer(start)
        if end is None:
            end = start + 1
        typechecks.require_integer(end)
        self.contig = normalize_chromosome(contig)
        self.start = start
        self.end = end

        # user might supply Ensembl release as an integer
        if isinstance(ensembl, int):
            ensembl = EnsemblRelease(release=ensembl)
        typechecks.require_instance(ensembl, EnsemblRelease, "ensembl")
        self.ensembl = ensembl

    def __hash__(self):
        return hash(
            (self.contig, self.start, self.end, self.ensembl.release))

    def __len__(self):
        return self.end - self.start

    def __str__(self):
        if len(self) == 1:
            template = "<Locus {contig}:{start} release={release}>"
        else:
            template = "<Locus {contig}:[{start}, {end}) release={release}>"
        return template.format(release=self.ensembl.release, **self.__dict__)

    def __repr__(self):
        return str(self)

    def __lt__(self, other):
        typechecks.require_instance(other, Locus, name="locus")
        if self.contig == other.contig:
            return self.start < other.start
        return self.contig < other.contig

    @property
    def base1_start(self):
        return self.start + 1

    @property
    def base1_end(self):
        if len(self.positions) == 0:
            raise ValueError(
                "Base 1 (inclusive) coordinates cannot represent "
                "zero-length locus.")
        return self.end  # change to base1 coordinates cancels out change to inclusive interval

    @property
    def positions(self):
        '''
        A Python range object giving the bases included in this locus.
        '''
        return range(self.start, self.end)

    @property
    def position(self):
        '''
        If this locus spans a single base, this property gives that position.
        Otherwise, raises a ValueError.
        '''
        if self.end != self.start + 1:
            raise ValueError("Not a single base: %s" % str(self))
        return self.start

    # Factory functions.
    @staticmethod
    def from_base1_coordinates(contig, start, end=None):
        '''
        Given coordinates in 1-based coordinates that are inclusive on start
        and end, return a Locus instance. Locus instances are always 0-based
        "interbase" coordinates.
        '''
        typechecks.require_string(contig)
        typechecks.require_integer(start)
        if end is None:
            end = start
        typechecks.require_integer(end)
        return Locus(contig, start - 1, end)

        
