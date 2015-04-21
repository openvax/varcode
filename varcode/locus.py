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

class Locus(namedtuple("Locus", "contig start end")):
    '''
    A genomic interval in 0-indexed interbase coordinates.

    See this blog post for a discussion on coordinate systems:
        http://alternateallele.blogspot.com/2012/03/genome-coordinate-conventions.html
    '''

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
    def from_inclusive_coordinates(contig, start, end=None):
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

    @staticmethod
    def from_interbase_coordinates(contig, start, end=None):
        '''
        Given coordinates in 0-based interbase coordinates, return a Locus
        instance.
        '''
        typechecks.require_string(contig)
        typechecks.require_integer(start)
        if end is None:
            end = start + 1
        typechecks.require_integer(end)
        return Locus(contig, start, end)
