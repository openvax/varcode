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

from collections import OrderedDict


class Pileup(object):
    '''
    A Pileup is a collection of PileupElement instances at a particular locus.

    Attributes
    ----------
    locus : Varcode.Locus
        The reference locus. Must be length 1, i.e. a single base.

    elements : OrderedDict of PileupElement instances
        This is logically and ordered set, which we implement as an OrderedDict
        with all values mapping to None.
    '''
    def __init__(self, locus, elements):
        '''
        Construct a new Pileup.

        Parameters
        ----------
        locus : Varcode.Locus
            The reference locus. Must be length 1, i.e. a single base.

        elements : iterable of PileupElement
            The pileup elements. The locus field of these instances must 
            match the locus parameter.
        '''
        self.locus = locus
        self.elements = OrderedDict((e, None) for e in elements)
        assert all(e.locus == self.locus for e in self.elements)

    def __iter__(self):
        return iter(self.elements)

    def __len__(self):
        return len(self.elements)

    def append(self, element):
        '''
        Append a PileupElement to this Pileup. If an identical PileupElement is
        already part of this Pileup, do nothing.
        '''
        assert element.locus == self.locus, (
            "Element locus (%s) != Pileup locus (%s)"
            % (element.locus, self.locus))
        self.elements[element] = None

    def update(self, other):
        '''
        Add all pileup elements from other into self.
        '''
        assert self.locus == other.locus
        self.elements.update(other.elements)

    def filter(self, filters):
        '''
        Apply filters to the pileup elements, and return a new Pileup with the
        filtered elements removed.

        Parameters
        ----------
        filters : list of PileupElement -> bool callables
            A PileupUp element is retained if all filters return True when
            called on it.
        '''
        new_elements = [
            e for e in self.elements
            if all(function(e) for function in filters)]
        return Pileup(self.locus, new_elements)
