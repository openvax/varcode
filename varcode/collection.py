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
import os.path
from collections import defaultdict

class Collection(object):
    """
    Methods shared by EffectCollection and VariantCollection.
    """

    def __init__(
            self,
            elements,
            path=None,
            distinct=False):
        self.distinct = distinct
        if distinct:
            elements = set(elements)
        self._elements = list(sorted(elements))
        self.path = path
        if path:
            # get the filename without any directory prefix
            self.filename = os.path.split(path)[1]
        else:
            self.filename = None

    def short_string(self):
        """
        Compact string representation which doesn't print any of the
        collection elements.
        """
        file_str = ""
        if self.filename:
            file_str = " from '%s'" % self.filename
        return "<%s%s with %d elements>" % (
            self.__class__.__name__,
            file_str,
            len(self))

    def to_string(self, limit=None):
        """
        Create a string representation of this collection, showing up to
        `limit` items.
        """
        header = self.short_string()
        if len(self) == 0:
            return header
        contents = ""
        element_lines = [
            "  -- %s" % (element,)
            for element in self._elements[:limit]
        ]
        contents = "\n".join(element_lines)

        if limit is not None and len(self._elements) > limit:
            contents += "\n  ... and %d more" % (len(self) - limit)
        return "%s\n%s" % (header, contents)

    def __str__(self):
        return self.to_string(limit=50)

    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self._elements)

    def __iter__(self):
        return iter(self._elements)

    def __hash__(self):
        return hash(len(self))

    def __getitem__(self, idx):
        return self._elements[idx]

    def __eq__(self, other):
        return (
            self.__class__ == other.__class__ and
            len(self) == len(other) and
            all(x == y for (x, y) in zip(self._elements, other._elements)))

    def clone_with_new_elements(self, new_elements):
        """
        Create another collection of the same class and  with same metadata
        but possibly different Variant or Effect entries.
        """
        return self.__class__(
            new_elements,
            path=self.path,
            distinct=self.distinct)

    def filter(self, filter_fn):
        return self.clone_with_new_elements([
            element
            for element in self._elements
            if filter_fn(element)])

    def groupby(self, key_fn):
        """
        Parameters
        ----------
        key_fn : function
            Takes an effect or variant, returns a grouping key.
        """
        result_dict = defaultdict(list)

        for x in self:
            result_dict[key_fn(x)].append(x)

        # convert result lists into same Collection type as this one
        return {
            k: self.clone_with_new_elements(elements)
            for (k, elements)
            in result_dict.items()
        }

    def multi_groupby(self, key_fn):
        """
        Like a groupby but expect the key_fn to return multiple keys for
        each element.
        """
        result_dict = defaultdict(list)

        for x in self:
            for key in key_fn(x):
                result_dict[key].append(x)

        # convert result lists into same Collection type as this one
        return {
            k: self.clone_with_new_elements(elements)
            for (k, elements)
            in result_dict.items()
        }
