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
            distinct=False,
            sort_key=None):
        self.distinct = distinct
        if distinct:
            elements = set(elements)
        self.elements = sorted(elements, key=sort_key)
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
            for element in self.elements[:limit]
        ]
        contents = "\n".join(element_lines)

        if limit is not None and len(self.elements) > limit:
            contents += "\n  ... and %d more" % (len(self) - limit)
        return "%s\n%s" % (header, contents)

    def __str__(self):
        return self.to_string(limit=50)

    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self.elements)

    def __iter__(self):
        return iter(self.elements)

    def __hash__(self):
        return hash(len(self))

    def __getitem__(self, idx):
        return self.elements[idx]

    def __eq__(self, other):
        if self is other:
            return True
        return (
            self.__class__ == other.__class__ and
            len(self) == len(other) and
            all(x == y for (x, y) in zip(self, other)))

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
            for element in self.elements
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

    def filter_above_threshold(
            self,
            key_fn,
            value_dict,
            threshold,
            default_value=0.0):
        """The code for filtering by gene or transcript expression was pretty
        much identical aside from which identifier you pull off an effect.
        So, factored out the common operations for filtering an effect
        collection into this helper method.

        Parameters
        ----------
        key_fn : callable
            Given an element of this collection, returns a key into `value_dict`

        value_dict : dict
            Dict from keys returned by `extract_key_fn` to float values

        threshold : float
            Only keep elements whose value in `value_dict` is above this
            threshold.

        default_value : float
            Value to use for elements whose key is not in `value_dict`
        """
        def filter_fn(x):
            key = key_fn(x)
            value = value_dict.get(key, default_value)
            return value > threshold
        return self.filter(filter_fn)

    def filter_any_above_threshold(
            self,
            multi_key_fn,
            value_dict,
            threshold,
            default_value=0.0):
        """Like filter_above_threshold but `multi_key_fn` returns multiple
        keys and the element is kept if any of them have a value above
        the given threshold.

        Parameters
        ----------
        multi_key_fn : callable
            Given an element of this collection, returns multiple keys
            into `value_dict`

        value_dict : dict
            Dict from keys returned by `extract_key_fn` to float values

        threshold : float
            Only keep elements whose value in `value_dict` is above this
            threshold.

        default_value : float
            Value to use for elements whose key is not in `value_dict`
        """
        def filter_fn(x):
            for key in multi_key_fn(x):
                value = value_dict.get(key, default_value)
                if value > threshold:
                    return True
            return False
        return self.filter(filter_fn)
