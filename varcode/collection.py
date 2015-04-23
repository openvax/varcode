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

from collections import OrderedDict

class Collection(object):
    """
    Methods shared by EffectCollection and VariantCollection.
    """

    def __init__(self, elements, filename=None, distinct=False):
        if distinct:
            elements = set(elements)
        self._elements = list(sorted(elements))
        self.filename = filename

    def __str__(self):
        file_str = ""
        if self.filename:
            file_str = " from '%s'" % self.filename
        count = len(self)
        contents = ""
        if 1 <= count <= 2:
            contents = ": %s" % ", ".join("%s" % x for x in self)
        elif count > 2:
            contents = ": %s, ..., %s" % (self[0], self[-1])
        return "<%s%s with %d elements%s>" % (
            self.__class__.__name__,
            file_str,
            count,
            contents)

    def __len__(self):
        return len(self._elements)

    def __iter__(self):
        return iter(self._elements)

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash(len(self))

    def __getitem__(self, idx):
        return self._elements[idx]

    def __eq__(self, other):
        return (
            self.__class__ == other.__class__ and
            len(self) == len(other) and
            all(x == y for (x, y) in zip(self._elements, other._elements)))

    def _clone_metadata(self, new_elements):
        """
        Create copy of VariantCollection with same metadata but possibly
        different Variant entries.
        """
        return self.__class__(
            elements=new_elements,
            filename=self.filename)

    def groupby(self, key_fn):
        """
        Parameters
        ----------
        key_fn : function
            Takes an effect or variant, returns a grouping key.
        """
        result_dict = OrderedDict()

        for x in self:
            key = key_fn(x)
            if key in result_dict:
                result_dict[key].append(x)
            else:
                result_dict[key] = [x]

        my_class = self.__class__

        # convert result lists into same Collection type as this one
        return OrderedDict(
            (k, my_class(elements))
            for (k, elements)
            in result_dict.items())
