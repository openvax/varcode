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

from collections import OrderedDict, Counter

from .common import memoize

class BaseCollection(object):
    """
    Methods shared by EffectCollection and VariantCollection
    """
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

    @memoize
    def reference_names(self):
        """
        All distinct reference names used by Variants in this
        collection.
        """
        return set(variant.reference_name for variant in self.variants)

    @memoize
    def gene_counts(self):
        """
        Count how many variants overlap each gene name.
        """
        counter = Counter()
        for variant in self.variants:
            for gene_name in variant.gene_names():
                counter[gene_name] += 1
        return counter
