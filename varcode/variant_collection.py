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
from collections import Counter

from .collection import Collection
from .effect_collection import EffectCollection
from .common import memoize

class VariantCollection(Collection):
    def __init__(self, variants, path=None, distinct=True):
        """Construct a VariantCollection from a list of Variant records.

        Parameters
        ----------
        variants : iterable
            Variant objects contained in this VariantCollection

        path : str, optional
            File path from which we loaded variants, though the current
            VariantCollection may only contain a subset of them.

        distinct : bool
            Don't keep repeated variants
        """
        Collection.__init__(
            self,
            elements=variants,
            path=path,
            distinct=distinct)

    @memoize
    def dataframe(self):
        """
        Construct a dataframe of the core Variant fields for all variants
        in this collection.
        """
        columns = [
            "reference_name"
            "contig",
            "start",
            "end",
            "ref",
            "alt",
            "gene_names",
        ]
        self._dataframe_from_columns(columns)

    @memoize
    def effects(self, raise_on_error=True):
        """
        Parameters
        ----------
        raise_on_error : bool, optional
            If exception is raised while determining effect of variant on a
            transcript, should it be raised? This default is True, meaning
            errors result in raised exceptions, otherwise they are only logged.

        """
        return EffectCollection([
            effect
            for variant in self
            for effect in variant.effects(raise_on_error=raise_on_error)
        ])

    @memoize
    def gene_counts(self):
        """
        Count how many variants overlap each gene name.
        """
        counter = Counter()
        for variant in self:
            for gene_name in variant.gene_names():
                counter[gene_name] += 1
        return counter

    @memoize
    def reference_names(self):
        """
        All distinct reference names used by Variants in this
        collection.
        """
        return set(variant.reference_name for variant in self)

    @memoize
    def multi_groupby_gene_name(self):
        """
        Group variants by the gene names they overlap, which may put each
        variant in multiple groups.
        """
        return self.multi_groupby(lambda x: x.gene_names())

    @memoize
    def detailed_string(self):
        lines = []
        gene_groups = self.multi_groupby_gene_name()
        for gene_name in sorted(gene_groups.keys()):
            lines.append("  %s:" % gene_name)
            for variant in gene_groups[gene_name]:
                lines.append("  -- %s" % variant)
        header = self.short_string()
        joined_lines = "\n".join(lines)
        return "%s\n%s" % (header, joined_lines)
