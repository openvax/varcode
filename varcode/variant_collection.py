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
    def __init__(self, variants, filename=None):
        """
        Construct a VariantCollection from a list of Variant records and
        the name of a reference genome.

        Parameters
        ----------
        variants : iterable
            Variant objects contained in this VariantCollection

        filename : str, optional
            File from which we loaded variants, though the current
            VariantCollection may only contain a subset of them.
        """
        Collection.__init__(
            self,
            elements=variants,
            filename=filename,
            distinct=True)

    @memoize
    def summary_string(self):
        """
        Returns a string indicating each variant in the collection.
        """
        fields = [
            ("n_variants", len(self)),
            ("reference", ", ".join(self.reference_names()))
        ]

        if self.filename:
            fields.append(("filename", self.filename))

        s = "VariantCollection(%s)" % (
            ", ".join(
                "%s=%s" % (k, v) for (k, v) in fields))
        for variant in self:
            gene_names = variant.gene_names()
            if len(gene_names):
                gene_names_string = " : %s" % ", ".join(gene_names)
            else:
                gene_names_string = ""
            s += "\n\t%s%s" % (variant, gene_names_string)
        return s

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
