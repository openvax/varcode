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
# used for local import of EffectCollection
import importlib

from typechecks import require_iterable_of

from .base_collection import BaseCollection
from .common import memoize
from .variant import Variant

class VariantCollection(BaseCollection):
    def __init__(self, variants, original_filename=None):
        """
        Construct a VariantCollection from a list of Variant records and
        the name of a reference genome.

        Parameters
        ----------

        variants : iterable
            Variant objects contained in this VariantCollection

        original_filename : str, optional
            File from which we loaded variants, though the current
            VariantCollection may only contain a subset of them.
        """
        require_iterable_of(variants, Variant, "variants")
        self.variants = list(sorted(set(variants)))
        self.original_filename = original_filename

    def __len__(self):
        return len(self.variants)

    def __iter__(self):
        return iter(self.variants)

    def __hash__(self):
        return hash(len(self.variants))

    def __eq__(self, other):
        return (
            isinstance(other, VariantCollection) and
            len(self.variants) == len(other.variants) and
            all(v1 == v2 for (v1, v2) in zip(self.variants, other.variants)))

    def summary_string(self):
        """
        Returns a string indicating each variant in the collection.
        """
        fields = [
            ("n_variants", len(self.variants)),
            ("reference", ", ".join(self.reference_names()))
        ]

        if self.original_filename:
            fields.append(("filename", self.original_filename))

        s = "VariantCollection(%s)" % (
            ", ".join(
                "%s=%s" % (k, v) for (k, v) in fields))
        for variant in self.variants:
            gene_names = variant.gene_names()
            if len(gene_names):
                gene_names_string = " : %s" % ", ".join(gene_names)
            else:
                gene_names_string = ""
            s += "\n\t%s%s" % (variant, gene_names_string)
        return s

    def __str__(self):
        suffix = ""
        if self.original_filename:
            suffix = ' from "%s"' % self.original_filename
        return ("<VariantCollection of %d variants%s>" %
            (len(self.variants), suffix))

    def __repr__(self):
        return str(self)

    def _clone_metadata(self, new_variants):
        """
        Create copy of VariantCollection with same metadata but possibly
        different Variant entries.
        """
        return VariantCollection(
            variants=new_variants,
            original_filename=self.original_filename)

    @memoize
    def effects(self, raise_on_error=True):
        """
        Generator of effects from all variants.

        Parameters
        ----------
        raise_on_error : bool, optional
            If exception is raised while determining effect of variant on a
            transcript, should it be raised? This default is True, meaning
            errors result in raised exceptions, otherwise they are only logged.

        min_priority_effect_class : type, optional
            Only return variants with priority at least as great the given type.
            Typical value  for coding effects is `effects.Substitution`

        only_gene_ids : set, optional
            If given, then only return effects for gene IDs present in this set

        only_transcript_ids : set, optional
            If given, then only return effects for transcript IDs present
            in this set.
        """
        # importing EffectCollection locally since Python won't
        # let us express a mutual dependency between two modules
        effect_collection = importlib.import_module(".effect_collection")

        effect_list = []
        for variant in self.variants:
            for effect in variant.effects(raise_on_error=raise_on_error):
                effect_list.append(effect)
        return effect_collection.EffectCollection(effect_list)
