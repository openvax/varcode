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
from collections import Counter, OrderedDict
from itertools import groupby

from typechecks import require_iterable_of

from .common import memoize
from .effects import Substitution
from .effect_ordering import effect_priority, transcript_effect_priority_dict
from .variant import Variant

class VariantCollection(object):

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
        self._cached_values = {}

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
    def high_priority_effects(self, *args, **kwargs):
        """Like VariantCollection.effects() but only returns effects whose
        priority is at least as high as a missense mutation
        (e.g. frameshifts, splice site mutations, &c).

        All arguments are passed on to Variant.top_effect(*args, **kwargs).
        """
        min_priority = transcript_effect_priority_dict[Substitution]
        results = OrderedDict()
        for variant in self.variants:
            effect = variant.top_effect(*args, **kwargs)
            priority = effect_priority(effect)
            if priority > min_priority:
                results[variant] = effect
        return results

    @memoize
    def effects(
            self,
            raise_on_error=True):
        """
        Returns an OrderedDict mapping each variant to list of its
        MutationEffect annotation objects.

        Parameters
        ----------
        raise_on_error : bool, optional
            If exception is raised while determining effect of variant on a
            transcript, should it be raised? This default is True, meaning
            errors result in raised exceptions, otherwise they are only logged.
        """
        all_effects = OrderedDict()
        for variant in self.variants:
            all_effects[variant] = variant.effects(
                raise_on_error=raise_on_error)
        return all_effects

    def full_effect_string(self, *args, **kwargs):
        """
        Create a long string with all transcript effects for each mutation,
        grouped by gene (if a mutation affects multiple genes).

        Arguments are passed on to self.variant_effects(*args, **kwargs).
        """
        lines = []
        effect_dict = self.variant_effects(*args, **kwargs)
        for variant, effects in effect_dict.items():
            lines.append("\n%s" % variant)

            gene_effects_groups = groupby(effects, lambda e: e.gene_id())
            for (gene_id, gene_effects) in gene_effects_groups:
                if gene_id:
                    gene_name = variant.ensembl.gene_name_of_gene_id(gene_id)
                    lines.append("  Gene: %s (%s)" % (gene_name, gene_id))
                # place transcript effects with more significant impact
                # on top (e.g. FrameShift should go before NoncodingTranscript)
                for effect in sorted(
                        gene_effects,
                        key=effect_priority,
                        reverse=True):
                    lines.append("  -- %s" % effect)
            # if we only printed one effect for this gene then
            # it's redundant to print it again as the highest priority effect
            if len(effects) > 1:
                # since summary effect also calls into Variant.effects,
                # give it the same arguments
                # (this is the downside of not having a VariantEffectCollection)
                best = effect.summary_effect(*args, **kwargs)
                lines.append("  Highest Priority Effect: %s" % best)
        return "\n".join(lines)

    @memoize
    def reference_names(self):
        """
        All distinct reference names used by Variants in this
        collection.
        """
        return {
            variant.reference_name
            for variant in self.variants
        }

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
