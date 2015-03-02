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

from memoized_property import memoized_property
from typechecks import require_iterable_of

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
        self.variants = set(variants)
        self.original_filename = original_filename

    def __len__(self):
        return len(self.variants)

    def __iter__(self):
        return iter(sorted(self.variants))

    def __eq__(self, other):
        return (
            isinstance(other, VariantCollection) and
            self.variants == other.variants)

    def variant_summary(self):
        """
        Returns a string indicating each variant in the collection.
        """
        fields = [
            ("n_variants", len(self.variants)),
            ("reference", ", ".join(self.reference_names))
        ]

        if self.original_filename:
            fields.append(("filename", self.original_filename))

        s = "VariantCollection(%s)" % (
            ", ".join(
                "%s=%s" % (k, v) for (k, v) in fields))
        for variant in self.variants:
            s += "\n\t%s" % variant
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

    def variant_effects(
            self,
            high_impact=False,
            only_coding_transcripts=False,
            raise_on_error=True):
        """
        Returns a list containing one VariantEffectCollection object for each
        variant in this VariantCollection.

        Parameters
        ----------
        high_impact : bool, optional
            Only show effects which make changes to coding transcripts that
            are predicted to be deleterious (amino acid sequence changes,
            changes to essential splice sites, frameshifts, &c)

        only_coding_transcripts : bool, optional
            Only annotate variant effects on coding transcripts.

        raise_on_error : bool, optional
            If exception is raised while determining effect of variant on a
            transcript, should it be raised? This default is True, meaning
            errors result in raised exceptions. If raise_on_error=False then
            exceptions are logged in VariantEffectCollection.errors.
        """
        results = []

        if high_impact:
            min_priority = transcript_effect_priority_dict[Substitution]
        else:
            min_priority = -1

        for variant in self.variants:
            variant_effect_collection = variant.effects(
                only_coding_transcripts=only_coding_transcripts,
                raise_on_error=raise_on_error)

            if only_coding_transcripts and len(variant_effect_collection) == 0:
                # if we only want coding transcripts, then skip all
                # intergenic and non-coding gene variants
                continue

            best_effect = variant_effect_collection.highest_priority_effect
            # either this variant is intergenic and there's no minimum
            # threshold for effect priority or the highest impact effect
            # is higher priority than the min_priority
            if ((best_effect is None and min_priority < 0) or
                    (effect_priority(best_effect) >= min_priority)):
                results.append(variant_effect_collection)
        return results

    def effect_summary(self, *args, **kwargs):
        """
        Create a long string with all transcript effects for each mutation,
        grouped by gene (if a mutation affects multiple genes).

        Arguments are passed on to self.variant_effects(*args, **kwargs).
        """
        lines = []
        for effect_collection in self.variant_effects(*args, **kwargs):
            variant = effect_collection.variant
            lines.append("\n%s" % variant)
            transcript_effect_lists = effect_collection.gene_effect_groups
            for (gene_id, effects) in transcript_effect_lists.items():
                gene_name = variant.ensembl.gene_name_of_gene_id(gene_id)
                lines.append("  Gene: %s (%s)" % (gene_name, gene_id))
                # place transcript effects with more significant impact
                # on top (e.g. FrameShift should go before NoncodingTranscript)
                for effect in sorted(
                        effects,
                        key=effect_priority,
                        reverse=True):
                    lines.append("  -- %s" % effect)
            # if we only printed one effect for this gene then
            # it's redundant to print it again as the highest priority effect
            if len(effect_collection) > 1:
                best = effect_collection.highest_priority_effect
                lines.append("  Highest Priority Effect: %s" % best)
        return "\n".join(lines)

    @memoized_property
    def reference_names(self):
        """
        All distinct reference names used by Variants in this
        collection.
        """
        return {variant.reference_name for variant in self.variants}

    @memoized_property
    def gene_counts(self, only_coding=False):
        """
        Count how many variants overlap each gene name.
        """
        counter = Counter()
        for variant in self.variants:
            for gene_name in variant.gene_names():
                counter[gene_name] += 1
        return counter
