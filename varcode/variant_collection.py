# Copyright (c) 2014. Mount Sinai School of Medicine
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

from .effects import Substitution
from .effect_ordering import effect_priority, transcript_effect_priority_dict

class VariantCollection(object):

    def __init__(self, variants, original_filename=None):
        """
        Construct a VariantCollection from a list of Variant records and
        the name of a reference genome.

        Parameters
        ----------

        variants : list
            Variant objects contained in this VariantCollection

        original_filename : str, optional
            File from which we loaded variants, though the current
            VariantCollection may only contain a subset of them.
        """
        self.variants = variants
        self.original_filename = original_filename

    def __len__(self):
        return len(self.variants)

    def __iter__(self):
        return iter(self.variants)

    def __eq__(self, other):
        return (
            isinstance(other, VariantCollection) and
            len(self.variants) == len(other.variants) and
            all(v1 == v2 for (v1, v2) in zip(self.variants, other.variants)))

    def __str__(self):
        fields = [
            ("n_variants", len(self.variants)),
            ("reference", self.reference_name)
        ]

        if self.original_filename:
            fields.append(("filename", self.original_filename))

        s = "VariantCollection(%s)" % (
            ", ".join(
                "%s=%s" % (k,v) for (k,v) in fields))
        for variant in self.variants:
            s += "\n\t%s" % variant
        return s

    def __repr__(self):
        return str(self)

    def clone_metadata(self, new_variants):
        """
        Create copy of VariantCollection with same metadata but possibly
        different Variant entries.
        """

        return VariantCollection(
            variants=new_variants,
            original_filename=self.original_filename)

    def clone(self):
        return self.clone_metadata(list(self.variants))

    def drop_duplicates(self):
        """
        Create a new VariantCollection without any duplicate variants.
        """
        return self.clone_metadata(set(self.variants))

    def variant_effects(
            self,
            min_effect_class=None,
            only_coding_transcripts=False,
            raise_on_error=True):
        """
        Returns a list containing one VariantEffectCollection object for each
        variant in this VariantCollection.

        Parameters
        ----------
        min_effect_class : TranscriptMutationEffect, optional
            Only return EffectCollections for variants whose highest priority
            effect is at least as significant as this effect class.

        only_coding_transcripts : bool, optional
            Only annotate variant effects on coding transcripts.

        raise_on_error : bool, optional
            If exception is raised while determining effect of variant on a
            transcript, should it be raised? This default is True, meaning
            errors result in raised exceptions. If raise_on_error=False then
            exceptions are logged in VariantEffectCollection.errors.
        """
        results = []

        if min_effect_class:
            min_priority = transcript_effect_priority_dict[min_effect_class]
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
                    (effect_priority(best_effect) > min_priority)):
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
            transcript_effect_lists = effect_collection.gene_transcript_effects
            for gene_id, effects in transcript_effect_lists.iteritems():
                gene_name = variant.ensembl.gene_name_of_gene_id(gene)
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
            if len(effect_collection.transcript_effects) > 1:
                best = effect_collection.highest_priority_effect
                lines.append("  Highest Priority Effect: %s" % best)
        return "\n".join(lines)

    def reference_names(self):
        """
        All distinct reference names used by Variants in this
        collection.
        """
        return { variant.reference_name for variant in self.variants }
