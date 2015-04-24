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

from collections import Counter, OrderedDict

from .collection import Collection
from .common import memoize
from .effect_ordering import (
    effect_priority,
    effect_sort_key,
    top_priority_effect
)

class EffectCollection(Collection):
    """
    Collection of MutationEffect objects and helpers for grouping or filtering
    them.
    """
    def __init__(self, effects, filename=None):
        Collection.__init__(
            self,
            elements=effects,
            filename=filename)

    def groupby_variant(self):
        return self.groupby(key_fn=lambda effect: effect.variant)

    def groupby_gene(self):
        return self.groupby(key_fn=lambda effect: effect.gene)

    def groupby_gene_name(self):
        return self.groupby(key_fn=lambda effect: effect.gene_name)

    def groupby_gene_id(self):
        return self.groupby(key_fn=lambda effect: effect.gene_id)

    def groupby_transcript(self):
        return self.groupby(key_fn=lambda effect: effect.transcript)

    def groupby_transcript_name(self):
        return self.groupby(key_fn=lambda effect: effect.transcript_name)

    def groupby_transcript_id(self):
        return self.groupby(key_fn=lambda effect: effect.transcript_id)

    def summary_string(self):
        """
        Create a long string with all transcript effects for each mutation,
        grouped by gene (if a mutation affects multiple genes).
        """
        lines = []
        for variant, variant_effects in self.groupby_variant():
            lines.append("\n%s" % variant)

            gene_effects_groups = variant_effects.groupby_gene_id()
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
            if len(variant_effects) > 1:
                # since summary effect also calls into Variant.effects,
                # give it the same arguments
                # (this is the downside of not having a VariantEffectCollection)
                best = effect.top_effect()
                lines.append("  Highest Priority Effect: %s" % best)
        return "\n".join(lines)

    @memoize
    def top_priority_effect(self):
        """Highest priority MutationEffect of all genes/transcripts overlapped
        by this variant. If this variant doesn't overlap anything, then this
        this method will return an Intergenic effect.

        If multiple effects have the same priority, then return the one
        which is associated with the longest transcript.
        """
        return top_priority_effect(self._elements)

    # TODO: find a way to express these kinds of methods without
    # duplicating every single groupby_* method
    @memoize
    def top_priority_effect_per_variant(self):
        """Highest priority effect for each unique variant"""
        return OrderedDict(
            (variant, top_priority_effect(variant_effects))
            for (variant, variant_effects)
            in self.groupby_variant().items())

    @memoize
    def top_priority_effect_per_transcript_id(self):
        """Highest priority effect for each unique transcript ID"""
        return OrderedDict(
            (transcript_id, top_priority_effect(variant_effects))
            for (transcript_id, variant_effects)
            in self.groupby_transcript_id().items())

    @memoize
    def top_priority_effect_per_gene_id(self):
        """Highest priority effect for each unique gene ID"""
        return OrderedDict(
            (gene_id, top_priority_effect(variant_effects))
            for (gene_id, variant_effects)
            in self.groupby_gene_id().items())

    def effect_expression(self, expression_levels):
        """
        Parameters
        ----------
        expression_levels : dict
            Dictionary mapping transcript IDs to length-normalized expression
            levels (either FPKM or TPM)

        Returns dictionary mapping each transcript effect to an expression
        quantity. Effects that don't have an associated transcript
        (e.g. Intergenic) will not be included.
        """
        return OrderedDict(
            (effect, expression_levels.get(effect.transcript.id, 0.0))
            for effect in self
            if effect.transcript is not None)

    def top_expression_effect(self, expression_levels):
        """
        Return effect whose transcript has the highest expression level.
        If none of the effects are expressed or have associated transcripts,
        then return None. In case of ties, add lexicographical sorting by
        effect priority and transcript length.
        """
        effect_expression_dict = self.effect_expression(expression_levels)
        if len(effect_expression_dict) == 0:
            return None
        return max(
            effect_expression_dict.items(),
            key=lambda effect, fpkm: (fpkm, effect_sort_key(effect)))

    @memoize
    def gene_counts(self):
        counter = Counter()
        for effect in self:
            counter[effect.gene_name] += 1
        return counter
