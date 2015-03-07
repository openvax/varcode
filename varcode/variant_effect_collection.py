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

from .effects import Intergenic, Intragenic, TranscriptMutationEffect
from .effect_ordering import top_priority_transcript_effect


class VariantEffects(object):
    """Collection of all the MutationEffect objects associated with
    each transcript overlapping a given variant, as well as some properties
    that attempt to summarize the impact of the mutation across all
    genes/transcripts.
    """
    def __init__(
            self,
            variant,
            gene_effect_groups,
            errors={}):
        """
        variant : Variant

        gene_effect_groups : dict
            Dictionary from gene ID to list of transcript variant effects

        errors : dict, optional
            Mapping from transcript objects to strings explaining error which
            was encountered while trying to annotate them.
        """
        self.variant = variant
        self.gene_effect_groups = gene_effect_groups

        # dictionary mapping from transcript IDs to transcript mutation effects
        # may only contain a subset of all effects due to those falling outside
        # of known transcripts
        self.transcript_effect_dict = {}
        for (_, effect_list) in self.gene_effect_groups.items():
            for effect in effect_list:
                if isinstance(effect, TranscriptMutationEffect):
                    self.transcript_effect_dict[effect.transcript.id] = effect

        # if our variant overlaps any genes, then choose the highest
        # priority transcript variant, otherwise call the variant "Intergenic"
        if len(self.transcript_effect_dict) > 0:
            self.summary_effect = top_priority_transcript_effect(
                self.transcript_effect_dict.values())
        # intragenic variant overlaps a gene but not any transcripts
        elif len(variant.gene_names()) > 0:
            self.summary_effect = Intragenic(variant)
        else:
            # variant does not fall within any known genes
            # TODO: annotate Upstream and Downstream variant
            # effects, since those might have functional significance
            self.summary_effect = Intergenic(variant)

        self.errors = errors

    def __str__(self):
        fields = [
            ("variant", self.variant.short_description()),
            ("genes", self.gene_names()),
            ("effects", self.effects())
        ]
        if self.errors:
            fields.append(("errors", self.errors))
        return "VariantEffects(%s)" % (
            ", ".join(["%s=%s" % (k, v) for (k, v) in fields]))

    def __repr__(self):
        return str(self)

    def __len__(self):
        """
        Length of a VariantEffectCollection is the number of effect objects
        it contains.
        """
        return len(self.transcript_effect_dict)
