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

from .effect_ordering import top_priority_transcript_effect

from pyensembl.biotypes import is_coding_biotype

class Annotation(object):
    """
    An Annotation object is a container for all the TranscriptEffects of
    a particular mutation, as well as some properties that attempt to
    summarize the impact of the mutation across all genes/transcripts.
    """

    def __init__(
            self,
            variant,
            genes,
            gene_transcript_effects,
            errors={}):
        """
        variant : Variant

        genes : list
            List of Gene objects

        gene_transcript_effects : dict
            Dictionary from gene ID to list of transcript variant effects

        errors : dict, optional
            Mapping from transcript objects to strings explaining error which
            was encountered while trying to annotate them.
        """
        self.variant = variant
        self.genes = genes
        self.gene_transcript_effects = gene_transcript_effects

        # dictionary mapping from transcript IDs to transcript mutation effects
        self.transcript_effects = {}
        for (_, transcript_effects) in self.gene_transcript_effects.items():
            for effect in transcript_effects:
                self.transcript_effects[effect.transcript.id] = effect

        # if our variant overlaps any genes, then choose the highest
        # priority transcript variant, otherwise call the variant "Intergenic"
        if len(self.transcript_effects) > 0:
            self.highest_priority_effect = top_priority_transcript_effect(
                self.transcript_effects.values())
            highest_priority_class = self.highest_priority_effect.__class__
            self.variant_summary = highest_priority_class.__name__
        else:
            self.highest_priority_effect = None
            self.variant_summary = "Intergenic"

        self.errors = errors

    def __str__(self):
        fields = [
            ("variant", self.variant.short_description()),
            ("genes", [gene.name for gene in self.genes]),
            ("transcript_effects", self.transcript_effects)
        ]
        if self.errors:
            fields.append( ("errors", self.errors) )
        return "VariantEffect(%s)" % (
            ", ".join(["%s=%s" % (k,v) for (k,v) in fields]))

    def __repr__(self):
        return str(self)
