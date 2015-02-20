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

from .effect_ordering import effect_priority

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

    def annotations(self, raise_on_error=True):
        """
        Determine the impact of each variant, return a list of
        Annotation objects.

        Parameters
        ----------

        raise_on_error : bool, optional
            Raise exception if error is encountered while annotating
            transcripts, otherwise track errors in Annotation.errors
            dictionary (default=True).

        """
        return [
           variant.annotate(
                raise_on_error=raise_on_error)
            for variant
            in self.variants
        ]

    def high_impact_variants(self, *args, **kwargs):
        """
        Returns a list of VariantEffect objects for variants predicted
        to have a significant impact on some transcript's function.

        All arguments are passed on to variant_effects(*args, **kwargs).
        """
        effects = self.variant_effects(*args, **kwargs)

    def effects_to_string(self, *args, **kwargs):
        """
        Create a long string with all transcript effects for each mutation,
        grouped by gene (if a mutation affects multiple genes).

        Arguments are passed on to variant_effects(*args, **kwargs).
        """
        lines = []
        for annotation in self.annotations(*args, **kwargs):
            transcript_effect_count = 0
            lines.append("\n%s" % annotation.variant)
            transcript_effect_lists = annotation.gene_transcript_effects
            for gene, transcript_effects in transcript_effect_lists.iteritems():
                gene_name = self.annot.ensembl.gene_name_of_gene_id(gene)
                lines.append("  Gene: %s (%s)" % (gene_name, gene))
                # print transcript effects with more significant impact
                # on top (e.g. FrameShift should go before NoncodingTranscript)
                for transcript_effect in sorted(
                        transcript_effects,
                        key=effect_priority,
                        reverse=True):
                    transcript_effect_count += 1
                    lines.append("  -- %s" % transcript_effect)
            # if we only printed one effect for this gene then
            # it's redundant to print it again as the highest priority effect
            if transcript_effect_count > 1:
                best = annotation.highest_priority_effect
                lines.append("  Highest Priority Effect: %s" % best)
        return "\n".join(lines)

    def print_effects(self, *args, **kwargs):
        """
        Print all variants and their transcript effects (grouped by gene).

        Arguments are passed on to effects_to_string(*args, **kwargs).
        """
        print(self.effects_to_string(*args, **kwargs))

    def reference_names(self):
        """
        All distinct reference names used by Variants in this
        collection.
        """
        return { variant.reference_name for variant in self.variants }
