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
from .reference_name import (
    infer_reference_name,
    ensembl_release_number_for_reference_name
)
from .variant_annotator import VariantAnnotator

class VariantCollection(object):

    def __init__(
            self,
            variants,
            reference_path=None,
            reference_name=None,
            ensembl_release=None,
            original_filename=None):
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

        reference_path : str, optional
            Path to reference FASTA file.

        reference_name : str, optional
            Name of reference genome (e.g. "GRCh37", "hg18"). If not given
            infer from reference path or from ensembl_release.

        ensembl_release : int, optional
            If not specified, infer Ensembl release from reference_name
        """
        self.variants = variants
        self.reference_path = reference_path
        if reference_name:
            # convert from e.g. "hg19" to "GRCh37"
            #
            # TODO: actually handle the differences between these references
            # instead of just treating them as interchangeable
            self.reference_name = infer_reference_name(reference_name)
        else:
            if reference_path:
                self.reference_name = infer_reference_name(reference_path)
            else:
                raise ValueError(
                    "Must specify one of reference_path or reference_name")

        if ensembl_release:
            self.ensembl_release = ensembl_release
        else:
            self.ensembl_release = ensembl_release_number_for_reference_name(
                self.reference_name)

        self.annot = VariantAnnotator(ensembl_release=self.ensembl_release)
        self.original_filename = original_filename

    def __len__(self):
        return len(self.variants)

    def __iter__(self):
        return iter(self.variants)

    def __str__(self):
        s = "VariantCollection(filename=%s, reference=%s)" % (
            self.filename, self.reference_name)
        for variant in self.variants:
            s += "\n\t%s" % variant
        return s

    def __repr__(self):
        return str(self)

    def clone(self, new_variants=None):
        """
        Create copy of VariantCollection with same metadata but possibly
        different Variant entries. If no variants provided, then just make a
        copy of self.variants.
        """
        if new_variants is None:
            new_variants = self.variants
        return VariantCollection(
            variants=list(new_variants),
            original_filename=self.original_filename,
            reference_path=self.reference_path,
            reference_name=self.reference_name,
            ensembl_release=self.ensembl_release)

    def drop_duplicates(self):
        """
        Create a new VariantCollection without any duplicate variants.
        """
        return self.clone(set(self.variants))

    def variant_effects(self, raise_on_error=True):
        """
        Determine the impact of each variant, return a list of
        VariantEffect objects.

        Parameters
        ----------

        raise_on_error : bool, optional
            Raise exception if error is encountered while annotating
            transcripts, otherwise track errors in VariantEffect.errors
            dictionary (default=True).

        """
        return [
            self.annot.effect(
                variant=variant,
                raise_on_error=raise_on_error)
            for variant
            in self.variants
        ]

    def effects_to_string(self, *args, **kwargs):
        """
        Create a long string with all transcript effects for each mutation,
        grouped by gene (if a mutation affects multiple genes).

        Arguments are passed on to variant_effects(*args, **kwargs).
        """
        lines = []
        for variant_effect in self.variant_effects(*args, **kwargs):
            transcript_effect_count = 0
            lines.append("\n%s" % variant_effect.variant)
            transcript_effect_lists = variant_effect.gene_transcript_effects
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
                best = variant_effect.highest_priority_effect
                lines.append("  Highest Priority Effect: %s" % best)
        return "\n".join(lines)

    def print_effects(self, *args, **kwargs):
        """
        Print all variants and their transcript effects (grouped by gene).

        Arguments are passed on to effects_to_string(*args, **kwargs).
        """
        print(self.effects_to_string(*args, **kwargs))
