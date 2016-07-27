# Copyright (c) 2016. Mount Sinai School of Medicine
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

import pandas as pd

from .collection import Collection
from .effect_collection import EffectCollection
from .common import memoize
from .variant import Variant, variant_ascending_position_sort_key


class VariantCollection(Collection):
    def __init__(
            self,
            variants,
            distinct=True,
            sort_key=variant_ascending_position_sort_key,
            source_to_metadata_dict={}):
        """
        Construct a VariantCollection from a list of Variant records.

        Parameters
        ----------
        variants : iterable
            Variant objects contained in this VariantCollection

        distinct : bool
            Don't keep repeated variants

        sort_key : callable

        source_to_metadata_dict : dict
            Dictionary mapping each source name (e.g. VCF path) to a dictionary
            from metadata attributes to values.
        """
        Collection.__init__(
            self,
            elements=variants,
            element_type=Variant,
            distinct=distinct,
            sort_key=sort_key,
            source_to_metadata_dict=source_to_metadata_dict)

    def to_dict(self):
        """
        Since Collection.to_dict() returns a state dictionary with an
        'elements' field we have to rename it to 'variants'.
        """
        state_dict = Collection.to_dict(self)
        state_dict["variants"] = state_dict.pop("elements")
        return state_dict

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

    def gene_counts(self):
        """
        Count how many variants overlap each gene name.
        """
        counter = Counter()
        for variant in self:
            for gene_name in variant.gene_names:
                counter[gene_name] += 1
        return counter

    @memoize
    def reference_names(self):
        """
        All distinct reference names used by Variants in this
        collection.
        """
        return set(variant.reference_name for variant in self)

    def groupby_gene_name(self):
        """
        Group variants by the gene names they overlap, which may put each
        variant in multiple groups.
        """
        return self.multi_groupby(lambda x: x.gene_names)

    def groupby_gene_id(self):
        return self.multi_groupby(lambda x: x.gene_ids)

    def detailed_string(self):
        lines = []
        gene_groups = self.groupby_gene_id()
        for gene_id, variant_group in sorted(gene_groups.items()):
            # in case we are combining variants with different Ensembl releases
            # use them all to gather gene names
            gene_names = {
                variant.ensembl.gene_name_of_gene_id(gene_id)
                for variant in variant_group
            }
            gene_name_str = ", ".join(sorted(gene_names))
            lines.append("  -- %s (%s):" % (gene_name_str, gene_id))
            for variant in variant_group:
                lines.append("  ---- %s" % variant)
        header = self.short_string()
        joined_lines = "\n".join(lines)
        return "%s\n%s" % (header, joined_lines)

    def clone_with_new_elements(self, new_elements):
        """
        Create another VariantCollection of the same class and with
        same state (including metadata) but possibly different entries.

        Warning: metadata is a dictionary keyed by variants. This method
        leaves that dictionary as-is, which may result in extraneous entries
        or missing entries.
        """
        return self.__class__(
            new_elements,
            source_to_metadata_dict=self.source_to_metadata_dict,
            distinct=self.distinct,
            sort_key=self.sort_key)

    def filter_by_transcript_expression(
            self,
            transcript_expression_dict,
            min_expression_value=0.0):
        """
        Filters variants down to those which have overlap a transcript whose
        expression value in the transcript_expression_dict argument is greater
        than min_expression_value.

        Parameters
        ----------
        transcript_expression_dict : dict
            Dictionary mapping Ensembl transcript IDs to expression estimates
            (either FPKM or TPM)

        min_expression_value : float
            Threshold above which we'll keep an effect in the result collection
        """
        return self.filter_any_above_threshold(
            multi_key_fn=lambda variant: variant.transcript_ids,
            value_dict=transcript_expression_dict,
            threshold=min_expression_value)

    def filter_by_gene_expression(
            self,
            gene_expression_dict,
            min_expression_value=0.0):
        """
        Filters variants down to those which have overlap a gene whose
        expression value in the transcript_expression_dict argument is greater
        than min_expression_value.

        Parameters
        ----------
        gene_expression_dict : dict
            Dictionary mapping Ensembl gene IDs to expression estimates
            (either FPKM or TPM)

        min_expression_value : float
            Threshold above which we'll keep an effect in the result collection
        """
        return self.filter_any_above_threshold(
            multi_key_fn=lambda effect: effect.gene_ids,
            value_dict=gene_expression_dict,
            threshold=min_expression_value)

    def exactly_equal(self, other):
        '''
        Comparison between VariantCollection instances that takes into account
        the info field of Variant instances.

        Returns
        ----------
        True if the variants in this collection equal the variants in the other
        collection. The Variant.info fields are included in the comparison.
        '''
        return (
            self.__class__ == other.__class__ and
            len(self) == len(other) and
            all(x.exactly_equal(y) for (x, y) in zip(self, other)))

    @classmethod
    def _merge_metadata_dictionaries(cls, dictionaries):
        """
        Helper function for combining variant collections: given multiple
        dictionaries mapping:
             source name -> (variant -> (attribute -> value))

        Returns dictionary with union of all variants and sources.
        """
        # three levels of nested dictionaries!
        #   {source name: {variant: {attribute: value}}}
        combined_dictionary = {}
        for source_to_metadata_dict in dictionaries:
            for source_name, variant_to_metadata_dict in source_to_metadata_dict.items():
                combined_dictionary.setdefault(source_name, {})
                combined_source_dict = combined_dictionary[source_name]
                for variant, metadata_dict in variant_to_metadata_dict.items():
                    combined_source_dict.setdefault(variant, {})
                    combined_source_dict[variant].update(metadata_dict)
        return combined_dictionary

    @classmethod
    def _combine_variant_collections(cls, combine_fn, variant_collections, kwargs):
        """
        Create a single VariantCollection from multiple different collections.

        Parameters
        ----------

        cls : class
            Should be VariantCollection

        combine_fn : function
            Function which takes any number of sets of variants and returns
            some combination of them (typically union or intersection).

        variant_collections : tuple of VariantCollection

        kwargs : dict
            Optional dictionary of keyword arguments to pass to the initializer
            for VariantCollection.
        """
        kwargs["variants"] = combine_fn(*[set(vc) for vc in variant_collections])
        kwargs["source_to_metadata_dict"] = cls._merge_metadata_dictionaries(
            [vc.source_to_metadata_dict for vc in variant_collections])
        for key, value in variant_collections[0].to_dict().items():
            # If some optional parameter isn't explicitly specified as an
            # argument to union() or intersection() then use the same value
            # as the first VariantCollection.
            #
            # I'm doing this so that the meaning of VariantCollection.union
            # and VariantCollection.intersection with a single argument is
            # the identity function (rather than setting optional parameters
            # to their default values.
            if key not in kwargs:
                kwargs[key] = value
        return cls(**kwargs)

    def union(self, *others, **kwargs):
        """
        Returns the union of variants in a several VariantCollection objects.
        """
        return self._combine_variant_collections(
            combine_fn=set.union,
            variant_collections=(self,) + others,
            kwargs=kwargs)

    def intersection(self, *others, **kwargs):
        """
        Returns the intersection of variants in several VariantCollection objects.
        """
        return self._combine_variant_collections(
            combine_fn=set.intersection,
            variant_collections=(self,) + others,
            kwargs=kwargs)

    def to_dataframe(self):
        """Build a DataFrame from this variant collection"""

        def row_from_variant(variant):
            return OrderedDict([
                ("chr", variant.contig),
                ("start", variant.original_start),
                ("ref", variant.original_ref),
                ("alt", variant.original_alt),
                ("gene_name", ";".join(variant.gene_names)),
                ("gene_id", ";".join(variant.gene_ids))
            ])
        rows = [row_from_variant(v) for v in self]
        if len(rows) == 0:
            # TODO: return a DataFrame with the appropriate columns
            return pd.DataFrame()
        return pd.DataFrame.from_records(rows, columns=rows[0].keys())
