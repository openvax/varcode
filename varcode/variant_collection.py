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
import json
from collections import Counter

from .collection import Collection
from .effect_collection import EffectCollection
from .common import memoize
from . import Variant

class VariantCollection(Collection):
    def __init__(
            self,
            variants,
            path=None,
            distinct=True,
            sort_key=None,
            metadata=None):
        """
        Construct a VariantCollection from a list of Variant records.

        Parameters
        ----------
        variants : iterable
            Variant objects contained in this VariantCollection

        path : str, optional
            File path from which we loaded variants, though the current
            VariantCollection may only contain a subset of them.

        distinct : bool
            Don't keep repeated variants

        sort_key : callable
        """
        if sort_key is None:
            # pylint: disable=function-redefined
            def sort_key(variant):
                return (variant.contig, variant.start)

        Collection.__init__(
            self,
            elements=variants,
            path=path,
            distinct=distinct,
            sort_key=sort_key)
        self.metadata = {} if metadata is None else metadata

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

    def to_json(self):
        """
        Returns a string giving this variant collection serialized in JSON
        format.
        """
        return json.dumps(self.__getstate__())

    def write_json_file(self, filename):
        """
        Serialize this VariantCollection to a JSON representation and write it
        out to a text file.
        """
        with open(filename, "w") as f:
            f.write(self.to_json())

    @classmethod
    def from_dict(cls, state):
        """
        Deserialize a VariantCollection encoded as a Python dict.
        """
        instance = cls.__new__(cls)
        instance.__setstate__(state)
        return instance

    @classmethod
    def from_json(cls, serialized):
        """
        Construct a VariantCollection from a JSON string.
        """
        return cls.from_dict(json.loads(serialized))

    @classmethod
    def read_json_file(cls, filename):
        """
        Construct a VariantCollection from a JSON file.
        """
        with open(filename, 'r') as f:
            json_string = f.read()
        return cls.from_json(json_string)

    def __getstate__(self):
        result = dict(self.__dict__)
        result['elements'] = [
            (e.__getstate__(), self.metadata.get(e)) for e in self]
        del result['metadata']
        return result

    def __setstate__(self, state):
        self.elements = []
        self.metadata = {}
        for (variant_data, meta_entry) in state.pop('elements'):
            variant = Variant.from_dict(variant_data)
            self.elements.append(variant)
            if meta_entry is not None:
                self.metadata[variant] = meta_entry
        self.__dict__.update(state)

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
