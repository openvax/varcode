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

from .cufflinks import (
    load_cufflinks_dataframe,
    load_cufflinks_fpkm_dict,
    remap_novel_gene_expression_onto_ensembl_ids
)


class ExpressionFilter(object):
    """
    Base class for callable expression filter objects
    """
    def __init__(
            self,
            fpkm_dict,
            min_fpkm=0.0):
        self.fpkm_dict = fpkm_dict
        self.min_fpkm = min_fpkm

    def any_expressed(self, feature_ids):
        return any(
            self.fpkm_dict.get(feature_id, 0.0) > self.min_fpkm
            for feature_id in feature_ids)

class VariantTranscriptExpressionFilter(ExpressionFilter):
    def __call__(self, variant):
        return self.any_expressed(variant.transcript_ids())

class VariantGeneExpressionFilter(ExpressionFilter):
    def __call__(self, variant):
        return self.any_expressed(variant.gene_ids())

class EffectTranscriptExpressionFilter(ExpressionFilter):
    def __call__(self, effect):
        return (
            effect.transcript_id is not None and
            self.fpkm_dict.get(effect.transcript_id, 0.0) > self.min_fpkm)

class EffectGeneExpressionFilter(ExpressionFilter):
    def __call__(self, effect):
        return (
            effect.gene_id is not None and
            self.fpkm_dict.get(effect.gene_id, 0.0) > self.min_fpkm)

def _load_gene_fpkm_dict(fpkm_path, remap_novel_genes):
    fpkm_df = load_cufflinks_dataframe(fpkm_path)
    # make a dictionary mapping gene IDs to FPKM either by pulling the
    # values directly from the dataframe or if we're remapping novel IDs
    # onto their overlapping Ensembl IDs, rely on
    # `remap_novel_gene_expression_onto_ensembl_ids` to return a dict
    if remap_novel_genes:
        return remap_novel_gene_expression_onto_ensembl_ids(fpkm_df)
    else:
        return {
            row.id: row.fpkm
            for (_, row) in fpkm_df.iterrows()
        }

def make_variant_gene_expression_filter(
        fpkm_path,
        min_fpkm=0.0,
        remap_novel_genes=False):
    return VariantGeneExpressionFilter(
        fpkm_dict=_load_gene_fpkm_dict(fpkm_path, remap_novel_genes),
        min_fpkm=min_fpkm)

def make_effect_gene_expression_filter(
        fpkm_path,
        min_fpkm=0.0,
        remap_novel_genes=False):
    return EffectGeneExpressionFilter(
        fpkm_dict=_load_gene_fpkm_dict(fpkm_path, remap_novel_genes),
        min_fpkm=min_fpkm)

def make_variant_transcript_expression_filter(fpkm_path, min_fpkm=0.0):
    fpkm_dict = load_cufflinks_fpkm_dict(fpkm_path)
    return VariantGeneExpressionFilter(fpkm_dict, min_fpkm)

def make_effect_transcript_expression_filter(fpkm_path, min_fpkm=0.0):
    fpkm_dict = load_cufflinks_fpkm_dict(fpkm_path)
    return EffectGeneExpressionFilter(fpkm_dict, min_fpkm)
