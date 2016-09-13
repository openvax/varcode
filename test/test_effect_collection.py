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

"""
Test properties of EffectCollection
"""

from nose.tools import eq_
from varcode.effects import IncompleteTranscript, Substitution
from .data import tcga_ov_variants, ov_wustle_variants
tcga_ov_effects = tcga_ov_variants.effects()
ov_wustle_effects = ov_wustle_variants.effects()

def test_to_dataframe():
    df = tcga_ov_effects.to_dataframe()
    eq_(len(tcga_ov_effects), len(df))

def test_effect_collection_gene_counts():
    # test that each gene is counted just once
    for gene, count in ov_wustle_effects.gene_counts().items():
        assert count > 1, \
            "Expected more than 1 effect for %s (got %d)" % (gene, count)

def test_effect_collection_groupby_gene():
    genes = ov_wustle_effects.groupby_gene().keys()
    # make sure that the IDs attached to Gene objects are the same as IDs
    # of groupby_gene_id
    gene_ids = set(ov_wustle_effects.groupby_gene_id().keys())
    eq_({gene.id for gene in genes}, gene_ids)

def test_effect_collection_groupby_gene_id():
    gene_ids = set(ov_wustle_effects.groupby_gene_id().keys())
    eq_(gene_ids, {
        'ENSG00000060718',
        'ENSG00000156876',
        'ENSG00000130939',
        'ENSG00000122477',
        'ENSG00000162688'
    })

def test_effect_collection_groupby_gene_name():
    gene_names = set(ov_wustle_effects.groupby_gene_name().keys())
    eq_(gene_names, {"AGL", "SASS6", "LRRC39", "UBE4B", "COL11A1"})

def test_effect_collection_groupby_variant():
    variants = set(ov_wustle_effects.groupby_variant().keys())
    # make sure that all the original variants are still present
    # in the group keys
    eq_(variants, set(ov_wustle_variants))

def test_effect_collection_filter_by_effect_priority():
    # every effect should be at least the same priority as "incomplete"
    eq_(
        tcga_ov_effects,
        tcga_ov_effects.filter_by_effect_priority(IncompleteTranscript))
    assert len(tcga_ov_effects) > len(
        tcga_ov_effects.filter_by_effect_priority(Substitution))

def test_effect_collection_drop_silent_and_noncoding():
    # some of the predicted effects are non-coding so should get dropped
    assert len(tcga_ov_effects) > len(tcga_ov_effects.drop_silent_and_noncoding())
