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

from nose.tools import eq_
from varcode import VariantCollection

from .data import (
    snp_rs4244285,
    snp_rs1537415
)

variants = VariantCollection([
    # gene ids: ['ENSG00000165841', 'ENSG00000276490']
    # transcript_ids : ['ENST00000371321', 'ENST00000464755']
    snp_rs4244285,
    # gene ids: ['ENSG00000204007']
    # transcript ids:  ['ENST00000371763', 'ENST00000613244']
    snp_rs1537415,
])

gene_fpkm_dict = {
    "ENSG00000165841": 10.0,
    "ENSG00000204007": 20.0,
    "ENSG00000276490": 30.0,
}

transcript_fpkm_dict = {
    "ENST00000371321": 10.0,
    "ENST00000464755": 20.0,
    "ENST00000371763": 30.0,
    "ENST00000613244": 40.0,
}

effects = variants.effects()

empty_variants = VariantCollection([])
empty_effects = empty_variants.effects()

def test_filter_variants():
    eq_(variants.filter(lambda _: True), variants)
    eq_(variants.filter(lambda _: False), empty_variants)

def test_filter_effects():
    eq_(effects.filter(lambda _: True), effects)
    eq_(effects.filter(lambda _: False), empty_effects)

def test_filter_variants_by_gene_expression():
    eq_(variants.filter_by_gene_expression(
        gene_fpkm_dict, 0.0), variants)
    eq_(variants.filter_by_gene_expression(
        gene_fpkm_dict, 100.0), empty_variants)

def test_filter_effects_by_gene_expression():
    eq_(effects.filter_by_gene_expression(
        gene_fpkm_dict, 0.0), effects)
    eq_(effects.filter_by_gene_expression(
        gene_fpkm_dict, 100.0), empty_effects)

def test_filter_variants_by_transcript_expression():
    expect_all = variants.filter_by_gene_expression(
        gene_fpkm_dict, 0.0)
    eq_(expect_all, variants)
    expect_none = variants.filter_by_gene_expression(
        gene_fpkm_dict, 100.0)
    eq_(expect_none, empty_variants)

def test_filter_effects_by_transcript_expression():
    expect_all = effects.filter_by_transcript_expression(
        transcript_fpkm_dict, 0.0)
    eq_(expect_all, effects)
    expect_none = effects.filter_by_transcript_expression(
        transcript_fpkm_dict, 100.0)
    eq_(expect_none, empty_effects)

def test_filter_silent_effects():
    # all dbSNP entries in the collection are silent
    assert len(effects.drop_silent_and_noncoding()) == 0
