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

"""
Test properties of VariantCollection objects other than effect annotations
"""
from collections import Counter
from nose.tools import eq_
from varcode import load_maf

from . import data_path

def test_reference_names():
    variants = load_maf(data_path("ov.wustle.subset5.maf"))
    eq_(variants.reference_names(), {'GRCh37'})

def test_summary_string():
    variants = load_maf(data_path("ov.wustle.subset5.maf"))
    summary_string = variants.summary_string()
    # expect one of the gene names from the MAF to be in the summary string
    assert "UBE4B" in summary_string, \
        "Expected gene name UBE4B in summary_string():\n%s" % summary_string
    assert "start=10238758, ref=G, alt=C" in summary_string, \
        "Expected variant g.10238758 G>C in summary_string():\n%s" % (
            summary_string,)

def test_gene_counts():
    variants = load_maf(data_path("tcga_ov.head.maf"))
    expected_coding_gene_counts = Counter()
    expected_coding_gene_counts["CDK11A"] = 1
    expected_coding_gene_counts["GNPAT"] = 1
    expected_coding_gene_counts["E2F2"] = 1
    expected_coding_gene_counts["VSIG2"] = 1
    all_gene_counts = variants.gene_counts()
    assert len(all_gene_counts) > len(expected_coding_gene_counts), \
        ("Gene counts for all genes must contain more elements than"
         " gene counts for only coding genes.")
    for (gene_name, count) in expected_coding_gene_counts.items():
        eq_(count, all_gene_counts[gene_name])

    # TODO: add `only_coding` parameter to gene_counts and then test
    # for exact equality between `coding_gene_counts` and
    # `expected_counts`
    #
    # coding_gene_counts = variants.gene_counts(only_coding=True)
    # eq_(coding_gene_counts, expected_counts)
