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
try:
    import cPickle as pickle
except ImportError:
    import pickle

from .data import ov_wustle_variants, tcga_ov_variants

from varcode import VariantCollection, Variant

def test_reference_names():
    eq_(ov_wustle_variants.reference_names(), {"GRCh37"})

def test_to_string():
    string_repr = str(ov_wustle_variants)
    assert "start=10238758, ref='G', alt='C'" in string_repr, \
        "Expected variant g.10238758 G>C in __str__:\n%s" % (
            string_repr,)

def test_detailed_string():
    detailed_string = ov_wustle_variants.detailed_string()
    # expect one of the gene names from the MAF to be in the summary string
    assert "UBE4B" in detailed_string, \
        "Expected gene name UBE4B in detailed_string():\n%s" % detailed_string
    assert "start=10238758, ref='G', alt='C'" in detailed_string, \
        "Expected variant g.10238758 G>C in detailed_string():\n%s" % (
            detailed_string,)

def test_gene_counts():
    expected_coding_gene_counts = Counter()
    expected_coding_gene_counts["CDK11A"] = 1
    expected_coding_gene_counts["GNPAT"] = 1
    expected_coding_gene_counts["E2F2"] = 1
    expected_coding_gene_counts["VSIG2"] = 1
    all_gene_counts = tcga_ov_variants.gene_counts()
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

def test_serialization():
    original = VariantCollection(
        [
            Variant(
                1, start=10, ref="AA", alt="AAT", ensembl=77),
            Variant(10, start=15, ref="A", alt="G"),
            Variant(20, start=150, ref="", alt="G"),
        ])
    original.metadata[original[0]] = {"a": "b"}
    original.metadata[original[2]] = {"bar": 2}

    # This causes the variants' ensembl objects to make a SQL connection,
    # which makes the ensembl object non-serializable. By calling this
    # method, we are checking that we don't attempt to directly serialize
    # the ensembl object.
    original.effects()

    # Test pickling.
    serialized = pickle.dumps(original)
    reconstituted = pickle.loads(serialized)
    eq_(original, reconstituted)
    eq_(reconstituted[0], original[0])
    eq_(reconstituted.metadata[original[0]], original.metadata[original[0]])

    # Test json.
    serialized = original.to_json()
    reconstituted = VariantCollection.from_json(serialized)
    eq_(original, reconstituted)
    eq_(reconstituted[0], original[0])
    eq_(reconstituted.metadata[original[0]], original.metadata[original[0]])
