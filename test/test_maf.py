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
from __future__ import absolute_import

from nose.tools import eq_
from pyensembl import ensembl_grch37 as ensembl

from varcode import Variant, load_maf, load_maf_dataframe

import pandas as pd

from .data import tcga_ov_variants, ov_wustle_variants, data_path

def test_maf():
    expected_tcga_ov_variants = [
        Variant(1, 1650797, "A", "G", ensembl),
        Variant(1, 23836447, "C", "A", ensembl),
        Variant(1, 231401797, "A", "C", ensembl),
        Variant(11, 124617502, "C", "G", ensembl),
    ]
    eq_(len(tcga_ov_variants), len(expected_tcga_ov_variants))
    for v_expect, v_maf in zip(expected_tcga_ov_variants, tcga_ov_variants):
        eq_(v_expect, v_maf)
        gene_name = tcga_ov_variants.metadata[v_maf]['Hugo_Symbol']
        assert any(gene.name == gene_name for gene in v_maf.genes), \
            "Expected gene name %s but got %s" % (gene_name, v_maf.genes)

def check_same_aa_change(variant, expected_aa_change):
    effect = variant.effects().top_priority_effect()
    change = effect.short_description
    eq_(
        change,
        expected_aa_change,
        "MAF file had annotation %s but Varcode gave %s" % (
            expected_aa_change, change))

def test_maf_aa_changes():
    # Parse a MAF file and make sure we're annotating the protein amino acid
    # changes in the same way.
    #
    # The data file used also contains spaces, which is good to test the parser
    # on.
    assert len(ov_wustle_variants) == 5

    expected_changes = {}
    # pylint: disable=no-member
    # pylint gets confused by read_csv
    maf_fields = pd.read_csv(
        ov_wustle_variants.path,
        sep="\t",
        comment="#")
    for _, row in maf_fields.iterrows():
        key = (str(row.Chromosome), row.Start_position)
        change = row.amino_acid_change
        # silent mutations just specificy which amino acid they affect via
        # e.g. "p.G384"
        if change[-1].isdigit():
            expected_changes[key] = "silent"
        else:
            expected_changes[key] = change

    for variant in ov_wustle_variants:
        key = (variant.contig, variant.start)
        expected = expected_changes[key]
        yield (check_same_aa_change, variant, expected)

def test_maf_number_entries_duplicates():
    # There are 3 duplicated mutations listed in the MAF
    path_to_maf_with_duplicates = data_path("duplicates.maf")
    variants = load_maf(path_to_maf_with_duplicates, distinct=True)
    assert len(variants) == 1
    variants = load_maf(path_to_maf_with_duplicates, distinct=False)
    assert len(variants) == 3

def test_load_maf():
    for raise_on_error in [True, False]:
        variants = load_maf(
            data_path("ov.wustle.subset5.maf"), raise_on_error=raise_on_error)
        eq_(len(variants), 5)


def test_load_maf_dataframe():
    for raise_on_error in [True, False]:
        variants_df = load_maf_dataframe(
            data_path("ov.wustle.subset5.maf"), raise_on_error=raise_on_error)
        eq_(len(variants_df), 5)


def test_xy_contigs():
    """
    Test MAFs with X and Y chromosomes rather than just numerical chromosomes.
    """
    for raise_on_error in [True, False]:
        variants = load_maf(
            data_path("tcga_ov.head.xychr.maf"), raise_on_error=True)
        eq_(len(variants), 4)


def test_load_utf8():
    """
    Test MAFs loaded with utf-8 encoding.
    """
    for raise_on_error in [True, False]:
        variants = load_maf(
            data_path("ov.wustle.subset5.maf"), raise_on_error=True, encoding="utf-8")
        eq_(len(variants), 5)
        # Make sure we avoid "TypeError: character mapping must return integer, None or unicode"
        # from Bio.Seq.
        _ = variants.effects()
