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

import pytest 

from pyensembl import ensembl_grch37 as ensembl

from varcode import SpliceOutcomeSet, Variant, load_maf, load_maf_dataframe

import pandas as pd

from .common import eq_
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


def generate_maf_aa_changes():
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
        # silent mutations in this MAF specify only the affected residue
        # (e.g. "p.G384") — varcode renders these in HGVS notation as
        # "p.G384=" per https://github.com/openvax/varcode/issues/217.
        if change[-1].isdigit():
            expected_changes[key] = change + "="
        else:
            expected_changes[key] = change

    for variant in ov_wustle_variants:
        key = (variant.contig, variant.start)
        expected = expected_changes[key]
        yield (variant, expected)

@pytest.mark.parametrize(['variant', 'expected_aa_change'], list(generate_maf_aa_changes()))
def test_maf_aa_changes(variant, expected_aa_change):
    effect = variant.effects().top_priority_effect()
    if isinstance(effect, SpliceOutcomeSet):
        # The MAF annotation describes the codon change under normal splicing.
        effect = effect.effect_if_splicing_unchanged
    change = effect.short_description
    eq_(
        change,
        expected_aa_change,
        "MAF file had annotation %s but Varcode gave %s" % (
            expected_aa_change, change))


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
        # Regression: annotating this utf-8 MAF once raised "TypeError: character mapping
        # must return integer, None or unicode" in the Biopython translation Varcode used then.
        _ = variants.effects()


def test_maf_column_case_normalization_preserves_order(tmp_path):
    from varcode.maf import MAF_COLUMN_NAMES

    path = tmp_path / "lowercase.maf"
    columns = [MAF_COLUMN_NAMES[0].lower()] + MAF_COLUMN_NAMES[1:]
    path.write_text("\t".join(columns) + "\n")
    assert load_maf_dataframe(str(path)).columns.tolist() == MAF_COLUMN_NAMES


def _write_maf(path, columns, rows):
    path.write_text("\n".join("\t".join(row) for row in [columns, *rows]) + "\n")
    return str(path)


def _maf_row(build, chromosome, position, ref, alt):
    return ["GENE", "0", ".", build, chromosome, position, position, "+",
            "Missense_Mutation", "SNP", ref, ref, alt, "novel", "", "T", "N", ref, ref]


def test_gdc_chr_prefixed_contigs_are_annotated_as_ensembl_contigs(tmp_path):
    # GDC GRCh38 MAFs write "chr7" and "chrM"; Ensembl names these "7" and "MT" (#680).
    from varcode.maf import MAF_COLUMN_NAMES

    maf = _write_maf(tmp_path / "gdc.maf", MAF_COLUMN_NAMES,
                     [_maf_row("GRCh38", "chr7", "140753336", "A", "T"),
                      _maf_row("GRCh38", "chrM", "6000", "C", "T")])
    braf, mt = load_maf(maf)
    assert (braf.contig, mt.contig) == ("7", "MT")
    assert braf == Variant("7", 140753336, "A", "T", "GRCh38")
    assert braf.effects().top_priority_effect().short_description == "p.V600E"
    assert mt.genes and mt.genes[0].contig == "MT"


def test_missing_required_columns_are_named_in_spec_order(tmp_path):
    # The eight columns someone might write by hand for one variant (#681).
    columns = ["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
               "Reference_Allele", "Tumor_Seq_Allele2", "NCBI_Build", "Variant_Classification"]
    maf = _write_maf(tmp_path / "in.maf", columns, [])
    message = (
        "is missing 11 of the 19 required columns: Entrez_Gene_Id, Center, Strand, "
        "Variant_Type, Tumor_Seq_Allele1, dbSNP_RS, dbSNP_Val_Status, Tumor_Sample_Barcode, "
        "Matched_Norm_Sample_Barcode, Match_Norm_Seq_Allele1, Match_Norm_Seq_Allele2 "
        "(found: %s)" % ", ".join(columns))
    # raise_on_error=False skips unreadable rows; without these columns none can be read (#708).
    for load in (load_maf_dataframe, load_maf):
        for raise_on_error in (True, False):
            with pytest.raises(ValueError) as caught:
                load(maf, raise_on_error=raise_on_error)
            assert str(caught.value) == "MAF file %s %s" % (maf, message)


@pytest.mark.parametrize("build", ["37", "GRCh37", "hg19"])
@pytest.mark.parametrize("contig", ["chrM", "M"])
def test_grch37_chrm_rows_keep_their_name_and_fail_annotation(tmp_path, build, contig):
    # hg19's chrM is not the rCRS that Ensembl's GRCh37 names MT, and its
    # positions are offset (hg19 chrM:12995 is rCRS MT:12994), #707.
    from varcode.errors import InvalidContigError
    from varcode.maf import MAF_COLUMN_NAMES

    maf = _write_maf(tmp_path / "hg19.maf", MAF_COLUMN_NAMES,
                     [_maf_row(build, contig, "12995", "A", "G"),
                      _maf_row(build, "MT", "12994", "A", "G")])
    hg19_mt, rcrs_mt = sorted(load_maf(maf), key=lambda v: v.contig != contig)
    assert (hg19_mt.contig, rcrs_mt.contig) == (contig, "MT")
    with pytest.raises(InvalidContigError, match="NC_001807"):
        hg19_mt.effects()


def test_empty_cells_are_none_and_the_collection_round_trips_through_json(tmp_path):
    # pandas reads an empty cell as NaN, which JSON can't hold (#704). GDC MAFs
    # leave dbSNP_Val_Status empty; some writers leave Tumor_Seq_Allele1 empty.
    from varcode import VariantCollection
    from varcode.maf import MAF_COLUMN_NAMES

    row = _maf_row("GRCh38", "chr7", "140753336", "A", "T")
    row[MAF_COLUMN_NAMES.index("dbSNP_RS")] = ""
    row[MAF_COLUMN_NAMES.index("Tumor_Seq_Allele1")] = ""
    maf = _write_maf(tmp_path / "gdc.maf", MAF_COLUMN_NAMES + ["t_alt_count", "vaf"],
                     [row + ["12", ""]])
    variants = load_maf(maf, optional_cols=["t_alt_count", "vaf"])
    variant, = variants
    assert variant.alt == "T"
    metadata = variants.metadata[variant]
    assert metadata["dbSNP_RS"] is None and metadata["dbSNP_Val_Status"] is None
    assert metadata["vaf"] is None and metadata["t_alt_count"] == 12
    restored = VariantCollection.from_json(variants.to_json())
    assert restored == variants
    assert restored.metadata[variant] == metadata
    # The bundled TCGA fixture has empty cells too.
    assert VariantCollection.from_json(tcga_ov_variants.to_json()) == tcga_ov_variants
