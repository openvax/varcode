"""Synthetic VCF fixtures for allele limits and source ALT provenance (#504)."""

import pandas as pd
import pytest

from varcode import load_vcf
from varcode.vcf import dataframes_to_variant_collection


@pytest.fixture
def mixed_vcf(tmp_path):
    """Construct the complete fixture here; no downloads or source cache."""
    path = tmp_path / "mixed.vcf"
    path.write_text(
        '##fileformat=VCFv4.2\n'
        '##reference=GRCh38\n'
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        '1\t50\tfiltered\tA\tC\t.\tLowQual\t.\n'
        '1\t60\tmissing\tA\t.\t.\tPASS\t.\n'
        '1\t70\tspanning\tA\t*\t.\tPASS\t.\n'
        '1\t100\tmixed\tA\tC,<DEL>,G\t.\tPASS\tEND=150\n'
        '1\t200\tdup\tA\t<DUP>\t.\tPASS\tEND=250\n'
        '1\t300\tlast\tA\tT\t.\tPASS\t.\n'
    )
    return str(path)


@pytest.mark.parametrize("chunk_size", [1, 2, 100])
@pytest.mark.parametrize("limit", [0, 1, 2, 3, 4, 5, 20, None])
def test_limit_counts_loaded_alleles_across_chunks(mixed_vcf, chunk_size, limit):
    variants = load_vcf(
        mixed_vcf, genome=81, max_variants=limit, chunk_size=chunk_size,
        parse_structural_variants=True, sort_key=None, distinct=False)
    expected = [("mixed", 0), ("mixed", 1), ("mixed", 2), ("dup", 0), ("last", 0)]
    if limit is not None:
        expected = expected[:limit]
    assert [(variants.metadata[v]["id"], variants.metadata[v]["alt_allele_index"])
            for v in variants] == expected
    assert len(variants.metadata) == len(expected)


@pytest.mark.parametrize("limit", [-1, 1.5, "1"])
def test_invalid_limit_is_rejected(mixed_vcf, limit):
    with pytest.raises(ValueError, match="max_variants"):
        load_vcf(mixed_vcf, genome=81, max_variants=limit)


def _frame(alt="C,G", info=False):
    columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
    row = ["1", 100, ".", "A", alt, ".", "PASS"]
    if info:
        columns.append("INFO")
        row.append(".")
    return pd.DataFrame([row], columns=columns)


def test_limit_does_not_consume_a_later_chunk():
    def chunks():
        yield _frame()
        pytest.fail("loader consumed a chunk after reaching its allele limit")

    variants = dataframes_to_variant_collection(
        chunks(), "synthetic.vcf", max_variants=1, variant_kwargs={"genome": 81})
    assert len(variants) == 1
    assert variants[0].alt == "C"


def test_zero_limit_does_not_consume_chunks():
    def chunks():
        pytest.fail("zero limit consumed input")
        yield _frame()

    assert len(dataframes_to_variant_collection(
        chunks(), "synthetic.vcf", max_variants=0)) == 0


def test_parser_stop_iteration_is_not_silently_swallowed():
    def broken_parser(info):
        raise StopIteration("parser failed")

    with pytest.raises(StopIteration, match="parser failed"):
        dataframes_to_variant_collection(
            [_frame(info=True)], "synthetic.vcf", info_parser=broken_parser,
            variant_kwargs={"genome": 81})


def test_invalid_dataframe_columns_raise_value_error():
    with pytest.raises(ValueError, match="columns"):
        dataframes_to_variant_collection([pd.DataFrame()], "synthetic.vcf")


def test_limit_is_applied_before_deduplication():
    variants = dataframes_to_variant_collection(
        [_frame("C,C,G")], "synthetic.vcf", max_variants=2,
        variant_kwargs={"genome": 81})
    assert len(variants) == 1
    assert variants[0].alt == "C"


def test_structural_alt_index_drives_sample_selection(tmp_path):
    path = tmp_path / "multiallelic-sv.vcf"
    path.write_text(
        '##fileformat=VCFv4.2\n'
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor\n'
        '1\t100\tmixed\tA\tC,<DEL>\t.\tPASS\tEND=150\tGT\t0/2\n'
    )
    variants = load_vcf(str(path), genome=81, parse_structural_variants=True)
    selected = variants.for_sample("tumor")
    assert len(selected) == 1
    assert selected[0].sv_type == "DEL"
