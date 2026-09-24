"""CLI input/diagnostic regressions; synthetic VCF construction lives here."""

import pandas as pd
import pytest
from pathlib import Path

from varcode import EffectCollection, Variant
from varcode.cli.effects_script import main as effects_main
from varcode.cli.genes_script import main as genes_main
from varcode.cli.variant_args import make_variants_parser, variant_collection_from_args
from varcode.effects import Failure


@pytest.fixture
def mixed_vcf(tmp_path):
    path = tmp_path / "mixed.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
        '##FILTER=<ID=REJECT,Description="Rejected">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr12\t25245350\tpass\tC\tT\t.\tPASS\t.\n"
        "chr12\t25245349\tfiltered\tC\tA,G\t.\tREJECT\t.\n"
        "chr12\t25245340\tsv\tN\t<DEL>\t.\t.\tEND=25245342\n")
    return str(path)


def test_filter_record_count_and_cli_sv_default(mixed_vcf):
    parser = make_variants_parser()
    args = ["--genome", "GRCh38", "--vcf", mixed_vcf]
    with pytest.warns(UserWarning, match=r"Skipped 1 VCF record.*--include-filtered"):
        variants = variant_collection_from_args(parser.parse_args(args))
    assert len(variants) == 2
    assert sum(getattr(v, "is_structural", False) for v in variants) == 1
    variants = variant_collection_from_args(parser.parse_args(args + ["--include-filtered"]))
    assert len(variants) == 4
    assert {v.contig for v in variants} == {"12"}


@pytest.mark.parametrize("main", [effects_main, genes_main])
def test_all_esvee_structural_records_reach_output(main, tmp_path):
    output = tmp_path / "sv.csv"
    main(["--genome", "GRCh38", "--vcf", str(Path(__file__).parent / "data/osteosarc_esvee_somatic.vcf"),
          "--output-csv", str(output)])
    table = pd.read_csv(output)
    assert len(table[["sv_type", "start", "end", "mate_contig", "mate_start"]].drop_duplicates()) == 8


@pytest.mark.parametrize("main", [effects_main, genes_main])
def test_skip_errors_retains_bad_contig_and_valid_chr_variant(main, tmp_path):
    output = tmp_path / "result.csv"
    args = ["--genome", "GRCh38", "--variant", "chr12", "25245350", "C", "T",
            "--variant", "99", "150", "A", "T", "--output-csv", str(output)]
    with pytest.raises(ValueError, match="Invalid contig"):
        main(args)
    if main is effects_main:
        args += ["--only-coding", "--one-per-variant"]
    main(args + ["--skip-errors"])
    table = pd.read_csv(output)
    status = "effect_type" if main is effects_main else "annotation_status"
    failed = table[table[status] == "Failure"]
    assert len(failed) == 1
    assert "Invalid contig" in failed.iloc[0].annotation_error
    assert "KRAS" in set(table.gene_name)
    if main is effects_main:
        assert "p.G12D" in set(table.effect)


def test_reference_mismatch_keeps_failures_with_coding_and_priority_filters(tmp_path):
    output = tmp_path / "mismatch.csv"
    effects_main([
        "--genome", "GRCh38", "--variant", "chr12", "25245350", "A", "T",
        "--only-coding", "--one-per-variant", "--skip-errors", "--output-csv", str(output)])
    table = pd.read_csv(output)
    failures = table[table.effect_type == "Failure"]
    assert len(failures) > 1  # keep each transcript's error, not just one
    assert failures.annotation_error.notna().all()


def test_transcript_free_failure_json_round_trip():
    variant = Variant("99", 150, "A", "T", genome=81)
    effects = EffectCollection([Failure(variant, error="unknown contig")])
    restored = EffectCollection.from_json(effects.to_json())
    assert restored[0].variant == variant
    assert restored[0].transcript is None
    assert restored[0].error == "unknown contig"
