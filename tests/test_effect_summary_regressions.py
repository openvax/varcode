"""Summaries preserve splice alternatives and omit reference-identical VCF ALTs."""

from pathlib import Path

import pandas as pd
import pytest

from varcode import EffectCollection, SpliceOutcomeSet, Variant, load_vcf
from varcode.cli.effects_script import main as effects_main
from varcode.effects import ExonicSpliceSite, effect_priority

# Legacy wrapper contract; the default composition path has separate regression tests.
pytestmark = pytest.mark.usefixtures("legacy_annotator")


@pytest.fixture
def splice_set():
    variant = Variant("1", 53513530, "A", "C", genome=75)
    transcript = variant.genome.transcript_by_id("ENST00000371514")
    return variant.effect_on_transcript(transcript)


def test_priority_selection_retains_splice_alternatives_without_materializing(splice_set):
    collection = EffectCollection([splice_set])
    assert splice_set._candidates_cache is None
    assert collection.top_priority_effect() is splice_set
    assert collection.top_priority_effect_per_variant()[splice_set.variant] is splice_set
    assert collection.top_priority_effect_per_transcript_id()[splice_set.transcript_id] is splice_set
    assert collection.top_priority_effect_per_gene_id()[splice_set.gene_id] is splice_set
    assert effect_priority(splice_set) == effect_priority(splice_set.effect_if_splicing_unchanged)
    assert splice_set._candidates_cache is None


def test_splice_dataframe_and_csv_show_every_candidate(splice_set, tmp_path):
    assert isinstance(splice_set, SpliceOutcomeSet)
    assert splice_set.disrupted_signal_class is ExonicSpliceSite
    collection = EffectCollection([splice_set])
    frame = collection.to_dataframe()
    path = tmp_path / "effects.csv"
    collection.to_csv(path)
    csv = pd.read_csv(path, comment="#")
    for table in (frame, csv):
        row = table.iloc[0]
        assert row.effect_type == "SpliceOutcomeSet"
        assert row.splice_signal == "ExonicSpliceSite"
        assert "p.D490A" in row.effect
        assert "p.D490fs" in row.effect
        for candidate in splice_set.candidates:
            assert candidate.effect.short_description in row.effect
    restored = EffectCollection.from_csv(path)
    assert restored[0].short_description == splice_set.short_description


def test_one_per_variant_cli_keeps_splice_alternatives(tmp_path):
    path = tmp_path / "one.csv"
    effects_main(["--genome", "GRCh37:75", "--variant", "1", "53513530", "A", "C",
                  "--one-per-variant", "--only-coding", "--output-csv", str(path)])
    table = pd.read_csv(path)
    assert len(table) == 1
    row = table.iloc[0]
    assert row.effect_type == "SpliceOutcomeSet"
    assert row.splice_signal == "ExonicSpliceSite"
    assert "normal-splicing (p.D490A)" in row.effect
    assert "p.D490fs" in row.effect
    assert "cryptic-donor (predicted)" in row.effect
    assert "intron-retention (donor side, predicted)" in row.effect


def test_reference_identical_alts_do_not_consume_limit_or_change_allele_indices(tmp_path):
    path = tmp_path / "no-change.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t100\t.\tAAA\tAAA\t.\tPASS\t.\n"
        "1\t200\t.\tA\ta,G\t.\tPASS\t.\n"
        "1\t300\t.\tT\tC\t.\tPASS\t.\n")
    with pytest.warns(UserWarning, match="Skipped 2 reference-identical ALT allele"):
        variants = load_vcf(str(path), genome=75, max_variants=1)
    variant, = variants
    assert (variant.start, variant.ref, variant.alt) == (200, "A", "G")
    metadata = variants.metadata[variant]
    assert metadata["alt_allele_index"] == 1
    assert metadata["vcf_alt_alleles"] == ("a", "G")


def test_no_change_vcf_record_is_not_exported_as_silent(tmp_path):
    path = Path(__file__).parent / "data/somatic_hg19_14muts.vcf"
    with pytest.warns(UserWarning, match="Skipped 1 reference-identical ALT allele"):
        variants = load_vcf(str(path))
    assert len(variants) == 13
    assert all(v.original_ref != v.original_alt for v in variants)
    table = variants.effects().to_dataframe()
    assert not ((table.contig == "11") & (table.start == 118244286)).any()
    assert "p.335=" not in set(table.effect)


@pytest.mark.parametrize("start,ref,alt,expected", [
    (100, "C", "", "chr1 g.100delC"),
    (100, "AC", "A", "chr1 g.101delC"),
    (100, "AC", "", "chr1 g.100_101delAC"),
    (100, "G", "G", "chr1 g.100="),
    (100, "AG", "AG", "chr1 g.100_101="),
])
def test_genomic_descriptions_use_single_positions_and_explicit_equality(start, ref, alt, expected):
    assert Variant("1", start, ref, alt, genome=75).short_description == expected
