"""Explicit annotation dates retain their assembly and never imply downloads."""

import re

import pandas as pd
from pyensembl import EnsemblRelease
import pytest

from varcode import Variant, VariantCollection
from varcode.cli.effects_script import arg_parser, main
from varcode.cli.variant_args import variant_collection_from_args
from varcode.reference import infer_genome, infer_genome_for_reference_name


@pytest.mark.parametrize("name,assembly,date,species,accession,ucsc", [
    ("GRCh38:2026_04", "GRCh38", "2026_04", "human", "GCA_000001405.29", False),
    (" GRCh38 : 2026_04 ", "GRCh38", "2026_04", "human", "GCA_000001405.29", False),
    ("GRCh38.2026_04", "GRCh38", "2026_04", "human", "GCA_000001405.29", False),
    ("GRCh38_2026_04", "GRCh38", "2026_04", "human", "GCA_000001405.29", False),
    ("GRCh38-2026_04", "GRCh38", "2026_04", "human", "GCA_000001405.29", False),
    ("GRCh38 2026_04", "GRCh38", "2026_04", "human", "GCA_000001405.29", False),
    ("GRCh37:2013_09", "GRCh37", "2013_09", "human", "GCA_000001405.14", False),
    ("B37:2013_09", "GRCh37", "2013_09", "human", "GCA_000001405.14", False),
    ("hg19:2013_09", "GRCh37", "2013_09", "human", "GCA_000001405.14", True),
    ("ROS_Cfam_1.0:2020_10", "ROS_Cfam_1.0", "2020_10", "dog", "GCA_014441545.1", False),
])
def test_dated_selection_preserves_assembly_and_upstream_cache(
        name, assembly, date, species, accession, ucsc):
    expected = EnsemblRelease.cached(date, species=species, assembly=assembly)
    genome, converted = infer_genome(name)
    assert genome.release == date
    assert genome.assembly == assembly
    assert genome.assembly_accession == accession
    assert converted is ucsc
    # Historical dated datasets use accession-based names/cache directories.
    # Do not relabel them or create a second cache under the requested alias.
    assert genome is expected


@pytest.mark.parametrize("name,message", [
    ("GRCh37:2026_04", "provides GRCh38, not GRCh37"),
    ("GRCh38:2013_09", "provides GRCh37, not GRCh38"),
    ("CanFam3.1:2020_10", "not CanFam3.1"),
])
def test_dated_assembly_mismatch_is_rejected(name, message):
    with pytest.raises(ValueError, match=re.escape(message)):
        infer_genome(name)


@pytest.mark.parametrize("name", [
    "GRCh38:2026_4", "GRCh38:2026-04", "GRCh38:latest", "GRCh38:",
    "GRCh38:2026_00", "GRCh38:2026_13",
])
def test_invalid_explicit_date_cannot_fall_back_to_another_release(name):
    with pytest.raises(ValueError, match="release|annotation date"):
        infer_genome(name)


@pytest.mark.parametrize("name", ["GRCh38:2026_04", "GRCh37:2013_09"])
def test_dated_inference_is_offline(name, monkeypatch):
    def unexpected(*args, **kwargs):
        pytest.fail("Reference selection attempted a download")

    monkeypatch.setattr("requests.sessions.Session.request", unexpected)
    monkeypatch.setattr("urllib.request.urlopen", unexpected)
    monkeypatch.setattr(EnsemblRelease, "download", unexpected)
    monkeypatch.setattr(EnsemblRelease, "download_genome_fasta", unexpected)
    # Bypass Varcode's memoization so the resolver actually runs.
    genome, _ = infer_genome_for_reference_name.__wrapped__(name)
    assert genome.release == name.split(":")[1]


@pytest.mark.parametrize("name", ["GRCh38:2026_04", "GRCh37:2013_09"])
def test_dated_variant_json_preserves_reference_identity(name):
    variant = Variant("7", 100, "A", "T", genome=name)
    restored = Variant.from_json(variant.to_json())
    assert restored.genome.release == name.split(":")[1]
    assert restored.genome.assembly == name.split(":")[0]
    assert restored.genome.assembly_accession == variant.genome.assembly_accession
    assert restored.genome.download_cache.cache_directory_path == (
        variant.genome.download_cache.cache_directory_path)


def test_dated_object_pass_through_retains_attached_dna(tmp_path):
    fasta = tmp_path / "custom.fa"
    fasta.write_text(">99\nACGT\n")
    genome = EnsemblRelease.cached(
        "2020_10", species="dog", assembly="ROS_Cfam_1.0", genome_fasta=fasta,
        cache_directory_path=str(tmp_path / "cache"))
    selected, converted = infer_genome(genome)
    assert selected is genome and not converted
    assert selected.sequence("99", 1, 4) == "ACGT"


@pytest.mark.parametrize("name,date", [
    ("GRCh38:2026_04", "2026_04"), ("GRCh37:2013_09", "2013_09"),
])
def test_cli_retains_the_requested_annotation_date(name, date):
    args = arg_parser.parse_args(["--genome", name, "--variant", "7", "100", "A", "T"])
    variants = variant_collection_from_args(args)
    assert variants[0].genome.release == date
    assert variants[0].genome.assembly == name.split(":")[0]


def test_cli_dated_dog_annotation_matches_native_ros_call(tmp_path):
    path = tmp_path / "effects.csv"
    main(["--genome", "ROS_Cfam_1.0:2020_10",
          "--variant", "37", "14828457", "G", "T", "--output-csv", str(path)])
    rows = pd.read_csv(path)
    assert "p.G1780V" in set(rows.loc[rows.gene_name == "ZDBF2", "effect"])


def test_cli_dated_mismatch_fails_before_annotation(monkeypatch, capsys):
    def unexpected(*args, **kwargs):
        pytest.fail("Assembly mismatch reached annotation")

    monkeypatch.setattr(VariantCollection, "effects", unexpected)
    with pytest.raises(SystemExit) as error:
        main(["--genome", "GRCh38:2013_09", "--variant", "7", "100", "A", "T"])
    assert error.value.code == 1
    assert "provides GRCh37, not GRCh38" in capsys.readouterr().err
