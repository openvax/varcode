"""Offline checks against the bundled, pinned public osteosarc snapshot.

See tests/README.md for acquisition and invocation. Source annotations are
provenance, not an independent oracle for protein consequences.
"""

import hashlib
import os
import shutil
import zipfile
from dataclasses import replace
from importlib.metadata import version
from pathlib import Path

import pytest
from pyensembl import cached_release

from varcode import Variant
from varcode.mutant_transcript import apply_variant_to_transcript

from .osteosarc_variants import SNAPSHOT_ID, read_fixture


SNAPSHOT_NAME = "2026-09-18t"
SNAPSHOT_ARCHIVE = Path(__file__).parent / "data" / "osteosarc_snapshot_2026-09-18t.zip"
SNAPSHOT_ARCHIVE_SHA256 = "cbb688ca5cbe775fc4e6a826124f861d6d605c65e34473582c5d86853a6018fa"


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    snapshot = os.environ.get("OSTEOSARC_TEST_SNAPSHOT")
    cache_root = os.environ.get("OSTEOSARC_TEST_CACHE")
    if not snapshot:
        pytest.importorskip("osteosarc", reason="install .[test-data]")
        # Hash before unpacking the trusted fixture into an isolated cache.
        # Dataset.open additionally verifies every source receipt and object.
        assert hashlib.sha256(SNAPSHOT_ARCHIVE.read_bytes()).hexdigest() == SNAPSHOT_ARCHIVE_SHA256
        cache_root = tmp_path_factory.mktemp("osteosarc-snapshot")
        with zipfile.ZipFile(SNAPSHOT_ARCHIVE) as archive:
            archive.extractall(cache_root)
        snapshot = SNAPSHOT_NAME
    # Once an external snapshot is requested, missing dependencies/cache
    # entries fail rather than skipping or falling back to the fixture.
    from osteosarc import Cache, Dataset

    assert version("osteosarc") == "0.7.0"
    cache = Cache(cache_root, offline=True)
    data = Dataset.open(snapshot, cache=cache, offline=True)
    assert data.id == SNAPSHOT_ID
    return data


def test_osteosarc_native_variants_preserve_alleles_and_provenance(dataset):
    entries = dataset.variants(status="ready")
    assert len(dataset.variants()) == 182
    assert len(entries) == 179
    genome = cached_release(81)
    native = entries.to_varcode(genome=genome)
    expected_ids = {entry.id for entry in entries}
    actual_ids = set()
    for variant in native:
        assert variant.genome is genome
        metadata = native.metadata[variant]
        assert metadata["source"]["snapshot_id"] == SNAPSHOT_ID
        assert metadata["source"]["receipts"]
        for entry in metadata["entries"]:
            actual_ids.add(entry["id"])
            contig, position, ref, alt = entries[entry["id"]].allele
            assert variant == Variant(
                contig, position, ref, alt, genome,
                convert_ucsc_contig_names=True)
    assert actual_ids == expected_ids


def test_osteosarc_corrected_map2_is_a_distinct_complex_allele(dataset):
    entry = dataset.variants()["MAP2-chr2-209694768"]
    assert entry.allele == (
        "chr2", 209694768, "CCTGGGCTACTGTGTGTTCAATAAGTACACAGT", "CAGGG")
    assert "allele-MAP2-chr2-209694768" in entry.annotations["corrections"]
    genome = cached_release(81)
    selected = dataset.variants().select(ids=entry.id).to_varcode(genome=genome)
    variant, = selected
    old_deletion = Variant("2", 209694768, "CCTGGGCTACTGTGTGTTCAATA", "C", genome)
    assert variant != old_deletion
    t = genome.transcript_by_id("ENST00000360351")
    model = apply_variant_to_transcript(variant, t)
    assert model is not None
    assert model.mutant_protein_sequence is not None
    assert model.mutant_protein_sequence != t.protein_sequence


def test_historical_fixture_alleles_survive_package_update(dataset):
    # The checked-in corpus remains the 0.1.0 export. New curation can resolve
    # additional entries or enrich provenance without rewriting that history.
    fixture = read_fixture()
    entries = dataset.variants()
    assert {entry.id for entry in entries} == {entry["id"] for entry in fixture["entries"]}
    for entry in fixture["entries"]:
        if entry["status"] == "ready":
            assert entries[entry["id"]].allele == tuple(entry["alleles"][0])


def test_mitochondrial_native_conversion_can_be_annotated(dataset):
    genome = cached_release(81)
    entries = dataset.variants(status="ready")
    entry, = [entry for entry in entries if entry.allele[0] == "chrM"]
    variants = entries.select(ids=entry.id).to_varcode(genome=genome)
    variant, = variants
    assert variant.contig == "MT"
    assert variant.original_contig == "chrM"
    effects = variant.effects(raise_on_error=True)
    assert any(effect.gene_name == "MT-ND5"
               and effect.short_description == "p.A220T" for effect in effects)


def test_shared_cache_offline_reuse_and_integrity(dataset, tmp_path):
    from osteosarc import Cache, Dataset, IntegrityError, OfflineError
    from osteosarc.cache import Receipt

    shared = Cache(dataset.cache.root, offline=True)
    reopened = Dataset.open(dataset.manifest["name"],
                            cache=shared, offline=True)
    assert reopened.id == dataset.id
    receipt = Receipt(**dataset.variants().source["receipts"]["source_variants"])
    source_path = shared.path(receipt)
    assert Cache(shared.root, offline=True).path(receipt) == source_path
    with pytest.raises(OfflineError):
        Cache(tmp_path / "empty", offline=True).path(receipt)

    # Modify a private copy, never the user's shared cache.
    isolated = Cache(tmp_path / "tampered", offline=True)
    isolated.objects.mkdir(parents=True)
    target = isolated.objects / source_path.name
    shutil.copyfile(source_path, target)
    assert isolated.path(receipt) == target
    with pytest.raises(IntegrityError):
        isolated.path(replace(receipt, size=receipt.size + 1))
    with target.open("r+b") as handle:
        byte = handle.read(1)
        handle.seek(0)
        handle.write(bytes([byte[0] ^ 1]))
    with pytest.raises(IntegrityError):
        isolated.path(receipt)
