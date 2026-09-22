"""Opt-in, offline checks against the pinned public osteosarc snapshot.

See tests/README.md for acquisition and invocation. Source annotations are
provenance, not an independent oracle for protein consequences.
"""

import os
import shutil
from dataclasses import replace
from importlib.metadata import version

import pytest
from pyensembl import cached_release

from varcode import Variant
from varcode.mutant_transcript import apply_variant_to_transcript

from .osteosarc_variants import SNAPSHOT_ID, read_fixture



@pytest.fixture(scope="module")
def dataset():
    snapshot = os.environ.get("OSTEOSARC_TEST_SNAPSHOT")
    if not snapshot:
        pytest.skip("set OSTEOSARC_TEST_SNAPSHOT to enable the offline corpus checks")
    # Once explicitly requested, missing dependencies/cache entries must fail.
    from osteosarc import Cache, Dataset

    assert version("osteosarc") == "0.1.4"
    cache = Cache(os.environ.get("OSTEOSARC_TEST_CACHE"), offline=True)
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

    shared = Cache(os.environ.get("OSTEOSARC_TEST_CACHE"), offline=True)
    reopened = Dataset.open(os.environ["OSTEOSARC_TEST_SNAPSHOT"],
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
