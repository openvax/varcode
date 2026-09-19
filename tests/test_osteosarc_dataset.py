"""Opt-in, offline checks against the pinned public osteosarc snapshot.

See tests/README.md for acquisition and invocation. Source annotations are
provenance, not an independent oracle for protein consequences.
"""

import os
from importlib.metadata import version

import pytest
from pyensembl import cached_release

from varcode import Variant
from varcode.mutant_transcript import apply_variant_to_transcript


SNAPSHOT_ID = "9b34ea0e13f9c1c35c3c88b7e646c0e608b86a143dee0e668bf3f74b909f815c"


@pytest.fixture(scope="module")
def dataset():
    snapshot = os.environ.get("OSTEOSARC_TEST_SNAPSHOT")
    if not snapshot:
        pytest.skip("set OSTEOSARC_TEST_SNAPSHOT to enable the offline corpus checks")
    # Once explicitly requested, missing dependencies/cache entries must fail.
    from osteosarc import Cache, Dataset

    assert version("osteosarc") == "0.1.0"
    cache = Cache(os.environ.get("OSTEOSARC_TEST_CACHE"), offline=True)
    data = Dataset.open(snapshot, cache=cache, offline=True)
    assert data.id == SNAPSHOT_ID
    return data


def test_osteosarc_native_variants_preserve_alleles_and_provenance(dataset):
    entries = dataset.variants(status="ready")
    assert len(dataset.variants()) == 182
    assert len(entries) == 177
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
            assert variant == Variant(contig.removeprefix("chr"), position, ref, alt, genome)
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
