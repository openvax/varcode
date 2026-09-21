"""Run real public site alleles through annotation without an external cache."""

from collections import Counter

import pytest
from pyensembl import cached_release

from varcode.effects import Failure, Unresolved

from .osteosarc_variants import (
    OSTEOSARC_VERSION,
    SNAPSHOT_ID,
    load_variants,
    read_fixture,
    variant_from_entry,
)


FIXTURE = read_fixture()
READY = [entry for entry in FIXTURE["entries"] if entry["status"] == "ready"]
UNRESOLVED = [entry for entry in FIXTURE["entries"] if entry["status"] != "ready"]


def test_collected_variants_keep_provenance_and_unresolved_entries():
    assert FIXTURE["source"]["snapshot_id"] == SNAPSHOT_ID
    assert FIXTURE["osteosarc_version"] == OSTEOSARC_VERSION
    assert FIXTURE["selection"] == {"set": "site", "corrections": True}
    assert len(READY) == 177
    assert len(UNRESOLVED) == 5
    assert Counter(entry["status"] for entry in UNRESOLVED) == {
        "missing_literal_allele": 3, "non_literal_allele": 2}
    assert len({entry["id"] for entry in FIXTURE["entries"]}) == 182
    assert all(len(entry["osteosarc_entry_sha256"]) == 64 for entry in FIXTURE["entries"])
    assert all(receipt["url"] and len(receipt["sha256"]) == 64
               for receipt in FIXTURE["source"]["receipts"].values())
    variants = load_variants()
    assert variants.source == "osteosarc:" + SNAPSHOT_ID
    assert {entry["id"] for metadata in variants.metadata.values()
            for entry in metadata["entries"]} == {entry["id"] for entry in READY}


@pytest.mark.parametrize("entry", READY, ids=lambda entry: entry["id"])
def test_collected_variant_can_be_annotated(entry, dual_annotator):
    variant = variant_from_entry(entry, cached_release(81))
    effects = variant.effects(raise_on_error=True, annotator=dual_annotator)
    assert effects
    assert not any(isinstance(effect, (Failure, Unresolved)) for effect in effects)


@pytest.mark.parametrize("entry", UNRESOLVED, ids=lambda entry: entry["id"])
def test_unresolved_entries_are_not_fabricated_as_native_variants(entry):
    with pytest.raises(ValueError, match="unique literal allele"):
        variant_from_entry(entry, cached_release(81))


def test_collected_variants_reject_wrong_assembly():
    with pytest.raises(ValueError, match="assembly"):
        load_variants(genome=cached_release(75))


def test_collected_mitochondrial_variant_preserves_original_contig():
    entry, = [entry for entry in READY if entry["alleles"][0][0] == "chrM"]
    variant = variant_from_entry(entry, cached_release(81))
    assert variant.original_contig == "chrM"
    assert variant.contig == "MT"
    assert variant.transcripts


def test_collected_map2_retains_corrected_allele_and_distinct_split_call():
    entries = {entry["id"]: entry for entry in READY}
    corrected = entries["MAP2-chr2-209694768"]
    assert corrected["alleles"] == [[
        "chr2", 209694768, "CCTGGGCTACTGTGTGTTCAATAAGTACACAGT", "CAGGG"]]
    assert "allele-MAP2-chr2-209694768" in corrected["corrections"]
    split = entries["MAP2-chr2-209694772"]
    assert "map2-split-representations" in split["corrections"]
    genome = cached_release(81)
    assert variant_from_entry(corrected, genome) != variant_from_entry(split, genome)
