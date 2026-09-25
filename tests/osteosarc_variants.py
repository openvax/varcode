"""Portable, provenance-preserving public variants for regression tests."""

import json
from pathlib import Path

from pyensembl import cached_release

from varcode import Variant, VariantCollection


FIXTURE_PATH = Path(__file__).parent / "data" / "osteosarc_variants.json"
SNAPSHOT_ID = "9b34ea0e13f9c1c35c3c88b7e646c0e608b86a143dee0e668bf3f74b909f815c"
OSTEOSARC_VERSION = "0.7.0"


def read_fixture(path=FIXTURE_PATH):
    with open(path) as handle:
        return json.load(handle)


def variant_from_entry(entry, genome):
    """Use the published assembly and allele, with explicit contig conversion."""
    if entry["assembly"] != genome.reference_name:
        raise ValueError("Fixture assembly does not match the supplied reference")
    if entry["status"] != "ready" or len(entry["alleles"]) != 1:
        raise ValueError("Fixture entry has no unique literal allele: %s" % entry["id"])
    contig, position, ref, alt = entry["alleles"][0]
    return Variant(contig, position, ref, alt, genome=genome,
                   convert_ucsc_contig_names=True)


def load_variants(path=FIXTURE_PATH, genome=None):
    """Load all ready entries, retaining every source ID after normalization.

    Unresolved entries stay in the JSON fixture for inspection. No reference
    files or source data are downloaded by this loader.
    """
    fixture = read_fixture(path)
    genome = cached_release(81) if genome is None else genome
    metadata = {}
    for entry in fixture["entries"]:
        if entry["status"] != "ready":
            continue
        variant = variant_from_entry(entry, genome)
        metadata.setdefault(variant, {"entries": []})["entries"].append(entry)
    source = "osteosarc:" + fixture["source"]["snapshot_id"]
    return VariantCollection(
        metadata, source_to_metadata_dict={source: metadata})
