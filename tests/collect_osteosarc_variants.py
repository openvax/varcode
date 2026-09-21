"""Export a compact test corpus from a verified offline osteosarc snapshot.

Run with ``python -m tests.collect_osteosarc_variants --help``. Acquisition
is separate; ordinary tests only read the checked-in JSON fixture.
"""

import argparse
import hashlib
import json
from importlib.metadata import version
from pathlib import Path

from .osteosarc_variants import FIXTURE_PATH, OSTEOSARC_VERSION, SNAPSHOT_ID


def collect_variants(dataset):
    """Keep all site entries and source receipts, including unresolved alleles."""
    entries = dataset.variants("site")
    records = []
    for record in entries.to_records():
        # Source records contain large nested count/peptide tables. Retain
        # their identity without copying them into an allele-only fixture.
        digest = hashlib.sha256(json.dumps(
            record, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        annotations = record.pop("annotations")
        record["corrections"] = annotations.get("corrections", ())
        record["count_corrections"] = annotations.get("count_corrections", ())
        record["osteosarc_entry_sha256"] = digest
        records.append(record)
    return {
        "schema_version": 1,
        "dataset": "osteosarc",
        "dataset_revision": "site-variants-v1",
        "osteosarc_version": OSTEOSARC_VERSION,
        "selection": {"set": "site", "corrections": True},
        "source": entries.source,
        "entries": sorted(records, key=lambda record: record["id"]),
    }


def render_fixture(fixture):
    return json.dumps(fixture, indent=2, sort_keys=True) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", required=True, help="Name of an existing offline snapshot")
    parser.add_argument("--cache", help="Shared cache root (otherwise osteosarc's default)")
    parser.add_argument("--expected-snapshot-id", default=SNAPSHOT_ID,
                        help="Required snapshot identity; override explicitly for a new corpus")
    parser.add_argument("--output", type=Path, default=FIXTURE_PATH)
    args = parser.parse_args()
    if version("osteosarc") != OSTEOSARC_VERSION:
        parser.error("Collection requires osteosarc==%s" % OSTEOSARC_VERSION)
    from osteosarc import Cache, Dataset

    dataset = Dataset.open(args.snapshot, cache=Cache(args.cache, offline=True),
                           offline=True, corrections=True)
    if dataset.id != args.expected_snapshot_id:
        parser.error("Snapshot identity does not match --expected-snapshot-id")
    fixture = collect_variants(dataset)
    args.output.write_text(render_fixture(fixture))
    ready = sum(record["status"] == "ready" for record in fixture["entries"])
    print("Collected %d ready and %d unresolved entries into %s" % (
        ready, len(fixture["entries"]) - ready, args.output))


if __name__ == "__main__":
    main()
