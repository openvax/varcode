"""Varcode's observed-junction fixtures match their reads in openvax-v1.

openvax-v1 is the OpenVax libraries' shared Sid test data, published by
osteosarc (iskandr/osteosarc#56). Each record in
``tests/data/osteosarc_observed_junctions.json`` keeps only the junction
window and its annotations; the openvax-v1 member named after it holds the ONT
read the window came from. The first run downloads the bundle (28 MB) into the
osteosarc cache; later runs reuse it offline.
"""

import hashlib
import json
from pathlib import Path

import pytest


FIXTURE = Path(__file__).parent / "data" / "osteosarc_observed_junctions.json"
MEMBER_PREFIX = "varcode/osteosarc_observed_junctions.json#"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def junction_records():
    return json.loads(FIXTURE.read_text())["records"]


@pytest.fixture(scope="module")
def shared_reads(tmp_path_factory):
    """SAM fields of each fixture's openvax-v1 member, by junction label."""
    osteosarc = pytest.importorskip("osteosarc", reason="install .[test-data]")
    bundle = osteosarc.fetch_bundle("openvax-v1")
    members = [MEMBER_PREFIX + row["label"] for row in junction_records()]
    paths = osteosarc.export_bundle(
        bundle, tmp_path_factory.mktemp("openvax-v1"), members=members, format="sam")
    return {
        member[len(MEMBER_PREFIX):]: [
            line.split("\t") for line in Path(path).read_text().splitlines()
            if not line.startswith("@")]
        for member, path in paths.items()}


def test_observed_junction_windows_come_from_their_shared_reads(shared_reads):
    records = junction_records()
    assert sorted(shared_reads) == sorted(row["label"] for row in records)
    for row in records:
        (fields,) = shared_reads[row["label"]]
        assert fields[0] == row["read_name"]
        # The audit hashed each read in one of its two orientations, and took
        # the junction window from that same orientation.
        sequence = fields[9]
        orientations = [sequence, sequence.translate(COMPLEMENT)[::-1]]
        (source,) = [
            candidate for candidate in orientations
            if hashlib.sha256(candidate.encode()).hexdigest()
            == row["source_sequence_sha256"]]
        assert row["sequence"] in source
