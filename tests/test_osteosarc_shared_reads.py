"""Varcode's observed-junction fixtures match their reads in openvax-v1.

openvax-v1 is the OpenVax libraries' shared Sid test data, published by
osteosarc (iskandr/osteosarc#56). Each record in
``tests/data/osteosarc_observed_junctions.json`` keeps only the junction
window and its annotations; the openvax-v1 member named after it holds the ONT
read the window came from. The first run downloads the bundle (28 MB) and
exports its members into the osteosarc cache; later runs reuse them offline.
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


def shared_read_fields(osteosarc, label):
    """SAM fields of each record in one junction's openvax-v1 member."""
    path = osteosarc.bundle_file("openvax-v1", MEMBER_PREFIX + label, format="sam")
    return [line.split("\t") for line in Path(path).read_text().splitlines()
            if not line.startswith("@")]


def test_observed_junction_windows_come_from_their_shared_reads():
    osteosarc = pytest.importorskip("osteosarc", reason="install .[test-data]")
    for row in junction_records():
        (fields,) = shared_read_fields(osteosarc, row["label"])
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
