"""Reacquire and check the primary sources of the curated ROS_Cfam fixture.

Usage: python check_ros_sources.py /path/to/sources [--download]
Without --download, all checks are offline. Reference annotations are installed
separately with PyEnsembl and their checksums are checked by test_canine_ros.py.
"""

import argparse
import hashlib
import json
from pathlib import Path
import re
from urllib.request import urlopen
import xml.etree.ElementTree as ET


def check_sources(directory, download=False):
    fixture = json.loads(Path(__file__).with_name("canine-ros-native.json").read_text())
    directory.mkdir(parents=True, exist_ok=True)
    articles = {}
    for key, source in fixture["sources"].items():
        path = directory / source["source_file"]
        if download:
            with urlopen(source["source_url"], timeout=60) as response:
                path.write_bytes(response.read())
        raw = path.read_bytes()
        assert hashlib.sha256(raw).hexdigest() == source["source_sha256"], path
        articles[key] = ET.fromstring(raw)

    for record in fixture["records"]:
        article = articles[record["source"]]
        table = article.find('.//table-wrap[@id="%s"]' % record["table"])
        rows = [[" ".join("".join(cell.itertext()).split()) for cell in row]
                for row in table.findall(".//tr")[1:]]
        if record["source"] == "melanoma":
            cells = rows[record["table_data_row"] - 1][4:]
            assert cells == ["SPTBN1", "56527053–56527054", "-", "Frameshift"]
            assert "Asn795fs" in re.sub(r"\s+", "", "".join(article.itertext()))
        else:
            # Only the gene column spans rows in these primary tables.
            gene = None
            expanded = []
            for row in rows:
                if len(row) == 9:
                    gene, *row = row
                assert len(row) == 8
                expanded.append([gene, *row])
            cells = expanded[record["table_data_row"] - 1]
            assert cells[:3] == [record["gene"], record["contig"], str(record["start"])]
            assert cells[4:] == [record["ref"], record["alt"], record["published_consequence"],
                                 record["published_vaf"], record["samples"]]
            text = re.sub(r"\s+", "", "".join(article.itertext()))
            assert record["published_cdna_change"] in text

        window = record["reference_window"]
        path = directory / ("varcode-ros-window-%s.json" % record["id"])
        if download:
            with urlopen(window["source_url"], timeout=60) as response:
                path.write_bytes(response.read())
        response = json.loads(path.read_bytes())
        # JSON key order is not stable across REST requests; the sequence and
        # versioned assembly/interval identity, rather than formatting, are pinned.
        assert response["id"] == window["id"], path
        assert response["seq"].upper() == window["sequence"], path
        assert hashlib.sha256(response["seq"].upper().encode()).hexdigest() == window["sequence_sha256"]
    print("Verified three publications and ten genomic reference windows")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_directory", type=Path)
    parser.add_argument("--download", action="store_true")
    args = parser.parse_args()
    check_sources(args.source_directory, download=args.download)
