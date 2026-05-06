# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""
Hard-coded expectations against the canonical VCF spec example fixtures
under ``tests/data/spec_examples/``.

The oracle tests in test_vcf_parsing.py already pin our parser to PyVCF3's
behaviour for *every* fixture. These tests instead check that our parser
returns exactly the values published in the VCF spec — so a regression
shows up as a diff against the spec, not just against PyVCF3.

Sources (samtools/hts-specs, MIT):
- VCFv4.2 §1.1 — three-sample trio example
- VCFv4.3 §5.4 — structural variant example
"""
from __future__ import annotations

import os

import pytest

from varcode.vcf_parsing import VCFHeader, FieldDef


SPEC_DIR = os.path.join(os.path.dirname(__file__), "data", "spec_examples")
TRIO_PATH = os.path.join(SPEC_DIR, "vcf42_spec_trio.vcf")
SV_PATH = os.path.join(SPEC_DIR, "vcf43_spec_sv.vcf")


# ---------------------------------------------------------------------------
# VCF 4.2 §1.1 — trio example
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def trio_header():
    return VCFHeader.from_path(TRIO_PATH)


class TestVCF42TrioHeader:
    """The header of the canonical 4.2 example exercises every metadata
    shape the spec defines: singular keys, multi-valued keys, INFO and
    FORMAT declarations across all five Type values, and a structured
    contig line with a quoted species field."""

    def test_fileformat(self, trio_header):
        assert trio_header.metadata["fileformat"] == "VCFv4.2"

    def test_filedate_is_singular_string(self, trio_header):
        # In SINGULAR_METADATA — stored as plain string, not list.
        assert trio_header.metadata["fileDate"] == "20090805"

    def test_reference_is_singular_string(self, trio_header):
        assert trio_header.metadata["reference"] == (
            "file:///seq/references/1000GenomesPilot-NCBI36.fasta")

    def test_non_singular_metadata_is_list(self, trio_header):
        # Non-SINGULAR keys (source, phasing, contig) collect into list[str].
        assert trio_header.metadata["source"] == ["myImputationProgramV3.1"]
        assert trio_header.metadata["phasing"] == ["partial"]
        # The contig line keeps its <...> body as the single string element.
        assert len(trio_header.metadata["contig"]) == 1
        assert 'species="Homo sapiens"' in trio_header.metadata["contig"][0]

    def test_samples_are_three_named_individuals(self, trio_header):
        assert trio_header.samples == ["NA00001", "NA00002", "NA00003"]

    def test_info_field_count_and_order(self, trio_header):
        # Order in the file: NS, DP, AF, AA, DB, H2.
        assert list(trio_header.info_fields) == ["NS", "DP", "AF", "AA", "DB", "H2"]

    def test_info_ns_declaration(self, trio_header):
        assert trio_header.info_fields["NS"] == FieldDef(
            id="NS", number=1, type="Integer",
            description="Number of Samples With Data")

    def test_info_af_is_per_alt(self, trio_header):
        af = trio_header.info_fields["AF"]
        assert af.number == -1  # "A" — one value per ALT
        assert af.type == "Float"

    def test_info_db_is_flag(self, trio_header):
        db = trio_header.info_fields["DB"]
        assert db.number == 0
        assert db.type == "Flag"

    def test_info_db_description_has_comma(self, trio_header):
        # "dbSNP membership, build 129" — the embedded comma is the canonical
        # case our quote-aware splitter has to handle.
        assert (trio_header.info_fields["DB"].description
                == "dbSNP membership, build 129")

    def test_format_field_count_and_order(self, trio_header):
        assert list(trio_header.format_fields) == ["GT", "GQ", "DP", "HQ"]

    def test_format_hq_is_pair(self, trio_header):
        hq = trio_header.format_fields["HQ"]
        assert hq.number == 2
        assert hq.type == "Integer"
        assert hq.description == "Haplotype Quality"


class TestVCF42TrioInfoParsing:
    """Each row in the spec example was hand-decoded from the published spec."""

    def test_row1_snp_with_flags(self, trio_header):
        # rs6054257: G>A, NS=3;DP=14;AF=0.5;DB;H2
        info = trio_header.parse_info("NS=3;DP=14;AF=0.5;DB;H2")
        assert info == {
            "NS": 3,         # Number=1 -> unwrapped scalar
            "DP": 14,
            "AF": [0.5],     # Number=A -> always list, even for one ALT
            "DB": True,      # Flag, present
            "H2": True,
        }

    def test_row2_low_quality(self, trio_header):
        info = trio_header.parse_info("NS=3;DP=11;AF=0.017")
        assert info == {"NS": 3, "DP": 11, "AF": [0.017]}

    def test_row3_multiallelic_per_alt_af(self, trio_header):
        # rs6040355: A>G,T — AF must be parallel to ALTs.
        info = trio_header.parse_info(
            "NS=2;DP=10;AF=0.333,0.667;AA=T;DB")
        assert info == {
            "NS": 2,
            "DP": 10,
            "AF": [0.333, 0.667],
            "AA": "T",
            "DB": True,
        }

    def test_row4_no_alt(self, trio_header):
        # ALT=. site — AF absent from INFO.
        info = trio_header.parse_info("NS=3;DP=13;AA=T")
        assert info == {"NS": 3, "DP": 13, "AA": "T"}

    def test_row5_microsat(self, trio_header):
        info = trio_header.parse_info("NS=3;DP=9;AA=G")
        assert info == {"NS": 3, "DP": 9, "AA": "G"}

    def test_dot_returns_empty(self, trio_header):
        assert trio_header.parse_info(".") == {}


class TestVCF42TrioSampleParsing:
    """Per-sample FORMAT decoding for each row's three samples."""

    def test_row1_phased_genotypes(self, trio_header):
        out = trio_header.parse_samples(
            ["0|0:48:1:51,51", "1|0:48:8:51,51", "1/1:43:5:.,."],
            "GT:GQ:DP:HQ")
        assert list(out) == ["NA00001", "NA00002", "NA00003"]
        assert out["NA00001"]["GT"] == "0|0"           # phased preserved
        assert out["NA00001"]["GQ"] == 48
        assert out["NA00001"]["DP"] == 1
        assert out["NA00001"]["HQ"] == [51, 51]
        # NA00003 has missing HQ values: ".,." -> [None, None]
        assert out["NA00003"]["HQ"] == [None, None]

    def test_row2_trailing_missing_hq(self, trio_header):
        out = trio_header.parse_samples(
            ["0|0:49:3:58,50", "0|1:3:5:65,3", "0/0:41:3"],
            "GT:GQ:DP:HQ")
        # NA00003 has only GT:GQ:DP populated; HQ trails missing.
        assert out["NA00003"]["GT"] == "0/0"
        assert out["NA00003"]["DP"] == 3
        assert out["NA00003"]["HQ"] is None

    def test_row3_mixed_phasing(self, trio_header):
        out = trio_header.parse_samples(
            ["1|2:21:6:23,27", "2|1:2:0:18,2", "2/2:35:4"],
            "GT:GQ:DP:HQ")
        assert out["NA00001"]["GT"] == "1|2"
        assert out["NA00002"]["GT"] == "2|1"
        assert out["NA00003"]["GT"] == "2/2"  # last sample unphased
        # Trailing HQ missing on NA00003.
        assert out["NA00003"]["HQ"] is None

    def test_row5_microsat_format_drops_hq(self, trio_header):
        # Only GT:GQ:DP for the microsat row — HQ not in FORMAT, so absent.
        out = trio_header.parse_samples(
            ["0/1:35:4", "0/2:17:2", "1/1:40:3"],
            "GT:GQ:DP")
        assert "HQ" not in out["NA00001"]
        assert out["NA00001"]["GT"] == "0/1"
        assert out["NA00002"]["GT"] == "0/2"
        assert out["NA00003"]["GT"] == "1/1"


# ---------------------------------------------------------------------------
# VCF 4.3 §5.4 — structural variant example
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def sv_header():
    return VCFHeader.from_path(SV_PATH)


class TestVCF43StructuralVariantHeader:
    def test_fileformat_is_4_3(self, sv_header):
        assert sv_header.metadata["fileformat"] == "VCFv4.3"

    def test_sv_specific_info_fields_declared(self, sv_header):
        assert "SVTYPE" in sv_header.info_fields
        assert "END" in sv_header.info_fields
        assert "CIPOS" in sv_header.info_fields
        assert "IMPRECISE" in sv_header.info_fields

    def test_cipos_is_pair_of_integers(self, sv_header):
        # The spec defines CIPOS as Number=2,Type=Integer.
        assert sv_header.info_fields["CIPOS"].number == 2
        assert sv_header.info_fields["CIPOS"].type == "Integer"

    def test_svlen_is_variable_count(self, sv_header):
        # SVLEN is Number=. (variable) -> encoded as None.
        assert sv_header.info_fields["SVLEN"].number is None

    def test_imprecise_is_flag(self, sv_header):
        f = sv_header.info_fields["IMPRECISE"]
        assert f.type == "Flag"
        assert f.number == 0


class TestVCF43StructuralVariantInfoParsing:
    def test_imprecise_deletion_with_confidence_intervals(self, sv_header):
        info = sv_header.parse_info(
            "IMPRECISE;SVTYPE=DEL;END=321887;SVLEN=-205;"
            "CIPOS=-56,20;CIEND=-10,62")
        assert info == {
            "IMPRECISE": True,
            "SVTYPE": "DEL",        # Number=1 -> scalar string
            "END": 321887,          # Number=1 -> scalar int
            "SVLEN": [-205],        # Number=. -> list
            "CIPOS": [-56, 20],     # Number=2 -> list
            "CIEND": [-10, 62],
        }

    def test_inversion_zero_length(self, sv_header):
        info = sv_header.parse_info(
            "SVTYPE=INV;END=9425916;SVLEN=0;CIPOS=-50,50;CIEND=-50,50")
        assert info["SVTYPE"] == "INV"
        assert info["SVLEN"] == [0]
        assert info["CIPOS"] == [-50, 50]


class TestVCF43StructuralVariantSampleParsing:
    def test_single_sample_genotype_quality(self, sv_header):
        assert sv_header.samples == ["SAMPLE"]
        out = sv_header.parse_samples(["0/1:12"], "GT:GQ")
        assert out["SAMPLE"]["GT"] == "0/1"
        assert out["SAMPLE"]["GQ"] == 12


# ---------------------------------------------------------------------------
# Documentation-style: full trip through every record in the trio fixture.
# ---------------------------------------------------------------------------

class TestTrioEndToEndExpectedShapes:
    """Walk the whole trio fixture and verify the structure of every parsed
    row. Useful as living documentation for what the parser returns end-to-end.
    """

    @pytest.fixture(scope="class")
    def parsed_rows(self):
        h = VCFHeader.from_path(TRIO_PATH)
        rows = []
        with open(TRIO_PATH) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                rows.append({
                    "chrom": cols[0],
                    "pos": int(cols[1]),
                    "id": cols[2],
                    "ref": cols[3],
                    "alt": cols[4],
                    "info": h.parse_info(cols[7]),
                    "format": cols[8],
                    "samples": h.parse_samples(cols[9:], cols[8]),
                })
        return rows

    def test_five_records(self, parsed_rows):
        assert len(parsed_rows) == 5

    def test_first_record_is_known_dbsnp_snp(self, parsed_rows):
        r = parsed_rows[0]
        assert r["id"] == "rs6054257"
        assert r["ref"] == "G"
        assert r["alt"] == "A"
        assert r["info"]["DB"] is True

    def test_third_record_is_multiallelic(self, parsed_rows):
        r = parsed_rows[2]
        assert r["alt"].split(",") == ["G", "T"]
        assert r["info"]["AF"] == [0.333, 0.667]

    def test_fourth_record_has_no_alt(self, parsed_rows):
        # ALT="." is preserved as-is at this layer; load_vcf filters such
        # rows out at a higher level.
        assert parsed_rows[3]["alt"] == "."

    def test_microsat_row_missing_hq(self, parsed_rows):
        r = parsed_rows[4]
        # FORMAT only has GT:GQ:DP — HQ not present in this row's samples.
        for sname in ("NA00001", "NA00002", "NA00003"):
            assert "HQ" not in r["samples"][sname]
            assert "GT" in r["samples"][sname]


class TestUsageExampleScript:
    """The standalone example in examples/vcf_header_usage.py exists to
    document the public surface; pin it with a smoke test so it doesn't
    silently rot."""

    def test_usage_example_runs_cleanly(self):
        import subprocess
        import sys
        repo_root = os.path.dirname(os.path.dirname(__file__))
        script = os.path.join(repo_root, "examples", "vcf_header_usage.py")
        result = subprocess.run(
            [sys.executable, script], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        # Spot-check a couple of expected lines so a silently-broken parser
        # would fail this test even if the script doesn't error out.
        assert "fileformat       : VCFv4.2" in result.stdout
        assert "samples          : ['NA00001', 'NA00002', 'NA00003']" in result.stdout
        assert "rs6054257" in result.stdout
