# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""
Hand-verified expectations against fixtures from real somatic callers.

These complement the parametrized PyVCF3 oracle tests by pinning specific
field values that *should* survive parsing — so a regression on, say, MuTect2
phasing or Strelka2 per-base tier counts shows up as a meaningful failure
against the published format, not just "PyVCF3 says X."

Fixtures live under ``tests/data/real_callers/``; see the README there for
provenance and the format quirks each one exercises.
"""
from __future__ import annotations

import os

import pytest

from varcode import load_vcf
from varcode.vcf_parsing import VCFHeader


CALLER_DIR = os.path.join(os.path.dirname(__file__), "data", "real_callers")
MUTECT2 = os.path.join(CALLER_DIR, "mutect2_example.vcf")
STRELKA2 = os.path.join(CALLER_DIR, "strelka2_somatic_snvs.vcf")
VEP = os.path.join(CALLER_DIR, "vep_annotated_csq.vcf")


# ---------------------------------------------------------------------------
# GATK4 MuTect2
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def mutect2_header():
    return VCFHeader.from_path(MUTECT2)


class TestMuTect2Header:
    def test_filter_vocabulary(self, mutect2_header):
        # MuTect2's distinctive filter set — not present in MuTect1 fixtures.
        # We don't store FILTER declarations in metadata (varcode reads
        # filters per-row), but the header should parse without error.
        assert mutect2_header.metadata["fileformat"] == "VCFv4.2"
        assert mutect2_header.metadata["source"] == ["Mutect2"]

    def test_tumor_normal_metadata(self, mutect2_header):
        # MuTect2 records which sample is tumor vs. normal in dedicated
        # metadata lines (non-singular -> list).
        assert mutect2_header.metadata["tumor_sample"] == ["TUMOR"]
        assert mutect2_header.metadata["normal_sample"] == ["NORMAL"]

    def test_per_alt_info_fields_have_number_a(self, mutect2_header):
        # TLOD, NLOD, NALOD, POPAF, MPOS — one value per ALT.
        for key in ("TLOD", "NLOD", "NALOD", "POPAF", "MPOS"):
            assert mutect2_header.info_fields[key].number == -1, (
                f"{key} should be Number=A (encoded as -1)")

    def test_per_allele_info_fields_have_number_r(self, mutect2_header):
        # MBQ, MFRL, MMQ — one value per (REF + each ALT).
        for key in ("MBQ", "MFRL", "MMQ"):
            assert mutect2_header.info_fields[key].number == -3, (
                f"{key} should be Number=R (encoded as -3)")

    def test_str_flag(self, mutect2_header):
        # STR is a Flag — Number=0, Type=Flag.
        assert mutect2_header.info_fields["STR"].type == "Flag"
        assert mutect2_header.info_fields["STR"].number == 0

    def test_rpa_is_variable_count(self, mutect2_header):
        # RPA repeats-per-allele has Number=. (variable).
        assert mutect2_header.info_fields["RPA"].number is None

    def test_phasing_format_fields(self, mutect2_header):
        # PGT and PID are MuTect2's physical phasing tags.
        assert "PGT" in mutect2_header.format_fields
        assert "PID" in mutect2_header.format_fields
        assert mutect2_header.format_fields["PGT"].type == "String"


class TestMuTect2InfoParsing:
    def test_simple_pass_record_per_alt_arrays(self, mutect2_header):
        # chr1:1041196 — single ALT, all per-alt fields are 1-element lists.
        info = mutect2_header.parse_info(
            "DP=87;ECNT=1;GERMQ=89;MBQ=20,30;MFRL=183,200;MMQ=60,60;"
            "MPOS=24;NALOD=1.4;NLOD=12.31;POPAF=6.00;ROQ=93;TLOD=53.78")
        assert info["TLOD"] == [53.78]
        assert info["NLOD"] == [12.31]
        assert info["MBQ"] == [20, 30]      # Number=R: ref + 1 alt
        assert info["MMQ"] == [60, 60]
        assert info["DP"] == 87             # Number=1 -> scalar
        assert info["ECNT"] == 1
        assert info["ROQ"] == 93.0          # Float scalar

    def test_multiallelic_per_alt_arrays(self, mutect2_header):
        # chr1:1102328 — two ALTs, A and T. Number=A fields have 2 values;
        # Number=R fields have 3 (ref + 2 alts).
        info = mutect2_header.parse_info(
            "DP=64;ECNT=2;GERMQ=72;MBQ=29,28,29;MFRL=193,178,201;"
            "MMQ=60,60,60;MPOS=12,30;NALOD=1.32,1.32;NLOD=6.32,6.32;"
            "POPAF=6.00,6.00;ROQ=93;TLOD=4.31,3.85")
        assert info["TLOD"] == [4.31, 3.85]
        assert info["MPOS"] == [12, 30]
        assert info["MBQ"] == [29, 28, 29]      # Number=R for multi-allelic
        assert info["MFRL"] == [193, 178, 201]

    def test_str_contraction_with_repeat_info(self, mutect2_header):
        # chr2:16092200 — STR flag set, RPA Number=. captures both ref+alt counts.
        info = mutect2_header.parse_info(
            "DP=58;ECNT=1;GERMQ=70;MBQ=29,29;MFRL=180,182;MMQ=60,60;"
            "MPOS=29;NALOD=1.31;NLOD=6.32;POPAF=4.30;ROQ=93;"
            "RPA=4,3;RU=CAA;STR;TLOD=18.27")
        assert info["STR"] is True              # Flag presence
        assert info["RU"] == "CAA"              # Number=1 -> scalar string
        assert info["RPA"] == [4, 3]            # Number=. -> list


class TestMuTect2SampleParsing:
    def test_normal_tumor_pair_simple(self, mutect2_header):
        # chr1:1041196 — normal hom-ref, tumor het.
        out = mutect2_header.parse_samples(
            ["0/0:42,0:0.023:42:18,0:24,0:99:0,120,1800",
             "0/1:30,15:0.336:45:14,7:16,8:99:451,0,925"],
            "GT:AD:AF:DP:F1R2:F2R1:GQ:PL")
        assert out["NORMAL"]["GT"] == "0/0"
        assert out["NORMAL"]["AD"] == [42, 0]
        assert out["NORMAL"]["AF"] == [0.023]   # Number=A even single-element
        assert out["TUMOR"]["GT"] == "0/1"
        assert out["TUMOR"]["AD"] == [30, 15]
        assert out["TUMOR"]["AF"] == [0.336]
        assert out["TUMOR"]["DP"] == 45

    def test_multiallelic_with_phasing_and_triallelic_gt(self, mutect2_header):
        # chr1:1102328 — TUMOR is genotype 0/1/2 (triallelic) with phased
        # PGT and a PID linking to the previous-record phasing group.
        out = mutect2_header.parse_samples(
            ["0/0:21,0,0:0.045,0.045:21:8,0,0:13,0,0:60:.:.:0,63,945,63,945,945",
             "0/1/2:14,5,4:0.245,0.214:23:8,2,2:6,3,2:54:0|1:1102328_G_A:54,0,389,141,182,498"],
            "GT:AD:AF:DP:F1R2:F2R1:GQ:PGT:PID:PL")

        # Triallelic GT preserved as raw string (not parsed structurally).
        assert out["TUMOR"]["GT"] == "0/1/2"
        # AD has 3 elements (ref + 2 alts).
        assert out["TUMOR"]["AD"] == [14, 5, 4]
        # AF has 2 elements (per-alt).
        assert out["TUMOR"]["AF"] == [0.245, 0.214]
        # PL has 6 elements: G(n_alleles=3) -> 6 genotype combos.
        assert out["TUMOR"]["PL"] == [54, 0, 389, 141, 182, 498]
        # PGT scalar string with phased pipe preserved.
        assert out["TUMOR"]["PGT"] == "0|1"
        assert out["TUMOR"]["PID"] == "1102328_G_A"

        # NORMAL has missing PGT/PID (".").
        assert out["NORMAL"]["PGT"] is None
        assert out["NORMAL"]["PID"] is None
        # NORMAL still has 3-element AD because ALT is multi-allelic.
        assert out["NORMAL"]["AD"] == [21, 0, 0]


class TestMuTect2EndToEnd:
    def test_load_vcf_filters_to_pass(self):
        vc = load_vcf(MUTECT2, genome="GRCh38")
        # 5 records in file; 2 are PASS (chr1:1041196, chr2:16092115).
        assert len(vc) == 2

    def test_load_vcf_keeps_filter_failures_with_only_passing_false(self):
        vc = load_vcf(MUTECT2, genome="GRCh38", only_passing=False)
        # 5 records, but the multi-allelic G→A,T splits into 2 variants,
        # so total Variant count is 6.
        assert len(vc) == 6


# ---------------------------------------------------------------------------
# Strelka2 somatic SNVs
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def strelka2_header():
    return VCFHeader.from_path(STRELKA2)


class TestStrelka2Header:
    def test_strelka2_metadata(self, strelka2_header):
        assert strelka2_header.metadata["source"] == ["strelka"]
        assert strelka2_header.metadata["source_version"] == ["2.9.10"]

    def test_cmdline_metadata_with_paths_and_flags(self, strelka2_header):
        # ##cmdline= contains spaces, slashes, and CLI flags. Stored verbatim
        # as a single string element in the list-valued metadata key.
        cmdline = strelka2_header.metadata["cmdline"]
        assert isinstance(cmdline, list) and len(cmdline) == 1
        assert "configureStrelkaSomaticWorkflow.py" in cmdline[0]
        assert "--normalBam normal.bam" in cmdline[0]

    def test_per_base_tier_format_fields_are_pairs(self, strelka2_header):
        # AU/CU/GU/TU are Number=2,Type=Integer (tier1, tier2 counts).
        for base in ("AU", "CU", "GU", "TU"):
            f = strelka2_header.format_fields[base]
            assert f.number == 2
            assert f.type == "Integer"

    def test_no_gt_in_format_fields(self, strelka2_header):
        # Strelka2 SNV records famously have no GT — header reflects that.
        assert "GT" not in strelka2_header.format_fields


class TestStrelka2InfoParsing:
    def test_sgt_arrow_notation(self, strelka2_header):
        # SGT="GG->AG" — the `>` is just a string char, no special handling.
        info = strelka2_header.parse_info(
            "SOMATIC;QSS=42;TQSS=1;NT=ref;QSS_NT=42;TQSS_NT=1;"
            "SGT=GG->AG;DP=58;MQ=60.00;MQ0=0;ReadPosRankSum=0.00;SNVSB=0.00")
        assert info["SGT"] == "GG->AG"
        assert info["SOMATIC"] is True
        assert info["QSS"] == 42
        assert info["MQ"] == 60.0
        assert info["NT"] == "ref"


class TestStrelka2SampleParsing:
    def test_per_base_tier_counts_for_g_to_a_snv(self, strelka2_header):
        # chr1:1234567 — G>A SNV.
        # NORMAL: 30 reads, all G (GU=30,30 in tier1,tier2).
        # TUMOR:  28 reads, 11 are A (AU=11,11), 17 are G (GU=17,17).
        out = strelka2_header.parse_samples(
            ["30:0:0:0:0,0:0,0:30,30:0,0",
             "28:0:0:0:11,11:0,0:17,17:0,0"],
            "DP:FDP:SDP:SUBDP:AU:CU:GU:TU")
        assert out["NORMAL"]["DP"] == 30
        assert out["NORMAL"]["GU"] == [30, 30]   # all reads are G
        assert out["NORMAL"]["AU"] == [0, 0]
        assert out["TUMOR"]["AU"] == [11, 11]    # 11 alt-allele reads in both tiers
        assert out["TUMOR"]["GU"] == [17, 17]    # 17 ref-allele reads
        # No GT field anywhere.
        assert "GT" not in out["NORMAL"]
        assert "GT" not in out["TUMOR"]


class TestStrelka2EndToEnd:
    def test_load_vcf_pass_only(self):
        vc = load_vcf(STRELKA2, genome="hg19")
        # 3 records, 2 are PASS.
        assert len(vc) == 2

    def test_no_gt_does_not_break_load_vcf(self):
        # The whole point of this fixture: load_vcf must not assume GT exists.
        vc = load_vcf(STRELKA2, genome="hg19", only_passing=False)
        # Pull sample_info from metadata and confirm GT is absent.
        meta = vc.source_to_metadata_dict[STRELKA2]
        for variant in vc:
            for sample_name, fields in meta[variant]["sample_info"].items():
                assert "GT" not in fields, (
                    f"Strelka2 SNV records should have no GT; got {sample_name}")
                # Per-base tier counts should be present.
                assert "AU" in fields and "CU" in fields
                assert "GU" in fields and "TU" in fields


# ---------------------------------------------------------------------------
# VEP-annotated CSQ
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def vep_header():
    return VCFHeader.from_path(VEP)


class TestVEPHeader:
    def test_csq_description_with_pipes_survives_quoted_split(self, vep_header):
        # The CSQ Description embeds the pipe-delimited subfield list. Our
        # quote-aware splitter must capture the whole thing intact — a naive
        # split-on-comma would merge it with neighbour fields.
        csq = vep_header.info_fields["CSQ"]
        assert csq.number is None      # Number=. (per-transcript list)
        assert csq.type == "String"
        assert "Format: Allele|Consequence" in csq.description

    def test_csq_subfield_count(self, vep_header):
        # 23 pipe-separated subfield names in this fixture's Format spec.
        csq = vep_header.info_fields["CSQ"]
        format_spec = csq.description.split("Format: ", 1)[1]
        subfields = format_spec.split("|")
        assert len(subfields) == 23
        assert subfields[0] == "Allele"
        assert subfields[-1] == "HGNC_ID"

    def test_vep_command_line_metadata_with_hyphenated_key(self, vep_header):
        # Metadata key "VEP-command-line" contains a hyphen; our _META_LINE
        # regex (`.+?` non-greedy on key) handles it correctly.
        assert "VEP-command-line" in vep_header.metadata
        cmd = vep_header.metadata["VEP-command-line"]
        assert isinstance(cmd, list) and "vep --input_file" in cmd[0]


class TestVEPCSQValueParsing:
    """The CSQ INFO value itself is a comma-separated list of per-transcript
    annotation blocks; each block is pipe-delimited. Our INFO parser splits
    on commas and hands back a list of strings — downstream code splits the
    pipes per record."""

    def test_single_transcript_record(self, vep_header):
        info = vep_header.parse_info(
            "AC=1;AF=0.5;AN=2;DP=42;"
            "CSQ=G|missense_variant|MODERATE|GENE1|ENSG00000000001|Transcript|"
            "ENST00000000001|protein_coding|2/9||ENST00000000001.1:c.34A>G|"
            "ENSP00000000001.1:p.Met12Val|34|34|12|M/V|Atg/Gtg|rs1||1||HGNC|1")
        assert len(info["CSQ"]) == 1
        record = info["CSQ"][0]
        # Pipe-split should produce 23 subfields matching the format spec.
        parts = record.split("|")
        assert len(parts) == 23
        assert parts[0] == "G"                  # Allele
        assert parts[1] == "missense_variant"   # Consequence
        assert parts[2] == "MODERATE"           # IMPACT
        assert parts[3] == "GENE1"              # SYMBOL
        # HGVSc preserves `>` and `:` and `.`.
        assert parts[10] == "ENST00000000001.1:c.34A>G"
        assert parts[15] == "M/V"               # Amino_acids with `/`

    def test_multi_transcript_comma_split(self, vep_header):
        # Two-transcript TP53 stop_gained — comma between the two records.
        info = vep_header.parse_info(
            "AC=1;AF=0.5;AN=2;DP=88;"
            "CSQ=T|stop_gained|HIGH|TP53|ENSG00000141510|Transcript|"
            "ENST00000269305|protein_coding|7/11||"
            "ENST00000269305.9:c.844C>T|ENSP00000269305.4:p.Arg282Ter|"
            "844|844|282|R/*|Cga/Tga|rs28934578||-1||HGNC|11998,"
            "T|stop_gained|HIGH|TP53|ENSG00000141510|Transcript|"
            "ENST00000420246|protein_coding|6/10||"
            "ENST00000420246.6:c.715C>T|ENSP00000391127.2:p.Arg239Ter|"
            "715|715|239|R/*|Cga/Tga|rs28934578||-1||HGNC|11998")
        assert len(info["CSQ"]) == 2
        # Each record has 23 pipe-separated subfields.
        for record in info["CSQ"]:
            assert len(record.split("|")) == 23
        # Different transcripts -> different ENST IDs.
        assert info["CSQ"][0].split("|")[6] == "ENST00000269305"
        assert info["CSQ"][1].split("|")[6] == "ENST00000420246"

    def test_intergenic_record_with_many_empty_subfields(self, vep_header):
        # Real VEP output emits exactly N-1 pipes for N format subfields,
        # even when most subfields are empty. With 23 subfields declared
        # in the format spec we expect 22 pipes / 23 fields after split,
        # of which only Allele/Consequence/IMPACT have content.
        info = vep_header.parse_info(
            "AC=1;AF=0.5;AN=2;DP=18;"
            "CSQ=T|intergenic_variant|MODIFIER||||||||||||||||||||")
        assert len(info["CSQ"]) == 1
        parts = info["CSQ"][0].split("|")
        assert len(parts) == 23
        assert parts[:3] == ["T", "intergenic_variant", "MODIFIER"]
        assert all(p == "" for p in parts[3:])

    def test_percent_encoded_value_is_preserved_verbatim(self, vep_header):
        # VCF 4.3 §1.2 specifies %3D for `=`. Both we and PyVCF3 keep it as
        # the literal string — see _split_struct_fields docstring.
        info = vep_header.parse_info(
            "AC=1;AF=0.25;AN=4;DP=60;"
            "CSQ=A|synonymous_variant|LOW|TP53|ENSG00000141510|Transcript|"
            "ENST00000269305|protein_coding|7/11||"
            "ENST00000269305.9:c.825G>A|ENSP00000269305.4:p.Lys275%3D|"
            "825|825|275|K|aaG/aaA|||-1||HGNC|11998")
        record = info["CSQ"][0]
        assert "p.Lys275%3D" in record   # not decoded — literal


class TestVEPEndToEnd:
    def test_load_vcf_splits_multiallelic_and_keeps_csq(self):
        vc = load_vcf(VEP, genome="GRCh38")
        # 4 lines, but chr17:7676200 has ALT="A,C" → splits to 2 variants.
        # So total = 5 Variant objects.
        assert len(vc) == 5
        # CSQ should be a list of strings on every variant's metadata.
        meta = vc.source_to_metadata_dict[VEP]
        for variant in vc:
            csq = meta[variant]["info"]["CSQ"]
            assert isinstance(csq, list)
            assert all(isinstance(r, str) for r in csq)
