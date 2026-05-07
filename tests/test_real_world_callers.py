# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""
Hand-verified expectations against fixtures from real somatic callers.

Expected values are pinned to specific field shapes that *should* survive
parsing — so a regression on, say, MuTect2 phasing or Strelka2 per-base
tier counts shows up as a meaningful failure against the published caller
format.

Fixtures live under ``tests/data/real_callers/``; see the README there for
provenance (which fixtures are verbatim from public corpora vs. which are
format-faithful reconstructions and why).
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
# GATK4 MuTect2 — verbatim from broadinstitute/gatk @ 3021e6924aeb
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def mutect2_header():
    return VCFHeader.from_path(MUTECT2)


class TestMuTect2Header:
    def test_source_is_mutect2(self, mutect2_header):
        assert mutect2_header.metadata["fileformat"] == "VCFv4.2"
        assert mutect2_header.metadata["source"] == ["Mutect2"]

    def test_metadata_key_with_internal_whitespace(self, mutect2_header):
        # ``##Mutect Version=2.1`` — the key contains a space. PyVCF3's
        # `.+?` non-greedy regex tolerates this; ours mirrors that.
        assert mutect2_header.metadata.get("Mutect Version") == ["2.1"]

    def test_filtering_status_metadata_with_long_value(self, mutect2_header):
        # ##filtering_status= contains a long human-readable warning string
        # with internal whitespace and punctuation.
        status = mutect2_header.metadata["filtering_status"]
        assert isinstance(status, list)
        assert "FilterMutectCalls" in status[0]

    def test_per_alt_info_fields(self, mutect2_header):
        # TLOD/POPAF/N_ART_LOD/NLOD/GERMQ are per-alt (Number=A → -1).
        for key in ("TLOD", "POPAF", "N_ART_LOD", "NLOD", "GERMQ"):
            assert mutect2_header.info_fields[key].number == -1, key

    def test_per_allele_format_fields(self, mutect2_header):
        # AD/F1R2/F2R1/MBQ/MFRL are Number=R (ref + each alt).
        for key in ("AD", "F1R2", "F2R1", "MBQ", "MFRL"):
            assert mutect2_header.format_fields[key].number == -3, key

    def test_str_flag_and_rpa_variable(self, mutect2_header):
        assert mutect2_header.info_fields["STR"].type == "Flag"
        assert mutect2_header.info_fields["STR"].number == 0
        assert mutect2_header.info_fields["RPA"].number is None  # Number=.

    def test_potential_polymorphic_numt_is_string(self, mutect2_header):
        # GATK4 MuTect2's mitochondrial mode emits a String "true"/"false"
        # tag for sites likely to be NuMTs (nuclear copies of mtDNA).
        f = mutect2_header.format_fields["POTENTIAL_POLYMORPHIC_NUMT"]
        assert f.type == "String"
        assert f.number == 1

    def test_single_sample(self, mutect2_header):
        assert mutect2_header.samples == ["NA12878"]


class TestMuTect2InfoParsing:
    def test_simple_snv_per_alt_arrays(self, mutect2_header):
        # chrM:152 T>C — single ALT, Number=A fields are 1-element lists.
        info = mutect2_header.parse_info(
            "DP=1582;ECNT=1;TLOD=5266.19;POPAF=5.000e-08;OCM=0")
        assert info["DP"] == 1582
        assert info["ECNT"] == 1
        assert info["TLOD"] == [5266.19]
        assert info["POPAF"] == [5e-08]
        assert info["OCM"] == 0

    def test_quadallelic_per_alt_arrays(self, mutect2_header):
        # chrM:302 A>AC,C,ACC — three ALTs, TLOD/POPAF have 3 values each.
        info = mutect2_header.parse_info(
            "DP=659;ECNT=4;TLOD=891.23,10.66,67.66;"
            "POPAF=5.000e-08,5.000e-08,5.000e-08;OCM=0")
        assert info["TLOD"] == [891.23, 10.66, 67.66]
        assert info["POPAF"] == [5e-08, 5e-08, 5e-08]

    def test_str_with_repeat_unit(self, mutect2_header):
        # chrM:310 T>TC — STR site with mononucleotide repeat (RU=C).
        info = mutect2_header.parse_info(
            "DP=705;ECNT=4;TLOD=1974.89;POPAF=5.000e-08;"
            "RPA=5,6;RU=C;STR;OCM=0")
        assert info["STR"] is True       # Flag presence
        assert info["RU"] == "C"         # Number=1 -> scalar string
        assert info["RPA"] == [5, 6]     # Number=. -> list (ref=5, alt=6)


class TestMuTect2SampleParsing:
    def test_simple_snv_with_potential_numt_tag(self, mutect2_header):
        # chrM:152 includes POTENTIAL_POLYMORPHIC_NUMT="true" as a string.
        out = mutect2_header.parse_samples(
            ["0/1:3,1556:0.998:2,777:1,779:30,30:16270,369:60,60:42:true"],
            "GT:AD:AF:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:POTENTIAL_POLYMORPHIC_NUMT")
        d = out["NA12878"]
        assert d["GT"] == "0/1"
        assert d["AD"] == [3, 1556]                              # Number=R
        assert d["AF"] == [0.998]                                # Number=A
        assert d["MBQ"] == [30, 30]
        assert d["MFRL"] == [16270, 369]
        assert d["MMQ"] == [60, 60]
        assert d["MPOS"] == [42]
        assert d["POTENTIAL_POLYMORPHIC_NUMT"] == "true"        # Type=String

    def test_quadallelic_four_way_genotype(self, mutect2_header):
        # chrM:302 — 4 alleles total, GT is 0/1/2/3 (real production output).
        out = mutect2_header.parse_samples(
            ["0/1/2/3:5,401,67,49:0.768,0.128,0.094:2,163,35,20:"
             "3,238,32,29:20,20,30,20:419,316,340,278:60,60,60,60:41,33,38"],
            "GT:AD:AF:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS")
        d = out["NA12878"]
        # Quad-allelic GT preserved as raw string.
        assert d["GT"] == "0/1/2/3"
        # AD has 4 elements (1 ref + 3 alts).
        assert d["AD"] == [5, 401, 67, 49]
        # AF has 3 elements (per ALT).
        assert d["AF"] == [0.768, 0.128, 0.094]
        # MPOS Number=A → 3 elements.
        assert d["MPOS"] == [41, 33, 38]
        # MBQ Number=R → 4 elements.
        assert d["MBQ"] == [20, 20, 30, 20]
        # MMQ is declared Number=A in the header but this 2020-era GATK4
        # build actually emits Number=R-style values (4 elements including
        # ref). PyVCF3 and we both just split-and-parse without validating
        # against the declared count, so both return 4 elements. Documenting
        # this as a real-world GATK quirk — varcode shouldn't reject such
        # output, and the test pins the lenient behaviour.
        assert d["MMQ"] == [60, 60, 60, 60]


class TestMuTect2EndToEnd:
    def test_load_vcf_keeps_all_unfiltered_records(self):
        # The fixture has no PASS filter on any record (FILTER="."), so
        # only_passing keeps everything: 7 input lines * (1 + #alts beyond 1)
        # alt-splits = 7 + 4 (chrM:302 has 3 alts -> 2 extra; chrM:400 too)
        # = 7 + 2 + 2 = 11.
        vc = load_vcf(MUTECT2, genome="GRCh37")
        assert len(vc) == 11

    def test_filter_preserves_dot_as_none(self):
        # Records with FILTER="." should have filter=None in metadata.
        vc = load_vcf(MUTECT2, genome="GRCh37")
        meta = vc.source_to_metadata_dict[MUTECT2]
        for variant in vc:
            assert meta[variant]["filter"] is None


# ---------------------------------------------------------------------------
# Strelka2 (format-faithful reconstruction — see README)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def strelka2_header():
    return VCFHeader.from_path(STRELKA2)


class TestStrelka2Header:
    def test_strelka2_metadata(self, strelka2_header):
        assert strelka2_header.metadata["source"] == ["strelka"]
        assert strelka2_header.metadata["source_version"] == ["2.9.10"]

    def test_cmdline_metadata_with_paths_and_flags(self, strelka2_header):
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
        assert out["TUMOR"]["AU"] == [11, 11]    # 11 alt-allele reads
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
        meta = vc.source_to_metadata_dict[STRELKA2]
        for variant in vc:
            for sample_name, fields in meta[variant]["sample_info"].items():
                assert "GT" not in fields, (
                    f"Strelka2 SNV records should have no GT; got {sample_name}")
                # Per-base tier counts should be present.
                assert "AU" in fields and "CU" in fields
                assert "GU" in fields and "TU" in fields


# ---------------------------------------------------------------------------
# VEP-annotated CSQ (format-faithful reconstruction — see README)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def vep_header():
    return VCFHeader.from_path(VEP)


class TestVEPHeader:
    def test_csq_description_with_pipes_survives_quoted_split(self, vep_header):
        # The CSQ Description embeds the pipe-delimited subfield list. Our
        # quote-aware splitter must capture the whole thing intact.
        csq = vep_header.info_fields["CSQ"]
        assert csq.number is None      # Number=. (per-transcript list)
        assert csq.type == "String"
        assert "Format: Allele|Consequence" in csq.description

    def test_csq_subfield_count(self, vep_header):
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
    """The CSQ INFO value is a comma-separated list of per-transcript blocks;
    each block is pipe-delimited. Our INFO parser splits on commas; downstream
    code splits the pipes per record."""

    def test_single_transcript_record(self, vep_header):
        info = vep_header.parse_info(
            "AC=1;AF=0.5;AN=2;DP=42;"
            "CSQ=G|missense_variant|MODERATE|GENE1|ENSG00000000001|Transcript|"
            "ENST00000000001|protein_coding|2/9||ENST00000000001.1:c.34A>G|"
            "ENSP00000000001.1:p.Met12Val|34|34|12|M/V|Atg/Gtg|rs1||1||HGNC|1")
        assert len(info["CSQ"]) == 1
        record = info["CSQ"][0]
        parts = record.split("|")
        assert len(parts) == 23
        assert parts[0] == "G"                  # Allele
        assert parts[1] == "missense_variant"   # Consequence
        assert parts[2] == "MODERATE"           # IMPACT
        # HGVSc preserves `>` and `:` and `.`.
        assert parts[10] == "ENST00000000001.1:c.34A>G"
        assert parts[15] == "M/V"               # Amino_acids with `/`

    def test_multi_transcript_comma_split(self, vep_header):
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
        for record in info["CSQ"]:
            assert len(record.split("|")) == 23
        # Different transcripts -> different ENST IDs.
        assert info["CSQ"][0].split("|")[6] == "ENST00000269305"
        assert info["CSQ"][1].split("|")[6] == "ENST00000420246"

    def test_intergenic_record_with_many_empty_subfields(self, vep_header):
        # Real VEP output emits N-1 pipes for N format subfields, even when
        # most subfields are empty.
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
        assert "p.Lys275%3D" in info["CSQ"][0]   # not decoded — literal


class TestVEPEndToEnd:
    def test_load_vcf_splits_multiallelic_and_keeps_csq(self):
        vc = load_vcf(VEP, genome="GRCh38")
        # 4 lines; chr17:7676200 ALT="A,C" splits to 2 variants. Total = 5.
        assert len(vc) == 5
        meta = vc.source_to_metadata_dict[VEP]
        for variant in vc:
            csq = meta[variant]["info"]["CSQ"]
            assert isinstance(csq, list)
            assert all(isinstance(r, str) for r in csq)
