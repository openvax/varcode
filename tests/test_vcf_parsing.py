# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""
Tests for ``varcode.vcf_parsing``.

The expected values in these tests are PyVCF3's outputs, captured once and
hardcoded as literals — PyVCF3 itself is not a test dependency. Per-fixture
oracle coverage (4.2/4.3 spec examples, MuTect2/Strelka2/VEP) lives in
``test_vcf_spec_examples.py`` and ``test_real_world_callers.py``; this file
covers parser internals and edge cases the fixtures don't reach.
"""
from __future__ import annotations

import os
from collections import OrderedDict

import pytest

from varcode.vcf_parsing import (
    FieldDef,
    VCFHeader,
    RESERVED_INFO,
    RESERVED_FORMAT,
    _split_struct_fields,
    _parse_field_decl,
)


DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


class TestStructFieldSplit:
    def test_simple_fields(self):
        d = _split_struct_fields("ID=AF,Number=A,Type=Float")
        assert d == OrderedDict([("ID", "AF"), ("Number", "A"), ("Type", "Float")])

    def test_quoted_value_with_commas(self):
        body = 'ID=X,Number=1,Type=Integer,Description="a, b, c"'
        d = _split_struct_fields(body)
        assert d["Description"] == "a, b, c"
        assert d["ID"] == "X"

    def test_optional_source_and_version(self):
        body = (
            'ID=AF,Number=A,Type=Float,Description="freq",'
            'Source="dbSNP",Version="3"')
        d = _split_struct_fields(body)
        assert d["Source"] == "dbSNP"
        assert d["Version"] == "3"

    def test_empty_description(self):
        body = 'ID=X,Number=1,Type=Integer,Description=""'
        d = _split_struct_fields(body)
        assert d["Description"] == ""

    def test_field_decl_dataclass(self):
        fd = _parse_field_decl(
            'ID=AF,Number=A,Type=Float,Description="freq, allele"')
        assert fd == FieldDef(
            id="AF", number=-1, type="Float",
            description="freq, allele")

    def test_empty_body_yields_empty_dict(self):
        assert _split_struct_fields("") == OrderedDict()

    def test_trailing_comma_is_tolerated(self):
        # Real VCF tools sometimes emit a trailing comma; the regex's
        # finditer just doesn't match past it. No spurious empty key.
        d = _split_struct_fields("ID=X,Number=1,")
        assert d == OrderedDict([("ID", "X"), ("Number", "1")])

    def test_whitespace_around_equals(self):
        # Permissive: PyVCF3's per-line regex doesn't allow this, but we do.
        # No real fixture exercises it; it's a robustness extra.
        d = _split_struct_fields("ID = X,Number = 1")
        assert d == OrderedDict([("ID", "X"), ("Number", "1")])

    def test_plain_value_with_embedded_equals(self):
        # The plain branch is `[^,]*`, so `=` inside an unquoted value is
        # captured as part of the value (the key/value split happens once
        # at the leading `=`).
        d = _split_struct_fields("ID=X,URL=http://example.com/?a=1&b=2")
        assert d["URL"] == "http://example.com/?a=1&b=2"

    def test_mixed_quoted_and_bare_preserves_order(self):
        body = (
            'ID=AF,Number=A,Type=Float,'
            'Description="freq, per ALT",'
            'Source="dbSNP",Version=3')
        d = _split_struct_fields(body)
        assert list(d) == [
            "ID", "Number", "Type", "Description", "Source", "Version"]
        assert d["Description"] == "freq, per ALT"
        assert d["Version"] == "3"

    def test_quoted_value_preserves_internal_whitespace(self):
        # A quoted value isn't stripped; bare values are. The trio fixture's
        # contig line relies on this for ``species="Homo sapiens"``.
        d = _split_struct_fields('ID=20,species="Homo sapiens"')
        assert d["species"] == "Homo sapiens"

    def test_bare_value_is_stripped(self):
        d = _split_struct_fields("ID=  X  ,Number= 1 ")
        assert d == OrderedDict([("ID", "X"), ("Number", "1")])


class TestNumberEncoding:
    @pytest.mark.parametrize("raw,expected", [
        (".", None),
        ("A", -1),
        ("G", -2),
        ("R", -3),
        ("0", 0),
        ("1", 1),
        ("3", 3),
    ])
    def test_number_specials(self, raw, expected):
        from varcode.vcf_parsing import _parse_number
        assert _parse_number(raw) == expected


def _parser(header_lines, samples=("S1",)):
    """Build a VCFHeader from in-memory header lines plus a column-header row."""
    all_lines = list(header_lines) + [
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
        "\t".join(samples) + "\n"
    ]
    return VCFHeader.from_lines(all_lines)


class TestParseInfoEdgeCases:
    """Edge cases for ``VCFHeader.parse_info``.

    Expected values were captured from PyVCF3 once and hardcoded here; PyVCF3
    is not imported at test time.
    """

    def test_info_dot_returns_empty_dict(self):
        h = _parser(['##fileformat=VCFv4.2\n'])
        assert h.parse_info(".") == {}

    def test_flag_present_with_no_value(self):
        h = _parser([
            '##fileformat=VCFv4.2\n',
            '##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP">\n',
        ])
        assert h.parse_info("DB") == {"DB": True}

    def test_reserved_flag_undeclared(self):
        # SOMATIC is in RESERVED_INFO. With no header decl, it's still
        # recognized as a Flag (matches PyVCF3's RESERVED_INFO fallback).
        h = _parser(['##fileformat=VCFv4.2\n'])
        assert h.parse_info("SOMATIC") == {"SOMATIC": True}

    def test_integer_number_one_unwraps(self):
        h = _parser([
            '##fileformat=VCFv4.2\n',
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">\n',
        ])
        assert h.parse_info("DP=42") == {"DP": 42}

    def test_integer_number_a_stays_list(self):
        h = _parser([
            '##fileformat=VCFv4.2\n',
            '##INFO=<ID=AC,Number=A,Type=Integer,Description="alt count">\n',
        ])
        assert h.parse_info("AC=3,5,7") == {"AC": [3, 5, 7]}

    def test_float_with_missing_value_in_list(self):
        h = _parser([
            '##fileformat=VCFv4.2\n',
            '##INFO=<ID=AF,Number=A,Type=Float,Description="freq">\n',
        ])
        assert h.parse_info("AF=0.5,.") == {"AF": [0.5, None]}

    def test_undeclared_string_field_stays_list(self):
        # Fully unknown field with a value -> typed as String, list-shaped
        # (no num=1 unwrap because there's no header decl).
        h = _parser(['##fileformat=VCFv4.2\n'])
        assert h.parse_info("UNKNOWN=foo") == {"UNKNOWN": ["foo"]}

    def test_reserved_integer_stays_list_when_undeclared(self):
        # DP is reserved Integer; without a header decl PyVCF3 doesn't unwrap.
        h = _parser(['##fileformat=VCFv4.2\n'])
        assert h.parse_info("DP=42") == {"DP": [42]}

    def test_multiple_fields(self):
        h = _parser([
            '##fileformat=VCFv4.2\n',
            '##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP">\n',
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">\n',
            '##INFO=<ID=AF,Number=A,Type=Float,Description="freq">\n',
        ])
        assert h.parse_info("DB;DP=100;AF=0.1,0.2") == {
            "DB": True,
            "DP": 100,
            "AF": [0.1, 0.2],
        }

    def test_bare_declared_non_flag_returns_presence_marker(self):
        """Deliberate divergence from PyVCF3: a key declared as Integer (or
        any non-Flag type) appearing bare (no `=`) returns ``{key: True}``
        — flag-like presence — instead of crashing with IndexError. This
        only surfaces on malformed input; the more useful answer is to
        keep parsing.
        """
        h = _parser([
            '##fileformat=VCFv4.2\n',
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">\n',
        ])
        assert h.parse_info("DP") == {"DP": True}


class TestParseSamplesEdgeCases:
    def _make(self, header_lines, samples=("S1", "S2")):
        return _parser(header_lines, samples=samples)

    def test_gt_is_passthrough(self):
        h = self._make([
            '##fileformat=VCFv4.2\n',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="genotype">\n',
        ])
        out = h.parse_samples(["0/1", "1|1"], "GT")
        assert out["S1"]["GT"] == "0/1"
        assert out["S2"]["GT"] == "1|1"  # phased pipe preserved

    def test_missing_cell_becomes_none(self):
        h = self._make([
            '##fileformat=VCFv4.2\n',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">\n',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="depth">\n',
        ])
        out = h.parse_samples(["0/1:.", "0/0:5"], "GT:DP")
        assert out["S1"]["DP"] is None
        assert out["S2"]["DP"] == 5

    def test_multivalue_integer_field(self):
        h = self._make([
            '##fileformat=VCFv4.2\n',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">\n',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="depths">\n',
        ])
        out = h.parse_samples(["0/1:10,20", "0/0:5,0"], "GT:AD")
        assert out["S1"]["AD"] == [10, 20]
        assert out["S2"]["AD"] == [5, 0]

    def test_trailing_missing_fields_default_to_none(self):
        h = self._make([
            '##fileformat=VCFv4.2\n',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">\n',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="depth">\n',
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="qual">\n',
        ])
        # Sample 1 has only GT:DP, GQ is missing.
        out = h.parse_samples(["0/1:10", "0/0:5:99"], "GT:DP:GQ")
        assert out["S1"]["GQ"] is None
        assert out["S2"]["GQ"] == 99

    def test_undeclared_format_field_treated_as_string(self):
        h = self._make([
            '##fileformat=VCFv4.2\n',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">\n',
        ])
        out = h.parse_samples(["0/1:foo", "0/0:bar"], "GT:UNKNOWN")
        # Undeclared -> String, num=None -> stays as a list (post-comma-split).
        assert out["S1"]["UNKNOWN"] == ["foo"]
        assert out["S2"]["UNKNOWN"] == ["bar"]

    def test_sample_count_mismatch_raises(self):
        h = self._make([
            '##fileformat=VCFv4.2\n',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">\n',
        ], samples=("S1", "S2"))
        with pytest.raises(ValueError):
            h.parse_samples(["0/1"], "GT")


class TestReservedTables:
    def test_reserved_info_includes_spec_keys(self):
        for k in ("AC", "AF", "AN", "DP", "DB", "SOMATIC"):
            assert k in RESERVED_INFO

    def test_reserved_format_includes_spec_keys(self):
        for k in ("GT", "DP", "GQ", "PS"):
            assert k in RESERVED_FORMAT


class TestGetMetadata:
    """Convenience accessor that smooths over the str-vs-list[str] asymmetry
    in the metadata dict. Lifted from the GATK MuTect2 fixture which has
    examples of both shapes."""

    @pytest.fixture(scope="class")
    def header(self):
        path = os.path.join(
            DATA_DIR, "real_callers", "mutect2_example.vcf")
        return VCFHeader.from_path(path)

    def test_singular_key_returns_string(self, header):
        # fileformat is in SINGULAR_METADATA -> stored as plain string.
        assert header.get_metadata("fileformat") == "VCFv4.2"

    def test_single_occurrence_list_unwraps(self, header):
        # ##source=Mutect2 -> stored as ["Mutect2"], unwrapped to "Mutect2".
        assert header.get_metadata("source") == "Mutect2"

    def test_key_with_internal_whitespace(self, header):
        # ##Mutect Version=2.1 — embedded space in key, unwrapped to scalar.
        assert header.get_metadata("Mutect Version") == "2.1"

    def test_missing_key_returns_default(self, header):
        assert header.get_metadata("does_not_exist") is None
        assert header.get_metadata("does_not_exist", "fallback") == "fallback"

    def test_multi_valued_key_raises(self):
        # Construct a header with a repeated non-singular key (real example:
        # multi-contig VCFs).
        h = VCFHeader.from_lines([
            "##fileformat=VCFv4.2\n",
            "##contig=<ID=1,length=10>\n",
            "##contig=<ID=2,length=20>\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        ])
        # The full list is still queryable via the underlying dict.
        assert len(h.metadata["contig"]) == 2
        # The accessor refuses to silently pick one.
        with pytest.raises(ValueError, match="2 values"):
            h.get_metadata("contig")


class TestNoRpy2OnImport:
    """`import varcode` must never pull rpy2 or PyVCF3 into sys.modules.

    This was the original motivation for #302. The earlier fix used a lazy
    `__getattr__` shim in `varcode/__init__.py` to defer the import; the
    runtime VCF parser rewrite removed PyVCF3 from the import graph entirely,
    so the invariant is now upheld unconditionally — but it's load-bearing
    enough that we still pin it.
    """

    def test_bare_import_does_not_pull_pyvcf3_or_rpy2(self):
        import subprocess
        import sys
        code = (
            "import sys; "
            "import varcode; "
            "bad = [m for m in sys.modules "
            "if m == 'vcf' or m.startswith('vcf.') or m.startswith('rpy2')]; "
            "assert not bad, bad; "
            "print('ok')"
        )
        result = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        assert "ok" in result.stdout
