# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""
Oracle tests for varcode.vcf_parsing.

PyVCF3 is the reference implementation: for each fixture VCF, we parse with
both PyVCF3 and our new VCFHeader, then assert the outputs match for every
header field, every row's INFO, and every row's per-sample data. Targeted
unit tests cover edge cases that the fixtures don't reach.
"""
from __future__ import annotations

import glob
import gzip
import os
from collections import OrderedDict

import pytest

import vcf as pyvcf  # PyVCF3 — kept as a test-only dependency oracle.

from varcode.vcf_parsing import (
    FieldDef,
    VCFHeader,
    RESERVED_INFO,
    RESERVED_FORMAT,
    _split_struct_fields,
    _parse_field_decl,
)


DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

# A fixture that's intentionally header-less. PyVCF3 mis-parses it (treats the
# first data row as the column-header row, producing garbage sample names).
# varcode's own load_vcf tests already skip it (see tests/test_vcf_output.py),
# so we exclude it from the oracle comparison rather than mirror PyVCF3's bug.
_KNOWN_BROKEN_IN_PYVCF3 = {"mutect-example-headerless.vcf"}

FIXTURES = [
    p for p in (
        sorted(glob.glob(os.path.join(DATA_DIR, "**", "*.vcf"), recursive=True)) +
        sorted(glob.glob(os.path.join(DATA_DIR, "**", "*.vcf.gz"), recursive=True)))
    if os.path.basename(p) not in _KNOWN_BROKEN_IN_PYVCF3]


def _pyvcf_reader(path):
    return pyvcf.Reader(filename=path, strict_whitespace=True)


def _pyvcf_call_to_dict(call):
    """Match the shape produced by varcode.vcf.pyvcf_calls_to_sample_info_list."""
    return OrderedDict(call.data._asdict())


# ---------------------------------------------------------------------------
# Oracle tests across every fixture.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("path", FIXTURES, ids=os.path.basename)
def test_header_samples_match_pyvcf(path):
    ours = VCFHeader.from_path(path)
    theirs = _pyvcf_reader(path)
    assert ours.samples == theirs.samples


@pytest.mark.parametrize("path", FIXTURES, ids=os.path.basename)
def test_header_singular_metadata_matches_pyvcf(path):
    ours = VCFHeader.from_path(path)
    theirs = _pyvcf_reader(path)
    # PyVCF3 stores SINGULAR_METADATA keys (fileformat / fileDate / reference)
    # as plain strings; everything else is a list. We mirror that.
    for key in ("fileformat", "fileDate", "reference"):
        if key in theirs.metadata:
            assert ours.metadata.get(key) == theirs.metadata[key], (
                "metadata[%r] mismatch in %s" % (key, path))


@pytest.mark.parametrize("path", FIXTURES, ids=os.path.basename)
def test_header_info_fields_match_pyvcf(path):
    ours = VCFHeader.from_path(path)
    theirs = _pyvcf_reader(path)
    assert list(ours.info_fields) == list(theirs.infos), (
        "INFO field order/keys differ in %s" % path)
    for key, our_def in ours.info_fields.items():
        their_def = theirs.infos[key]
        assert our_def.id == their_def.id
        assert our_def.number == their_def.num, (
            "INFO[%s].number: %r vs %r" % (key, our_def.number, their_def.num))
        assert our_def.type == their_def.type, (
            "INFO[%s].type: %r vs %r" % (key, our_def.type, their_def.type))


@pytest.mark.parametrize("path", FIXTURES, ids=os.path.basename)
def test_header_format_fields_match_pyvcf(path):
    ours = VCFHeader.from_path(path)
    theirs = _pyvcf_reader(path)
    assert list(ours.format_fields) == list(theirs.formats), (
        "FORMAT field order/keys differ in %s" % path)
    for key, our_def in ours.format_fields.items():
        their_def = theirs.formats[key]
        assert our_def.id == their_def.id
        assert our_def.number == their_def.num, (
            "FORMAT[%s].number: %r vs %r" % (key, our_def.number, their_def.num))
        assert our_def.type == their_def.type, (
            "FORMAT[%s].type: %r vs %r" % (key, our_def.type, their_def.type))


def _iter_data_rows(path):
    """Yield non-header rows of a VCF as tab-split column lists."""
    open_fn = gzip.open if path.endswith(".gz") else open
    with open_fn(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            yield line.rstrip("\n").rstrip("\r").split("\t")


@pytest.mark.parametrize("path", FIXTURES, ids=os.path.basename)
def test_parse_info_matches_pyvcf_for_every_row(path):
    ours = VCFHeader.from_path(path)
    raw_info_strings = [cols[7] for cols in _iter_data_rows(path) if len(cols) >= 8]
    expected = list(_pyvcf_reader(path))
    assert len(raw_info_strings) == len(expected)
    for raw, record in zip(raw_info_strings, expected):
        our_info = ours.parse_info(raw)
        assert our_info == record.INFO, (
            "INFO mismatch on %s row %s: ours=%r theirs=%r" % (
                path, raw, our_info, record.INFO))


@pytest.mark.parametrize("path", FIXTURES, ids=os.path.basename)
def test_parse_samples_matches_pyvcf_for_every_row(path):
    ours = VCFHeader.from_path(path)
    if not ours.samples:
        pytest.skip("no samples in %s" % path)
    raw_rows = [
        (cols[8], cols[9:]) for cols in _iter_data_rows(path) if len(cols) >= 10
    ]
    expected = list(_pyvcf_reader(path))
    for (fmt, sample_cells), record in zip(raw_rows, expected):
        if fmt == "." or fmt is None:
            continue
        ours_parsed = ours.parse_samples(sample_cells, fmt)
        theirs_parsed = OrderedDict(
            (call.sample, _pyvcf_call_to_dict(call))
            for call in record.samples)
        assert list(ours_parsed) == list(theirs_parsed)
        for sample_name in theirs_parsed:
            our_row = ours_parsed[sample_name]
            their_row = theirs_parsed[sample_name]
            assert list(our_row) == list(their_row), (
                "FORMAT field order mismatch for %s in %s" % (sample_name, path))
            for k in their_row:
                assert our_row[k] == their_row[k], (
                    "%s sample=%s field=%s ours=%r theirs=%r" % (
                        path, sample_name, k, our_row[k], their_row[k]))


# ---------------------------------------------------------------------------
# Targeted unit tests.
# ---------------------------------------------------------------------------

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


class TestParseInfoEdgeCases:
    """Cases the fixtures may not exercise; both parsers tested explicitly."""

    def _both(self, header_lines, info_str, fixture_samples=None):
        """
        Build a tiny VCFHeader + a PyVCF3 Reader from the same in-memory header,
        return (ours, theirs) parsed dicts.
        """
        if fixture_samples is None:
            fixture_samples = ["S1"]
        all_lines = list(header_lines) + [
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
            "\t".join(fixture_samples) + "\n"
        ]
        ours = VCFHeader.from_lines(all_lines)
        # PyVCF3 needs a real stream of header + at least one record? Let's give
        # it the header + a dummy record so it can construct the Reader.
        record_line = "1\t1\t.\tA\tT\t.\tPASS\t" + info_str + "\tGT\t0/0\n"
        import io
        stream = io.StringIO("".join(all_lines + [record_line]))
        theirs = pyvcf.Reader(fsock=stream, strict_whitespace=True)
        their_info = theirs._parse_info(info_str)
        return ours.parse_info(info_str), their_info

    def test_info_dot_returns_empty_dict(self):
        ours, theirs = self._both(
            ['##fileformat=VCFv4.2\n'], ".")
        assert ours == {}
        assert ours == theirs

    def test_flag_present_with_no_value(self):
        header = [
            '##fileformat=VCFv4.2\n',
            '##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP">\n',
        ]
        ours, theirs = self._both(header, "DB")
        assert ours == theirs == {"DB": True}

    def test_reserved_flag_undeclared(self):
        # SOMATIC is in RESERVED_INFO. With no header decl, both parsers
        # should still recognize it as a Flag.
        ours, theirs = self._both(['##fileformat=VCFv4.2\n'], "SOMATIC")
        assert ours == theirs == {"SOMATIC": True}

    def test_integer_number_one_unwraps(self):
        header = [
            '##fileformat=VCFv4.2\n',
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">\n',
        ]
        ours, theirs = self._both(header, "DP=42")
        assert ours == theirs == {"DP": 42}

    def test_integer_number_a_stays_list(self):
        header = [
            '##fileformat=VCFv4.2\n',
            '##INFO=<ID=AC,Number=A,Type=Integer,Description="alt count">\n',
        ]
        ours, theirs = self._both(header, "AC=3,5,7")
        assert ours == theirs == {"AC": [3, 5, 7]}

    def test_float_with_missing_value_in_list(self):
        header = [
            '##fileformat=VCFv4.2\n',
            '##INFO=<ID=AF,Number=A,Type=Float,Description="freq">\n',
        ]
        ours, theirs = self._both(header, "AF=0.5,.")
        assert ours == theirs == {"AF": [0.5, None]}

    def test_undeclared_string_field_stays_list(self):
        # Fully unknown field with a value -> typed as String, list-shaped
        # (no num=1 unwrap because there's no header decl).
        ours, theirs = self._both(['##fileformat=VCFv4.2\n'], "UNKNOWN=foo")
        assert ours == theirs == {"UNKNOWN": ["foo"]}

    def test_reserved_integer_stays_list_when_undeclared(self):
        # DP is reserved Integer; without a header decl PyVCF3 doesn't unwrap.
        ours, theirs = self._both(['##fileformat=VCFv4.2\n'], "DP=42")
        assert ours == theirs == {"DP": [42]}

    def test_multiple_fields(self):
        header = [
            '##fileformat=VCFv4.2\n',
            '##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP">\n',
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">\n',
            '##INFO=<ID=AF,Number=A,Type=Float,Description="freq">\n',
        ]
        ours, theirs = self._both(header, "DB;DP=100;AF=0.1,0.2")
        assert ours == theirs

    def test_bare_declared_integer_diverges_from_pyvcf3(self):
        """Deliberate divergence: a key declared as Integer (or any non-Flag
        type) that appears bare (no `=`) crashes PyVCF3 with IndexError but
        we return ``{key: True}`` — flag-like presence. This only surfaces
        on malformed input; the more useful answer is to keep parsing.
        """
        header = [
            '##fileformat=VCFv4.2\n',
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">\n',
        ]
        all_lines = list(header) + [
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"]
        ours = VCFHeader.from_lines(all_lines)
        record_line = "1\t1\t.\tA\tT\t.\tPASS\tDP\tGT\t0/0\n"
        import io
        stream = io.StringIO("".join(all_lines + [record_line]))
        theirs = pyvcf.Reader(fsock=stream, strict_whitespace=True)

        assert ours.parse_info("DP") == {"DP": True}
        with pytest.raises(IndexError):
            theirs._parse_info("DP")


class TestParseSamplesEdgeCases:
    def _make(self, header_lines, samples=("S1", "S2")):
        all_lines = list(header_lines) + [
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
            "\t".join(samples) + "\n"
        ]
        return VCFHeader.from_lines(all_lines)

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
