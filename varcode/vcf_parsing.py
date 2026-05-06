# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
VCF header / INFO / FORMAT parsing.

Replaces the runtime PyVCF3 dependency that previously dragged rpy2/R into
sys.modules just for header metadata and per-cell INFO/FORMAT decoding (#302).

PyVCF3's behaviour is the reference implementation; tests in
``tests/test_vcf_parsing.py`` pin this module's output to PyVCF3's for the
existing fixture VCFs, so future regressions surface as test failures.
"""
from __future__ import annotations

import gzip
import re
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Iterable, Optional, Union


# Reserved INFO/FORMAT keys per VCF spec. Used as a fallback when a header
# doesn't declare a key. Tables mirror PyVCF3's RESERVED_INFO / RESERVED_FORMAT
# verbatim so we behave identically on minimally-headered VCFs.
RESERVED_INFO = {
    "AA": "String", "AC": "Integer", "AF": "Float", "AN": "Integer",
    "BQ": "Float", "CIGAR": "String", "DB": "Flag", "DP": "Integer",
    "END": "Integer", "H2": "Flag", "H3": "Flag", "MQ": "Float",
    "MQ0": "Integer", "NS": "Integer", "SB": "String", "SOMATIC": "Flag",
    "VALIDATED": "Flag", "1000G": "Flag",
    "IMPRECISE": "Flag", "NOVEL": "Flag", "SVTYPE": "String",
    "SVLEN": "Integer", "CIPOS": "Integer", "CIEND": "Integer",
    "HOMLEN": "Integer", "HOMSEQ": "String", "BKPTID": "String",
    "MEINFO": "String", "METRANS": "String", "DGVID": "String",
    "DBVARID": "String", "DBRIPID": "String", "MATEID": "String",
    "PARID": "String", "EVENT": "String", "CILEN": "Integer",
    "DPADJ": "Integer", "CN": "Integer", "CNADJ": "Integer",
    "CICN": "Integer", "CICNADJ": "Integer",
}

RESERVED_FORMAT = {
    "GT": "String", "DP": "Integer", "FT": "String", "GL": "Float",
    "GLE": "String", "PL": "Integer", "GP": "Float", "GQ": "Integer",
    "HQ": "Integer", "PS": "Integer", "PQ": "Integer", "EC": "Integer",
    "MQ": "Integer",
    "CN": "Integer", "CNQ": "Float", "CNL": "Float", "NQ": "Integer",
    "HAP": "Integer", "AHAP": "Integer",
}

# Header keys whose value is treated as a single string (overwrite on repeat).
# Other ##KEY=VAL lines accumulate into a list[str]. Mirrors PyVCF3's
# SINGULAR_METADATA.
SINGULAR_METADATA = ("fileformat", "fileDate", "reference")

# "Missing" sentinels used by PyVCF3's _map: any of these as a *whole* split
# element becomes None. We replicate the exact set so our typed lists agree
# with PyVCF3's element-by-element.
_MISSING_VALUES = frozenset((".", "", "NA"))

# Number=. / A / G / R encoded the same way PyVCF3 does (so a header round-trip
# yields the same .number value).
_NUMBER_SPECIALS = {".": None, "A": -1, "G": -2, "R": -3}


@dataclass(frozen=True)
class FieldDef:
    """Header declaration for an INFO or FORMAT field."""

    id: str
    number: Optional[int]   # int >=0 for fixed counts; None for ".", -1/-2/-3 for A/G/R
    type: str               # "Integer" | "Float" | "Flag" | "Character" | "String"
    description: str = ""


_INFO_LINE = re.compile(r"^##INFO=<(.+)>\s*$")
_FORMAT_LINE = re.compile(r"^##FORMAT=<(.+)>\s*$")
_META_LINE = re.compile(r"^##(?P<key>.+?)=(?P<val>.*)$")


def _split_struct_fields(body: str) -> "OrderedDict[str, str]":
    """
    Split a ``key=val,key="quoted, with commas",...`` body into a dict.

    Comma-splitting is quote-aware: commas inside double-quoted values don't
    terminate a field. ``Description`` strings frequently embed commas, so a
    naïve ``.split(",")`` would corrupt the parse.
    """
    fields: "OrderedDict[str, str]" = OrderedDict()
    i = 0
    n = len(body)
    while i < n:
        eq = body.find("=", i)
        if eq == -1:
            break
        key = body[i:eq].strip()
        i = eq + 1
        if i < n and body[i] == '"':
            i += 1
            start = i
            while i < n and body[i] != '"':
                if body[i] == "\\" and i + 1 < n:
                    i += 2
                else:
                    i += 1
            val = body[start:i]
            if i < n and body[i] == '"':
                i += 1
            while i < n and body[i] != ",":
                i += 1
        else:
            start = i
            while i < n and body[i] != ",":
                i += 1
            val = body[start:i].strip()
        if i < n and body[i] == ",":
            i += 1
        fields[key] = val
    return fields


def _parse_number(s: str) -> Optional[int]:
    if s in _NUMBER_SPECIALS:
        return _NUMBER_SPECIALS[s]
    return int(s)


def _parse_field_decl(body: str) -> FieldDef:
    f = _split_struct_fields(body)
    return FieldDef(
        id=f["ID"],
        number=_parse_number(f["Number"]),
        type=f["Type"],
        description=f.get("Description", ""),
    )


def _parse_filter(filt_str: str) -> Optional[list]:
    """Mirror PyVCF3's `_parse_filter`: '.' -> None, 'PASS' -> [], else split on ';'."""
    if filt_str == ".":
        return None
    if filt_str == "PASS":
        return []
    return filt_str.split(";")


def _iter_header_lines(fh) -> Iterable[str]:
    """Yield lines until and including the #CHROM column header row."""
    for line in fh:
        yield line
        if line.startswith("#CHROM"):
            break


@dataclass
class VCFHeader:
    """Parsed VCF header: metadata, INFO/FORMAT field declarations, sample names.

    Drop-in replacement for the subset of PyVCF3's `Reader` that varcode
    actually consumed: ``metadata`` dict, ``samples`` list, INFO and FORMAT
    cell parsing.
    """

    metadata: "OrderedDict[str, Union[str, list]]" = field(default_factory=OrderedDict)
    info_fields: "OrderedDict[str, FieldDef]" = field(default_factory=OrderedDict)
    format_fields: "OrderedDict[str, FieldDef]" = field(default_factory=OrderedDict)
    samples: list = field(default_factory=list)

    @classmethod
    def from_path(cls, path: str) -> "VCFHeader":
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            return cls.from_lines(_iter_header_lines(fh))

    @classmethod
    def from_lines(cls, lines: Iterable[str]) -> "VCFHeader":
        h = cls()
        for line in lines:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith("##"):
                m = _INFO_LINE.match(line)
                if m:
                    fd = _parse_field_decl(m.group(1))
                    h.info_fields[fd.id] = fd
                    continue
                m = _FORMAT_LINE.match(line)
                if m:
                    fd = _parse_field_decl(m.group(1))
                    h.format_fields[fd.id] = fd
                    continue
                m = _META_LINE.match(line)
                if m:
                    key = m.group("key")
                    val = m.group("val")
                    if key in SINGULAR_METADATA:
                        h.metadata[key] = val
                    else:
                        bucket = h.metadata.setdefault(key, [])
                        if isinstance(bucket, list):
                            bucket.append(val)
                continue
            if line.startswith("#"):
                # #CHROM column-header row; tab-split and grab samples.
                fields = line[1:].split("\t")
                h.samples = fields[9:] if len(fields) > 9 else []
                break
        return h

    # ---- INFO parsing ------------------------------------------------------

    def parse_info(self, info_str: str) -> dict:
        """
        Parse an INFO column string like ``"AF=0.5;DB;DP=100"`` into a typed dict.

        Pinned to PyVCF3's `_parse_info` semantics; see oracle tests.
        """
        if info_str == ".":
            return {}
        out: dict = {}
        for entry in info_str.split(";"):
            kv = entry.split("=", 1)
            key = kv[0]
            has_val = len(kv) > 1
            decl = self.info_fields.get(key)
            if decl is not None:
                ftype = decl.type
            elif key in RESERVED_INFO:
                ftype = RESERVED_INFO[key]
            else:
                ftype = "String" if has_val else "Flag"

            if ftype == "Flag" or not has_val:
                out[key] = True
                continue

            vals = kv[1].split(",")
            if ftype == "Integer":
                try:
                    val = [int(x) if x not in _MISSING_VALUES else None for x in vals]
                except ValueError:
                    val = [float(x) if x not in _MISSING_VALUES else None for x in vals]
            elif ftype == "Float":
                val = [float(x) if x not in _MISSING_VALUES else None for x in vals]
            else:  # String / Character
                val = [x if x not in _MISSING_VALUES else None for x in vals]

            # PyVCF3 unwraps single-element list iff the *header declared* num == 1,
            # i.e. not for reserved-but-undeclared keys.
            if decl is not None and decl.number == 1:
                val = val[0]
            out[key] = val
        return out

    # ---- FORMAT/sample parsing --------------------------------------------

    def parse_samples(
        self,
        sample_strings: list,
        format_string: str,
    ) -> "OrderedDict[str, OrderedDict[str, object]]":
        """
        Parse colon-delimited per-sample strings using ``format_string``.

        Returns ``OrderedDict[sample_name -> OrderedDict[field_name -> value]]``,
        matching what `pyvcf_calls_to_sample_info_list` produced from PyVCF3's
        ``_parse_samples`` output.
        """
        if len(sample_strings) != len(self.samples):
            raise ValueError(
                "Got %d sample columns but header declares %d samples" % (
                    len(sample_strings), len(self.samples)))

        fields = format_string.split(":")
        # Resolve each field's (type, num) from header → reserved → String/None.
        field_specs = []
        for f in fields:
            decl = self.format_fields.get(f)
            if decl is not None:
                field_specs.append((f, decl.type, decl.number))
            elif f in RESERVED_FORMAT:
                field_specs.append((f, RESERVED_FORMAT[f], None))
            else:
                field_specs.append((f, "String", None))

        result: "OrderedDict[str, OrderedDict[str, object]]" = OrderedDict()
        for sample_name, sample_str in zip(self.samples, sample_strings):
            cells = sample_str.split(":")
            row: "OrderedDict[str, object]" = OrderedDict()
            for i, (fname, ftype, fnum) in enumerate(field_specs):
                cell = cells[i] if i < len(cells) else None
                if fname == "GT":
                    row[fname] = cell
                    continue
                if fname == "FT":
                    row[fname] = _parse_filter(cell) if cell is not None else None
                    continue
                if cell is None or cell == "" or cell == ".":
                    row[fname] = None
                    continue
                if fnum == 1:
                    if ftype == "Integer":
                        try:
                            row[fname] = int(cell)
                        except ValueError:
                            row[fname] = float(cell)
                    elif ftype == "Float":
                        row[fname] = float(cell)
                    else:
                        row[fname] = cell
                    continue
                vals = cell.split(",")
                if ftype == "Integer":
                    try:
                        row[fname] = [int(x) if x not in _MISSING_VALUES else None for x in vals]
                    except ValueError:
                        row[fname] = [float(x) if x not in _MISSING_VALUES else None for x in vals]
                elif ftype == "Float":
                    row[fname] = [float(x) if x not in _MISSING_VALUES else None for x in vals]
                else:
                    # PyVCF3 leaves string lists un-None-mapped on this branch.
                    row[fname] = vals
            result[sample_name] = row
        return result
