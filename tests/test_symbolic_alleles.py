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
Regression tests for https://github.com/openvax/varcode/issues/88

VCFs (e.g. from 1000 Genomes) contain symbolic alleles for structural
variants such as `<CN0>` or `<INS:ME:ALU>`, and breakend notation like
`G]17:198982]`. Varcode previously crashed on these; it should instead
skip them (optionally with a warning).
"""

import os
import tempfile

import pytest

from varcode import load_vcf


VCF_HEADER = """##fileformat=VCFv4.1
##reference=GRCh37
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def _write_vcf(body):
    fd, path = tempfile.mkstemp(suffix=".vcf")
    with os.fdopen(fd, "w") as f:
        f.write(VCF_HEADER)
        f.write(body)
    return path


def test_symbolic_cnv_allele_is_skipped():
    path = _write_vcf(
        "22\t51179178\tBI_GS_DEL1\tA\t<CN0>\t100\tPASS\t.\n"
        "22\t51179100\tnormal_snv\tA\tG\t100\tPASS\t.\n"
    )
    try:
        vc = load_vcf(path, genome="GRCh37")
    finally:
        os.unlink(path)
    # The normal SNV should load, the symbolic allele should be skipped.
    assert len(vc) == 1, \
        "Expected 1 variant after skipping symbolic allele, got %d" % len(vc)
    assert vc[0].ref == "A" and vc[0].alt == "G"


def test_symbolic_insertion_allele_is_skipped():
    path = _write_vcf(
        "22\t51068654\tALU_test\tG\t<INS:ME:ALU>\t100\tPASS\t.\n"
    )
    try:
        vc = load_vcf(path, genome="GRCh37")
    finally:
        os.unlink(path)
    assert len(vc) == 0, \
        "Expected 0 variants (symbolic allele skipped), got %d" % len(vc)


def test_breakend_allele_is_skipped():
    # VCF 4.1 breakend notation, e.g. G]17:198982]
    path = _write_vcf(
        "17\t198982\tbnd_W\tG\tG]17:321681]\t100\tPASS\t.\n"
        "17\t198983\tnormal\tC\tT\t100\tPASS\t.\n"
    )
    try:
        vc = load_vcf(path, genome="GRCh37")
    finally:
        os.unlink(path)
    assert len(vc) == 1, \
        "Expected 1 variant after skipping breakend, got %d" % len(vc)
    assert vc[0].ref == "C" and vc[0].alt == "T"


def test_multiallelic_with_one_symbolic_and_one_normal_keeps_normal():
    # When a row has multiple alts, only the symbolic ones should be skipped.
    path = _write_vcf(
        "22\t51179178\tmixed\tA\tG,<CN0>\t100\tPASS\t.\n"
    )
    try:
        vc = load_vcf(path, genome="GRCh37")
    finally:
        os.unlink(path)
    assert len(vc) == 1, \
        "Expected 1 variant (normal SNV kept), got %d" % len(vc)
    assert vc[0].ref == "A" and vc[0].alt == "G"


def test_mixed_vcf_with_many_symbolic_alleles_does_not_crash():
    # Previously, a VCF with symbolic alleles caused the loader to crash
    # with AttributeError: '_SV' object has no attribute 'sequence'
    # (or ValueError on invalid nucleotides). The loader should be robust.
    path = _write_vcf(
        "22\t51179178\ta\tA\t<CN0>\t100\tPASS\t.\n"
        "22\t51068654\tb\tG\t<INS:ME:ALU>\t100\tPASS\t.\n"
        "22\t51068700\tc\tG\t<DEL>\t100\tPASS\t.\n"
        "22\t51179100\tnormal\tA\tG\t100\tPASS\t.\n"
    )
    try:
        vc = load_vcf(path, genome="GRCh37")
    finally:
        os.unlink(path)
    assert len(vc) == 1, "Expected 1 variant after filtering symbolic alleles"


def test_user_is_warned_when_symbolic_alleles_are_skipped():
    # The user should be informed via warnings.warn (not just DEBUG log)
    # when alleles are skipped, so silent data loss doesn't go unnoticed.
    import warnings
    path = _write_vcf(
        "22\t51179178\ta\tA\t<CN0>\t100\tPASS\t.\n"
        "22\t51068654\tb\tG\t<INS:ME:ALU>\t100\tPASS\t.\n"
        "22\t51179100\tnormal\tA\tG\t100\tPASS\t.\n"
    )
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            load_vcf(path, genome="GRCh37")
    finally:
        os.unlink(path)
    matching = [
        w for w in caught
        if "symbolic" in str(w.message).lower()
        or "breakend" in str(w.message).lower()
    ]
    assert len(matching) >= 1, \
        "Expected at least one warning about skipped symbolic alleles"
    assert "2" in str(matching[0].message), \
        "Warning should report the count (2) of skipped alleles: %r" % (
            str(matching[0].message),)


def test_no_warning_when_no_symbolic_alleles_skipped():
    # If the VCF has only normal variants, no skip warning should fire.
    import warnings
    path = _write_vcf(
        "22\t51179100\tnormal\tA\tG\t100\tPASS\t.\n"
    )
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            load_vcf(path, genome="GRCh37")
    finally:
        os.unlink(path)
    symbolic_warnings = [
        w for w in caught
        if "symbolic" in str(w.message).lower()
        or "breakend" in str(w.message).lower()
    ]
    assert len(symbolic_warnings) == 0, \
        "Should not emit symbolic-allele warnings when nothing was skipped"
