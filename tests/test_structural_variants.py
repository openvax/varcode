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

"""Tests for :class:`StructuralVariant` and the VCF symbolic-allele
parser (PR 8; tracked in #252 / #264).

Coverage:

* ``StructuralVariant`` construction — named fields, defaults, bad-type
  rejection, ``isinstance(v, Variant)`` still holds.
* ``parse_symbolic_alt`` — canonical symbolic forms (``<DEL>``,
  ``<INS:ME:ALU>``, ``<CN0>``), breakend grammar, unknown tokens fall
  back gracefully, ``*`` returns None.
* VCF loader integration — ``parse_structural_variants=True`` yields
  ``StructuralVariant`` rows; default behaviour (off) still filters
  them for back-compat.
"""

import os
import tempfile

import pytest

from varcode import (
    StructuralVariant,
    SV_TYPES,
    Variant,
    load_vcf,
    parse_symbolic_alt,
)


# --------------------------------------------------------------------
# StructuralVariant construction
# --------------------------------------------------------------------


def test_sv_deletion_basic_fields():
    sv = StructuralVariant(
        contig="7",
        start=117_480_000,
        end=117_580_000,
        sv_type="DEL",
        alt="<DEL>",
        genome="GRCh38",
    )
    assert sv.sv_type == "DEL"
    assert sv.start == 117_480_000
    assert sv.end == 117_580_000
    assert sv.length == 100_001
    assert sv.symbolic_alt == "<DEL>"
    assert sv.is_structural is True


def test_sv_is_a_variant_subclass():
    """Downstream code that checks ``isinstance(x, Variant)`` sees
    SVs too. Critical for generic variant handling pipelines."""
    sv = StructuralVariant(
        contig="7", start=1000, sv_type="INV", end=2000, genome="GRCh38")
    assert isinstance(sv, Variant)


def test_sv_unknown_type_rejected():
    with pytest.raises(ValueError):
        StructuralVariant(
            contig="7", start=1, sv_type="NOT_A_TYPE", genome="GRCh38")


def test_sv_breakend_fields():
    sv = StructuralVariant(
        contig="17",
        start=198_982,
        sv_type="BND",
        alt="G]17:198982]",
        mate_contig="17",
        mate_start=198_982,
        mate_orientation="]]",
        genome="GRCh38",
    )
    assert sv.mate_contig == "17"
    assert sv.mate_start == 198_982
    assert sv.mate_orientation == "]]"
    assert sv.length is None  # undefined for breakends


def test_sv_insertion_length_is_zero_on_reference():
    sv = StructuralVariant(
        contig="1", start=100, sv_type="INS",
        alt="<INS>", genome="GRCh38")
    # INS is zero-width on the reference; the inserted sequence's
    # length lives in alt_assembly / info, not in ref coords.
    assert sv.length == 0


def test_sv_carries_confidence_intervals():
    sv = StructuralVariant(
        contig="1", start=100, end=200, sv_type="DEL",
        ci_start=(-10, 10), ci_end=(-5, 5),
        genome="GRCh38")
    assert sv.ci_start == (-10, 10)
    assert sv.ci_end == (-5, 5)


def test_sv_alt_assembly_hook_preserved():
    """Long-read / targeted-assembly pipelines can attach the
    resolved rearranged allele. varcode stashes it verbatim; the
    SV annotator will prefer it when present."""
    sv = StructuralVariant(
        contig="1", start=100, end=200, sv_type="DEL",
        alt_assembly="ACGTACGTACGT", genome="GRCh38")
    assert sv.alt_assembly == "ACGTACGTACGT"


def test_sv_info_bag_round_trips():
    sv = StructuralVariant(
        contig="1", start=100, end=200, sv_type="DEL",
        info={"SVMETHOD": "manta", "HOMLEN": 3},
        genome="GRCh38")
    assert sv.info["SVMETHOD"] == "manta"
    assert sv.info["HOMLEN"] == 3


def test_sv_types_coverage():
    """Every SV_TYPES member should construct cleanly."""
    for t in SV_TYPES:
        StructuralVariant(
            contig="1", start=100, end=200, sv_type=t, genome="GRCh38")


# --------------------------------------------------------------------
# Symbolic-allele parser
# --------------------------------------------------------------------


def test_parse_symbolic_del():
    sv = parse_symbolic_alt(
        contig="1", start=100, ref="N", alt="<DEL>",
        info={"END": 500}, genome="GRCh38")
    assert sv is not None
    assert sv.sv_type == "DEL"
    assert sv.end == 500


def test_parse_symbolic_ins_me_alu():
    """Subtype notation ``<INS:ME:ALU>`` collapses to the top-level
    ``INS`` type with the subtype stashed in ``info``."""
    sv = parse_symbolic_alt(
        contig="1", start=100, ref="N", alt="<INS:ME:ALU>",
        info={}, genome="GRCh38")
    assert sv is not None
    assert sv.sv_type == "INS"
    assert sv.info["symbolic_subtype"] == "ME:ALU"


def test_parse_symbolic_cn0_maps_to_cnv():
    sv = parse_symbolic_alt(
        contig="1", start=100, ref="N", alt="<CN0>",
        info={"END": 500}, genome="GRCh38")
    assert sv is not None
    assert sv.sv_type == "CNV"
    # The integer is preserved so consumers can distinguish a
    # zero-copy deletion from a higher-copy duplication.
    assert sv.info["copy_number"] == 0


def test_parse_symbolic_cn3_distinct_from_cn0():
    """``<CN3>`` and ``<CN0>`` both map to ``sv_type="CNV"`` but
    must remain distinguishable — one's a deletion, one's a
    duplication. The integer rides along in ``info["copy_number"]``."""
    cn0 = parse_symbolic_alt(
        contig="1", start=100, ref="N", alt="<CN0>",
        info={"END": 500}, genome="GRCh38")
    cn3 = parse_symbolic_alt(
        contig="1", start=100, ref="N", alt="<CN3>",
        info={"END": 500}, genome="GRCh38")
    assert cn0.info["copy_number"] == 0
    assert cn3.info["copy_number"] == 3


def test_parse_symbolic_cnv_no_number_omits_copy_number():
    """``<CNV>`` (no integer) means "copy number altered, count
    unspecified" — the parser shouldn't fabricate a number."""
    sv = parse_symbolic_alt(
        contig="1", start=100, ref="N", alt="<CNV>",
        info={"END": 500}, genome="GRCh38")
    assert sv.sv_type == "CNV"
    assert "copy_number" not in sv.info


def test_parse_symbolic_ins_extracts_insseq_into_alt_assembly():
    """Long-read callers (Manta) ship the inserted sequence in
    INFO/INSSEQ; the parser routes it into ``alt_assembly`` so the
    SV annotator's assembly path picks it up."""
    sv = parse_symbolic_alt(
        contig="1", start=100, ref="A", alt="<INS>",
        info={"INSSEQ": "ACGTACGT"}, genome="GRCh38")
    assert sv is not None
    assert sv.sv_type == "INS"
    assert sv.alt_assembly == "ACGTACGT"


def test_parse_symbolic_ins_extracts_svinsseq_alias():
    """Some Sniffles/PBSV variants use SVINSSEQ instead of INSSEQ.
    Both are accepted."""
    sv = parse_symbolic_alt(
        contig="1", start=100, ref="A", alt="<INS>",
        info={"SVINSSEQ": "TTTTGGGG"}, genome="GRCh38")
    assert sv.alt_assembly == "TTTTGGGG"


def test_parse_symbolic_insseq_unwraps_single_element_list():
    """Undeclared INFO Strings come back from the header parser as
    1-element lists (no ``Number=1`` to unwrap on). The SV parser
    unwraps those rather than passing the list through to
    ``alt_assembly``."""
    sv = parse_symbolic_alt(
        contig="1", start=100, ref="A", alt="<INS>",
        info={"INSSEQ": ["ACGT"]}, genome="GRCh38")
    assert sv.alt_assembly == "ACGT"


def test_parse_symbolic_non_ins_ignores_insseq():
    """``INSSEQ`` only applies to INS — a DEL row that happens to
    carry the field shouldn't have it land in ``alt_assembly``."""
    sv = parse_symbolic_alt(
        contig="1", start=100, ref="N", alt="<DEL>",
        info={"INSSEQ": "ACGT", "END": 500}, genome="GRCh38")
    assert sv.alt_assembly is None


def test_parse_symbolic_svtype_info_overrides_alt():
    """When INFO/SVTYPE is present, it takes precedence over the
    ALT token. Useful for VCFs that use a non-canonical ALT but
    a canonical SVTYPE."""
    sv = parse_symbolic_alt(
        contig="1", start=100, ref="N", alt="<MYTOOL>",
        info={"SVTYPE": "DEL", "END": 200}, genome="GRCh38")
    assert sv is not None
    assert sv.sv_type == "DEL"


def test_parse_breakend_orientation_suffix():
    """``G]17:198982]`` — prefix-anchored breakend (joined-after)."""
    sv = parse_symbolic_alt(
        contig="2", start=321_681, ref="G",
        alt="G]17:198982]", info={}, genome="GRCh38")
    assert sv is not None
    assert sv.sv_type == "BND"
    assert sv.mate_contig == "17"
    assert sv.mate_start == 198_982
    assert sv.mate_orientation == "]]"


def test_parse_breakend_orientation_prefix():
    """``]13:123456]T`` — suffix-anchored breakend (joined-before)."""
    sv = parse_symbolic_alt(
        contig="2", start=100, ref="T",
        alt="]13:123456]T", info={}, genome="GRCh38")
    assert sv is not None
    assert sv.mate_contig == "13"
    assert sv.mate_orientation == "]]"


def test_parse_breakend_plus_strand():
    """``[13:123456[T`` — opposite strand orientation."""
    sv = parse_symbolic_alt(
        contig="2", start=100, ref="T",
        alt="[13:123456[T", info={}, genome="GRCh38")
    assert sv is not None
    assert sv.mate_orientation == "[["


def test_parse_spanning_deletion_star_returns_none():
    """``*`` is a cross-row reference, not an SV to annotate."""
    sv = parse_symbolic_alt(
        contig="1", start=100, ref="A", alt="*",
        info={}, genome="GRCh38")
    assert sv is None


def test_parse_non_symbolic_returns_none():
    """Plain nucleotide ALTs aren't SVs — parser declines to handle."""
    sv = parse_symbolic_alt(
        contig="1", start=100, ref="A", alt="G",
        info={}, genome="GRCh38")
    assert sv is None


def test_parse_unknown_symbolic_falls_back_to_bnd():
    """Future VCF extensions (``<MYTOOL>``, etc.) don't crash the
    loader — they come back as ``BND`` with the original token
    preserved in ``info["symbolic_subtype"]``."""
    sv = parse_symbolic_alt(
        contig="1", start=100, ref="N", alt="<UNKNOWN_SV>",
        info={}, genome="GRCh38")
    assert sv is not None
    assert sv.sv_type == "BND"
    assert sv.info["symbolic_subtype"] == "UNKNOWN_SV"


# --------------------------------------------------------------------
# VCF loader integration
# --------------------------------------------------------------------


def _write_vcf(body: str) -> str:
    fd, path = tempfile.mkstemp(suffix=".vcf")
    with os.fdopen(fd, "w") as f:
        f.write(body)
    return path


def test_vcf_loader_skips_svs_by_default():
    """Default behaviour (``parse_structural_variants=False``)
    preserves the 2.x contract — symbolic alleles dropped with a
    warning."""
    body = (
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38\n"
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"type\">\n"
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"end\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "22\t51179178\t.\tA\t<DEL>\t100\tPASS\tSVTYPE=DEL;END=51179500\n"
        "22\t51179700\t.\tC\tT\t100\tPASS\t.\n"
    )
    path = _write_vcf(body)
    try:
        vc = load_vcf(path, genome="GRCh38")
    finally:
        os.unlink(path)
    # Only the normal SNV survives.
    assert len(vc) == 1
    assert vc[0].alt == "T"


def test_vcf_loader_parses_svs_when_opted_in():
    """``parse_structural_variants=True`` yields
    :class:`StructuralVariant` rows alongside normal variants."""
    body = (
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38\n"
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"type\">\n"
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"end\">\n"
        "##INFO=<ID=MATEID,Number=1,Type=String,Description=\"mate\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "22\t51179178\tsv1\tA\t<DEL>\t100\tPASS\tSVTYPE=DEL;END=51179500\n"
        "22\t51179700\tsnv1\tC\tT\t100\tPASS\t.\n"
        "22\t51180000\tbnd1\tG\tG]17:198982]\t100\tPASS\tMATEID=bnd2\n"
    )
    path = _write_vcf(body)
    try:
        vc = load_vcf(path, genome="GRCh38", parse_structural_variants=True)
    finally:
        os.unlink(path)
    # All three variants survive.
    assert len(vc) == 3
    kinds = {type(v).__name__ for v in vc}
    assert "StructuralVariant" in kinds
    assert "Variant" in kinds
    sv_rows = [v for v in vc if isinstance(v, StructuralVariant)]
    assert len(sv_rows) == 2
    sv_types = {sv.sv_type for sv in sv_rows}
    assert sv_types == {"DEL", "BND"}


def test_vcf_loader_spanning_star_always_skipped():
    """``*`` is skipped regardless of the ``parse_structural_variants``
    flag — it's not an SV to parse."""
    body = (
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "22\t51179200\t.\tA\t*\t100\tPASS\t.\n"
        "22\t51179201\t.\tC\tT\t100\tPASS\t.\n"
    )
    path = _write_vcf(body)
    try:
        vc = load_vcf(path, genome="GRCh38", parse_structural_variants=True)
    finally:
        os.unlink(path)
    assert len(vc) == 1
    assert vc[0].ref == "C" and vc[0].alt == "T"


def test_vcf_loader_default_warning_mentions_flag():
    """When SVs are dropped because ``parse_structural_variants=False``,
    the warning should tell the user about the flag — not just point
    at an issue tracker."""
    import warnings
    body = (
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38\n"
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"type\">\n"
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"end\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "22\t51179178\t.\tA\t<DEL>\t100\tPASS\tSVTYPE=DEL;END=51179500\n"
        "22\t51179700\t.\tC\tT\t100\tPASS\t.\n"
    )
    path = _write_vcf(body)
    try:
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter("always")
            load_vcf(path, genome="GRCh38")
    finally:
        os.unlink(path)
    skip_warnings = [
        w for w in captured if "symbolic/breakend" in str(w.message)]
    assert len(skip_warnings) == 1
    assert "parse_structural_variants=True" in str(skip_warnings[0].message)


def test_vcf_loader_sv_carries_sample_info():
    """SV rows from somatic callers (Manta, DELLY, GRIDSS) ship
    FORMAT/sample data — the loader must parse it just like the
    simple-ALT path, not stash ``None``."""
    body = (
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38\n"
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"type\">\n"
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"end\">\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"genotype\">\n"
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n"
        "##FORMAT=<ID=PR,Number=2,Type=Integer,Description=\"paired-read support\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n"
        "22\t51179178\tsv1\tA\t<DEL>\t100\tPASS\tSVTYPE=DEL;END=51179500"
        "\tGT:DP:PR\t0/0:50:48,2\t0/1:55:30,25\n"
    )
    path = _write_vcf(body)
    try:
        vc = load_vcf(path, genome="GRCh38", parse_structural_variants=True)
    finally:
        os.unlink(path)
    assert len(vc) == 1
    sv = vc[0]
    assert isinstance(sv, StructuralVariant)
    sample_info = vc.metadata[sv]["sample_info"]
    assert sample_info is not None, (
        "SV rows must carry per-sample FORMAT data through to metadata; "
        "previously sample_info was always None on the SV branch.")
    assert list(sample_info) == ["NORMAL", "TUMOR"]
    assert sample_info["TUMOR"]["GT"] == "0/1"
    assert sample_info["TUMOR"]["DP"] == 55
    assert sample_info["TUMOR"]["PR"] == [30, 25]


def test_vcf_loader_sv_yields_structural_variant_effect():
    """End-to-end: a VCF with an SV row, loaded with the flag on,
    produces a typed :class:`StructuralVariantEffect` when annotated.

    This is the gate that protects #257 / #252 / #261: if this test
    works, downstream issues can build SV-typed pipelines on top of
    ``load_vcf`` instead of constructing :class:`StructuralVariant`
    objects by hand.
    """
    from pyensembl import cached_release
    from varcode.effects.effect_classes import StructuralVariantEffect
    ensembl_grch38 = cached_release(81)
    cftr = ensembl_grch38.transcript_by_id("ENST00000003084")
    body = (
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38\n"
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"type\">\n"
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"end\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "7\t%d\t.\tA\t<DEL>\t100\tPASS\tSVTYPE=DEL;END=%d\n"
    ) % (cftr.start + 100, cftr.start + 50_000)
    path = _write_vcf(body)
    try:
        vc = load_vcf(
            path, genome=ensembl_grch38, parse_structural_variants=True)
        effects = vc.effects()
    finally:
        os.unlink(path)
    sv_effects = [e for e in effects if isinstance(e, StructuralVariantEffect)]
    assert len(sv_effects) >= 1, (
        "Expected at least one StructuralVariantEffect from the loaded "
        "<DEL> row; got effect types %s" % [type(e).__name__ for e in effects])
