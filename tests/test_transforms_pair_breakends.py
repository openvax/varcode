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

"""Tests for :func:`varcode.transforms.pair_breakends` (#364).

Covers the test matrix from PR #367:

* Paired BND collapse across Manta / DELLY / GRIDSS-style inputs.
* Pass-through cases — non-BND, single-row TRA, single-ended BND.
* Edge cases — mate-missing-from-VC, N-way collision, asymmetric
  MATEID, disagreeing genotypes, disagreeing alt_assembly.
* Idempotence, mixed VC, empty VC, source_variants provenance.
"""

import os
import tempfile
import warnings

import pytest

from varcode import StructuralVariant, Variant, load_vcf
from varcode.transforms import pair_breakends


# --------------------------------------------------------------------
# Fixture helpers
# --------------------------------------------------------------------


def _write_vcf(body: str) -> str:
    fd, path = tempfile.mkstemp(suffix=".vcf")
    with os.fdopen(fd, "w") as f:
        f.write(body)
    return path


_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##reference=GRCh38\n"
    "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"type\">\n"
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"end\">\n"
    "##INFO=<ID=MATEID,Number=1,Type=String,Description=\"mate\">\n"
    "##INFO=<ID=PARID,Number=1,Type=String,Description=\"partner\">\n"
    "##INFO=<ID=INSSEQ,Number=1,Type=String,Description=\"inserted seq\">\n"
    "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"mate chrom\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def _load(body: str, only_passing: bool = True):
    path = _write_vcf(_HEADER + body)
    try:
        return load_vcf(
            path, genome="GRCh38",
            parse_structural_variants=True,
            only_passing=only_passing)
    finally:
        os.unlink(path)


def _bnds(vc):
    return [v for v in vc if getattr(v, "sv_type", None) == "BND"]


# --------------------------------------------------------------------
# 1. Manta paired BND
# --------------------------------------------------------------------


def test_manta_paired_bnd_collapses_to_single_combined():
    """Two symmetric BND rows merge into one combined StructuralVariant
    with source_variants pointing at both originals."""
    body = (
        "19\t15250000\tbnd1\tG\tG]15:34350000]\t100\tPASS\tMATEID=bnd2\n"
        "15\t34350000\tbnd2\tA\tA[19:15250000[\t100\tPASS\tMATEID=bnd1\n"
    )
    vc = _load(body)
    assert len(_bnds(vc)) == 2

    out = pair_breakends(vc)
    bnds = _bnds(out)
    assert len(bnds) == 1
    combined = bnds[0]
    # Deterministic primary: lex-earlier source row ID is "bnd1" -> chr19.
    assert combined.contig == "19"
    assert combined.start == 15_250_000
    assert combined.mate_contig == "15"
    assert combined.mate_start == 34_350_000
    assert len(combined.source_variants) == 2
    source_contigs = sorted(v.contig for v in combined.source_variants)
    assert source_contigs == ["15", "19"]


# --------------------------------------------------------------------
# 2. DELLY-style paired BND (different ID format)
# --------------------------------------------------------------------


def test_delly_paired_bnd_with_long_ids():
    """DELLY-style IDs are longer (``DEL00000001`` etc.) but the
    MATEID mechanism is identical."""
    body = (
        "19\t15250000\tDEL00000001_1\tG\tG]15:34350000]\t100\tPASS\t"
        "MATEID=DEL00000001_2\n"
        "15\t34350000\tDEL00000001_2\tA\tA[19:15250000[\t100\tPASS\t"
        "MATEID=DEL00000001_1\n"
    )
    vc = _load(body)
    out = pair_breakends(vc)
    bnds = _bnds(out)
    assert len(bnds) == 1
    assert len(bnds[0].source_variants) == 2


# --------------------------------------------------------------------
# 3. GRIDSS paired BND with PARID
# --------------------------------------------------------------------


def test_gridss_paired_bnd_with_parid_alias():
    """Older GRIDSS uses PARID instead of MATEID. The transform treats
    them as aliases."""
    body = (
        "19\t15250000\tgridss1o\tG\tG]15:34350000]\t100\tPASS\t"
        "PARID=gridss1h\n"
        "15\t34350000\tgridss1h\tA\tA[19:15250000[\t100\tPASS\t"
        "PARID=gridss1o\n"
    )
    vc = _load(body)
    out = pair_breakends(vc)
    bnds = _bnds(out)
    assert len(bnds) == 1
    assert len(bnds[0].source_variants) == 2


# --------------------------------------------------------------------
# 4. BreakDancer-style single-row TRA (no MATEID)
# --------------------------------------------------------------------


def test_single_row_translocation_passes_through():
    """A BreakDancer/CREST-style single-row TRA (symbolic ALT, no
    MATEID) passes through unchanged."""
    body = (
        "19\t15250000\ttra1\tN\t<TRA>\t100\tPASS\tSVTYPE=TRA;CHR2=15;END=34350000\n"
    )
    vc = _load(body)
    # The symbolic <TRA> currently falls back to sv_type="BND" via the
    # parser's unknown-token path, which is the right shape for pairing
    # to skip it (no MATEID).
    out = pair_breakends(vc)
    assert len(out) == 1
    assert out[0].source_variants == ()


# --------------------------------------------------------------------
# 5. Single-ended BND (no MATEID)
# --------------------------------------------------------------------


def test_single_ended_bnd_without_mateid_passes_through():
    """A BND without MATEID — e.g. GRIDSS's unresolved single-end
    breakend representation — passes through."""
    body = (
        "19\t15250000\tbnd1\tG\tG]15:34350000]\t100\tPASS\t.\n"
    )
    vc = _load(body)
    out = pair_breakends(vc)
    bnds = _bnds(out)
    assert len(bnds) == 1
    assert bnds[0].source_variants == ()


# --------------------------------------------------------------------
# 6. Mate-missing-from-VC
# --------------------------------------------------------------------


def test_mate_missing_from_vc_warns_and_passes_through():
    """When MATEID references an ID not in this collection, the
    half passes through with a warning."""
    body = (
        "19\t15250000\tbnd1\tG\tG]15:34350000]\t100\tPASS\tMATEID=bnd2\n"
    )
    vc = _load(body)
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        out = pair_breakends(vc)
    msgs = [str(w.message) for w in captured]
    assert any("mate" in m.lower() and "bnd2" in m for m in msgs), msgs
    bnds = _bnds(out)
    assert len(bnds) == 1
    assert bnds[0].source_variants == ()


# --------------------------------------------------------------------
# 7. N-way MATEID collision
# --------------------------------------------------------------------


def test_three_way_mateid_collision_leaves_group_unpaired():
    """Three rows referencing the same MATEID group are left
    unpaired with a warning. (Each row's MATEID points at the
    same id, creating a degenerate pairing group.)"""
    body = (
        "19\t15250000\tbnd1\tG\tG]15:34350000]\t100\tPASS\tMATEID=bnd2\n"
        "15\t34350000\tbnd2\tA\tA[19:15250000[\t100\tPASS\tMATEID=bnd1\n"
        "19\t15260000\tbnd3\tG\tG]15:34350000]\t100\tPASS\tMATEID=bnd2\n"
    )
    vc = _load(body)
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        out = pair_breakends(vc)
    msgs = [str(w.message) for w in captured]
    assert any("ambiguous" in m.lower() or "members" in m.lower()
               for m in msgs), msgs
    # All three pass through unpaired.
    bnds = _bnds(out)
    assert len(bnds) == 3
    assert all(b.source_variants == () for b in bnds)


# --------------------------------------------------------------------
# 8. Asymmetric MATEID
# --------------------------------------------------------------------


def test_asymmetric_mateid_leaves_pair_unpaired():
    """A.mateid == B.id but B.mateid != A.id is asymmetric. The pair
    is left unpaired with a warning."""
    body = (
        "19\t15250000\tbnd1\tG\tG]15:34350000]\t100\tPASS\tMATEID=bnd2\n"
        "15\t34350000\tbnd2\tA\tA[19:15250000[\t100\tPASS\tMATEID=bndX\n"
    )
    vc = _load(body)
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        out = pair_breakends(vc)
    msgs = [str(w.message) for w in captured]
    assert any("asymmetric" in m.lower() or "mate" in m.lower()
               for m in msgs), msgs
    bnds = _bnds(out)
    assert len(bnds) == 2
    assert all(b.source_variants == () for b in bnds)


# --------------------------------------------------------------------
# 9. Disagreeing genotypes raise
# --------------------------------------------------------------------


def test_disagreeing_genotypes_raise():
    """Both halves of a paired BND describe the same biological event;
    disagreeing per-sample GT indicates a caller bug. The transform
    raises rather than silently coercing."""
    header_with_sample = (
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38\n"
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"t\">\n"
        "##INFO=<ID=MATEID,Number=1,Type=String,Description=\"m\">\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n"
        "19\t15250000\tbnd1\tG\tG]15:34350000]\t100\tPASS\tMATEID=bnd2\tGT\t0/1\n"
        "15\t34350000\tbnd2\tA\tA[19:15250000[\t100\tPASS\tMATEID=bnd1\tGT\t1/1\n"
    )
    path = _write_vcf(header_with_sample)
    try:
        vc = load_vcf(path, genome="GRCh38", parse_structural_variants=True)
    finally:
        os.unlink(path)
    with pytest.raises(ValueError, match="disagree on genotype"):
        pair_breakends(vc)


# --------------------------------------------------------------------
# 10. Disagreeing alt_assembly warns, A wins
# --------------------------------------------------------------------


def test_disagreeing_alt_assembly_warns_a_wins():
    """When the two halves carry different INSSEQ assemblies, the
    transform warns and keeps A's (lex-earlier ID). Both originals
    remain reachable via source_variants."""
    # INS-style assembly carried on BND rows is non-standard but
    # representative — Manta's BND output sometimes ships partial
    # INSSEQ when local assembly partially resolved the join.
    a = StructuralVariant(
        contig="19", start=15_250_000, sv_type="BND",
        alt="G]15:34350000]", mate_contig="15",
        mate_start=34_350_000, alt_assembly="ACGTA",
        info={"mateid": "bnd2"},
        genome="GRCh38",
    )
    b = StructuralVariant(
        contig="15", start=34_350_000, sv_type="BND",
        alt="A[19:15250000[", mate_contig="19",
        mate_start=15_250_000, alt_assembly="TTTTT",
        info={"mateid": "bnd1"},
        genome="GRCh38",
    )
    from varcode import VariantCollection
    metadata = {
        "synth": {
            a: {"id": "bnd1", "info": {"MATEID": "bnd2"}},
            b: {"id": "bnd2", "info": {"MATEID": "bnd1"}},
        }
    }
    vc = VariantCollection(
        variants=[a, b],
        sources={"synth"},
        source_to_metadata_dict=metadata,
    )
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        out = pair_breakends(vc)
    msgs = [str(w.message) for w in captured]
    assert any("alt_assembly" in m for m in msgs), msgs
    bnds = _bnds(out)
    assert len(bnds) == 1
    # A's id is "bnd1" -> lex-earlier; combined inherits its assembly.
    assert bnds[0].alt_assembly == "ACGTA"


# --------------------------------------------------------------------
# 11. Mixed VC: SNVs + paired BNDs + DEL + single-row TRA
# --------------------------------------------------------------------


def test_mixed_vc_only_paired_bnds_collapse():
    body = (
        "1\t1000\t.\tA\tT\t100\tPASS\t.\n"
        "1\t2000\t.\tG\tC\t100\tPASS\t.\n"
        "22\t51179178\tdel1\tA\t<DEL>\t100\tPASS\tSVTYPE=DEL;END=51179500\n"
        "19\t15250000\tbnd1\tG\tG]15:34350000]\t100\tPASS\tMATEID=bnd2\n"
        "15\t34350000\tbnd2\tA\tA[19:15250000[\t100\tPASS\tMATEID=bnd1\n"
    )
    vc = _load(body)
    assert len(vc) == 5  # 2 SNVs + DEL + 2 BNDs
    out = pair_breakends(vc)
    assert len(out) == 4  # 2 SNVs + DEL + 1 combined BND

    bnds = _bnds(out)
    assert len(bnds) == 1
    assert len(bnds[0].source_variants) == 2

    # Non-BND variants pass through with empty source_variants.
    non_bnd = [v for v in out if getattr(v, "sv_type", None) != "BND"]
    assert all(v.source_variants == () for v in non_bnd)


# --------------------------------------------------------------------
# 12. Idempotence
# --------------------------------------------------------------------


def test_idempotent_on_already_paired_input():
    """Running pair_breakends twice produces the same VC the second
    time — the combined variant's source_variants is non-empty, so it
    skips re-pairing."""
    body = (
        "19\t15250000\tbnd1\tG\tG]15:34350000]\t100\tPASS\tMATEID=bnd2\n"
        "15\t34350000\tbnd2\tA\tA[19:15250000[\t100\tPASS\tMATEID=bnd1\n"
    )
    vc = _load(body)
    once = pair_breakends(vc)
    twice = pair_breakends(once)
    assert len(once) == len(twice) == 1
    # Same combined variant (object identity not required; semantic
    # equality is enough).
    assert once[0].contig == twice[0].contig
    assert once[0].start == twice[0].start
    assert once[0].mate_contig == twice[0].mate_contig
    assert once[0].mate_start == twice[0].mate_start
    # source_variants preserved across the second call.
    assert len(twice[0].source_variants) == 2


# --------------------------------------------------------------------
# 13. Empty VC
# --------------------------------------------------------------------


def test_empty_vc_passes_through():
    from varcode import VariantCollection
    vc = VariantCollection(variants=[], source_to_metadata_dict={})
    out = pair_breakends(vc)
    assert len(out) == 0


# --------------------------------------------------------------------
# 14. source_variants is not part of hash/equality
# --------------------------------------------------------------------


def test_source_variants_not_in_hash_or_equality():
    """``source_variants`` is provenance, not identity. Two variants
    with the same (contig, start, ref, alt, reference_name) but
    different source_variants compare equal and hash identically."""
    v_plain = Variant("1", 100, "A", "T", "GRCh38")
    v_with_source = Variant("1", 100, "A", "T", "GRCh38")
    v_with_source.source_variants = (Variant("1", 99, "A", "T", "GRCh38"),)
    assert v_plain == v_with_source
    assert hash(v_plain) == hash(v_with_source)


# --------------------------------------------------------------------
# 15. Metadata round-trip — combined variant has its own metadata
# --------------------------------------------------------------------


def test_combined_variant_has_merged_metadata_entry():
    """After pairing, the combined variant has a single metadata
    entry per source path with the merged qual/filter/info."""
    body = (
        "19\t15250000\tbnd1\tG\tG]15:34350000]\t50\tPASS\tMATEID=bnd2\n"
        "15\t34350000\tbnd2\tA\tA[19:15250000[\t30\tLowQual\tMATEID=bnd1\n"
    )
    vc = _load(body, only_passing=False)
    out = pair_breakends(vc)
    combined = _bnds(out)[0]
    # Single source path, one metadata entry for the combined variant.
    source_path = next(iter(out.source_to_metadata_dict))
    meta = out.source_to_metadata_dict[source_path][combined]
    # qual: min(50, 30) = 30.
    assert meta["qual"] == 30
    # filter: union; LowQual wins over PASS.
    assert "LowQual" in meta["filter"]
    assert "PASS" not in meta["filter"]
    # info: paired_with populated to B's ID.
    assert meta["info"]["paired_with"] == "bnd2"


# --------------------------------------------------------------------
# 16. End-to-end: combined fusion produces single GeneFusion class
# --------------------------------------------------------------------


def test_paired_bnd_in_genes_produces_single_event_after_pairing():
    """Reciprocal BND between BRD4 (chr19) and NUTM1 (chr15) gives N
    GeneFusion effects per BRD4 transcript pre-pair (one fusion
    direction per half × transcripts), and only the primary-side
    direction × transcripts after pair_breakends — the second
    direction is reachable via combined.source_variants."""
    body = (
        "19\t15250000\tbnd1\tG\tG]15:34350000]\t100\tPASS\tMATEID=bnd2\n"
        "15\t34350000\tbnd2\tA\tA[19:15250000[\t100\tPASS\tMATEID=bnd1\n"
    )
    vc = _load(body)
    out = pair_breakends(vc)
    effects_post = out.effects(raise_on_error=False)
    fusion_post = [e for e in effects_post
                   if type(e).__name__ == "GeneFusion"]
    # All post-pair fusion effects share the same combined variant.
    assert fusion_post, "expected at least one GeneFusion post-pair"
    variants_seen = {id(e.variant) for e in fusion_post}
    assert len(variants_seen) == 1
    # That variant has both halves in source_variants.
    assert len(fusion_post[0].variant.source_variants) == 2
