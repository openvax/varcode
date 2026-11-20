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

"""Tests for :func:`varcode.transforms.combine_cis_snvs` (#368).

Exercises all scientific and provenance invariants on real GRCh38 (Ensembl 81)
transcripts:

* Explicit phase/sample scope requirement and validation.
* Both strands (CFTR ``+`` strand with 5' UTR = 132; TP53 ``-`` strand with
  5' UTR = 190, not divisible by 3, catching raw ``spliced_offset // 3`` bugs).
* Codon positions 1+2, 1+3 (reference-validated gap fill), and 1+2+3.
* Codon boundaries (adjacent bases in consecutive codons do not combine).
* Split codons across exon-exon junctions (never absorb intervening introns).
* Differing overlapping transcript CDS frames (CDKN2A p16INK4a vs p14ARF:
  declined without lex-first tie-breaking).
* Trans, unknown, contradictory, and mixed hom+het (``1/1`` + ``1|0``) phasing.
* Multi-sample phase conflicts vs compatible somatic/germline multi-sample VCFs.
* Same-position multi-allelic alternatives preserved unchanged.
* Reference-validated gap base refusal on ref mismatch or unavailable bases.
* Full metadata provenance, idempotence, and consequence equivalence against
  the haplotype annotation path.
"""

from __future__ import annotations

import warnings

import pytest
from pyensembl import cached_release

from varcode import (
    MolecularPhaseResolver,
    SampleNotFoundError,
    VCFPhaseResolver,
    Variant,
    VariantCollection,
    load_vcf,
)
from varcode.transforms import combine_cis_snvs, left_align_indels, pair_breakends


ENSEMBL81 = cached_release(81)


def _write_phased_vcf(tmp_path, name, samples, rows):
    """Write a minimal GRCh38 VCF with GT and PS FORMAT fields."""
    header = [
        "##fileformat=VCFv4.2",
        "##reference=GRCh38",
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]
    path = tmp_path / name
    path.write_text("\n".join(header + rows) + "\n")
    return str(path)


class _ExplicitCisSource:
    """Minimal ReadPhasingSource/in_cis test stub for MolecularPhaseResolver."""

    def __init__(self, cis_pairs=(), trans_pairs=()):
        self._cis = {frozenset(p) for p in cis_pairs}
        self._trans = {frozenset(p) for p in trans_pairs}
        self._covered = {v for p in cis_pairs for v in p}

    def has_evidence(self, variant) -> bool:
        return variant in self._covered

    def partners_in_cis(self, variant):
        out = []
        for pair in self._cis:
            if variant in pair:
                out.extend(v for v in pair if v != variant)
        return tuple(out)

    def in_cis(self, v1, v2, transcript=None):
        key = frozenset((v1, v2))
        if key in self._cis:
            return True
        if key in self._trans:
            return False
        return None


def _assert_consequence_equivalent_to_haplotype(vc, combined_vc, resolver, transcript_id):
    """Verify combined MNV produces the exact same protein consequence as
    haplotype annotation on ``transcript_id``."""
    hap_effects = vc.effects(phase_resolver=resolver)
    joint_effects = [
        e for e in hap_effects
        if e.transcript_id == transcript_id and len(getattr(e, "variants", ())) >= 2
    ]
    assert len(joint_effects) == 1, (
        "Expected one joint haplotype effect on %s, got %r" % (transcript_id, joint_effects)
    )
    mnv_effects = [
        e for e in combined_vc.effects()
        if e.transcript_id == transcript_id
    ]
    assert len(mnv_effects) == 1, (
        "Expected one combined MNV effect on %s, got %r" % (transcript_id, mnv_effects)
    )
    joint = joint_effects[0]
    mnv = mnv_effects[0]
    assert type(mnv) is type(joint)
    assert mnv.short_description == joint.short_description
    assert mnv.mutant_protein_sequence == joint.mutant_protein_sequence
    assert mnv.original_protein_sequence == joint.original_protein_sequence


# --------------------------------------------------------------------
# Explicit scope validation
# --------------------------------------------------------------------


def test_requires_explicit_phase_resolver_or_sample():
    vc = VariantCollection([Variant("7", 117559590, "A", "G", genome=ENSEMBL81)])
    with pytest.raises(ValueError, match="phase_resolver or sample"):
        combine_cis_snvs(vc)


def test_rejects_both_phase_resolver_and_sample():
    vc = VariantCollection([Variant("7", 117559590, "A", "G", genome=ENSEMBL81)])
    resolver = VCFPhaseResolver(vc, "tumor")
    with pytest.raises(ValueError, match="not both"):
        combine_cis_snvs(vc, phase_resolver=resolver, sample="tumor")


def test_unknown_sample_raises_sample_not_found(tmp_path):
    path = _write_phased_vcf(
        tmp_path,
        "single.vcf",
        ["tumor"],
        [
            "7\t117559590\tsnv1\tA\tG\t50\tPASS\tDP=40\tGT:PS\t1|0:100",
            "7\t117559591\tsnv2\tT\tC\t60\tPASS\tDP=42\tGT:PS\t1|0:100",
        ],
    )
    vc = load_vcf(path, genome=ENSEMBL81)
    with pytest.raises(SampleNotFoundError):
        combine_cis_snvs(vc, sample="missing_sample")


# --------------------------------------------------------------------
# Plus strand (CFTR, 5' UTR = 132) — 1+2, 1+3 gap fill, 1+2+3, hom-alt
# --------------------------------------------------------------------


def test_plus_strand_adjacent_cis_snvs_combine_and_match_haplotype(tmp_path):
    # CFTR ENST00000003084 (+) codon 506: 117559590..117559592 = ATC (Ile)
    # 117559590 A>G + 117559591 T>C -> GCC (Ala)
    path = _write_phased_vcf(
        tmp_path,
        "cftr_12.vcf",
        ["tumor"],
        [
            "7\t117559590\tsnv1\tA\tG\t50\tPASS\tDP=40\tGT:PS\t1|0:1000",
            "7\t117559591\tsnv2\tT\tC\t35\tPASS\tDP=45\tGT:PS\t1|0:1000",
        ],
    )
    vc = load_vcf(path, genome=ENSEMBL81)
    out = combine_cis_snvs(vc, sample="tumor")

    assert len(out) == 1
    mnv = out[0]
    assert (mnv.contig, mnv.start, mnv.ref, mnv.alt) == ("7", 117559590, "AT", "GC")
    assert mnv.source_variants == (vc[0], vc[1])

    meta = out.metadata[mnv]
    assert meta["id"] == "snv1"
    assert meta["qual"] == 35.0
    assert meta["filter"] == []
    assert meta["info"]["combined_from_cis_snvs"] is True
    assert meta["info"]["original_starts"] == (117559590, 117559591)
    assert meta["info"]["source_ids"] == ("snv1", "snv2")
    assert len(meta["info"]["source_metadata"]) == 2
    assert meta["sample_info"]["tumor"]["GT"] == "1|0"
    assert meta["sample_info"]["tumor"]["PS"] == 1000

    # Equivalence against haplotype annotation path
    _assert_consequence_equivalent_to_haplotype(
        vc, out, VCFPhaseResolver(vc, "tumor"), "ENST00000003084")

    # Idempotence
    out2 = combine_cis_snvs(out, sample="tumor")
    assert list(out2) == list(out)
    assert out2[0].source_variants == mnv.source_variants


def test_plus_strand_codon_positions_1_and_3_fills_reference_gap(tmp_path):
    # CFTR codon 506: 117559590..117559592 = ATC
    # Pos 1 (A>G) + Pos 3 (C>A) -> middle base T filled from reference -> GTA (Val)
    path = _write_phased_vcf(
        tmp_path,
        "cftr_13.vcf",
        ["tumor"],
        [
            "7\t117559590\tsnv1\tA\tG\t50\tPASS\tDP=30\tGT:PS\t1|0:500",
            "7\t117559592\tsnv3\tC\tA\t40\tPASS\tDP=32\tGT:PS\t1|0:500",
        ],
    )
    vc = load_vcf(path, genome=ENSEMBL81)
    out = combine_cis_snvs(vc, sample="tumor")

    assert len(out) == 1
    mnv = out[0]
    assert (mnv.contig, mnv.start, mnv.ref, mnv.alt) == ("7", 117559590, "ATC", "GTA")
    assert mnv.source_variants == (vc[0], vc[1])
    _assert_consequence_equivalent_to_haplotype(
        vc, out, VCFPhaseResolver(vc, "tumor"), "ENST00000003084")


def test_three_cis_snvs_same_codon_combine_into_one_mnv(tmp_path):
    # CFTR codon 506: ATC -> TGG (Trp) via all 3 homozygous-alt SNVs
    path = _write_phased_vcf(
        tmp_path,
        "cftr_123.vcf",
        ["tumor"],
        [
            "7\t117559590\ts1\tA\tT\t50\tPASS\t.\tGT:PS\t1/1:.",
            "7\t117559591\ts2\tT\tG\t50\tPASS\t.\tGT:PS\t1/1:.",
            "7\t117559592\ts3\tC\tG\t50\tPASS\t.\tGT:PS\t1/1:.",
        ],
    )
    vc = load_vcf(path, genome=ENSEMBL81)
    out = combine_cis_snvs(vc, sample="tumor")

    assert len(out) == 1
    mnv = out[0]
    assert (mnv.contig, mnv.start, mnv.ref, mnv.alt) == ("7", 117559590, "ATC", "TGG")
    assert mnv.source_variants == (vc[0], vc[1], vc[2])
    assert out.metadata[mnv]["sample_info"]["tumor"]["GT"] == "1/1"
    _assert_consequence_equivalent_to_haplotype(
        vc, out, VCFPhaseResolver(vc, "tumor"), "ENST00000003084")


# --------------------------------------------------------------------
# Minus strand with nonzero 5' UTR not divisible by 3 (TP53: 5' UTR = 190)
# --------------------------------------------------------------------


def test_minus_strand_nonzero_utr_cds_relative_codon_and_boundary(tmp_path):
    # TP53 ENST00000269305 (- strand, 5' UTR = 190; 190 % 3 == 1!):
    # Codon 174 (R175): genomic 7675087=G (pos 2), 7675088=C (pos 1), 7675089=G (pos 0).
    # Raw spliced_offset // 3 would wrongly separate 7675087 from 7675088/7675089.
    path = _write_phased_vcf(
        tmp_path,
        "tp53_r175.vcf",
        ["tumor"],
        [
            "17\t7675087\tp2\tG\tA\t60\tPASS\t.\tGT:PS\t0|1:900",
            "17\t7675088\tp1\tC\tT\t55\tPASS\t.\tGT:PS\t0|1:900",
        ],
    )
    vc = load_vcf(path, genome=ENSEMBL81)
    out = combine_cis_snvs(vc, sample="tumor")

    assert len(out) == 1
    mnv = out[0]
    assert (mnv.contig, mnv.start, mnv.ref, mnv.alt) == ("17", 7675087, "GC", "AT")
    assert mnv.source_variants == (vc[0], vc[1])
    assert out.metadata[mnv]["sample_info"]["tumor"]["GT"] == "0|1"
    _assert_consequence_equivalent_to_haplotype(
        vc, out, VCFPhaseResolver(vc, "tumor"), "ENST00000269305")


def test_codon_boundary_cis_snvs_do_not_combine_both_strands(tmp_path):
    # Plus strand CFTR: 117559592 (codon 506 pos 2) + 117559593 (codon 507 pos 0)
    plus_path = _write_phased_vcf(
        tmp_path,
        "cftr_boundary.vcf",
        ["tumor"],
        [
            "7\t117559592\tc506\tC\tA\t50\tPASS\t.\tGT:PS\t1|0:100",
            "7\t117559593\tc507\tT\tC\t50\tPASS\t.\tGT:PS\t1|0:100",
        ],
    )
    vc_plus = load_vcf(plus_path, genome=ENSEMBL81)
    assert list(combine_cis_snvs(vc_plus, sample="tumor")) == list(vc_plus)

    # Minus strand TP53: 7675086 (codon 175 pos 0) + 7675087 (codon 174 pos 2)
    # Note: raw spliced_offset // 3 would have WRONGLY merged these two!
    minus_path = _write_phased_vcf(
        tmp_path,
        "tp53_boundary.vcf",
        ["tumor"],
        [
            "17\t7675086\tc175\tA\tG\t50\tPASS\t.\tGT:PS\t1|0:200",
            "17\t7675087\tc174\tG\tA\t50\tPASS\t.\tGT:PS\t1|0:200",
        ],
    )
    vc_minus = load_vcf(minus_path, genome=ENSEMBL81)
    assert list(combine_cis_snvs(vc_minus, sample="tumor")) == list(vc_minus)


# --------------------------------------------------------------------
# Split codons across exon-exon junctions never absorb introns
# --------------------------------------------------------------------


def test_split_codon_across_exon_junction_is_never_merged_both_strands(tmp_path):
    # CFTR (+) exon 1 -> exon 2 split codon 17 (24 kb intron between 117480147 and 117504253)
    cftr_split = _write_phased_vcf(
        tmp_path,
        "cftr_split.vcf",
        ["tumor"],
        [
            "7\t117480147\tex1\tG\tA\t50\tPASS\t.\tGT:PS\t1|0:10",
            "7\t117504253\tex2\tC\tT\t50\tPASS\t.\tGT:PS\t1|0:10",
        ],
    )
    vc_cftr = load_vcf(cftr_split, genome=ENSEMBL81)
    out_cftr = combine_cis_snvs(vc_cftr, sample="tumor")
    assert list(out_cftr) == list(vc_cftr)
    assert all(v.source_variants == () for v in out_cftr)

    # TP53 (-) exon 2 -> exon 3 split codon 24 (118 bp intron between 7676521 and 7676403)
    tp53_split = _write_phased_vcf(
        tmp_path,
        "tp53_split.vcf",
        ["tumor"],
        [
            "17\t7676403\tex3\tT\tC\t50\tPASS\t.\tGT:PS\t1|0:20",
            "17\t7676521\tex2\tA\tG\t50\tPASS\t.\tGT:PS\t1|0:20",
        ],
    )
    vc_tp53 = load_vcf(tp53_split, genome=ENSEMBL81)
    out_tp53 = combine_cis_snvs(vc_tp53, sample="tumor")
    assert list(out_tp53) == list(vc_tp53)


# --------------------------------------------------------------------
# Differing overlapping transcript CDS frames (CDKN2A p16INK4a vs p14ARF)
# --------------------------------------------------------------------


def test_differing_transcript_cds_frames_declines_reduction_without_lex_first(tmp_path):
    # CDKN2A exon 2 (chr9, - strand):
    # 21971119 (T): p16 codon 79 pos 2, ARF codon 94 pos 0
    # 21971120 (C): p16 codon 79 pos 1, ARF codon 93 pos 2
    # 21971121 (G): p16 codon 79 pos 0, ARF codon 93 pos 1
    # Two SNVs (1119 + 1120): same codon in p16, DIFFERENT codons in ARF -> declined!
    path_two = _write_phased_vcf(
        tmp_path,
        "cdkn2a_two.vcf",
        ["tumor"],
        [
            "9\t21971119\ts1\tT\tA\t50\tPASS\t.\tGT:PS\t1|0:300",
            "9\t21971120\ts2\tC\tA\t50\tPASS\t.\tGT:PS\t1|0:300",
        ],
    )
    vc_two = load_vcf(path_two, genome=ENSEMBL81)
    with pytest.warns(UserWarning, match="incompatible CDS/codon placements"):
        out_two = combine_cis_snvs(vc_two, sample="tumor")
    assert list(out_two) == list(vc_two)

    # Three SNVs (1119 + 1120 + 1121): p16 groups {1119,1120,1121}, ARF groups {1120,1121}
    # -> overlapping incompatible proposals; declined without picking lex-first transcript!
    path_three = _write_phased_vcf(
        tmp_path,
        "cdkn2a_three.vcf",
        ["tumor"],
        [
            "9\t21971119\ts1\tT\tA\t50\tPASS\t.\tGT:PS\t1|0:300",
            "9\t21971120\ts2\tC\tA\t50\tPASS\t.\tGT:PS\t1|0:300",
            "9\t21971121\ts3\tG\tT\t50\tPASS\t.\tGT:PS\t1|0:300",
        ],
    )
    vc_three = load_vcf(path_three, genome=ENSEMBL81)
    with pytest.warns(UserWarning, match="incompatible transcript interpretations"):
        out_three = combine_cis_snvs(vc_three, sample="tumor")
    assert list(out_three) == list(vc_three)


# --------------------------------------------------------------------
# Trans, unknown, contradictory, mixed hom+het, multi-sample conflicts
# --------------------------------------------------------------------


def test_trans_and_unknown_phase_and_same_ps_trans_remain_unchanged(tmp_path):
    # Same PS=100, but 1|0 vs 0|1 (TRANS!) — raw PS equality must NOT merge them!
    trans_path = _write_phased_vcf(
        tmp_path,
        "trans_same_ps.vcf",
        ["tumor"],
        [
            "7\t117559590\ts1\tA\tG\t50\tPASS\t.\tGT:PS\t1|0:100",
            "7\t117559591\ts2\tT\tC\t50\tPASS\t.\tGT:PS\t0|1:100",
        ],
    )
    vc_trans = load_vcf(trans_path, genome=ENSEMBL81)
    assert list(combine_cis_snvs(vc_trans, sample="tumor")) == list(vc_trans)

    # Unphased 0/1 with same PS tag -> unknown phase -> unchanged
    unphased_path = _write_phased_vcf(
        tmp_path,
        "unphased.vcf",
        ["tumor"],
        [
            "7\t117559590\ts1\tA\tG\t50\tPASS\t.\tGT:PS\t0/1:100",
            "7\t117559591\ts2\tT\tC\t50\tPASS\t.\tGT:PS\t0/1:100",
        ],
    )
    vc_unphased = load_vcf(unphased_path, genome=ENSEMBL81)
    assert list(combine_cis_snvs(vc_unphased, sample="tumor")) == list(vc_unphased)


def test_mixed_hom_and_het_declines_to_avoid_losing_second_haplotype(tmp_path):
    # v1 is 1/1 (on both haplotypes), v2 is 1|0 (only hap 0) -> hap 1 still carries v1 alone!
    path = _write_phased_vcf(
        tmp_path,
        "hom_plus_het.vcf",
        ["tumor"],
        [
            "7\t117559590\ts1\tA\tG\t50\tPASS\t.\tGT:PS\t1/1:100",
            "7\t117559591\ts2\tT\tC\t50\tPASS\t.\tGT:PS\t1|0:100",
        ],
    )
    vc = load_vcf(path, genome=ENSEMBL81)
    with pytest.warns(UserWarning, match="incompatible genotypes"):
        out = combine_cis_snvs(vc, sample="tumor")
    assert list(out) == list(vc)


def test_contradictory_resolver_evidence_is_declined():
    v1 = Variant("7", 117559590, "A", "G", genome=ENSEMBL81)
    v2 = Variant("7", 117559591, "T", "C", genome=ENSEMBL81)
    v3 = Variant("7", 117559592, "C", "A", genome=ENSEMBL81)
    vc = VariantCollection([v1, v2, v3])
    # Contradictory: v1~v2 cis, v2~v3 cis, v1~v3 trans!
    resolver = MolecularPhaseResolver(
        _ExplicitCisSource(cis_pairs=[(v1, v2), (v2, v3)], trans_pairs=[(v1, v3)]))
    with pytest.warns(UserWarning, match="contradictory"):
        out = combine_cis_snvs(vc, phase_resolver=resolver)
    assert list(out) == list(vc)


def test_multi_sample_somatic_cis_combines_but_conflict_declines(tmp_path):
    # Somatic MNV: tumor has 1|0 + 1|0 (same PS), normal is 0/0 + 0/0 -> compatible!
    somatic_path = _write_phased_vcf(
        tmp_path,
        "somatic_pair.vcf",
        ["tumor", "normal"],
        [
            "7\t117559590\ts1\tA\tG\t50\tPASS\t.\tGT:PS\t1|0:100\t0/0:.",
            "7\t117559591\ts2\tT\tC\t50\tPASS\t.\tGT:PS\t1|0:100\t0/0:.",
        ],
    )
    vc_som = load_vcf(somatic_path, genome=ENSEMBL81)
    out_som = combine_cis_snvs(vc_som, sample="tumor")
    assert len(out_som) == 1
    assert out_som.metadata[out_som[0]]["sample_info"]["tumor"]["GT"] == "1|0"
    assert out_som.metadata[out_som[0]]["sample_info"]["normal"]["GT"] == "0/0"

    # Multi-sample conflict: tumor has 1|0 + 1|0 (cis), normal has 1|0 + 0|1 (trans!)
    conflict_path = _write_phased_vcf(
        tmp_path,
        "multisample_conflict.vcf",
        ["tumor", "normal"],
        [
            "7\t117559590\ts1\tA\tG\t50\tPASS\t.\tGT:PS\t1|0:100\t1|0:200",
            "7\t117559591\ts2\tT\tC\t50\tPASS\t.\tGT:PS\t1|0:100\t0|1:200",
        ],
    )
    vc_conf = load_vcf(conflict_path, genome=ENSEMBL81)
    with pytest.warns(UserWarning, match="normal"):
        out_conf = combine_cis_snvs(vc_conf, sample="tumor")
    assert list(out_conf) == list(vc_conf)


def test_molecular_phase_resolver_combines_unphased_dna_vcf(tmp_path):
    # Unphased 0/1 somatic calls phased externally by RNA reads (MolecularPhaseResolver)
    path = _write_phased_vcf(
        tmp_path,
        "unphased_dna.vcf",
        ["tumor", "normal"],
        [
            "7\t117559590\ts1\tA\tG\t50\tPASS\t.\tGT:PS\t0/1:.\t0/0:.",
            "7\t117559591\ts2\tT\tC\t50\tPASS\t.\tGT:PS\t0/1:.\t0/0:.",
        ],
    )
    vc = load_vcf(path, genome=ENSEMBL81)
    resolver = MolecularPhaseResolver(_ExplicitCisSource(cis_pairs=[(vc[0], vc[1])]))
    out = combine_cis_snvs(vc, phase_resolver=resolver)
    assert len(out) == 1
    assert out[0].ref == "AT" and out[0].alt == "GC"
    assert out.metadata[out[0]]["sample_info"]["tumor"]["GT"] == "0/1"
    assert out.metadata[out[0]]["sample_info"]["normal"]["GT"] == "0/0"


# --------------------------------------------------------------------
# Same-position alternatives, reference validation, and composition
# --------------------------------------------------------------------


def test_same_position_multiallelic_alternatives_preserved_unchanged(tmp_path):
    path = _write_phased_vcf(
        tmp_path,
        "multiallelic.vcf",
        ["tumor"],
        [
            "7\t117559590\ts1\tA\tG,C\t50\tPASS\t.\tGT:PS\t1|2:100",
            "7\t117559591\ts2\tT\tC\t50\tPASS\t.\tGT:PS\t1|0:100",
        ],
    )
    vc = load_vcf(path, genome=ENSEMBL81)
    out = combine_cis_snvs(vc, sample="tumor")
    assert list(out) == list(vc)


def test_reference_mismatch_declines_combination():
    # True reference at CFTR 117559590 is A, not C -> ref validation must decline!
    v_bad = Variant("7", 117559590, "C", "G", genome=ENSEMBL81)
    v_ok = Variant("7", 117559591, "T", "C", genome=ENSEMBL81)
    vc = VariantCollection([v_bad, v_ok])
    resolver = MolecularPhaseResolver(_ExplicitCisSource(cis_pairs=[(v_bad, v_ok)]))
    with pytest.warns(UserWarning, match="disagrees with reference"):
        out = combine_cis_snvs(vc, phase_resolver=resolver)
    assert list(out) == list(vc)


def test_composition_with_left_align_and_pair_breakends_preserves_non_snvs(tmp_path):
    path = _write_phased_vcf(
        tmp_path,
        "mixed.vcf",
        ["tumor"],
        [
            "7\t117559590\ts1\tA\tG\t50\tPASS\t.\tGT:PS\t1|0:100",
            "7\t117559591\ts2\tT\tC\t50\tPASS\t.\tGT:PS\t1|0:100",
            "7\t117559600\tind\tGA\tG\t50\tPASS\t.\tGT:PS\t0/1:.",
        ],
    )
    vc = load_vcf(path, genome=ENSEMBL81)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        out1 = combine_cis_snvs(left_align_indels(pair_breakends(vc)), sample="tumor")
        out2 = left_align_indels(pair_breakends(combine_cis_snvs(vc, sample="tumor")))
    assert list(out1) == list(out2)
    assert len(out1) == 2
