"""Tests for the default ``dna`` annotator and the minimal 3-noun API
(``MutationEffect``, ``EffectCandidate``, ``MutantTranscript``).
"""

from __future__ import annotations

import json
from pathlib import Path
import time

import pytest
from pyensembl import cached_release

from varcode import (
    DNAEffectAnnotator,
    EffectCandidate,
    EffectCollection,
    GermlineContext,
    MolecularPhaseResolver,
    MutantTranscript,
    ReferenceMismatchError,
    StructuralVariant,
    Variant,
    VariantCollection,
    classification_priority,
    discover_created_orfs,
    effect_priority,
    get_annotator,
    get_default_annotator,
    iter_candidate_effects,
    iter_candidate_paths,
    parse_symbolic_alt,
    registered_annotators,
    to_proteome,
)
from varcode.effects.effect_classes import (
    ComplexSubstitution,
    Deletion,
    FrameShift,
    GeneFusion,
    Insertion,
    Substitution,
    TranslocationToIntergenic,
)
from varcode.nucleotides import reverse_complement


GENOME = cached_release(81)


def _protein(effect):
    mt = getattr(effect, "mutant_transcript", None)
    if mt is not None and mt.mutant_protein_sequence is not None:
        return mt.mutant_protein_sequence
    return getattr(effect, "mutant_protein_sequence", None)


class _CisPair:
    source = "test_dna_phase"

    def __init__(self, v1, v2, answer=True):
        self._pair = (v1, v2)
        self._answer = answer

    def has_evidence(self, variant) -> bool:
        return variant in self._pair

    def partners_in_cis(self, variant):
        if self._answer is True and variant in self._pair:
            return tuple(v for v in self._pair if v != variant)
        return ()

    def in_cis(self, v1, v2, transcript=None):
        if {v1, v2} == set(self._pair):
            return self._answer
        return None


def test_dna_annotator_is_default_with_slim_api_and_deprecations():
    assert registered_annotators() == ("dna", "fast", "transcript_model")
    assert registered_annotators(include_deprecated=True) == (
        "dna",
        "fast",
        "protein_diff",
        "transcript_model",
    )
    assert isinstance(get_annotator("dna"), DNAEffectAnnotator)
    assert get_default_annotator().name in ("dna", "protein_diff")
    assert get_annotator("fast").name == "fast"
    assert get_annotator("transcript_model").name == "transcript_model"

    with pytest.warns(DeprecationWarning, match="deprecated"):
        pdiff = get_annotator("protein_diff")
    assert pdiff.name == "protein_diff"


def test_dna_flat_effect_candidates_and_uncertainty_axes():
    braf = GENOME.transcript_by_id("ENST00000288602")  # - strand
    met = GENOME.transcript_by_id("ENST00000397752")  # + strand

    # 1. Deterministic interior exonic SNV (BRAF V600E): fast-path, no uncertainty.
    v_braf = Variant("7", 140753336, "A", "T", genome=GENOME)
    eff_braf = v_braf.effect_on_transcript(braf, annotator="dna")
    assert isinstance(eff_braf, Substitution)
    assert eff_braf.short_description == "p.V600E"
    assert not hasattr(eff_braf, "realized_candidates")
    assert len(eff_braf.candidates) == 1
    cand = eff_braf.candidates[0]
    assert isinstance(cand, EffectCandidate)
    assert cand.uncertainty == ()
    assert cand.mechanisms == ("normal",)
    assert cand.effect.mutant_protein_sequence == eff_braf.mutant_protein_sequence

    # 2. Pure splice-site variant (MET exon 14 donor +1 G>A): uncertainty=('splicing',).
    ex14 = met.exons[13]
    v_met_don = Variant("7", ex14.end + 1, "G", "A", genome=GENOME)
    eff_met = v_met_don.effect_on_transcript(met, annotator="dna")
    assert all(isinstance(c, EffectCandidate) for c in eff_met.candidates)
    assert all(c.uncertainty == ("splicing",) for c in eff_met.candidates)
    assert any("exon_skip" in c.mechanisms for c in eff_met.candidates)

    # 3. Unphased germline variant in same codon: uncertainty=('phasing',).
    g_braf = Variant("7", 140753335, "C", "T", genome=GENOME)
    ctx = GermlineContext.from_variants([g_braf])
    eff_phase = v_braf.effect_on_transcript(braf, annotator="dna", germline=ctx)
    assert len(eff_phase.candidates) == 2
    assert all(isinstance(c, EffectCandidate) for c in eff_phase.candidates)
    assert all(c.uncertainty == ("phasing",) for c in eff_phase.candidates)
    states = {state for c in eff_phase.candidates for _, state in c.phase}
    assert states == {"cis", "trans"}

    # 4. Joint phasing + splicing uncertainty on MET exon 14 last base + unphased germline.
    ex14_last_ref = met.sequence[met.spliced_offset(ex14.end)]
    ex14_last_alt = "A" if ex14_last_ref != "A" else "C"
    v_met_ex = Variant("7", ex14.end, ex14_last_ref, ex14_last_alt, genome=GENOME)
    g_met_pos = ex14.start + 30
    g_met_ref = met.sequence[met.spliced_offset(g_met_pos)]
    g_met_alt = "A" if g_met_ref != "A" else "G"
    g_met = Variant("7", g_met_pos, g_met_ref, g_met_alt, genome=GENOME)
    met_ctx = GermlineContext.from_variants([g_met])
    eff_both = v_met_ex.effect_on_transcript(met, annotator="dna", germline=met_ctx)
    assert all(isinstance(c, EffectCandidate) for c in eff_both.candidates)
    assert all(c.uncertainty == ("phasing", "splicing") for c in eff_both.candidates)
    assert list(iter_candidate_effects(eff_both))
    assert all(
        all(isinstance(item, EffectCandidate) for item in path)
        for _, path in iter_candidate_paths(eff_both)
    )


def test_mutant_transcript_novel_junctions_are_ascending_on_both_strands():
    """MutantTranscript.novel_junctions always has intron_start <= intron_end on + and - strands."""
    # + strand: MET exon 14 skipping joins exon 13 end -> exon 15 start
    met = GENOME.transcript_by_id("ENST00000397752")
    ex13, ex14, ex15 = met.exons[12], met.exons[13], met.exons[14]
    eff_met = Variant("7", ex14.end + 1, "G", "A", genome=GENOME).effect_on_transcript(
        met, annotator="dna"
    )
    skip_met = next(c for c in eff_met.candidates if "exon_skip" in c.mechanisms)
    mt_met = skip_met.effect.mutant_transcript
    assert mt_met is not None and mt_met.is_complete
    assert len(mt_met.novel_junctions) == 1
    j_met = mt_met.novel_junctions[0]
    assert j_met == ("7", ex13.end + 1, ex15.start - 1, "+")
    assert j_met[1] <= j_met[2]

    # - strand: TP53 exon 5 skipping joins exon 4 start -> exon 6 end
    tp53 = GENOME.transcript_by_id("ENST00000269305")
    ex4, ex5, ex6 = tp53.exons[3], tp53.exons[4], tp53.exons[5]
    # - strand donor +1 is ex5.start - 1 (G on - strand = C on + strand)
    eff_tp53 = Variant("17", ex5.start - 1, "C", "T", genome=GENOME).effect_on_transcript(
        tp53, annotator="dna"
    )
    skip_tp53 = next(c for c in eff_tp53.candidates if "exon_skip" in c.mechanisms)
    mt_tp53 = skip_tp53.effect.mutant_transcript
    assert mt_tp53 is not None and mt_tp53.is_complete
    assert len(mt_tp53.novel_junctions) == 1
    j_tp53 = mt_tp53.novel_junctions[0]
    # Ascending genomic intron coordinates even on - strand!
    assert j_tp53 == ("17", ex6.end + 1, ex4.start - 1, "-")
    assert j_tp53[1] <= j_tp53[2]


def test_dna_matches_transcript_model_across_cancer_hard_cases():
    """Default dna annotator matches transcript_model on EML4-ALK, KIAA1549-BRAF, FLT3-ITD, rescue."""
    # EML4--ALK inversion fusion
    eml4 = GENOME.transcript_by_id("ENST00000318522")
    alk = GENOME.transcript_by_id("ENST00000389048")
    sv_ea = parse_symbolic_alt(
        "2", eml4.exons[12].end + 500, "N", f"N]{alk.contig}:{alk.exons[19].end + 500}]",
        genome=GENOME,
    )
    tm_ea = sv_ea.effect_on_transcript(eml4, annotator="transcript_model")
    u_ea = sv_ea.effect_on_transcript(eml4, annotator="dna")
    assert isinstance(u_ea, GeneFusion)
    assert _protein(u_ea) == _protein(tm_ea)
    assert all(isinstance(c, EffectCandidate) for c in u_ea.candidates)

    # KIAA1549--BRAF tandem duplication fusion
    kiaa = GENOME.transcript_by_id("ENST00000440172")
    braf = GENOME.transcript_by_id("ENST00000288602")
    sv_kb = StructuralVariant(
        "7", kiaa.exons[15].start - 500, "DUP", end=braf.exons[8].end + 500, genome=GENOME
    )
    tm_kb = sv_kb.effect_on_transcript(kiaa, annotator="transcript_model")
    u_kb = sv_kb.effect_on_transcript(kiaa, annotator="dna")
    assert isinstance(u_kb, GeneFusion)
    assert _protein(u_kb) == _protein(tm_kb)

    # FLT3-ITD intra-exonic DUP
    flt3 = GENOME.transcript_by_id("ENST00000241453")
    ex14 = flt3.exons[13]
    mid = (ex14.start + ex14.end) // 2
    sv_flt3 = StructuralVariant("13", mid, "DUP", end=mid + 20, genome=GENOME)
    tm_flt3 = sv_flt3.effect_on_transcript(flt3, annotator="transcript_model")
    u_flt3 = sv_flt3.effect_on_transcript(flt3, annotator="dna")
    assert isinstance(u_flt3, Insertion)
    assert _protein(u_flt3) == _protein(tm_flt3)

    # Compound frameshift rescue in cis (EGFR exon 19)
    egfr = GENOME.transcript_by_id("ENST00000275493")
    ex19 = egfr.exons[18]
    p1 = ex19.start + 30
    p2 = p1 + 9
    o1 = egfr.spliced_offset(p1)
    o2 = egfr.spliced_offset(p2)
    r1 = Variant("7", p1, egfr.sequence[o1 : o1 + 2], egfr.sequence[o1], genome=GENOME)
    r2 = Variant("7", p2, egfr.sequence[o2], egfr.sequence[o2] + "G", genome=GENOME)
    resolver = MolecularPhaseResolver(_CisPair(r1, r2, answer=True))
    tm_effs = VariantCollection([r1, r2]).effects(
        annotator="transcript_model", phase_resolver=resolver
    )
    u_effs = VariantCollection([r1, r2]).effects(annotator="dna", phase_resolver=resolver)
    tm_joint = next(e for e in tm_effs if e.transcript_id == egfr.id and len(getattr(e, "variants", ())) == 2)
    u_joint = next(e for e in u_effs if e.transcript_id == egfr.id and len(getattr(e, "variants", ())) == 2)
    assert isinstance(u_joint, (Substitution, ComplexSubstitution))
    assert type(u_joint) is type(tm_joint)
    assert _protein(u_joint) == _protein(tm_joint)


def test_open_mutant_transcript_template_for_big_intronic_insertion():
    """Big intronic insertion (>=50 bp) attaches open MutantTranscript template without guessing splicing."""
    met = GENOME.transcript_by_id("ENST00000397752")  # + strand
    ex13, ex14, ex15 = met.exons[12], met.exons[13], met.exons[14]
    ins_pos = ex14.end + 400
    assert ins_pos < ex15.start - 50
    big_insert = "GCTAGC" * 80  # 480 bp
    sv_ins = Variant("7", ins_pos, "A", "A" + big_insert, genome=GENOME)

    eff = sv_ins.effect_on_transcript(met, annotator="dna")
    mt = eff.mutant_transcript
    assert isinstance(mt, MutantTranscript)
    ev = mt.evidence
    assert ev["last_intact_exon"] == 14
    assert ev["skip_fallback_exon"] == 13
    assert ev["three_prime_reentry_exon"] == 15
    assert ev["breakpoint_in_intron"] is True
    assert mt.mate_region is None
    assert ev["inserted_sequence"] == big_insert
    assert mt.exit_frame in (0, 1, 2)
    assert ev["three_prime_entry_frame"] in (0, 1, 2)
    assert mt.exit_frame == ev["three_prime_entry_frame"]

    wt_j = mt.wildtype_junctions
    assert ("7", ex13.end + 1, ex14.start - 1, "+") in wt_j
    assert ("7", ex14.end + 1, ex15.start - 1, "+") in wt_j

    assert [seg[0] for seg in ev["template_segments"]] == [
        "five_prime_host",
        "inserted_sequence",
        "three_prime_host",
    ]

    def fake_provider(contig, start, end):
        assert contig == "7" and start <= end
        return "A" * (end - start + 1)

    template = mt.template_sequence(fake_provider)
    assert big_insert in template
    seg5, _, seg3 = ev["template_segments"]
    expected_len = (seg5[3] - seg5[2] + 1) + len(big_insert) + (seg3[3] - seg3[2] + 1)
    assert len(template) == expected_len


def test_open_mutant_transcript_template_for_gene_to_intergenic_both_strands():
    """Gene->intergenic translocations attach open MutantTranscript templates on + and - strands."""
    # + strand: MET intron 14 -> intergenic chr6:44218908 moving '+'
    met = GENOME.transcript_by_id("ENST00000397752")
    ex14 = met.exons[13]
    sv_plus = parse_symbolic_alt(
        "7", ex14.end + 300, "N", "N[6:44218908[", genome=GENOME
    )
    eff_plus = sv_plus.effect_on_transcript(met, annotator="dna")
    assert isinstance(eff_plus, TranslocationToIntergenic)
    mt_p = eff_plus.mutant_transcript
    assert mt_p is not None and not mt_p.is_complete
    ev_p = mt_p.evidence
    assert ev_p["sequence_status"] == "retained_reference_fragment"
    assert ev_p["last_intact_exon"] == 14
    assert ev_p["three_prime_reentry_exon"] is None
    assert mt_p.mate_region == ("6", 44218908, 44218908 + 25_000 - 1, "+")
    assert [seg[0] for seg in ev_p["template_segments"]] == [
        "five_prime_host",
        "mate_search_window",
    ]

    # - strand: TP53 intron 5 (81 bp wide) -> intergenic chr10:40000000 moving '-'
    tp53 = GENOME.transcript_by_id("ENST00000269305")
    ex5 = tp53.exons[4]
    intergenic_pos = 40_000_000
    assert not GENOME.genes_at_locus("10", intergenic_pos)
    sv_minus = parse_symbolic_alt(
        "17", ex5.start - 40, "N", f"]10:{intergenic_pos}]N", genome=GENOME
    )
    eff_minus = sv_minus.effect_on_transcript(tp53, annotator="dna")
    assert isinstance(eff_minus, TranslocationToIntergenic)
    mt_m = eff_minus.mutant_transcript
    assert mt_m is not None and not mt_m.is_complete
    ev_m = mt_m.evidence
    assert ev_m["last_intact_exon"] == 5
    assert ev_m["skip_fallback_exon"] == 4
    assert ev_m["breakpoint_in_intron"] is True
    assert mt_m.mate_region == (
        "10",
        intergenic_pos - 25_000 + 1,
        intergenic_pos,
        "-",
    )
    assert mt_m.mate_region[1] <= mt_m.mate_region[2]

    def strand_provider(contig, start, end):
        return ("A" * (end - start)) + "C"

    tmpl = mt_m.template_sequence(strand_provider)
    assert tmpl.startswith("G")  # RC of trailing 'C' on '-' five_prime_host

    # 5'UTR-only TMPRSS2--ERG honestly reports exit_frame=None.
    tmprss2 = GENOME.transcript_by_id("ENST00000332149")
    erg = GENOME.transcript_by_id("ENST00000288319")
    sv_te = parse_symbolic_alt(
        "21", tmprss2.exons[0].start - 500, "N", f"]21:{erg.exons[3].end + 500}]N", genome=GENOME
    )
    eff_te = sv_te.effect_on_transcript(tmprss2, annotator="dna")
    assert eff_te.mutant_transcript.exit_frame is None


def test_fast_path_safety_guards_and_reference_mismatch():
    """Fast-path must never swallow reference mismatches, splice sites, start/stop, mito, or Sec."""
    braf = GENOME.transcript_by_id("ENST00000288602")  # - strand
    bad = Variant("7", 140753336, "C", "T", genome=GENOME)
    with pytest.raises(ReferenceMismatchError):
        bad.effect_on_transcript(braf, annotator="dna")

    start_pos = max(braf.start_codon_positions)
    v_start = Variant("7", start_pos, "T", "G", genome=GENOME)
    u_start = v_start.effect_on_transcript(braf, annotator="dna")
    tm_start = v_start.effect_on_transcript(braf, annotator="transcript_model")
    assert type(u_start) is type(tm_start)
    assert u_start.short_description == tm_start.short_description

    mt_v = Variant("MT", 12994, "G", "A", genome=GENOME)
    mt_u = mt_v.effects(annotator="dna")
    mt_tm = mt_v.effects(annotator="transcript_model")
    assert [(type(e).__name__, e.short_description, e.mutant_protein_sequence) for e in mt_u] == [
        (type(e).__name__, e.short_description, e.mutant_protein_sequence) for e in mt_tm
    ]

    gpx4 = GENOME.transcript_by_id("ENST00000354171")
    ex = gpx4.exons[2]
    mid = (ex.start + ex.end) // 2
    ref_b = gpx4.sequence[gpx4.spliced_offset(mid)]
    alt_b = "A" if ref_b != "A" else "G"
    v_sec = Variant(gpx4.contig, mid, ref_b, alt_b, genome=GENOME)
    u_sec = v_sec.effect_on_transcript(gpx4, annotator="dna")
    tm_sec = v_sec.effect_on_transcript(gpx4, annotator="transcript_model")
    assert type(u_sec) is type(tm_sec)
    assert u_sec.mutant_protein_sequence == tm_sec.mutant_protein_sequence


def test_dna_json_roundtrip_preserves_candidates_junctions_and_templates():
    """DNA effects with EffectCandidate, novel_junctions, and open MutantTranscript round-trip."""
    met = GENOME.transcript_by_id("ENST00000397752")
    ex14 = met.exons[13]
    e_splice = Variant("7", ex14.end + 1, "G", "A", genome=GENOME).effect_on_transcript(
        met, annotator="dna"
    )
    sv_trans = parse_symbolic_alt(
        "7", ex14.end + 300, "N", "N[6:44218908[", genome=GENOME
    )
    e_trans = sv_trans.effect_on_transcript(met, annotator="dna")

    col = EffectCollection([e_splice, e_trans])
    restored = EffectCollection.from_json(col.to_json())
    by_var = {e.variant: e for e in restored}

    r_splice = by_var[e_splice.variant]
    assert all(isinstance(c, EffectCandidate) for c in r_splice.candidates)
    r_skip = next(c for c in r_splice.candidates if "exon_skip" in c.mechanisms)
    assert len(r_skip.effect.mutant_transcript.novel_junctions) == 1
    j = r_skip.effect.mutant_transcript.novel_junctions[0]
    assert j[1] <= j[2]

    r_trans = by_var[e_trans.variant]
    assert isinstance(r_trans.mutant_transcript, MutantTranscript)
    assert r_trans.mutant_transcript.exit_frame == e_trans.mutant_transcript.exit_frame
    assert r_trans.mutant_transcript.mate_region == e_trans.mutant_transcript.mate_region
    assert r_trans.mutant_transcript.wildtype_junctions == e_trans.mutant_transcript.wildtype_junctions
    tmpl = r_trans.mutant_transcript.template_sequence(lambda c, s, e: "A" * (e - s + 1))
    assert len(tmpl) > 25_000


def test_dna_coding_snv_batch_parity_and_speedup_over_transcript_model():
    """On interior exonic SNVs, dna matches transcript_model byte-for-byte and runs faster."""
    met = GENOME.transcript_by_id("ENST00000397752")
    ex2 = met.exons[1]
    cds_start = min(met.start_codon_spliced_offsets)
    variants = []
    pos = ex2.start + 30
    while len(variants) < 60 and pos < ex2.end - 10:
        off = met.spliced_offset(pos)
        if off >= cds_start + 6:
            ref = met.sequence[off]
            alt = "A" if ref != "A" else "G"
            variants.append(Variant("7", pos, ref, alt, genome=GENOME))
        pos += 15
    assert len(variants) == 60

    variants[0].effect_on_transcript(met, annotator="transcript_model")
    variants[0].effect_on_transcript(met, annotator="dna")

    t0 = time.perf_counter()
    tm_effects = [v.effect_on_transcript(met, annotator="transcript_model") for v in variants]
    t_tm = time.perf_counter() - t0

    t1 = time.perf_counter()
    u_effects = [v.effect_on_transcript(met, annotator="dna") for v in variants]
    t_u = time.perf_counter() - t1

    for u_eff, tm_eff in zip(u_effects, tm_effects):
        assert type(u_eff) is type(tm_eff)
        assert u_eff.short_description == tm_eff.short_description
        assert u_eff.mutant_protein_sequence == tm_eff.mutant_protein_sequence

    assert t_u < t_tm, f"Expected dna ({t_u:.4f}s) faster than transcript_model ({t_tm:.4f}s)"


def test_dna_osteosarc_full_cohort_parity_with_transcript_model():
    """All 179 real clinical osteosarc variants match transcript_model byte-for-byte."""
    from .osteosarc_variants import load_variants

    variants = list(load_variants(genome=GENOME))
    assert len(variants) == 179
    for variant in variants:
        tm_effects = list(variant.effects(annotator="transcript_model", raise_on_error=True))
        u_effects = list(variant.effects(annotator="dna", raise_on_error=True))
        assert len(u_effects) == len(tm_effects), f"Effect count mismatch on {variant}"
        for u_eff, tm_eff in zip(u_effects, tm_effects):
            assert getattr(u_eff.transcript, "id", None) == getattr(tm_eff.transcript, "id", None)
            assert type(u_eff) is type(tm_eff), f"Class mismatch on {variant}"
            assert u_eff.short_description == tm_eff.short_description
            assert u_eff.mutant_protein_sequence == tm_eff.mutant_protein_sequence
            assert u_eff.modifies_protein_sequence == tm_eff.modifies_protein_sequence


def test_dna_interior_mnv_and_indel_fast_path_parity_both_strands():
    """Interior coding MNVs, in-frame indels, frameshifts, and Silent match transcript_model."""
    met = GENOME.transcript_by_id("ENST00000397752")
    brca1 = GENOME.transcript_by_id("ENST00000357654")
    comp = str.maketrans("ACGT", "TGCA")

    test_cases = []
    for tx, ex in ((met, met.exons[1]), (brca1, brca1.exons[9])):
        cds_start = min(tx.start_codon_spliced_offsets)
        cds_end = max(tx.stop_codon_spliced_offsets)
        step = max(40, (ex.end - ex.start) // 12)
        for pos in range(ex.start + 30, ex.end - 30, step):
            off = tx.spliced_offset(pos)
            if not (cds_start + 12 <= off <= cds_end - 20):
                continue
            if tx.strand == "+":
                ref4 = tx.sequence[off : off + 4]
            else:
                off_end = tx.spliced_offset(pos + 3)
                ref4 = tx.sequence[off_end : off + 1][::-1].translate(comp)
            if len(ref4) != 4 or any(b not in "ACGT" for b in ref4):
                continue
            b0 = ref4[0]
            alt0 = "A" if b0 != "A" else "G"
            alt1 = "C" if ref4[1] != "C" else "T"
            test_cases.extend([
                (tx, Variant(tx.contig, pos, b0, alt0, genome=GENOME)),
                (tx, Variant(tx.contig, pos, ref4[:2], alt0 + alt1, genome=GENOME)),
                (tx, Variant(tx.contig, pos, ref4, b0, genome=GENOME)),
                (tx, Variant(tx.contig, pos, b0, b0 + "GCA", genome=GENOME)),
                (tx, Variant(tx.contig, pos, b0, b0 + "TAA", genome=GENOME)),
                (tx, Variant(tx.contig, pos, ref4[:2], b0, genome=GENOME)),
                (tx, Variant(tx.contig, pos, b0, b0 + "G", genome=GENOME)),
            ])
    assert len(test_cases) >= 70

    for tx, variant in test_cases:
        tm_eff = variant.effect_on_transcript(tx, annotator="transcript_model")
        u_eff = variant.effect_on_transcript(tx, annotator="dna")
        assert type(u_eff) is type(tm_eff), f"Class mismatch on {variant} ({tx.id})"
        assert u_eff.short_description == tm_eff.short_description
        assert u_eff.mutant_protein_sequence == tm_eff.mutant_protein_sequence
        assert u_eff.modifies_protein_sequence == tm_eff.modifies_protein_sequence

    ex2 = met.exons[1]
    pos = ex2.start + 60
    real_seq = met.sequence[met.spliced_offset(pos) : met.spliced_offset(pos) + 4]
    wrong4 = "AAAA" if real_seq != "AAAA" else "CCCC"
    bad_del = Variant("7", pos, wrong4, wrong4[0], genome=GENOME)
    with pytest.raises(ReferenceMismatchError):
        bad_del.effect_on_transcript(met, annotator="dna")


def test_dna_interoperates_with_priority_proteome_orf_discovery_and_mt_polya():
    """DNA effects work with effect_priority, to_proteome, discover_created_orfs, and MT poly-A."""
    cftr = GENOME.transcript_by_id("ENST00000003084")
    weak_donor = Variant("7", 117531117, "A", "C", genome=GENOME)
    u_splice = weak_donor.effect_on_transcript(cftr, annotator="dna")
    del_alt = next(c.effect for c in u_splice.candidates if isinstance(c.effect, Deletion))
    assert effect_priority(u_splice) >= classification_priority(del_alt)
    assert list(EffectCollection([u_splice]).filter_by_effect_priority(Deletion)) == [u_splice]

    braf = GENOME.transcript_by_id("ENST00000288602")
    u_v600e = Variant("7", 140753336, "A", "T", genome=GENOME).effect_on_transcript(
        braf, annotator="dna"
    )
    germline_ctx = GermlineContext.from_variants([Variant("7", 140753335, "C", "T", genome=GENOME)])
    u_phased = Variant("7", 140753336, "A", "T", genome=GENOME).effect_on_transcript(
        braf, germline=germline_ctx, annotator="dna"
    )
    proteome = to_proteome([u_v600e, u_splice, u_phased])
    assert proteome.proteins
    assert all(p.producer.tool == "varcode" for p in proteome.proteins)

    fixture = json.loads((Path(__file__).parent / "data/created_orfs/CDKN2A.json").read_text())
    native = fixture["native"]
    g95 = cached_release(native["annotation_version"])
    cdkn2a_tx = g95.transcript_by_id(native["transcript_id"])
    allele = fixture["genomic"]
    v_start = Variant(allele["contig"], allele["start"], allele["ref"], allele["alt"], genome=g95)
    u_utr = v_start.effect_on_transcript(cdkn2a_tx, annotator="dna")
    report = discover_created_orfs(u_utr)
    assert report.complete
    assert len(report.effects) == 1
    assert report.effects[0].mutant_protein_sequence == fixture["refseq"]["expected_protein"]

    nd1 = GENOME.transcript_by_id("ENST00000361390")
    mt_dup = StructuralVariant("MT", 3500, "DUP", end=3560, genome=GENOME)
    u_mt = mt_dup.effect_on_transcript(nd1, annotator="dna")
    assert isinstance(u_mt.mutant_transcript, MutantTranscript)
    assert u_mt.mutant_transcript.exit_frame in (0, 1, 2)


def test_dna_fast_path_and_rna_refinement_preserve_observations_and_support():
    """#776: fast-path dna effects refine cleanly and preserve rna_observation / rna_support."""
    from types import SimpleNamespace
    from varcode import (
        RNAEvidence,
        apply_phase_resolver_to_effects,
        apply_rna_evidence_to_effects,
        make_rna_outcome,
    )
    from varcode.effects import Unresolved

    braf = GENOME.transcript_by_id("ENST00000288602")
    v_v600e = Variant("7", 140753336, "A", "T", genome=GENOME)
    fast_effect = v_v600e.effect_on_transcript(braf, annotator="dna")
    assert isinstance(fast_effect, Substitution)
    # Fast-path SNV does not carry paired transcript-model outcomes initially.
    assert fast_effect.candidates[0].outcomes == ()

    # 1. Partial RNA window supplied via phase_resolver must NOT overwrite Substitution.mutant_transcript;
    #    it must append an explicit Unresolved candidate with original_protein_sequence=None and rna_observation.
    partial_window = MutantTranscript.from_sequence(
        "ATGGCTTAA",
        reference_transcript=braf,
        mutant_protein_sequence="MA",
        annotator_name="isovar",
        evidence={"sequence_status": "junction_fragment", "read_count": 5},
    )
    resolver = SimpleNamespace(
        mutant_transcript=lambda variant, transcript: (
            partial_window if variant == v_v600e and transcript.id == braf.id else None
        )
    )
    refined_phase = apply_phase_resolver_to_effects(EffectCollection([fast_effect]), resolver)[0]
    assert isinstance(refined_phase, Substitution)
    assert refined_phase.mutant_protein_sequence == fast_effect.mutant_protein_sequence
    obs_cands = [
        c for c in refined_phase.candidates
        if isinstance(c.effect, Unresolved)
        and getattr(c.effect.mutant_transcript, "cdna_sequence", None) == "ATGGCTTAA"
    ]
    assert len(obs_cands) == 1
    assert obs_cands[0].effect.original_protein_sequence is None
    assert obs_cands[0].effect.modifies_protein_sequence is None
    assert obs_cands[0].evidence["rna_observation"]["read_count"] == 5

    # 2. Full predicted cDNA matched via apply_rna_evidence_to_effects preserves rna_support on EffectCandidate.
    full_cdna = refined_phase.mutant_transcript.cdna_sequence
    assert full_cdna is not None
    full_obs = Unresolved(v_v600e, braf, "observed_rna")
    full_obs.mutant_transcript = MutantTranscript.from_sequence(
        full_cdna,
        reference_transcript=braf,
        mutant_protein_sequence=refined_phase.mutant_protein_sequence,
        annotator_name="reads",
    )
    rna_evidence = RNAEvidence(
        [make_rna_outcome(full_obs, source="reads", transcript_model_id="full", read_count=12)]
    )
    refined_rna = apply_rna_evidence_to_effects(EffectCollection([fast_effect]), rna_evidence)[0]
    sub_cand = next(c for c in refined_rna.candidates if isinstance(c.effect, Substitution))
    assert sub_cand.evidence["rna_support"][0]["read_count"] == 12
    assert sub_cand.evidence["rna_support"][0]["transcript_model_id"] == "full"


def test_dna_effect_to_dict_and_json_roundtrip_preserves_candidates_and_phase():
    """#778: plain dna SNVs and phased/splice dna effects round-trip candidates losslessly."""
    from varcode import MutationEffect, get_effect_candidates, iter_candidate_paths

    assert not hasattr(MutationEffect, "_effect_graph_fields")

    braf = GENOME.transcript_by_id("ENST00000288602")
    v_v600e = Variant("7", 140753336, "A", "T", genome=GENOME)

    # Plain fast-path dna SNV (BRAF V600E) must preserve candidates on to_dict/from_dict & JSON.
    eff = v_v600e.effect_on_transcript(braf)
    assert isinstance(eff, Substitution)
    assert len(get_effect_candidates(eff)) == 1
    assert eff.candidates[0].source == "varcode.dna"

    restored_dict = type(eff).from_dict(eff.to_dict())
    assert len(get_effect_candidates(restored_dict)) == 1
    assert restored_dict.candidates[0].source == "varcode.dna"
    assert restored_dict.candidates[0].effect is restored_dict
    assert restored_dict.most_likely_candidate is restored_dict.candidates[0]
    assert restored_dict.most_likely_effect is restored_dict
    assert restored_dict.highest_priority_candidate is restored_dict.candidates[0]
    assert restored_dict.highest_priority_effect is restored_dict
    assert restored_dict.effects == (restored_dict,)
    assert "most_likely_candidate" not in restored_dict.__dict__
    assert "highest_priority_candidate" not in restored_dict.__dict__
    assert len(list(iter_candidate_paths(restored_dict))) == 1
    assert hash(restored_dict) == hash(eff)
    assert restored_dict.to_dict() == eff.to_dict()

    restored_json = type(eff).from_json(eff.to_json())
    assert len(get_effect_candidates(restored_json)) == 1
    assert restored_json.to_dict() == eff.to_dict()

    restored_coll = EffectCollection.from_json(EffectCollection([eff]).to_json())
    assert len(get_effect_candidates(restored_coll[0])) == 1
    assert restored_coll[0].to_dict() == eff.to_dict()

    # Legacy fast annotator plain SNV retains flat dict without effect_graph_schema.
    legacy = v_v600e.effect_on_transcript(braf, annotator="fast")
    assert "effect_graph_schema" not in legacy.to_dict()

    # Phased germline dna effect round-trips EffectCandidate.phase Variant references cleanly.
    germline_v = Variant("7", 140753335, "C", "T", genome=GENOME)
    phased = v_v600e.effect_on_transcript(
        braf, germline=GermlineContext.from_variants([germline_v])
    )
    restored_phased = type(phased).from_json(phased.to_json())
    assert {
        (v.contig, v.start, v.ref, v.alt, state)
        for c in restored_phased.candidates
        for v, state in c.phase
    } == {
        (v.contig, v.start, v.ref, v.alt, state)
        for c in phased.candidates
        for v, state in c.phase
    }
    assert restored_phased.to_dict() == phased.to_dict()

    # Splice-site dna effect preserves multi-candidate graph and novel_junctions.
    cftr = GENOME.transcript_by_id("ENST00000003084")
    splice_eff = Variant("7", 117531117, "A", "C", genome=GENOME).effect_on_transcript(cftr)
    restored_splice = type(splice_eff).from_json(splice_eff.to_json())
    assert len(restored_splice.candidates) == len(splice_eff.candidates)
    assert restored_splice.to_dict() == splice_eff.to_dict()


def test_uniform_mutation_effect_candidate_api_and_outcomes_roundtrip():
    """Uniform public candidate API on MutationEffect and lossless EffectCandidate.outcomes."""
    from varcode import (
        EffectCandidate,
        MultiOutcomeEffect,
        MutationEffect,
        RNAEvidence,
        apply_rna_evidence_to_effects,
        make_rna_outcome,
    )
    from varcode.effects import StructuralVariantEffect, Unresolved

    # Weird private-member plumbing is gone.
    assert not hasattr(MutationEffect, "_attached")
    assert not hasattr(MultiOutcomeEffect, "_combine_with_extra_candidates")
    assert not hasattr(StructuralVariantEffect, "_attach_primary_effects")
    assert not hasattr(StructuralVariantEffect, "_attach_cryptic_candidates")
    assert not hasattr(StructuralVariantEffect, "_attach_splice_outcomes")

    braf = GENOME.transcript_by_id("ENST00000288602")
    v_v600e = Variant("7", 140753336, "A", "T", genome=GENOME)

    # Bare fast-annotator effects expose the same uniform candidate interface.
    bare = v_v600e.effect_on_transcript(braf, annotator="fast")
    assert bare.candidates == ()
    assert bare.effects == ()
    assert bare.most_likely_candidate is None
    assert bare.most_likely_effect is None
    assert bare.highest_priority_candidate is None
    assert bare.highest_priority_effect is None
    extra = EffectCandidate(bare, source="custom", evidence={"note": "manual"})
    bare.add_candidates((extra,))
    assert bare.candidates == (extra,)
    assert bare.most_likely_candidate is extra
    assert bare.most_likely_effect is bare
    assert bare.highest_priority_candidate is extra
    assert bare.highest_priority_effect is bare

    # Phased dna effects carry public EffectCandidate.outcomes (no _paired_outcomes)
    # and preserve them across JSON round-trip so downstream RNA refinement works.
    germline_v = Variant("7", 140753335, "C", "T", genome=GENOME)
    phased = v_v600e.effect_on_transcript(
        braf, germline=GermlineContext.from_variants([germline_v])
    )
    assert any(c.outcomes for c in phased.candidates)
    assert all("_paired_outcomes" not in c.__dict__ for c in phased.candidates)

    restored = type(phased).from_json(phased.to_json())
    assert [len(c.outcomes) for c in restored.candidates] == [
        len(c.outcomes) for c in phased.candidates
    ]
    assert restored.most_likely_candidate is restored.candidates[0]

    # RNA refinement works on the deserialized effect because .outcomes survived JSON.
    observed = Unresolved(v_v600e, braf, "rna_transcript_comparison")
    observed.mutant_transcript = MutantTranscript.from_sequence(
        "ATGGCTTAA",
        reference_transcript=braf,
        mutant_protein_sequence="MA",
        annotator_name="isovar",
    )
    refined_coll = apply_rna_evidence_to_effects(
        EffectCollection([restored]),
        RNAEvidence([make_rna_outcome(observed, source="isovar", read_count=4)]),
    )
    assert any(c.evidence.get("rna_observation") for c in refined_coll[0].candidates)
