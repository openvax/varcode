"""Tests for the opt-in ``unified`` DNA-only annotator, flat hypotheses, and search anchors."""

from __future__ import annotations

import time

import pytest
from pyensembl import cached_release

from varcode import (
    DNAHypothesis,
    DNASearchAnchor,
    EffectCollection,
    GenomicSpliceJunction,
    GermlineContext,
    MolecularPhaseResolver,
    ReferenceMismatchError,
    StructuralVariant,
    UnifiedDNAEffectAnnotator,
    Variant,
    VariantCollection,
    build_dna_search_anchor,
    get_annotator,
    iter_candidate_effects,
    iter_candidate_paths,
    parse_symbolic_alt,
    registered_annotators,
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
        if self._answer is not True or variant not in self._pair:
            return ()
        return tuple(v for v in self._pair if v != variant)

    def in_cis(self, v1, v2, transcript=None):
        if {v1, v2} == set(self._pair):
            return self._answer
        return None


def test_unified_annotator_is_registered_opt_in_without_changing_default():
    assert "unified" in registered_annotators()
    assert isinstance(get_annotator("unified"), UnifiedDNAEffectAnnotator)
    assert get_annotator("fast").name == "fast"
    assert get_annotator("transcript_model").name == "transcript_model"


def test_unified_flat_dna_hypotheses_and_uncertainty_axes():
    braf = GENOME.transcript_by_id("ENST00000288602")  # - strand
    met = GENOME.transcript_by_id("ENST00000397752")  # + strand

    # 1. Deterministic interior exonic SNV (BRAF V600E): fast-path, no uncertainty.
    v_braf = Variant("7", 140753336, "A", "T", genome=GENOME)
    eff_braf = v_braf.effect_on_transcript(braf, annotator="unified")
    assert isinstance(eff_braf, Substitution)
    assert eff_braf.short_description == "p.V600E"
    assert not hasattr(eff_braf, "realized_candidates")
    assert len(eff_braf.candidates) == 1
    hyp = eff_braf.candidates[0]
    assert isinstance(hyp, DNAHypothesis)
    assert hyp.uncertainty == ()
    assert hyp.splice_mechanisms == ("normal",)
    assert hyp.mutant_protein == eff_braf.mutant_protein_sequence

    # 2. Splice-only uncertainty (MET exon 14 donor +1).
    ex14 = met.exons[13]
    v_met_don = Variant("7", ex14.end + 1, "G", "A", genome=GENOME)
    eff_met = v_met_don.effect_on_transcript(met, annotator="unified")
    assert all(isinstance(c, DNAHypothesis) for c in eff_met.candidates)
    assert eff_met.candidates[0].uncertainty == ("splicing",)
    mechs = {m for c in eff_met.candidates for m in c.splice_mechanisms}
    assert {"exon_skip", "intron_retention", "cryptic", "normal"} <= mechs

    # 3. Phasing-only uncertainty (BRAF V600E + unphased germline V600M).
    v_germ = Variant("7", 140753337, "C", "T", genome=GENOME)
    ctx = GermlineContext.from_variants([v_germ], reference_name="GRCh38")
    eff_phase = v_braf.effect_on_transcript(braf, annotator="unified", germline=ctx)
    assert len(eff_phase.candidates) == 2
    assert all(c.uncertainty == ("phasing",) for c in eff_phase.candidates)
    assert {
        (c.baseline_protein[599], c.mutant_protein[599], c.short_description)
        for c in eff_phase.candidates
    } == {("V", "E", "p.V600E"), ("M", "K", "p.M600K")}

    # 4. Combined splicing + phasing uncertainty (MET ex14 last base + unphased ex14 germline SNP).
    v_met_ex = Variant("7", ex14.end, "G", "A", genome=GENOME)
    g_pos = ex14.start + 30
    g_ref = met.sequence[met.spliced_offset(g_pos)]
    g_alt = "T" if g_ref != "T" else "C"
    met_ctx = GermlineContext.from_variants(
        [Variant("7", g_pos, g_ref, g_alt, genome=GENOME)], reference_name="GRCh38"
    )
    eff_both = v_met_ex.effect_on_transcript(met, annotator="unified", germline=met_ctx)
    assert eff_both.candidates[0].uncertainty == ("phasing", "splicing")

    # Works cleanly through shared traversal helpers.
    paths = list(iter_candidate_paths(eff_met))
    assert len(paths) == len(eff_met.candidates)
    assert all(isinstance(path[-1], DNAHypothesis) for _, path in paths)
    assert len(list(iter_candidate_effects(eff_met))) == len(eff_met.candidates)


def test_genomic_splice_junctions_are_ascending_on_both_strands():
    """Novel junctions always satisfy intron_start <= intron_end (BAM CIGAR N convention)."""
    # + strand: MET exon 14 skip joins exon 13 -> exon 15
    met = GENOME.transcript_by_id("ENST00000397752")
    ex13, ex14, ex15 = met.exons[12], met.exons[13], met.exons[14]
    eff_met = Variant("7", ex14.end + 1, "G", "A", genome=GENOME).effect_on_transcript(
        met, annotator="unified"
    )
    met_skip = next(c for c in eff_met.candidates if "exon_skip" in c.splice_mechanisms)
    (j_met,) = met_skip.novel_junctions
    assert j_met == GenomicSpliceJunction(
        contig="7",
        intron_start=ex13.end + 1,
        intron_end=ex15.start - 1,
        strand="+",
        donor_base=ex13.end,
        acceptor_base=ex15.start,
    )
    assert j_met.intron_start <= j_met.intron_end

    # - strand: TP53 exon 5 skip joins exon 4 -> exon 6
    tp53 = GENOME.transcript_by_id("ENST00000269305")
    ex4, ex5, ex6 = tp53.exons[3], tp53.exons[4], tp53.exons[5]
    eff_tp53 = Variant("17", ex5.start - 1, "C", "T", genome=GENOME).effect_on_transcript(
        tp53, annotator="unified"
    )
    tp53_skip = next(c for c in eff_tp53.candidates if "exon_skip" in c.splice_mechanisms)
    (j_tp53,) = tp53_skip.novel_junctions
    assert j_tp53 == GenomicSpliceJunction(
        contig="17",
        intron_start=ex6.end + 1,
        intron_end=ex4.start - 1,
        strand="-",
        donor_base=ex4.start,
        acceptor_base=ex6.end,
    )
    assert j_tp53.intron_start <= j_tp53.intron_end


def test_unified_matches_transcript_model_across_cancer_hard_cases():
    """Unified matches transcript_model across fusions, splice drivers, ITDs, and rescue."""
    # EML4--ALK (+ / - inversion)
    eml4 = GENOME.transcript_by_id("ENST00000318522")
    alk = GENOME.transcript_by_id("ENST00000389048")
    sv_ea = parse_symbolic_alt(
        "2", eml4.exons[12].end + 500, "N", f"N]2:{alk.exons[19].end + 500}]", genome=GENOME
    )
    u_ea = sv_ea.effect_on_transcript(eml4, annotator="unified")
    tm_ea = sv_ea.effect_on_transcript(eml4, annotator="transcript_model")
    assert isinstance(u_ea, GeneFusion)
    assert _protein(u_ea) == _protein(tm_ea)
    assert len(_protein(u_ea)) == 1059
    assert isinstance(u_ea.dna_search_anchor, DNASearchAnchor)
    assert u_ea.dna_search_anchor.mate_search_region[3] == "-"

    # KIAA1549--BRAF (- / - DUP)
    kiaa = GENOME.transcript_by_id("ENST00000440172")
    braf = GENOME.transcript_by_id("ENST00000288602")
    sv_kb = StructuralVariant(
        "7", kiaa.exons[15].start - 500, "DUP", end=braf.exons[8].end + 500, genome=GENOME
    )
    u_kb = sv_kb.effect_on_transcript(kiaa, annotator="unified")
    tm_kb = sv_kb.effect_on_transcript(kiaa, annotator="transcript_model")
    assert isinstance(u_kb, GeneFusion)
    assert _protein(u_kb) == _protein(tm_kb)
    assert len(_protein(u_kb)) == 2135

    # FLT3-ITD (- strand intra-exonic DUP)
    flt3 = GENOME.transcript_by_id("ENST00000241453")
    sv_flt3 = StructuralVariant("13", 28034110, "DUP", end=28034130, genome=GENOME)
    u_flt3 = sv_flt3.effect_on_transcript(flt3, annotator="unified")
    tm_flt3 = sv_flt3.effect_on_transcript(flt3, annotator="transcript_model")
    assert isinstance(u_flt3, Insertion)
    assert u_flt3.mutant_protein_sequence == tm_flt3.mutant_protein_sequence
    assert len(u_flt3.mutant_protein_sequence) == 1000
    assert isinstance(u_flt3.dna_search_anchor, DNASearchAnchor)

    # Compound frameshift rescue in EGFR exon 19
    egfr = GENOME.transcript_by_id("ENST00000275493")
    ex19 = egfr.exons[18]
    cds_start = min(egfr.start_codon_spliced_offsets)
    s_off = cds_start + 745 * 3 + 6
    pos_map = {egfr.spliced_offset(p): p for p in range(ex19.start, ex19.end + 1)}
    p0 = pos_map[s_off]
    r1 = Variant("7", p0, egfr.sequence[s_off], "", genome=GENOME)
    r2 = Variant("7", p0 + 11, "", "C", genome=GENOME)
    resolver = MolecularPhaseResolver(_CisPair(r1, r2, answer=True))
    u_effs = VariantCollection([r1, r2]).effects(annotator="unified", phase_resolver=resolver)
    tm_effs = VariantCollection([r1, r2]).effects(annotator="transcript_model", phase_resolver=resolver)
    u_joint = [
        e for e in u_effs
        if getattr(getattr(e, "transcript", None), "id", None) == egfr.id
        and len(getattr(e, "variants", ())) > 1
    ][0]
    tm_joint = [
        e for e in tm_effs
        if getattr(getattr(e, "transcript", None), "id", None) == egfr.id
        and len(getattr(e, "variants", ())) > 1
    ][0]
    assert isinstance(u_joint, ComplexSubstitution)
    assert u_joint.mutant_protein_sequence == tm_joint.mutant_protein_sequence
    assert len(u_joint.mutant_protein_sequence) == 1210


def test_dna_search_anchor_for_big_intronic_insertion_and_template_materialization():
    """Big intronic insertion gives Isovar exact 5'/3' anchors, frame, junctions, and template."""
    met = GENOME.transcript_by_id("ENST00000397752")  # + strand
    ex13, ex14, ex15 = met.exons[12], met.exons[13], met.exons[14]
    inserted = "ACGT" * 500  # 2 kb novel intronic insertion between exon 14 and exon 15
    bp = ex14.end + 250
    sv_ins = StructuralVariant(
        "7", bp, "INS", end=bp, alt_assembly=inserted, genome=GENOME
    )

    eff = sv_ins.effect_on_transcript(met, annotator="unified")
    anchor = eff.dna_search_anchor
    assert isinstance(anchor, DNASearchAnchor)
    assert anchor.five_prime_transcript_id == met.id
    assert anchor.five_prime_strand == "+"
    assert anchor.last_intact_exon_number == 14
    assert anchor.last_intact_exon_span == (ex14.start, ex14.end)
    assert anchor.skip_fallback_exon_number == 13
    assert anchor.skip_fallback_exon_span == (ex13.start, ex13.end)
    assert anchor.breakpoint_in_intron is True
    assert anchor.last_intact_donor_base == ex14.end
    assert anchor.five_prime_cds_prefix is not None
    assert anchor.exit_frame == len(anchor.five_prime_cds_prefix) % 3
    assert anchor.three_prime_reentry_exon_number == 15
    assert anchor.three_prime_reentry_exon_span == (ex15.start, ex15.end)
    assert anchor.three_prime_reentry_acceptor_base == ex15.start
    assert anchor.three_prime_entry_frame == anchor.exit_frame
    assert anchor.inserted_sequence == inserted
    assert anchor.mate_search_region is None

    # Wild-type junction exon 14 -> exon 15 is present for distinguishing mutant splices.
    wt_introns = {(j.intron_start, j.intron_end) for j in anchor.wildtype_junctions}
    assert (ex14.end + 1, ex15.start - 1) in wt_introns

    # Materialize mutant pre-mRNA template with a mock genomic provider.
    def provider(contig, start, end):
        return "A" * (end - start + 1)

    template = anchor.materialize_template(provider)
    assert inserted in template
    labels = [seg[0] for seg in anchor.template_segments]
    assert labels == ["five_prime_host", "inserted_sequence", "three_prime_host"]


def test_dna_search_anchor_for_gene_to_intergenic_translocations_both_strands():
    """Gene->intergenic translocations give oriented mate search window on + and - strands."""
    # + strand: MET intron 14 -> intergenic chr6:44218908 moving '+'
    met = GENOME.transcript_by_id("ENST00000397752")
    ex14 = met.exons[13]
    sv_plus = parse_symbolic_alt(
        "7", ex14.end + 300, "N", "N[6:44218908[", genome=GENOME
    )
    eff_plus = sv_plus.effect_on_transcript(met, annotator="unified")
    assert isinstance(eff_plus, TranslocationToIntergenic)
    anchor_p = eff_plus.dna_search_anchor
    assert anchor_p.five_prime_strand == "+"
    assert anchor_p.last_intact_exon_number == 14
    assert anchor_p.three_prime_reentry_exon_number is None
    assert anchor_p.mate_search_region == ("6", 44218908, 44218908 + 25_000 - 1, "+")
    assert [seg[0] for seg in anchor_p.template_segments] == [
        "five_prime_host",
        "mate_search_window",
    ]

    # - strand: TP53 intron 5 -> intergenic chr10:40000000 moving '-'
    tp53 = GENOME.transcript_by_id("ENST00000269305")
    ex5 = tp53.exons[4]
    intergenic_pos = 40_000_000
    assert not GENOME.genes_at_locus("10", intergenic_pos)
    sv_minus = parse_symbolic_alt(
        "17", ex5.start - 40, "N", f"]10:{intergenic_pos}]N", genome=GENOME
    )
    eff_minus = sv_minus.effect_on_transcript(tp53, annotator="unified")
    assert isinstance(eff_minus, TranslocationToIntergenic)
    anchor_m = eff_minus.dna_search_anchor
    assert anchor_m.five_prime_strand == "-"
    assert anchor_m.last_intact_exon_number == 5
    assert anchor_m.skip_fallback_exon_number == 4
    assert anchor_m.breakpoint_in_intron is True
    assert anchor_m.last_intact_donor_base == ex5.start
    assert anchor_m.mate_search_region == (
        "10",
        intergenic_pos - 25_000 + 1,
        intergenic_pos,
        "-",
    )
    assert anchor_m.mate_search_region[1] <= anchor_m.mate_search_region[2]

    # Verify reverse-complement handling in materialize_template for '-' segments.
    def strand_provider(contig, start, end):
        return ("A" * (end - start)) + "C"

    tmpl = anchor_m.materialize_template(strand_provider)
    assert tmpl.startswith("G")  # RC of trailing 'C' on '-' five_prime_host

    # 5'UTR-only TMPRSS2--ERG honestly reports no 5' CDS prefix / exit frame.
    tmprss2 = GENOME.transcript_by_id("ENST00000332149")
    erg = GENOME.transcript_by_id("ENST00000288319")
    sv_te = parse_symbolic_alt(
        "21", tmprss2.exons[0].start - 500, "N", f"]21:{erg.exons[3].end + 500}]N", genome=GENOME
    )
    anchor_te = build_dna_search_anchor(sv_te, tmprss2)
    assert anchor_te.last_intact_exon_number == 1
    assert anchor_te.five_prime_cds_prefix is None
    assert anchor_te.exit_frame is None


def test_fast_path_safety_guards_and_reference_mismatch():
    """Fast-path must never swallow reference mismatches, splice sites, start/stop, mito, or Sec."""
    braf = GENOME.transcript_by_id("ENST00000288602")  # - strand
    # Reference mismatch at BRAF V600 (true ref on + strand is A, pass C).
    bad = Variant("7", 140753336, "C", "T", genome=GENOME)
    with pytest.raises(ReferenceMismatchError):
        bad.effect_on_transcript(braf, annotator="unified")

    # Start codon SNV must still classify as StartLoss / AlternateStartCodon (not fast-path).
    start_pos = max(braf.start_codon_positions)  # - strand first base of ATG (ref T on +)
    v_start = Variant("7", start_pos, "T", "G", genome=GENOME)
    u_start = v_start.effect_on_transcript(braf, annotator="unified")
    tm_start = v_start.effect_on_transcript(braf, annotator="transcript_model")
    assert type(u_start) is type(tm_start)
    assert u_start.short_description == tm_start.short_description

    # Mitochondrial variant (non-standard genetic code) must match transcript_model.
    mt_v = Variant("MT", 12994, "G", "A", genome=GENOME)
    mt_u = mt_v.effects(annotator="unified")
    mt_tm = mt_v.effects(annotator="transcript_model")
    assert [(type(e).__name__, e.short_description, e.mutant_protein_sequence) for e in mt_u] == [
        (type(e).__name__, e.short_description, e.mutant_protein_sequence) for e in mt_tm
    ]

    # Selenocysteine transcript (GPX4) must match transcript_model.
    gpx4 = GENOME.transcript_by_id("ENST00000354171")
    ex = gpx4.exons[2]
    mid = (ex.start + ex.end) // 2
    ref_b = gpx4.sequence[gpx4.spliced_offset(mid)]
    alt_b = "A" if ref_b != "A" else "G"
    v_sec = Variant(gpx4.contig, mid, ref_b, alt_b, genome=GENOME)
    u_sec = v_sec.effect_on_transcript(gpx4, annotator="unified")
    tm_sec = v_sec.effect_on_transcript(gpx4, annotator="transcript_model")
    assert type(u_sec) is type(tm_sec)
    assert u_sec.mutant_protein_sequence == tm_sec.mutant_protein_sequence


def test_unified_json_roundtrip_preserves_hypotheses_junctions_and_anchors():
    """Unified effects with DNAHypothesis, GenomicSpliceJunction, and DNASearchAnchor round-trip."""
    met = GENOME.transcript_by_id("ENST00000397752")
    ex14 = met.exons[13]
    e_splice = Variant("7", ex14.end + 1, "G", "A", genome=GENOME).effect_on_transcript(
        met, annotator="unified"
    )
    sv_trans = parse_symbolic_alt(
        "7", ex14.end + 300, "N", "N[6:44218908[", genome=GENOME
    )
    e_trans = sv_trans.effect_on_transcript(met, annotator="unified")

    col = EffectCollection([e_splice, e_trans])
    restored = EffectCollection.from_json(col.to_json())
    by_var = {e.variant: e for e in restored}

    r_splice = by_var[e_splice.variant]
    assert all(isinstance(c, DNAHypothesis) for c in r_splice.candidates)
    r_skip = next(c for c in r_splice.candidates if "exon_skip" in c.splice_mechanisms)
    assert len(r_skip.novel_junctions) == 1
    assert isinstance(r_skip.novel_junctions[0], GenomicSpliceJunction)
    assert r_skip.novel_junctions[0].intron_start <= r_skip.novel_junctions[0].intron_end

    r_trans = by_var[e_trans.variant]
    assert isinstance(r_trans.dna_search_anchor, DNASearchAnchor)
    assert r_trans.dna_search_anchor == e_trans.dna_search_anchor


def test_unified_coding_snv_batch_parity_and_speedup_over_transcript_model():
    """On interior exonic SNVs, unified matches transcript_model byte-for-byte and runs faster."""
    met = GENOME.transcript_by_id("ENST00000397752")
    # Build 60 real interior coding SNVs across MET exon 2 (large coding exon, + strand).
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

    # Warmup
    variants[0].effect_on_transcript(met, annotator="transcript_model")
    variants[0].effect_on_transcript(met, annotator="unified")

    t0 = time.perf_counter()
    tm_effects = [v.effect_on_transcript(met, annotator="transcript_model") for v in variants]
    t_tm = time.perf_counter() - t0

    t1 = time.perf_counter()
    u_effects = [v.effect_on_transcript(met, annotator="unified") for v in variants]
    t_u = time.perf_counter() - t1

    for u_eff, tm_eff in zip(u_effects, tm_effects):
        assert type(u_eff) is type(tm_eff)
        assert u_eff.short_description == tm_eff.short_description
        assert u_eff.mutant_protein_sequence == tm_eff.mutant_protein_sequence

    assert t_u < t_tm, f"Expected unified ({t_u:.4f}s) faster than transcript_model ({t_tm:.4f}s)"


def test_unified_osteosarc_full_cohort_parity_with_transcript_model():
    """All 179 real clinical osteosarc variants match transcript_model byte-for-byte."""
    from .osteosarc_variants import load_variants

    variants = list(load_variants(genome=GENOME))
    assert len(variants) == 179
    for variant in variants:
        tm_effects = list(variant.effects(annotator="transcript_model", raise_on_error=True))
        u_effects = list(variant.effects(annotator="unified", raise_on_error=True))
        assert len(u_effects) == len(tm_effects), f"Effect count mismatch on {variant}"
        for u_eff, tm_eff in zip(u_effects, tm_effects):
            assert getattr(u_eff.transcript, "id", None) == getattr(tm_eff.transcript, "id", None)
            assert type(u_eff) is type(tm_eff), f"Class mismatch on {variant}"
            assert u_eff.short_description == tm_eff.short_description, (
                f"short_description mismatch on {variant}: {u_eff} vs {tm_eff}"
            )
            assert u_eff.mutant_protein_sequence == tm_eff.mutant_protein_sequence, (
                f"mutant_protein_sequence mismatch on {variant} ({type(u_eff).__name__})"
            )
            assert u_eff.modifies_protein_sequence == tm_eff.modifies_protein_sequence


def test_unified_interior_mnv_and_indel_fast_path_parity_both_strands():
    """Interior coding MNVs, in-frame indels, frameshifts, and Silent match transcript_model."""
    met = GENOME.transcript_by_id("ENST00000397752")  # + strand
    brca1 = GENOME.transcript_by_id("ENST00000357654")  # - strand
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
            # SNV, MNV, in-frame del, in-frame ins, pure-ins stop, frameshift del, frameshift ins
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
        u_eff = variant.effect_on_transcript(tx, annotator="unified")
        assert type(u_eff) is type(tm_eff), f"Class mismatch on {variant} ({tx.id})"
        assert u_eff.short_description == tm_eff.short_description, (
            f"short_description mismatch on {variant} ({tx.id}): {u_eff} vs {tm_eff}"
        )
        assert u_eff.mutant_protein_sequence == tm_eff.mutant_protein_sequence
        assert u_eff.modifies_protein_sequence == tm_eff.modifies_protein_sequence

    # Reference mismatch on a coding deletion must still raise ReferenceMismatchError.
    ex2 = met.exons[1]
    pos = ex2.start + 60
    real_seq = met.sequence[met.spliced_offset(pos) : met.spliced_offset(pos) + 4]
    wrong4 = "AAAA" if real_seq != "AAAA" else "CCCC"
    bad_del = Variant("7", pos, wrong4, wrong4[0], genome=GENOME)
    with pytest.raises(ReferenceMismatchError):
        bad_del.effect_on_transcript(met, annotator="unified")


def test_unified_interoperates_with_priority_proteome_orf_discovery_and_mt_polya():
    """Unified effects work with effect_priority, to_proteome, discover_created_orfs, and MT poly-A."""
    import json
    from pathlib import Path

    from varcode import (
        StructuralVariant,
        classification_priority,
        discover_created_orfs,
        effect_priority,
        to_proteome,
    )

    # 1. effect_priority & filter_by_effect_priority see coding alternatives on unified splice effects.
    cftr = GENOME.transcript_by_id("ENST00000003084")
    weak_donor = Variant("7", 117531117, "A", "C", genome=GENOME)
    u_splice = weak_donor.effect_on_transcript(cftr, annotator="unified")
    del_alt = next(c.effect for c in u_splice.candidates if isinstance(c.effect, Deletion))
    assert effect_priority(u_splice) >= classification_priority(del_alt)
    assert list(EffectCollection([u_splice]).filter_by_effect_priority(Deletion)) == [u_splice]

    # 2. to_proteome exports unified coding + splice + germline effects with tool='varcode'.
    braf = GENOME.transcript_by_id("ENST00000288602")
    u_v600e = Variant("7", 140753336, "A", "T", genome=GENOME).effect_on_transcript(
        braf, annotator="unified"
    )
    germline_ctx = GermlineContext.from_variants([Variant("7", 140753335, "C", "T", genome=GENOME)])
    u_phased = Variant("7", 140753336, "A", "T", genome=GENOME).effect_on_transcript(
        braf, germline=germline_ctx, annotator="unified"
    )
    proteome = to_proteome([u_v600e, u_splice, u_phased])
    assert proteome.proteins
    assert all(p.producer.tool == "varcode" for p in proteome.proteins)

    # 3. discover_created_orfs works on unified 5'UTR start-gain + unphased germline.
    fixture = json.loads((Path(__file__).parent / "data/created_orfs/CDKN2A.json").read_text())
    native = fixture["native"]
    g95 = cached_release(native["annotation_version"])
    cdkn2a_tx = g95.transcript_by_id(native["transcript_id"])
    allele = fixture["genomic"]
    v_start = Variant(allele["contig"], allele["start"], allele["ref"], allele["alt"], genome=g95)
    u_utr = v_start.effect_on_transcript(cdkn2a_tx, annotator="unified")
    report = discover_created_orfs(u_utr)
    assert report.complete
    assert len(report.effects) == 1
    assert report.effects[0].mutant_protein_sequence == fixture["refseq"]["expected_protein"]

    # 4. Mitochondrial transcript with poly-A stop codon (MT-ND1) builds DNASearchAnchor cleanly.
    nd1 = GENOME.transcript_by_id("ENST00000361390")
    mt_dup = StructuralVariant("MT", 3500, "DUP", end=3560, genome=GENOME)
    u_mt = mt_dup.effect_on_transcript(nd1, annotator="unified")
    assert isinstance(u_mt.dna_search_anchor, DNASearchAnchor)
    assert u_mt.dna_search_anchor.five_prime_cds_prefix is not None
    assert u_mt.dna_search_anchor.exit_frame in (0, 1, 2)
