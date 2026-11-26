"""Canonical cancer driver, fusion oncogene, and adversarial hard-case suite.

This Phase-0 golden suite pins independent biological expectations across every
major somatic variant class and scale (SNVs, MNVs, in-frame & frameshift indels,
intra-exonic tandem duplications, compound frameshift rescue, splice-altering
drivers on both strands, and fusion oncogenes across translocation, inversion,
deletion, and tandem-duplication topologies) so future annotator simplification
and fast-path work can be verified against both ``transcript_model`` and ``fast``
without silent regressions.

All expected mutant proteins are translated independently from Ensembl 81
(GRCh38) transcript nucleotide sequences using NCBI genetic code 1 rather than
comparing one Varcode path against another.
"""

from __future__ import annotations

import pytest
from pyensembl import cached_release

from varcode import (
    EffectCollection,
    GermlineContext,
    MolecularPhaseResolver,
    StructuralVariant,
    Variant,
    VariantCollection,
    iter_candidate_effects,
    parse_symbolic_alt,
)
from varcode.effects.effect_classes import (
    ComplexSubstitution,
    Deletion,
    FrameShift,
    GeneFusion,
    HaplotypeEffect,
    Insertion,
    Substitution,
    Unresolved,
)
from varcode.nucleotides import reverse_complement


GENOME = cached_release(81)

# Independent NCBI standard genetic code (table 1) so expected proteins do not
# depend on Varcode's own translator.
_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def _translate_independent(cds: str) -> str:
    """Translate a coding sequence to the first stop codon independently."""
    aa = []
    for i in range(0, len(cds) - 2, 3):
        residue = _CODON_TABLE[cds[i:i + 3].upper()]
        if residue == "*":
            break
        aa.append(residue)
    return "".join(aa)


def _fusion_protein(effect):
    mt = getattr(effect, "mutant_transcript", None)
    if mt is not None and mt.mutant_protein_sequence is not None:
        return mt.mutant_protein_sequence
    return getattr(effect, "mutant_protein_sequence", None)


class _StaticPairPhaseSource:
    """Minimal ReadPhasingSource for cis or trans pairs."""

    source = "test_dna_phase"

    def __init__(self, variants, *, in_cis_answer: bool | None):
        self._variants = tuple(variants)
        self._in_cis_answer = in_cis_answer

    def has_evidence(self, variant) -> bool:
        return variant in self._variants

    def partners_in_cis(self, variant):
        if self._in_cis_answer is not True or variant not in self._variants:
            return ()
        return tuple(v for v in self._variants if v != variant)

    def in_cis(self, v1, v2, transcript=None):
        if {v1, v2} == set(self._variants):
            return self._in_cis_answer
        return None


# ---------------------------------------------------------------------------
# 1. Fusion oncogenes across all 4 SV / strand topologies
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("annotator", ["transcript_model", "fast"])
def test_eml4_alk_inversion_fusion_v1_e13a20(annotator):
    """inv(2)(p21p23) (+ / -): EML4 exon 13 joined to ALK exon 20 (1,059 aa)."""
    eml4 = GENOME.transcript_by_id("ENST00000318522")  # + strand
    alk = GENOME.transcript_by_id("ENST00000389048")  # - strand

    e13_cds_end = (
        eml4.spliced_offset(eml4.exons[12].end)
        - min(eml4.start_codon_spliced_offsets)
        + 1
    )
    a20_cds_start = (
        alk.spliced_offset(alk.exons[19].end)
        - min(alk.start_codon_spliced_offsets)
    )
    expected_protein = _translate_independent(
        eml4.coding_sequence[:e13_cds_end] + alk.coding_sequence[a20_cds_start:]
    )
    assert len(expected_protein) == 1059
    assert "PGKGPKVYRRKHQE" in expected_protein

    # Keep LEFT of EML4 intron 13 (5' on +) joined to LEFT of ALK intron 19 (3' on -).
    bp1 = eml4.exons[12].end + 500
    bp2 = alk.exons[19].end + 500
    sv = parse_symbolic_alt("2", bp1, "N", f"N]2:{bp2}]", genome=GENOME)

    effect = sv.effect_on_transcript(eml4, annotator=annotator)
    assert isinstance(effect, GeneFusion)
    assert effect.five_prime_transcript.gene_name == "EML4"
    assert effect.three_prime_transcript.gene_name == "ALK"
    assert _fusion_protein(effect) == expected_protein

    restored = EffectCollection.from_json(EffectCollection([effect]).to_json())
    assert len(restored) == 1
    assert _fusion_protein(restored[0]) == expected_protein


@pytest.mark.parametrize("annotator", ["transcript_model", "fast"])
@pytest.mark.parametrize("representation", ["DUP", "BND"])
def test_kiaa1549_braf_tandem_duplication_fusion(annotator, representation):
    """dup(7)(q34) (- / -): KIAA1549 exon 16 joined to BRAF exon 9 (2,135 aa)."""
    kiaa = GENOME.transcript_by_id("ENST00000440172")  # - strand
    braf = GENOME.transcript_by_id("ENST00000288602")  # - strand

    kiaa_cds_end = (
        kiaa.spliced_offset(kiaa.exons[15].start)
        - min(kiaa.start_codon_spliced_offsets)
        + 1
    )
    braf_cds_start = (
        braf.spliced_offset(braf.exons[8].end)
        - min(braf.start_codon_spliced_offsets)
    )
    expected_protein = _translate_independent(
        kiaa.coding_sequence[:kiaa_cds_end] + braf.coding_sequence[braf_cds_start:]
    )
    assert len(expected_protein) == 2135
    assert "ANNPCSDLIRDQGF" in expected_protein

    kiaa_bp = kiaa.exons[15].start - 500
    braf_bp = braf.exons[8].end + 500
    if representation == "DUP":
        sv = StructuralVariant("7", kiaa_bp, "DUP", end=braf_bp, genome=GENOME)
    else:
        sv = parse_symbolic_alt("7", kiaa_bp, "N", f"]7:{braf_bp}]N", genome=GENOME)

    effect = sv.effect_on_transcript(kiaa, annotator=annotator)
    assert isinstance(effect, GeneFusion)
    assert effect.five_prime_transcript.gene_name == "KIAA1549"
    assert effect.three_prime_transcript.gene_name == "BRAF"
    assert _fusion_protein(effect) == expected_protein

    restored = EffectCollection.from_json(EffectCollection([effect]).to_json())
    assert _fusion_protein(restored[0]) == expected_protein


@pytest.mark.parametrize("annotator", ["transcript_model", "fast"])
@pytest.mark.parametrize("representation", ["DEL", "BND"])
def test_tmprss2_erg_5utr_deletion_fusion_does_not_fabricate_start(
    annotator, representation
):
    """del(21)(q22.2q22.3) (- / -): TMPRSS2 5'UTR exon 1 joined to ERG exon 4.

    TMPRSS2's annotated start codon is in exon 2, so the fusion retains the
    full spliced cDNA (4,507 nt) while leaving mutant_protein_sequence as None
    rather than guessing an unannotated internal ERG initiation site.
    """
    tmprss2 = GENOME.transcript_by_id("ENST00000332149")  # - strand
    erg = GENOME.transcript_by_id("ENST00000288319")  # - strand

    expected_cdna = (
        tmprss2.sequence[:tmprss2.spliced_offset(tmprss2.exons[0].start) + 1]
        + erg.sequence[erg.spliced_offset(erg.exons[3].end):]
    )
    assert len(expected_cdna) == 4507

    tmprss2_bp = tmprss2.exons[0].start - 500
    erg_bp = erg.exons[3].end + 500
    if representation == "DEL":
        sv = StructuralVariant("21", erg_bp, "DEL", end=tmprss2_bp, genome=GENOME)
    else:
        sv = parse_symbolic_alt("21", tmprss2_bp, "N", f"]21:{erg_bp}]N", genome=GENOME)

    effect = sv.effect_on_transcript(tmprss2, annotator=annotator)
    assert isinstance(effect, GeneFusion)
    assert effect.five_prime_transcript.gene_name == "TMPRSS2"
    assert effect.three_prime_transcript.gene_name == "ERG"
    assert effect.mutant_transcript is not None
    assert effect.mutant_transcript.cdna_sequence == expected_cdna
    assert _fusion_protein(effect) is None


@pytest.mark.parametrize("annotator", ["transcript_model", "fast"])
def test_bcr_abl1_p190_intronic_vs_exonic_breakpoint_insertions(annotator):
    """t(9;22) (+ / +) e1a2 p190: intronic insertions splice out; exonic ones stay."""
    bcr = GENOME.transcript_by_id("ENST00000305877")  # + strand
    abl1 = GENOME.transcript_by_id("ENST00000318560")  # + strand

    e1_cds_end = (
        bcr.spliced_offset(bcr.exons[0].end)
        - min(bcr.start_codon_spliced_offsets)
        + 1
    )
    a2_cds_start = (
        abl1.spliced_offset(abl1.exons[1].start)
        - min(abl1.start_codon_spliced_offsets)
    )
    expected_p190 = _translate_independent(
        bcr.coding_sequence[:e1_cds_end] + abl1.coding_sequence[a2_cds_start:]
    )
    assert len(expected_p190) == 1530
    assert "AFHGDAEALQRPVASDFEPQ" in expected_p190

    # Intronic breakpoints: +3bp and +1bp breakpoint insertions are spliced out.
    bcr_intron_bp = bcr.exons[0].end + 500
    abl_intron_bp = abl1.exons[1].start - 500
    for alt in (
        f"N[9:{abl_intron_bp}[",
        f"NACG[9:{abl_intron_bp}[",
        f"NA[9:{abl_intron_bp}[",
    ):
        sv = parse_symbolic_alt("22", bcr_intron_bp, "N", alt, genome=GENOME)
        eff = sv.effect_on_transcript(bcr, annotator=annotator)
        assert isinstance(eff, GeneFusion)
        assert _fusion_protein(eff) == expected_p190

    # Exonic breakpoints: +3bp inserts 1 aa, +1bp shifts reading frame.
    bcr_ex_bp = bcr.exons[0].end - 10
    abl_ex_bp = abl1.exons[1].start + 10
    clean = parse_symbolic_alt(
        "22", bcr_ex_bp, "N", f"N[9:{abl_ex_bp}[", genome=GENOME
    ).effect_on_transcript(bcr, annotator=annotator)
    ins3 = parse_symbolic_alt(
        "22", bcr_ex_bp, "N", f"NGCG[9:{abl_ex_bp}[", genome=GENOME
    ).effect_on_transcript(bcr, annotator=annotator)
    ins1 = parse_symbolic_alt(
        "22", bcr_ex_bp, "N", f"NG[9:{abl_ex_bp}[", genome=GENOME
    ).effect_on_transcript(bcr, annotator=annotator)

    assert len(ins3.mutant_transcript.cdna_sequence) == (
        len(clean.mutant_transcript.cdna_sequence) + 3
    )
    assert len(ins1.mutant_transcript.cdna_sequence) == (
        len(clean.mutant_transcript.cdna_sequence) + 1
    )
    assert len(_fusion_protein(ins3)) == len(_fusion_protein(clean)) + 1
    assert _fusion_protein(ins1) != _fusion_protein(clean)


# ---------------------------------------------------------------------------
# 2. Canonical splice-altering cancer drivers (+ and - strand)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("annotator", ["transcript_model", "fast"])
def test_met_exon14_skipping_donor_acceptor_and_exonic_last_base(annotator):
    """MET exon 14 (141 nt = 47 aa in-frame skip) across donor, acceptor, and c.3028G>A."""
    met = GENOME.transcript_by_id("ENST00000397752")  # + strand, 1390 aa
    ex14 = met.exons[13]
    cds_start = min(met.start_codon_spliced_offsets)
    ex14_cds_start = met.spliced_offset(ex14.start) - cds_start
    ex14_cds_end = met.spliced_offset(ex14.end) - cds_start + 1
    expected_skip_protein = _translate_independent(
        met.coding_sequence[:ex14_cds_start] + met.coding_sequence[ex14_cds_end:]
    )
    assert len(expected_skip_protein) == len(met.protein_sequence) - 47 == 1343

    # Canonical donor +1 and acceptor -1 both realize the 47-aa exon-14 deletion.
    for variant in (
        Variant("7", ex14.end + 1, "G", "A", genome=GENOME),
        Variant("7", ex14.start - 1, "G", "A", genome=GENOME),
    ):
        effect = variant.effect_on_transcript(met, annotator=annotator)
        terminal = [e for e in iter_candidate_effects(effect) if e is not None]
        proteins = {e.mutant_protein_sequence for e in terminal}
        assert expected_skip_protein in proteins
        if annotator == "transcript_model":
            # Without chromosome FASTA, intron retention & cryptic stay explicit Unresolved.
            unresolved = {
                e.short_description for e in terminal if isinstance(e, Unresolved)
            }
            assert "unresolved-intron-retention" in unresolved
            assert "unresolved-cryptic" in unresolved

    # Exonic last base c.3028G>A must retain BOTH normal-splicing p.D1010N AND exon-14 skip.
    exonic_donor = Variant("7", ex14.end, "G", "A", genome=GENOME)
    effect = exonic_donor.effect_on_transcript(met, annotator=annotator)
    terminal = [e for e in iter_candidate_effects(effect) if e is not None]
    proteins = {e.mutant_protein_sequence for e in terminal}
    mutated_cds = (
        met.coding_sequence[:ex14_cds_end - 1]
        + "A"
        + met.coding_sequence[ex14_cds_end:]
    )
    expected_d1010n = _translate_independent(mutated_cds)
    assert expected_d1010n in proteins
    assert expected_skip_protein in proteins

    restored = EffectCollection.from_json(EffectCollection([effect]).to_json())
    restored_proteins = {
        e.mutant_protein_sequence
        for e in iter_candidate_effects(restored[0])
        if e is not None
    }
    assert expected_d1010n in restored_proteins
    assert expected_skip_protein in restored_proteins


@pytest.mark.parametrize("annotator", ["transcript_model", "fast"])
def test_tp53_minus_strand_splice_donor_skips_exon5(annotator):
    """TP53 (- strand) exon 5 donor +1 yields exon-5 skip frameshift p.Y126fs."""
    tp53 = GENOME.transcript_by_id("ENST00000269305")  # - strand
    ex4, ex5, ex6 = tp53.exons[3], tp53.exons[4], tp53.exons[5]
    cds_start = min(tp53.start_codon_spliced_offsets)
    ex5_start = tp53.spliced_offset(ex5.end) - cds_start
    ex5_end = tp53.spliced_offset(ex5.start) - cds_start + 1
    expected_skip_protein = _translate_independent(
        tp53.coding_sequence[:ex5_start] + tp53.coding_sequence[ex5_end:]
    )

    # On - strand, intron 5 donor +1 is at ex5.start - 1 (ref C on + strand).
    donor = Variant("17", ex5.start - 1, "C", "T", genome=GENOME)
    effect = donor.effect_on_transcript(tp53, annotator=annotator)
    terminal = [e for e in iter_candidate_effects(effect) if e is not None]
    skip_effects = [
        e for e in terminal if e.mutant_protein_sequence == expected_skip_protein
    ]
    assert skip_effects
    if annotator == "transcript_model":
        assert isinstance(effect, FrameShift)
        assert effect.short_description == "p.Y126fs"
        # Verify exon 4 -> exon 6 non-reference junction coordinates on - strand.
        skip_cand = effect.realized_candidates[0]
        ((contig1, left_pos, contig2, right_pos),) = (
            skip_cand.outcomes[0].mutant.junction_signature
        )
        assert (contig1, contig2) == ("17", "17")
        assert (left_pos, right_pos) == (ex4.start, ex6.end)
        # Ascending genomic intron interval matches BAM CIGAR N bounds.
        assert (min(left_pos, right_pos) + 1, max(left_pos, right_pos) - 1) == (
            ex6.end + 1,
            ex4.start - 1,
        )


# ---------------------------------------------------------------------------
# 3. Complex coding drivers, MNVs, ITDs, and compound frameshift rescue
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("annotator", ["transcript_model", "fast"])
def test_braf_v600k_mnv_matches_cis_snvs_and_diverges_from_trans(annotator):
    """BRAF (- strand) V600K MNV == cis(V600E + V600M), whereas trans keeps them separate."""
    braf = GENOME.transcript_by_id("ENST00000288602")  # - strand
    v_mnv = Variant("7", 140753336, "AC", "TT", genome=GENOME)
    v_e = Variant("7", 140753336, "A", "T", genome=GENOME)  # p.V600E
    v_m = Variant("7", 140753337, "C", "T", genome=GENOME)  # p.V600M

    expected_v600k = braf.protein_sequence[:599] + "K" + braf.protein_sequence[600:]
    expected_v600e = braf.protein_sequence[:599] + "E" + braf.protein_sequence[600:]
    expected_v600m = braf.protein_sequence[:599] + "M" + braf.protein_sequence[600:]

    mnv_eff = v_mnv.effect_on_transcript(braf, annotator=annotator)
    assert isinstance(mnv_eff, Substitution)
    assert mnv_eff.short_description == "p.V600K"
    assert mnv_eff.mutant_protein_sequence == expected_v600k

    # Cis SNVs produce the exact same V600K protein alongside the individual effects.
    cis_resolver = MolecularPhaseResolver(
        _StaticPairPhaseSource((v_e, v_m), in_cis_answer=True)
    )
    cis_effects = VariantCollection([v_e, v_m]).effects(
        annotator=annotator, phase_resolver=cis_resolver
    )
    joint = [
        e
        for e in cis_effects
        if getattr(getattr(e, "transcript", None), "id", None) == braf.id
        and len(getattr(e, "variants", ())) > 1
    ]
    assert len(joint) == 1
    assert isinstance(
        joint[0], Substitution if annotator == "transcript_model" else HaplotypeEffect
    )
    assert joint[0].mutant_protein_sequence == expected_v600k

    # Trans SNVs emit no joint effect and keep V600E and V600M distinct.
    trans_resolver = MolecularPhaseResolver(
        _StaticPairPhaseSource((v_e, v_m), in_cis_answer=False)
    )
    trans_effects = VariantCollection([v_e, v_m]).effects(
        annotator=annotator, phase_resolver=trans_resolver
    )
    assert not any(len(getattr(e, "variants", ())) > 1 for e in trans_effects)
    braf_indiv = {
        e.variant: e.mutant_protein_sequence
        for e in trans_effects
        if getattr(getattr(e, "transcript", None), "id", None) == braf.id
    }
    assert braf_indiv == {v_e: expected_v600e, v_m: expected_v600m}


@pytest.mark.parametrize("annotator", ["transcript_model", "fast"])
def test_egfr_exon19_deletion_and_npm1_type_a_frameshift(annotator):
    """EGFR p.E746_A750del (in-frame) and NPM1 type A (+4 bp TCTG C-terminal neo-tail)."""
    egfr = GENOME.transcript_by_id("ENST00000275493")  # + strand
    ex19 = egfr.exons[18]
    cds_start = min(egfr.start_codon_spliced_offsets)
    s_off = cds_start + 745 * 3
    pos_map = {egfr.spliced_offset(p): p for p in range(ex19.start, ex19.end + 1)}
    v_egfr = Variant("7", pos_map[s_off], egfr.sequence[s_off:s_off + 15], "", genome=GENOME)
    expected_egfr = egfr.protein_sequence[:745] + egfr.protein_sequence[750:]

    egfr_eff = v_egfr.effect_on_transcript(egfr, annotator=annotator)
    assert isinstance(egfr_eff, Deletion)
    assert egfr_eff.short_description == "p.ELREA746del"
    assert egfr_eff.mutant_protein_sequence == expected_egfr

    npm1 = GENOME.transcript_by_id("ENST00000296930")  # + strand
    ex11 = npm1.exons[10]
    npm1_cds_start = min(npm1.start_codon_spliced_offsets)
    npm1_pos_map = {npm1.spliced_offset(p): p for p in range(ex11.start, ex11.end + 1)}
    ins_off = npm1_cds_start + 863
    v_npm1 = Variant("5", npm1_pos_map[ins_off - 1], "", "TCTG", genome=GENOME)
    mutant_npm1_rna = npm1.sequence[:ins_off] + "TCTG" + npm1.sequence[ins_off:]
    expected_npm1 = _translate_independent(mutant_npm1_rna[npm1_cds_start:])
    assert expected_npm1.endswith("QDLCLAVEEVSLRK")

    npm1_eff = v_npm1.effect_on_transcript(npm1, annotator=annotator)
    assert isinstance(npm1_eff, FrameShift)
    assert npm1_eff.short_description == "p.W288fs"
    assert npm1_eff.mutant_protein_sequence == expected_npm1

    restored = EffectCollection.from_json(
        EffectCollection([egfr_eff, npm1_eff]).to_json()
    )
    assert {e.variant: e.mutant_protein_sequence for e in restored} == {
        v_egfr: expected_egfr,
        v_npm1: expected_npm1,
    }


@pytest.mark.parametrize("annotator", ["transcript_model", "fast"])
def test_flt3_itd_intra_exonic_dup_matches_explicit_insertion(annotator):
    """Regression for #769: intra-exonic FLT3-ITD DUP matches explicit 21-bp insertion."""
    flt3 = GENOME.transcript_by_id("ENST00000241453")  # - strand, 993 aa
    dup_start, dup_end = 28034110, 28034130
    s_lo = flt3.spliced_offset(dup_end)
    s_hi = flt3.spliced_offset(dup_start)
    dup_transcript_seq = flt3.sequence[s_lo:s_hi + 1]
    dup_plus = reverse_complement(dup_transcript_seq)

    cds_start = min(flt3.start_codon_spliced_offsets)
    expected_rna = (
        flt3.sequence[:s_lo]
        + dup_transcript_seq
        + flt3.sequence[s_lo:]
    )
    expected_protein = _translate_independent(expected_rna[cds_start:])
    assert len(expected_protein) == len(flt3.protein_sequence) + 7 == 1000

    v_ins = Variant("13", dup_end, "", dup_plus, genome=GENOME)
    v_dup = StructuralVariant("13", dup_start, "DUP", end=dup_end, genome=GENOME)

    eff_ins = v_ins.effect_on_transcript(flt3, annotator=annotator)
    eff_dup = v_dup.effect_on_transcript(flt3, annotator=annotator)

    assert eff_ins.mutant_protein_sequence == expected_protein
    assert _fusion_protein(eff_dup) == expected_protein
    if annotator == "transcript_model":
        assert isinstance(eff_dup, Insertion)
        assert eff_dup.short_description == eff_ins.short_description == "p.603insYEYDLKW"
        assert len(eff_dup.candidates) == 1

    restored = EffectCollection.from_json(EffectCollection([eff_dup]).to_json())
    assert _fusion_protein(restored[0]) == expected_protein


@pytest.mark.parametrize("annotator", ["dna", "transcript_model", "fast"])
def test_egfr_exon20_itd_and_kinase_domain_duplication_kdd(annotator):
    """EGFR exon 20 intra-exonic ITD (p.770insSVD) and EGFR-KDD (exons 18-25 tandem DUP)."""
    egfr = GENOME.transcript_by_id("ENST00000275493")  # + strand, 28 exons, 1210 aa
    cds_start = min(egfr.start_codon_spliced_offsets)

    # 1. EGFR exon 20 intra-exonic 9-bp ITD (S768_D770dup -> p.770insSVD, 1,213 aa)
    ex20 = egfr.exons[19]
    s_off = cds_start + 767 * 3  # S768-V769-D770 (0-based 767..769)
    pos_map = {egfr.spliced_offset(p): p for p in range(ex20.start, ex20.end + 1)}
    dup_start, dup_end = pos_map[s_off], pos_map[s_off + 8]
    dup_seq = egfr.sequence[s_off:s_off + 9]
    expected_ex20_protein = _translate_independent(
        (egfr.sequence[:s_off] + dup_seq + egfr.sequence[s_off:])[cds_start:]
    )
    assert len(expected_ex20_protein) == len(egfr.protein_sequence) + 3 == 1213

    eff_ins = Variant("7", dup_end, "", dup_seq, genome=GENOME).effect_on_transcript(
        egfr, annotator=annotator
    )
    eff_dup = StructuralVariant(
        "7", dup_start, "DUP", end=dup_end, genome=GENOME
    ).effect_on_transcript(egfr, annotator=annotator)
    assert eff_ins.mutant_protein_sequence == expected_ex20_protein
    assert _fusion_protein(eff_dup) == expected_ex20_protein

    # 2. EGFR-KDD: multi-exon tandem duplication of exons 18..25 (intron 17 -> intron 25)
    ex18, ex25 = egfr.exons[17], egfr.exons[24]
    ex18_cds = egfr.spliced_offset(ex18.start) - cds_start
    ex25_cds = egfr.spliced_offset(ex25.end) - cds_start + 1
    kdd_cds = (
        egfr.coding_sequence[:ex25_cds]
        + egfr.coding_sequence[ex18_cds:ex25_cds]
        + egfr.coding_sequence[ex25_cds:]
    )
    expected_kdd_protein = _translate_independent(kdd_cds)
    assert len(expected_kdd_protein) == len(egfr.protein_sequence) + 351 == 1561

    # Intron 25 is 379 bp wide (ex25.end=55201355, ex26.start=55201735); +150 sits inside intron 25.
    sv_kdd = StructuralVariant(
        "7", ex18.start - 500, "DUP", end=ex25.end + 150, genome=GENOME
    )
    eff_kdd = sv_kdd.effect_on_transcript(egfr, annotator=annotator)
    assert _fusion_protein(eff_kdd) == expected_kdd_protein
    if annotator == "dna":
        assert isinstance(eff_kdd, Insertion)
        assert len(eff_kdd.candidates) == 1
        # The single diagnostic non-reference join is the EGFR-KDD back-junction exon 25 -> exon 18!
        assert eff_kdd.mutant_transcript.novel_junctions == ()
        (kdd_back,) = eff_kdd.mutant_transcript.back_junctions
        assert kdd_back == ("7", ex18.start + 1, ex25.end - 1, "+")

    restored = EffectCollection.from_json(EffectCollection([eff_kdd]).to_json())
    assert _fusion_protein(restored[0]) == expected_kdd_protein


@pytest.mark.parametrize("annotator", ["transcript_model", "fast"])
def test_compound_frameshift_rescue_restores_reading_frame(annotator):
    """Two cis frameshifting indels (-1 bp + +1 bp) restore frame as ComplexSubstitution."""
    egfr = GENOME.transcript_by_id("ENST00000275493")  # + strand, 1210 aa
    ex19 = egfr.exons[18]
    cds_start = min(egfr.start_codon_spliced_offsets)
    s_off = cds_start + 745 * 3 + 6
    pos_map = {egfr.spliced_offset(p): p for p in range(ex19.start, ex19.end + 1)}
    p0 = pos_map[s_off]
    r1 = Variant("7", p0, egfr.sequence[s_off], "", genome=GENOME)  # p.R748fs
    r2 = Variant("7", p0 + 11, "", "C", genome=GENOME)  # p.S752fs

    # Independently apply both edits to the transcript sequence.
    mutated_rna = (
        egfr.sequence[:s_off]
        + egfr.sequence[s_off + 1:s_off + 12]
        + "C"
        + egfr.sequence[s_off + 12:]
    )
    expected_protein = _translate_independent(mutated_rna[cds_start:])
    assert len(expected_protein) == len(egfr.protein_sequence) == 1210

    # Individually both are frameshifts.
    assert isinstance(r1.effect_on_transcript(egfr, annotator=annotator), FrameShift)
    assert isinstance(r2.effect_on_transcript(egfr, annotator=annotator), FrameShift)

    # In cis they rescue the reading frame.
    resolver = MolecularPhaseResolver(
        _StaticPairPhaseSource((r1, r2), in_cis_answer=True)
    )
    effects = VariantCollection([r1, r2]).effects(
        annotator=annotator, phase_resolver=resolver
    )
    joint = [
        e
        for e in effects
        if getattr(getattr(e, "transcript", None), "id", None) == egfr.id
        and len(getattr(e, "variants", ())) > 1
    ]
    assert len(joint) == 1
    assert joint[0].mutant_protein_sequence == expected_protein
    if annotator == "transcript_model":
        assert isinstance(joint[0], ComplexSubstitution)
        assert joint[0].short_description == "p.REAT748EKQH"


# ---------------------------------------------------------------------------
# 4. Adversarial fast-path & baseline-integrity boundaries
# ---------------------------------------------------------------------------


def test_exonic_splice_boundary_snvs_never_drop_splice_candidates_on_either_strand():
    """Exonic donor (last 3 bases) and acceptor (first base) keep multi-outcome sets."""
    met = GENOME.transcript_by_id("ENST00000397752")  # + strand
    tp53 = GENOME.transcript_by_id("ENST00000269305")  # - strand

    # MET exon 14 acceptor first base (+ strand) and TP53 exon 5 donor last base (- strand)
    met_acc_pos = met.exons[13].start
    met_acc_ref = met.sequence[met.spliced_offset(met_acc_pos)]
    met_acc_alt = "A" if met_acc_ref != "A" else "C"
    v_met_acc = Variant("7", met_acc_pos, met_acc_ref, met_acc_alt, genome=GENOME)

    tp53_don_pos = tp53.exons[4].start
    tp53_don_ref = reverse_complement(tp53.sequence[tp53.spliced_offset(tp53_don_pos)])
    tp53_don_alt = "A" if tp53_don_ref != "A" else "G"
    v_tp53_don = Variant("17", tp53_don_pos, tp53_don_ref, tp53_don_alt, genome=GENOME)

    for variant, transcript in ((v_met_acc, met), (v_tp53_don, tp53)):
        eff = variant.effect_on_transcript(transcript, annotator="transcript_model")
        assert len(eff.candidates) > 1
        mechanisms = {
            option.mechanism
            for cand in eff.realized_candidates
            for outcome in cand.outcomes
            if outcome.hypothesis.splice_plan is not None
            for _, option in outcome.hypothesis.splice_plan.choices
        }
        assert "normal" in mechanisms
        assert "exon_skip" in mechanisms


def test_unphased_same_codon_germline_preserves_distinct_patient_baselines():
    """BRAF V600E with unphased germline V600M retains both cis (M600K) and trans (V600E)."""
    braf = GENOME.transcript_by_id("ENST00000288602")  # - strand
    somatic = Variant("7", 140753336, "A", "T", genome=GENOME)  # V600E on ref
    germline = Variant("7", 140753337, "C", "T", genome=GENOME)  # V600M on ref
    ctx = GermlineContext.from_variants([germline], reference_name="GRCh38")

    effect = somatic.effect_on_transcript(
        braf, annotator="transcript_model", germline=ctx
    )
    assert len(effect.realized_candidates) == 2
    pairs = {
        (
            cand.outcomes[0].baseline.protein_sequence[599],
            cand.outcomes[0].mutant.protein_sequence[599],
            cand.effect.short_description,
        )
        for cand in effect.realized_candidates
    }
    assert pairs == {
        ("V", "E", "p.V600E"),  # trans: baseline is reference V600
        ("M", "K", "p.M600K"),  # cis: baseline carries germline M600
    }
