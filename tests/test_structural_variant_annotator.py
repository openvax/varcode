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

"""Tests for :class:`StructuralVariantAnnotator` (PR 10; #252).

These cover the minimal dispatch shape — SV type + overlap pattern
routed to the right effect class — and the ``MultiOutcomeEffect``
harmonized interface (``effect.outcomes`` tuple, each carrying
``Outcome`` with ``source="varcode"``).

Downstream integrations (external splice predictors, RNA evidence,
long-read assembly) layer additional outcomes on top; those aren't
tested here because varcode doesn't ship the integrations themselves
— only the shape they plug into.
"""

import pytest
from pyensembl import cached_release

from varcode import (
    Outcome,
    StructuralVariant,
    StructuralVariantAnnotator,
    get_annotator,
)
from varcode.effects import (
    GeneFusion,
    Intronic,
    Inversion,
    LargeDeletion,
    LargeDuplication,
    MultiOutcomeEffect,
    TranslocationToIntergenic,
)

ensembl_grch38 = cached_release(81)
_ANNOTATOR = StructuralVariantAnnotator()

CFTR_ID = "ENST00000003084"   # chr7, + strand, 1480 aa
BRCA1_ID = "ENST00000357654"  # chr17, - strand


def _cftr():
    return ensembl_grch38.transcript_by_id(CFTR_ID)


# --------------------------------------------------------------------
# Registry wiring
# --------------------------------------------------------------------


def test_sv_annotator_registered_by_name():
    annotator = get_annotator("structural_variant")
    assert isinstance(annotator, StructuralVariantAnnotator)


def test_sv_annotator_supports_sv_types():
    assert "DEL" in _ANNOTATOR.supports
    assert "BND" in _ANNOTATOR.supports
    assert "INV" in _ANNOTATOR.supports
    assert "DUP" in _ANNOTATOR.supports
    # Doesn't claim SNV/indel/mnv — those are for fast/protein_diff.
    assert "snv" not in _ANNOTATOR.supports


# --------------------------------------------------------------------
# Deletion
# --------------------------------------------------------------------


def test_sv_large_deletion_spans_cftr_exons():
    transcript = _cftr()
    # Deletion spanning several CFTR exons (coordinates from the CFTR
    # transcript envelope).
    sv = StructuralVariant(
        contig="7",
        start=transcript.start + 100,
        end=transcript.start + 50_000,
        sv_type="DEL",
        alt="<DEL>",
        genome=ensembl_grch38,
    )
    effect = _ANNOTATOR.annotate_on_transcript(sv, transcript)
    assert isinstance(effect, LargeDeletion)
    assert isinstance(effect, MultiOutcomeEffect)
    assert len(effect.affected_exons) >= 1
    assert effect.short_description == "sv-deletion"


def test_sv_large_deletion_outcomes_shape():
    """The harmonized interface: ``effect.outcomes`` is a tuple of
    :class:`Outcome` entries with ``source="varcode"``."""
    transcript = _cftr()
    sv = StructuralVariant(
        contig="7",
        start=transcript.start + 100,
        end=transcript.start + 50_000,
        sv_type="DEL",
        genome=ensembl_grch38,
    )
    effect = _ANNOTATOR.annotate_on_transcript(sv, transcript)
    outcomes = effect.outcomes
    assert isinstance(outcomes, tuple)
    assert len(outcomes) >= 1
    for o in outcomes:
        assert isinstance(o, Outcome)
        assert o.source == "varcode"
        assert o.probability is None  # varcode doesn't assign; integrations do


def test_sv_deletion_purely_intronic_returns_intronic():
    """A small DEL that lands entirely inside a CFTR intron should
    not be classified as LargeDeletion (no exons overlapped)."""
    transcript = _cftr()
    # CFTR intron 1 is roughly 117480148–117530899 (between exon 1
    # end and exon 2 start). Pick a small interval well inside that.
    sv = StructuralVariant(
        contig="7",
        start=117_490_000,
        end=117_490_100,
        sv_type="DEL",
        genome=ensembl_grch38,
    )
    effect = _ANNOTATOR.annotate_on_transcript(sv, transcript)
    assert isinstance(effect, Intronic)


# --------------------------------------------------------------------
# Duplication / inversion
# --------------------------------------------------------------------


def test_sv_duplication_on_cftr_exons():
    transcript = _cftr()
    sv = StructuralVariant(
        contig="7",
        start=transcript.start + 100,
        end=transcript.start + 10_000,
        sv_type="DUP",
        genome=ensembl_grch38,
    )
    effect = _ANNOTATOR.annotate_on_transcript(sv, transcript)
    assert isinstance(effect, LargeDuplication)


def test_sv_inversion_on_cftr_exons():
    transcript = _cftr()
    sv = StructuralVariant(
        contig="7",
        start=transcript.start + 100,
        end=transcript.start + 10_000,
        sv_type="INV",
        genome=ensembl_grch38,
    )
    effect = _ANNOTATOR.annotate_on_transcript(sv, transcript)
    assert isinstance(effect, Inversion)


# --------------------------------------------------------------------
# Breakend → fusion / translocation
# --------------------------------------------------------------------


def test_sv_breakend_without_mate_is_translocation_intergenic():
    """BND with no mate info → TranslocationToIntergenic."""
    transcript = _cftr()
    sv = StructuralVariant(
        contig="7",
        start=transcript.start + 500,
        sv_type="BND",
        alt="N[?:?[",
        genome=ensembl_grch38,
    )
    effect = _ANNOTATOR.annotate_on_transcript(sv, transcript)
    assert isinstance(effect, TranslocationToIntergenic)


def test_sv_breakend_with_mate_in_different_gene_is_fusion():
    """BND whose mate lands on another protein-coding transcript is
    classified as GeneFusion."""
    cftr = _cftr()
    brca1 = ensembl_grch38.transcript_by_id(BRCA1_ID)
    # Pick an explicit position on BRCA1.
    mate_pos = brca1.start + 1_000
    sv = StructuralVariant(
        contig="7",
        start=cftr.start + 500,
        sv_type="BND",
        alt="N]17:%d]" % mate_pos,
        mate_contig="17",
        mate_start=mate_pos,
        mate_orientation="]]",
        genome=ensembl_grch38,
    )
    effect = _ANNOTATOR.annotate_on_transcript(sv, cftr)
    assert isinstance(effect, GeneFusion)
    assert effect.partner_transcript is not None
    # At least that it's *some* protein-coding transcript overlapping
    # the mate position. We don't pin the exact ID because Ensembl 81
    # may have multiple overlapping transcripts.
    assert effect.partner_transcript.is_protein_coding is True


def test_sv_breakend_mate_intergenic_is_translocation():
    """BND whose mate lands in intergenic space is
    TranslocationToIntergenic."""
    cftr = _cftr()
    # chr22:15,500,000 is a verified intergenic position in the
    # chr22 p-arm heterochromatin — no overlapping transcripts in
    # Ensembl 81.
    sv = StructuralVariant(
        contig="7",
        start=cftr.start + 500,
        sv_type="BND",
        alt="N]22:15500000]",
        mate_contig="22",
        mate_start=15_500_000,
        mate_orientation="]]",
        genome=ensembl_grch38,
    )
    effect = _ANNOTATOR.annotate_on_transcript(sv, cftr)
    assert isinstance(effect, TranslocationToIntergenic)


# --------------------------------------------------------------------
# Outcome extensibility
# --------------------------------------------------------------------


def test_external_integrations_can_attach_scored_outcomes():
    """Smoke test: an external tool wraps a varcode-nominated
    effect in an Outcome with source/probability/evidence, and
    that Outcome works as the interchange format PR 7 defined."""
    transcript = _cftr()
    sv = StructuralVariant(
        contig="7",
        start=transcript.start + 100,
        end=transcript.start + 5_000,
        sv_type="DEL",
        genome=ensembl_grch38,
    )
    effect = _ANNOTATOR.annotate_on_transcript(sv, transcript)

    # An RNA-evidence tool (Isovar-style) would score the varcode
    # outcome with junction read support.
    rna_scored = Outcome(
        effect=effect,
        probability=0.94,
        source="isovar",
        evidence={"junction_reads": 87, "split_reads": 15})
    assert rna_scored.source == "isovar"
    assert rna_scored.probability == 0.94
    assert rna_scored.evidence["junction_reads"] == 87

    # A long-read assembly tool would attach its own resolution.
    lr_scored = Outcome(
        effect=effect,
        probability=0.99,
        source="longread_assembly",
        evidence={"assembled_by": "hifiasm",
                  "breakpoint_confidence": "exact"})
    assert lr_scored.source == "longread_assembly"


def test_sv_annotator_honors_alt_assembly_hook():
    """Variants carrying a long-read-derived ``alt_assembly`` flow
    into :attr:`StructuralVariantEffect.mutant_transcript` as a
    single-segment MutantTranscript whose ``cdna_sequence`` equals
    the assembly, bypassing the inference-from-breakpoints path
    (#338)."""
    transcript = _cftr()
    assembled_allele = "ACGTACGTACGTACGT"
    sv = StructuralVariant(
        contig="7",
        start=transcript.start + 100,
        end=transcript.start + 5_000,
        sv_type="DEL",
        alt_assembly=assembled_allele,
        genome=ensembl_grch38,
    )
    effect = _ANNOTATOR.annotate_on_transcript(sv, transcript)
    assert isinstance(effect, LargeDeletion)
    mt = effect.mutant_transcript
    assert mt is not None
    assert mt.cdna_sequence == assembled_allele
    assert len(mt.reference_segments) == 1
    segment = mt.reference_segments[0]
    assert segment.label == "alt_assembly"
    assert segment.source.sequence == assembled_allele
    assert mt.evidence == {"source": "alt_assembly"}


def test_sv_annotator_attaches_cryptic_outcomes_from_assembly():
    """When the SV's ``alt_assembly`` contains plausible donor +
    acceptor motif pairs, :class:`CrypticExonCandidate` entries
    attach to the effect's ``outcomes`` as ``source="varcode_motif"``
    (#337). This exercises the wiring independently of whether the
    ambient genome exposes a FASTA API."""
    from varcode.effects.effect_classes import CrypticExonCandidate
    transcript = _cftr()
    # Synthetic assembly with an exact donor/acceptor consensus pair
    # flanking a plausible exon: 3' acceptor (CAG|G) ... 5' donor
    # (AAG|GTAAGT) with an intron-length spacer that satisfies the
    # default 30-10000 range.
    assembly = ("A" * 20 + "CAGG" + "C" * 100 + "AAGGTAAGT" + "A" * 20)
    sv = StructuralVariant(
        contig="7",
        start=transcript.start + 100,
        end=transcript.start + 5_000,
        sv_type="DEL",
        alt_assembly=assembly,
        genome=ensembl_grch38)
    effect = _ANNOTATOR.annotate_on_transcript(sv, transcript)
    outcomes = effect.outcomes
    # Primary outcome + at least one cryptic candidate.
    assert len(outcomes) >= 2
    cryptic = [o for o in outcomes if o.source == "varcode_motif"]
    assert len(cryptic) >= 1
    for o in cryptic:
        assert isinstance(o.effect, CrypticExonCandidate)
        assert 0.0 <= o.probability <= 1.0
        assert "donor_score" in o.evidence
        assert "acceptor_score" in o.evidence


def test_sv_without_cryptic_motifs_has_single_outcome():
    """An SV whose flanking / alt_assembly has no plausible motifs
    should still produce a single primary outcome — the enumerator
    invocation must not error on empty results."""
    transcript = _cftr()
    sv = StructuralVariant(
        contig="7",
        start=transcript.start + 100,
        end=transcript.start + 5_000,
        sv_type="DEL",
        alt_assembly="A" * 200,  # uniform A — no donor/acceptor motifs
        genome=ensembl_grch38)
    effect = _ANNOTATOR.annotate_on_transcript(sv, transcript)
    outcomes = effect.outcomes
    # Only the primary varcode classification.
    assert all(o.source == "varcode" for o in outcomes)


def test_alt_assembly_also_applies_to_fusion_and_inversion():
    transcript = _cftr()
    assembly = "AACCGGTT" * 4
    for sv_type in ("DUP", "INV"):
        sv = StructuralVariant(
            contig="7",
            start=transcript.start + 100,
            end=transcript.start + 5_000,
            sv_type=sv_type,
            alt_assembly=assembly,
            genome=ensembl_grch38)
        effect = _ANNOTATOR.annotate_on_transcript(sv, transcript)
        mt = effect.mutant_transcript
        assert mt is not None, "alt_assembly ignored for %s" % sv_type
        assert mt.cdna_sequence == assembly
        assert len(mt.reference_segments) == 1


# --------------------------------------------------------------------
# MutantTranscript attachment (#335).
#
# Every SV effect carries a MutantTranscript with reference_segments
# populated. For DEL/DUP, cdna_sequence is derivable from cached
# pyensembl cDNA alone and is populated. For INV, the segment shape
# lands but cdna_sequence assembly is deferred. For BND (fusion /
# translocation) the segment shape lands; fused-protein math is #336.
# --------------------------------------------------------------------


def test_large_deletion_attaches_mutant_transcript_with_derivable_cdna():
    transcript = _cftr()
    sv = StructuralVariant(
        contig="7",
        start=transcript.start + 300,
        end=transcript.start + 30_000,
        sv_type="DEL",
        genome=ensembl_grch38)
    effect = _ANNOTATOR.annotate_on_transcript(sv, transcript)
    mt = effect.mutant_transcript
    assert mt is not None
    assert mt.is_structural
    assert len(mt.reference_segments) >= 2
    # cdna is derivable purely from cached cDNA; must be strictly
    # shorter than the reference (we deleted something).
    assert mt.cdna_sequence is not None
    assert 0 < len(mt.cdna_sequence) < len(str(transcript.sequence))
    # Assembled cDNA matches segment concatenation.
    assembled = "".join(
        str(transcript.sequence)[s.start:s.end]
        for s in mt.reference_segments)
    assert assembled == mt.cdna_sequence


def test_large_duplication_attaches_mutant_transcript_with_longer_cdna():
    transcript = _cftr()
    sv = StructuralVariant(
        contig="7",
        start=transcript.start + 300,
        end=transcript.start + 30_000,
        sv_type="DUP",
        genome=ensembl_grch38)
    effect = _ANNOTATOR.annotate_on_transcript(sv, transcript)
    mt = effect.mutant_transcript
    assert mt is not None
    assert mt.is_structural
    # DUP layout is pre + body(1+) + post — at least 3 segments.
    assert len(mt.reference_segments) >= 3
    # cdna is derivable; must be strictly longer than the reference.
    assert mt.cdna_sequence is not None
    assert len(mt.cdna_sequence) > len(str(transcript.sequence))


def test_inversion_attaches_mutant_transcript_with_strand_flipped_body():
    transcript = _cftr()
    sv = StructuralVariant(
        contig="7",
        start=transcript.start + 300,
        end=transcript.start + 30_000,
        sv_type="INV",
        genome=ensembl_grch38)
    effect = _ANNOTATOR.annotate_on_transcript(sv, transcript)
    mt = effect.mutant_transcript
    assert mt is not None
    # At least one segment is strand='-' — the inverted middle.
    strands = {s.strand for s in mt.reference_segments}
    assert "-" in strands
    # cdna is intentionally not yet assembled here (#335 scope
    # limits inversion to shape only; sequence assembly is a
    # follow-up that requires careful junction handling).
    assert mt.cdna_sequence is None


def test_gene_fusion_attaches_two_segment_mutant_transcript():
    cftr = _cftr()
    brca1 = ensembl_grch38.transcript_by_id(BRCA1_ID)
    mate_pos = brca1.start + 1_000
    sv = StructuralVariant(
        contig="7",
        start=cftr.start + 500,
        sv_type="BND",
        alt="N]17:%d]" % mate_pos,
        mate_contig="17",
        mate_start=mate_pos,
        mate_orientation="]]",
        genome=ensembl_grch38)
    effect = _ANNOTATOR.annotate_on_transcript(sv, cftr)
    assert isinstance(effect, GeneFusion)
    mt = effect.mutant_transcript
    assert mt is not None
    assert len(mt.reference_segments) == 2
    labels = {s.label for s in mt.reference_segments}
    assert "5p_partner" in labels
    assert "3p_partner" in labels


# --------------------------------------------------------------------
# Real oncogenic fusions (#336).
#
# Intron:intron breakpoints that produce in-frame exon:exon fusions.
# Each test pins the retained-cDNA boundaries to the expected exon
# cumulative offsets and asserts the fused protein carries the 5'
# partner's N-terminus.
# --------------------------------------------------------------------


def _sum_exon_lengths(transcript, indices):
    """Sum the lengths of a subset of exons (1-indexed)."""
    return sum(
        t.end - t.start + 1
        for i, t in enumerate(transcript.exons, start=1)
        if i in indices)


def test_brd4_nutm1_fusion_nut_midline_carcinoma():
    """BRD4 (chr19, reverse strand) intron 11-12 joined to NUTM1
    (chr15, forward strand) intron 1-2. This is the canonical NUT
    midline carcinoma fusion — BRD4 exons 1-11 fused to NUTM1 exons
    2-7 (in the NUTM1-002 transcript picked by the annotator)."""
    brd4 = ensembl_grch38.transcript_by_id("ENST00000263377")
    brd4_breakpoint = 15_250_000   # intron 11-12 (reverse strand)
    nutm1_breakpoint = 34_347_000  # intron 1-2 (forward strand, on NUTM1-002)
    sv = StructuralVariant(
        contig="19",
        start=brd4_breakpoint,
        sv_type="BND",
        alt="N]15:%d]" % nutm1_breakpoint,
        mate_contig="15",
        mate_start=nutm1_breakpoint,
        mate_orientation="]]",
        genome=ensembl_grch38)
    effect = _ANNOTATOR.annotate_on_transcript(sv, brd4)
    assert isinstance(effect, GeneFusion)
    mt = effect.mutant_transcript
    assert mt is not None
    assert len(mt.reference_segments) == 2

    # 5p BRD4 retains exons 1-11 (reverse-strand exons whose genomic
    # start is upstream of the breakpoint in transcript order).
    expected_5p_len = _sum_exon_lengths(brd4, set(range(1, 12)))
    assert mt.reference_segments[0].end == expected_5p_len

    # Fused protein starts with BRD4's N-terminus and extends
    # through the junction into NUTM1's body.
    assert mt.mutant_protein_sequence is not None
    brd4_ref = str(brd4.protein_sequence)
    assert mt.mutant_protein_sequence.startswith(brd4_ref[:40])
    # The fusion should not be shorter than the retained 5p
    # N-terminal portion — translation must continue through the
    # junction without hitting an immediate stop.
    retained_n_term_aa = (expected_5p_len - min(
        brd4.start_codon_spliced_offsets)) // 3
    assert len(mt.mutant_protein_sequence) > retained_n_term_aa // 2


def test_bcr_abl1_fusion_philadelphia_chromosome():
    """BCR (chr22, forward) intron 13-14 joined to ABL1 (chr9,
    forward) intron 1-2 — Philadelphia chromosome, chronic myeloid
    leukemia. Retains BCR exons 1-13 and ABL1 exons 2-11."""
    bcr = ensembl_grch38.transcript_by_id("ENST00000305877")
    bcr_breakpoint = 23_290_000      # intron 13-14
    abl1_breakpoint = 130_840_000    # intron 1-2
    sv = StructuralVariant(
        contig="22",
        start=bcr_breakpoint,
        sv_type="BND",
        alt="N]9:%d]" % abl1_breakpoint,
        mate_contig="9",
        mate_start=abl1_breakpoint,
        mate_orientation="]]",
        genome=ensembl_grch38)
    effect = _ANNOTATOR.annotate_on_transcript(sv, bcr)
    assert isinstance(effect, GeneFusion)
    mt = effect.mutant_transcript

    expected_5p_len = _sum_exon_lengths(bcr, set(range(1, 14)))
    assert mt.reference_segments[0].end == expected_5p_len

    # 3p starts from ABL1 exon 2 onward; exon 1 of ABL1-001 is 460 bp.
    abl1 = ensembl_grch38.transcript_by_id("ENST00000318560")
    expected_3p_offset = _sum_exon_lengths(abl1, {1})
    assert mt.reference_segments[1].start == expected_3p_offset

    assert mt.mutant_protein_sequence is not None
    bcr_ref = str(bcr.protein_sequence)
    assert mt.mutant_protein_sequence.startswith(bcr_ref[:40])


def test_ewsr1_fli1_fusion_ewing_sarcoma_type1():
    """EWSR1 (chr22, forward) intron 7-8 joined to FLI1 (chr11,
    forward) intron 5-6 — Ewing sarcoma type-1 fusion. Retains
    EWSR1 exons 1-7 and FLI1 exons 6-10 (in FLI1-201, the
    transcript the annotator auto-picks at the breakpoint)."""
    ewsr1 = ensembl_grch38.transcript_by_id("ENST00000414183")
    ewsr1_breakpoint = 29_285_000    # intron 7-8
    fli1_breakpoint = 128_800_000    # intron 5-6 in FLI1-201 coords
    sv = StructuralVariant(
        contig="22",
        start=ewsr1_breakpoint,
        sv_type="BND",
        alt="N]11:%d]" % fli1_breakpoint,
        mate_contig="11",
        mate_start=fli1_breakpoint,
        mate_orientation="]]",
        genome=ensembl_grch38)
    effect = _ANNOTATOR.annotate_on_transcript(sv, ewsr1)
    assert isinstance(effect, GeneFusion)
    mt = effect.mutant_transcript

    expected_5p_len = _sum_exon_lengths(ewsr1, set(range(1, 8)))
    assert mt.reference_segments[0].end == expected_5p_len

    assert mt.mutant_protein_sequence is not None
    ewsr1_ref = str(ewsr1.protein_sequence)
    assert mt.mutant_protein_sequence.startswith(ewsr1_ref[:40])


def test_fusion_cdna_equals_segment_concatenation():
    """Assembled cDNA must match the concatenation of the two
    ReferenceSegment slices — the segment layout and the sequence
    can't drift."""
    brd4 = ensembl_grch38.transcript_by_id("ENST00000263377")
    sv = StructuralVariant(
        contig="19",
        start=15_250_000,
        sv_type="BND",
        alt="N]15:34347000]",
        mate_contig="15",
        mate_start=34_347_000,
        mate_orientation="]]",
        genome=ensembl_grch38)
    effect = _ANNOTATOR.annotate_on_transcript(sv, brd4)
    mt = effect.mutant_transcript
    partner = effect.partner_transcript
    fivep = str(brd4.sequence)
    threep = str(partner.sequence)
    s0, s1 = mt.reference_segments
    assembled = fivep[s0.start:s0.end] + threep[s1.start:s1.end]
    assert assembled == mt.cdna_sequence


def test_fusion_protein_terminates_within_or_at_3p_partner_end():
    """Translation stops at the first stop codon — which for a true
    in-frame fusion lies inside the retained 3p cDNA. Assert we
    don't overshoot the fused cDNA (the translator would return
    whatever bases fit, without extending)."""
    brd4 = ensembl_grch38.transcript_by_id("ENST00000263377")
    sv = StructuralVariant(
        contig="19",
        start=15_250_000,
        sv_type="BND",
        alt="N]15:34347000]",
        mate_contig="15",
        mate_start=34_347_000,
        mate_orientation="]]",
        genome=ensembl_grch38)
    effect = _ANNOTATOR.annotate_on_transcript(sv, brd4)
    mt = effect.mutant_transcript
    # Codon-level bound: protein aa count * 3 <= len(cdna) - cds_start.
    cds_start = min(brd4.start_codon_spliced_offsets)
    assert len(mt.mutant_protein_sequence) * 3 <= (
        len(mt.cdna_sequence) - cds_start)


def test_fusion_with_reverse_strand_3p_partner():
    """Exercises the reverse-strand branch of
    :func:`_cdna_offset_at_3p_breakpoint` — CFTR (forward strand, 5p)
    fused to BRCA1 (reverse strand, 3p). Not a biologically observed
    fusion; the purpose is to pin the arithmetic for a forward-5p +
    reverse-3p pair."""
    cftr = _cftr()
    brca1 = ensembl_grch38.transcript_by_id(BRCA1_ID)
    # CFTR intron 1 (between exon 1 end 117480148 and exon 2 start
    # 117504253). Forward strand.
    cftr_breakpoint = 117_485_000
    # BRCA1 intron 2-3 (between exon 3 end 43115779 and exon 2 start
    # 43124017) on reverse strand. Pick 43120000, strictly in the
    # intron.
    brca1_breakpoint = 43_120_000
    sv = StructuralVariant(
        contig="7",
        start=cftr_breakpoint,
        sv_type="BND",
        alt="N]17:%d]" % brca1_breakpoint,
        mate_contig="17",
        mate_start=brca1_breakpoint,
        mate_orientation="]]",
        genome=ensembl_grch38)
    effect = _ANNOTATOR.annotate_on_transcript(sv, cftr)
    assert isinstance(effect, GeneFusion)
    mt = effect.mutant_transcript

    # CFTR retains exon 1 only (185 bp).
    expected_5p_len = _sum_exon_lengths(cftr, {1})
    assert mt.reference_segments[0].end == expected_5p_len

    # BRCA1 3p retains exons 3 onward; exons 1+2 sum to 199 bp in
    # transcript order (reverse-strand).
    assert mt.reference_segments[1].start == _sum_exon_lengths(
        effect.partner_transcript, {1, 2})


def test_fusion_protein_is_none_when_5p_cds_past_breakpoint():
    """When the 5p partner's CDS start offset lies beyond the
    retained cDNA, :func:`_translate_fused_cdna` returns None so the
    effect carries no fused-protein sequence (documented branch of
    #336)."""
    from varcode.annotators.structural_variant import _translate_fused_cdna
    cftr = _cftr()
    cds_start = min(cftr.start_codon_spliced_offsets)
    # Simulate a fused cDNA where the retained 5p portion is shorter
    # than the CDS start.
    fake_fused = "A" * (cds_start * 2)  # plenty of length overall
    result = _translate_fused_cdna(
        fake_fused, cftr, five_prime_len=cds_start - 1)
    assert result is None


def test_deletion_with_partial_exon_overlap():
    """A DEL that cuts mid-exon must produce segments whose cDNA
    concatenation equals the reference cDNA minus the deleted
    exonic bases. Exercises the partial-exon branches of
    :func:`_cdna_ranges_kept_after_deletion`."""
    transcript = _cftr()
    exon1 = transcript.exons[0]
    # Deletion that starts 50 bp into CFTR exon 1 and extends 10 kb
    # into the intron — partial exon overlap on the right side.
    del_start = exon1.start + 50
    del_end = exon1.end + 10_000
    sv = StructuralVariant(
        contig="7",
        start=del_start,
        end=del_end,
        sv_type="DEL",
        genome=ensembl_grch38)
    effect = _ANNOTATOR.annotate_on_transcript(sv, transcript)
    mt = effect.mutant_transcript
    assert mt is not None
    # 50 bp of exon 1 survive (the 5' portion before the cut).
    assert mt.reference_segments[0].start == 0
    assert mt.reference_segments[0].end == 50
    # Subsequent segments start at exon 2's cDNA offset (exon 1 is
    # 185 bp, all after cut position 50 is deleted).
    later_starts = [s.start for s in mt.reference_segments[1:]]
    assert all(s >= 185 for s in later_starts), (
        "Expected post-deletion segments at or after cDNA offset 185, "
        "got %r" % later_starts)
    # Assembled cDNA matches segment concatenation.
    cdna = str(transcript.sequence)
    assembled = "".join(cdna[s.start:s.end] for s in mt.reference_segments)
    assert assembled == mt.cdna_sequence


def test_warns_on_reverse_complement_mate_orientation():
    """Reverse-complement BND orientations (``[]`` / ``][``) trip a
    warning so a caller doesn't silently get canonical-direction
    output (#336)."""
    import warnings
    cftr = _cftr()
    brca1 = ensembl_grch38.transcript_by_id(BRCA1_ID)
    sv = StructuralVariant(
        contig="7",
        start=cftr.start + 500,
        sv_type="BND",
        alt="N[17:%d[" % (brca1.start + 1000),
        mate_contig="17",
        mate_start=brca1.start + 1000,
        mate_orientation="[]",
        genome=ensembl_grch38)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _ANNOTATOR.annotate_on_transcript(sv, cftr)
    msgs = [str(w.message) for w in caught]
    assert any("reverse-complement" in m for m in msgs), (
        "Expected reverse-complement orientation warning, got %r" % msgs)


def test_canonical_mate_orientation_does_not_warn():
    """Canonical ``]]`` / ``[[`` pairings — or a missing orientation
    — pass silently. Guards against the reverse-complement warning
    becoming noisy for the common case."""
    import warnings
    cftr = _cftr()
    brca1 = ensembl_grch38.transcript_by_id(BRCA1_ID)
    sv = StructuralVariant(
        contig="7",
        start=cftr.start + 500,
        sv_type="BND",
        alt="N]17:%d]" % (brca1.start + 1000),
        mate_contig="17",
        mate_start=brca1.start + 1000,
        mate_orientation="]]",
        genome=ensembl_grch38)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _ANNOTATOR.annotate_on_transcript(sv, cftr)
    assert not any(
        "reverse-complement" in str(w.message) for w in caught)


def test_translocation_to_intergenic_attaches_single_segment_mutant_transcript():
    cftr = _cftr()
    sv = StructuralVariant(
        contig="7",
        start=cftr.start + 500,
        sv_type="BND",
        alt="N]22:15500000]",
        mate_contig="22",
        mate_start=15_500_000,
        mate_orientation="]]",
        genome=ensembl_grch38)
    effect = _ANNOTATOR.annotate_on_transcript(sv, cftr)
    assert isinstance(effect, TranslocationToIntergenic)
    mt = effect.mutant_transcript
    assert mt is not None
    assert len(mt.reference_segments) == 1
