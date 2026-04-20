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
    """Variants carrying a long-read-derived ``alt_assembly`` still
    classify cleanly. The annotator doesn't *use* the assembly yet
    (that's a future refinement for computing fused-protein
    sequences) but must not crash on it."""
    transcript = _cftr()
    sv = StructuralVariant(
        contig="7",
        start=transcript.start + 100,
        end=transcript.start + 5_000,
        sv_type="DEL",
        alt_assembly="ACGTACGTACGTACGT",  # dummy assembled sequence
        genome=ensembl_grch38,
    )
    effect = _ANNOTATOR.annotate_on_transcript(sv, transcript)
    assert isinstance(effect, LargeDeletion)


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
    # Protein / cdna assembly is #336 — not populated here.
    assert mt.cdna_sequence is None


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
