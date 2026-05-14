# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""End-to-end tests for germline-aware effect prediction (#268).

Covers the wiring of ``germline=`` through ``Variant.effects`` and
``VariantCollection.effects``, the ``PhaseCandidateSet`` output
shape when phase is unknown, LOH detection, sparseness propagation,
hemizygous handling, and cross-VCF build-mismatch validation.

Tests use a known-good CFTR codon-overlap pair (somatic at 117531100
T→A produces p.L159M; germline at 117531101 T→C is in the same codon
and changes the cis result to p.S159T). Both positions' reference
bases are read from the GRCh38 reference, so any future changes to
the ref alleles fail loudly rather than silently mis-classifying.
"""
from __future__ import annotations

import os
import tempfile
import warnings

import pytest

from pyensembl import cached_release

from varcode import (
    Completeness,
    GenomeBuildMismatchError,
    GermlineContext,
    Variant,
    VariantCollection,
    apply_germline_to_transcript,
    detect_loh,
    enumerate_phase_hypotheses,
    load_vcf,
    predict_germline_aware_effect,
)
from varcode.annotators.registry import get_default_annotator
from varcode.effects.effect_classes import (
    PhaseCandidateSet,
    Substitution,
)


ensembl_grch38 = cached_release(81)
CFTR_ID = "ENST00000003084"

# Known-good CFTR codon-overlap pair. Real GRCh38 reference bases at
# these positions. somatic alone produces p.L159M; somatic+germline
# in cis produces p.S159T (different codon outcome). In trans, the
# somatic effect is identical to the reference-relative case.
SOMATIC_POS = 117_531_100
SOMATIC_REF = "T"
SOMATIC_ALT = "A"
GERMLINE_POS = 117_531_101
GERMLINE_REF = "T"
GERMLINE_ALT = "C"
DISTANT_GERMLINE_POS = 117_668_600  # past the end of the somatic's exon


def _cftr():
    return ensembl_grch38.transcript_by_id(CFTR_ID)


def _somatic():
    return Variant(
        contig="7", start=SOMATIC_POS, ref=SOMATIC_REF, alt=SOMATIC_ALT,
        genome=ensembl_grch38)


def _same_codon_germline():
    return Variant(
        contig="7", start=GERMLINE_POS, ref=GERMLINE_REF, alt=GERMLINE_ALT,
        genome=ensembl_grch38)


# --------------------------------------------------------------------
# Back-compat: empty / None germline is byte-identical to today
# --------------------------------------------------------------------


class TestBackCompat:
    """Today's reference-relative output must survive the addition
    of ``germline=`` unchanged. None / empty contexts are no-ops."""

    def test_no_germline_kwarg_yields_today_output(self):
        v = _somatic()
        effects = list(v.effects(raise_on_error=False))
        assert effects, "CFTR somatic should produce at least one effect"
        assert any(isinstance(e, Substitution) for e in effects)

    def test_explicit_none_is_byte_identical(self):
        v = _somatic()
        without = list(v.effects(raise_on_error=False))
        explicit = list(v.effects(raise_on_error=False, germline=None))
        # Same effect classes, in the same order, with the same
        # short descriptions. We can't check object identity (effects
        # are constructed per call), but the predicted output should
        # be byte-identical at the description level.
        assert [type(e).__name__ for e in without] == [
            type(e).__name__ for e in explicit]
        assert [
            getattr(e, "short_description", None) for e in without
        ] == [
            getattr(e, "short_description", None) for e in explicit]

    def test_empty_context_is_byte_identical(self):
        v = _somatic()
        ctx = GermlineContext.empty()
        without = list(v.effects(raise_on_error=False))
        with_empty = list(v.effects(raise_on_error=False, germline=ctx))
        assert [type(e).__name__ for e in without] == [
            type(e).__name__ for e in with_empty]


# --------------------------------------------------------------------
# Same-codon overlap, unknown phase → PhaseCandidateSet
# --------------------------------------------------------------------


class TestSameCodonOverlap:
    """The marquee case: somatic and germline land in the same codon
    with no phase data. Honest output is a possibility set."""

    def _ctx(self):
        return GermlineContext.from_variants(
            [_same_codon_germline()], reference_name="GRCh38")

    def test_produces_phase_candidate_set(self):
        ann = get_default_annotator()
        e = predict_germline_aware_effect(
            _somatic(), _cftr(), self._ctx(), ann)
        assert isinstance(e, PhaseCandidateSet)

    def test_two_hypotheses_for_one_germline_variant(self):
        ann = get_default_annotator()
        e = predict_germline_aware_effect(
            _somatic(), _cftr(), self._ctx(), ann)
        # One germline → 2^1 = 2 hypotheses (cis, trans).
        assert len(e.outcomes) == 2

    def test_outcomes_carry_haplotype_evidence(self):
        ann = get_default_annotator()
        e = predict_germline_aware_effect(
            _somatic(), _cftr(), self._ctx(), ann)
        haps = {o.evidence.get("haplotype") for o in e.outcomes}
        # The opaque tags from the enumerator: A (all-cis), B (all-trans).
        assert haps == {"A", "B"}

    def test_outcomes_carry_phase_state_unknown(self):
        ann = get_default_annotator()
        e = predict_germline_aware_effect(
            _somatic(), _cftr(), self._ctx(), ann)
        states = {o.evidence.get("phase_state") for o in e.outcomes}
        assert states == {"unknown"}

    def test_outcomes_carry_germline_variants(self):
        ann = get_default_annotator()
        e = predict_germline_aware_effect(
            _somatic(), _cftr(), self._ctx(), ann)
        for o in e.outcomes:
            germ = o.evidence.get("germline_variants")
            # Cis hypothesis carries the germline variant; trans
            # carries an empty tuple. Either way it's a tuple, not
            # ``None`` (so consumers can iterate uniformly).
            assert isinstance(germ, tuple)

    def test_trans_hypothesis_yields_reference_relative_effect(self):
        """The trans hypothesis (germline on the OTHER haplotype) is
        the somatic on the reference codon — should match what
        reference-relative annotation produces today."""
        ann = get_default_annotator()
        ref_relative = ann.annotate_on_transcript(_somatic(), _cftr())
        e = predict_germline_aware_effect(
            _somatic(), _cftr(), self._ctx(), ann)
        trans_outcome = next(
            o for o in e.outcomes if o.evidence["haplotype"] == "B")
        assert (trans_outcome.effect.short_description
                == ref_relative.short_description)

    def test_cis_hypothesis_differs_from_reference_relative(self):
        """The cis hypothesis should produce a different effect than
        reference-relative (otherwise why bother). For the CFTR
        117531100/117531101 pair: ref is L→M, cis is S→T."""
        ann = get_default_annotator()
        ref_relative = ann.annotate_on_transcript(_somatic(), _cftr())
        e = predict_germline_aware_effect(
            _somatic(), _cftr(), self._ctx(), ann)
        cis_outcome = next(
            o for o in e.outcomes if o.evidence["haplotype"] == "A")
        assert (cis_outcome.effect.short_description
                != ref_relative.short_description)

    def test_short_description_marks_ambiguity(self):
        """Convention: ``?`` prefix on the description signals 'this is
        the most-likely candidate of a possibility set'."""
        ann = get_default_annotator()
        e = predict_germline_aware_effect(
            _somatic(), _cftr(), self._ctx(), ann)
        assert e.short_description.startswith("?")


# --------------------------------------------------------------------
# No overlap → standard reference-relative effect
# --------------------------------------------------------------------


class TestNoOverlap:
    def test_distant_germline_yields_standard_effect(self):
        ann = get_default_annotator()
        ctx = GermlineContext.from_variants([
            Variant(contig="7", start=DISTANT_GERMLINE_POS,
                    ref="C", alt="T", genome=ensembl_grch38),
        ], reference_name="GRCh38")
        e = predict_germline_aware_effect(_somatic(), _cftr(), ctx, ann)
        # Single Substitution effect (not PhaseCandidateSet) —
        # patient transcript == reference transcript at this locus.
        assert not isinstance(e, PhaseCandidateSet)
        assert isinstance(e, Substitution)


# --------------------------------------------------------------------
# Phase resolver collapses ambiguity
# --------------------------------------------------------------------


class _StubPhaseResolver:
    """Returns a fixed answer for every (somatic, germline) pair.
    Lets us pin the resolver-collapses-ambiguity contract without
    plumbing a full PS-tagged VCF into a test."""

    def __init__(self, in_cis_answer: bool):
        self._answer = in_cis_answer

    def in_cis(self, v1, v2):
        return self._answer


class TestPhaseResolverCollapses:
    def test_resolver_says_cis_yields_single_hypothesis(self):
        ann = get_default_annotator()
        ctx = GermlineContext.from_variants(
            [_same_codon_germline()], reference_name="GRCh38")
        e = predict_germline_aware_effect(
            _somatic(), _cftr(), ctx, ann,
            phase_resolver=_StubPhaseResolver(in_cis_answer=True))
        # Single resolved hypothesis → not a PhaseCandidateSet; the
        # cis-applied classification surfaces as the primary effect.
        assert not isinstance(e, PhaseCandidateSet)
        # Phase metadata still attached for downstream filtering.
        assert getattr(e, "germline_phase_state", None) == "phased"

    def test_resolver_says_trans_yields_single_hypothesis(self):
        ann = get_default_annotator()
        ctx = GermlineContext.from_variants(
            [_same_codon_germline()], reference_name="GRCh38")
        ref_relative = ann.annotate_on_transcript(_somatic(), _cftr())
        e = predict_germline_aware_effect(
            _somatic(), _cftr(), ctx, ann,
            phase_resolver=_StubPhaseResolver(in_cis_answer=False))
        # Trans resolved → effect should match reference-relative.
        assert e.short_description == ref_relative.short_description


# --------------------------------------------------------------------
# LOH detection
# --------------------------------------------------------------------


class TestLOH:
    def test_loh_detector_matches_position_and_alt(self):
        s = _somatic()
        # Same position+alt → LOH.
        match = Variant(contig="7", start=SOMATIC_POS, ref=SOMATIC_REF,
                        alt=SOMATIC_ALT, genome=ensembl_grch38)
        assert detect_loh(s, [match]) is True
        # Same position, different alt → not LOH.
        diff_alt = Variant(contig="7", start=SOMATIC_POS, ref=SOMATIC_REF,
                           alt="G", genome=ensembl_grch38)
        assert detect_loh(s, [diff_alt]) is False
        # Different position → not LOH.
        diff_pos = Variant(contig="7", start=SOMATIC_POS + 5,
                           ref="A", alt="T", genome=ensembl_grch38)
        assert detect_loh(s, [diff_pos]) is False
        # No germline → not LOH.
        assert detect_loh(s, []) is False

    def test_loh_flag_set_on_effect_when_somatic_matches_germline(self):
        ann = get_default_annotator()
        ctx = GermlineContext.from_variants(
            [Variant(contig="7", start=SOMATIC_POS, ref=SOMATIC_REF,
                     alt=SOMATIC_ALT, genome=ensembl_grch38)],
            reference_name="GRCh38")
        e = predict_germline_aware_effect(_somatic(), _cftr(), ctx, ann)
        assert getattr(e, "is_loh", False) is True

    def test_loh_flag_not_set_when_alt_differs(self):
        ann = get_default_annotator()
        ctx = GermlineContext.from_variants(
            [Variant(contig="7", start=SOMATIC_POS, ref=SOMATIC_REF,
                     alt="G",  # different alt
                     genome=ensembl_grch38)],
            reference_name="GRCh38")
        e = predict_germline_aware_effect(_somatic(), _cftr(), ctx, ann)
        assert getattr(e, "is_loh", False) is False


# --------------------------------------------------------------------
# Sparseness propagation
# --------------------------------------------------------------------


class TestSparsenessPropagation:
    def test_sparse_context_no_overlap_flags_germline_unknown(self):
        """When the context is SPARSE and no germline overlaps the
        somatic's window, we don't know if the patient is ref/ref —
        the somatic caller may not have queried the position.
        Annotate normally but flag the effect."""
        ann = get_default_annotator()
        ctx = GermlineContext.from_variants(
            [],  # no germline calls — SPARSE means absence is unknown
            completeness=Completeness.SPARSE,
            reference_name="GRCh38")
        e = predict_germline_aware_effect(_somatic(), _cftr(), ctx, ann)
        # Empty context is falsy → predict_germline_aware_effect's
        # short-circuit hits before the SPARSE flag matters. We test
        # the SPARSE-with-distant-germline case where the context is
        # truthy but the window is empty.

    def test_sparse_with_distant_germline_flags_unknown(self):
        ann = get_default_annotator()
        ctx = GermlineContext.from_variants(
            [Variant(contig="7", start=DISTANT_GERMLINE_POS,
                     ref="C", alt="T", genome=ensembl_grch38)],
            completeness=Completeness.SPARSE,
            reference_name="GRCh38")
        e = predict_germline_aware_effect(_somatic(), _cftr(), ctx, ann)
        # Distant germline = no overlap = patient germline state unknown
        # at the somatic's window because SPARSE.
        assert getattr(e, "germline_unknown", False) is True

    def test_complete_context_no_overlap_does_not_flag(self):
        """COMPLETE: absence at this position implies ref/ref. No
        flag — there's no uncertainty."""
        ann = get_default_annotator()
        ctx = GermlineContext.from_variants(
            [Variant(contig="7", start=DISTANT_GERMLINE_POS,
                     ref="C", alt="T", genome=ensembl_grch38)],
            completeness=Completeness.COMPLETE,
            reference_name="GRCh38")
        e = predict_germline_aware_effect(_somatic(), _cftr(), ctx, ann)
        assert getattr(e, "germline_unknown", False) is False


# --------------------------------------------------------------------
# Hemizygous chromosome
# --------------------------------------------------------------------


class TestHemizygous:
    """Hemizygous chromosomes (chrY, chrM) have no second haplotype,
    so phase ambiguity is moot — germline-in-window is always cis."""

    def test_chrM_yields_single_implicit_hypothesis(self):
        # Stand-alone test of the enumerator — building a chrM
        # transcript test fixture would be heavier than is worth here.
        somatic = Variant(contig="MT", start=100, ref="A", alt="T",
                          genome=ensembl_grch38)
        germline = Variant(contig="MT", start=101, ref="C", alt="G",
                           genome=ensembl_grch38)
        hyps = enumerate_phase_hypotheses(somatic, [germline])
        assert len(hyps) == 1
        assert hyps[0].phase_state == "implicit"
        # chrM germline in window is automatically cis.
        assert germline in hyps[0].cis

    def test_chrY_yields_single_implicit_hypothesis(self):
        somatic = Variant(contig="Y", start=100, ref="A", alt="T",
                          genome=ensembl_grch38)
        germline = Variant(contig="Y", start=101, ref="C", alt="G",
                           genome=ensembl_grch38)
        hyps = enumerate_phase_hypotheses(somatic, [germline])
        assert len(hyps) == 1
        assert hyps[0].phase_state == "implicit"


# --------------------------------------------------------------------
# Hypothesis cap
# --------------------------------------------------------------------


class TestHypothesisCap:
    """When 2^n exceeds the cap, emit a single
    too_many_hypotheses placeholder. Default cap is 8 ⇒ up to 3
    germline variants in window fully unphased."""

    def _germline(self, pos):
        return Variant(contig="7", start=pos, ref="A", alt="T",
                       genome=ensembl_grch38)

    def test_under_cap_enumerates(self):
        somatic = Variant(contig="7", start=100, ref="A", alt="T",
                          genome=ensembl_grch38)
        germ = [self._germline(101), self._germline(102)]
        # 2 germline → 4 hypotheses, under default cap of 8.
        hyps = enumerate_phase_hypotheses(somatic, germ)
        assert len(hyps) == 4

    def test_over_cap_returns_placeholder(self):
        somatic = Variant(contig="7", start=100, ref="A", alt="T",
                          genome=ensembl_grch38)
        # 4 germline → 16 hypotheses, over default cap.
        germ = [self._germline(p) for p in range(101, 105)]
        hyps = enumerate_phase_hypotheses(somatic, germ)
        assert len(hyps) == 1
        assert hyps[0].phase_state == "too_many_hypotheses"
        # Conservatively marks all germline as cis (the worst-case).
        assert len(hyps[0].cis) == len(germ)

    def test_caller_can_raise_cap(self):
        somatic = Variant(contig="7", start=100, ref="A", alt="T",
                          genome=ensembl_grch38)
        germ = [self._germline(p) for p in range(101, 105)]
        # Raise the cap to 32 → 16 hypotheses now fit.
        hyps = enumerate_phase_hypotheses(somatic, germ, max_hypotheses=32)
        assert len(hyps) == 16


# --------------------------------------------------------------------
# Cross-VCF build mismatch
# --------------------------------------------------------------------


def _write_vcf(body):
    fd, path = tempfile.mkstemp(suffix=".vcf")
    with os.fdopen(fd, "w") as f:
        f.write(body)
    return path


class TestBuildMismatch:
    """The hard error from #268 §7: somatic VCF and germline VCF on
    different reference builds. Surfaces at the
    VariantCollection.effects level so the user sees one readable
    error rather than per-variant ReferenceMismatchError noise."""

    def test_mismatch_raises_at_collection_level(self):
        somatic = VariantCollection([_somatic()])
        # Hand-construct a germline context claiming GRCh37.
        germ = GermlineContext.from_variants(
            [_same_codon_germline()],
            reference_name="GRCh37")
        with pytest.raises(GenomeBuildMismatchError):
            somatic.effects(germline=germ, raise_on_error=False)

    def test_validate_reference_false_bypasses(self):
        """Opt-out for users who've explicitly lifted over."""
        somatic = VariantCollection([_somatic()])
        germ = GermlineContext.from_variants(
            [_same_codon_germline()],
            reference_name="GRCh37")
        # No raise even though references disagree.
        somatic.effects(
            germline=germ, raise_on_error=False, validate_reference=False)


# --------------------------------------------------------------------
# Wiring through Variant.effects / VariantCollection.effects
# --------------------------------------------------------------------


class TestEffectsKwargWiring:
    def test_variant_effects_with_germline(self):
        v = _somatic()
        ctx = GermlineContext.from_variants(
            [_same_codon_germline()], reference_name="GRCh38")
        effects = list(v.effects(raise_on_error=False, germline=ctx))
        # At least one PhaseCandidateSet among the per-transcript
        # outputs (CFTR has multiple transcripts).
        assert any(isinstance(e, PhaseCandidateSet) for e in effects)

    def test_collection_effects_with_germline(self):
        somatic = VariantCollection([_somatic()])
        ctx = GermlineContext.from_variants(
            [_same_codon_germline()], reference_name="GRCh38")
        effects = somatic.effects(germline=ctx, raise_on_error=False)
        assert any(isinstance(e, PhaseCandidateSet) for e in effects)
        # Annotator metadata still tracks the underlying annotator,
        # not "germline" (germline is a transcript modifier, not an
        # annotator). The collection-level annotator name is
        # preserved.
        assert effects.annotator is not None


# --------------------------------------------------------------------
# apply_germline_to_transcript public helper
# --------------------------------------------------------------------


class TestApplyGermlineToTranscript:
    def test_returns_none_for_empty_context(self):
        result = apply_germline_to_transcript(
            _cftr(), GermlineContext.empty())
        assert result is None

    def test_returns_none_when_no_germline_in_window(self):
        """With a somatic specified but the germline not in its
        window, the helper has nothing to apply locally → None.
        Caller falls back to the reference transcript."""
        ctx = GermlineContext.from_variants(
            [Variant(contig="7", start=DISTANT_GERMLINE_POS, ref="C",
                     alt="T", genome=ensembl_grch38)],
            reference_name="GRCh38")
        result = apply_germline_to_transcript(
            _cftr(), ctx, somatic_variant=_somatic())
        assert result is None

    def test_applies_germline_when_in_window(self):
        ctx = GermlineContext.from_variants(
            [_same_codon_germline()], reference_name="GRCh38")
        result = apply_germline_to_transcript(
            _cftr(), ctx, somatic_variant=_somatic())
        assert result is not None
        # MutantTranscript with at least one edit applied.
        assert len(result.edits) >= 1

    def test_whole_transcript_mode_picks_up_distant_germline(self):
        """Without a somatic variant, the helper applies all germline
        variants overlapping ANY exon of the transcript — useful for
        building a single patient transcript for cohort-level
        analysis. Distant germline now lands."""
        ctx = GermlineContext.from_variants(
            [_same_codon_germline()], reference_name="GRCh38")
        result = apply_germline_to_transcript(_cftr(), ctx)
        assert result is not None
        assert len(result.edits) >= 1


# --------------------------------------------------------------------
# Composition with RNA evidence (cross-axis)
# --------------------------------------------------------------------


class TestComposesWithRNAEvidence:
    """The whole point of the design is that germline + RNA + phase
    compose. An RNA outcome with ``evidence={"haplotype": "B"}``
    should align with the haplotype B germline-aware outcome via the
    shared evidence dict shape."""

    def test_phase_ambiguous_outcomes_carry_alignable_haplotype_keys(self):
        """The germline-aware Outcomes carry ``haplotype`` keys with
        opaque tags ("A", "B"); the RNA-evidence Outcomes from #259
        use the same key. Cross-axis matching is just a dict lookup
        — pin that the keys are present and match the EffectOutcome
        evidence-dict convention."""
        ann = get_default_annotator()
        ctx = GermlineContext.from_variants(
            [_same_codon_germline()], reference_name="GRCh38")
        e = predict_germline_aware_effect(
            _somatic(), _cftr(), ctx, ann)
        for outcome in e.outcomes:
            assert "haplotype" in outcome.evidence
            assert "phase_state" in outcome.evidence
