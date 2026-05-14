# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Tests for :mod:`varcode.rna_evidence` (#259).

The protocol + apply hook should:

* Be a no-op when no resolver is supplied (existing pipelines unchanged).
* Append RNA-observed outcomes to ``MultiOutcomeEffect.outcomes`` in
  the agreed-upon ordering (DNA-predicted first, RNA-observed after)
  and via the back-compat ``_with_extra_outcomes`` hook so subclass
  overrides of ``outcomes`` (SpliceOutcomeSet, StructuralVariantEffect)
  pick up evidence uniformly.
* Leave non-multi-outcome effects untouched even when the resolver
  has evidence — the contract is "refinement of a possibility set",
  not "replacement of a deterministic call".
* Compose with subsequent evidence calls so a pipeline can attach
  Isovar evidence first, then add Exacto evidence later.
"""
from __future__ import annotations

import os
import tempfile

import pytest

from pyensembl import cached_release

from varcode import (
    NullRNAEvidenceResolver,
    RNAEvidenceResolver,
    StructuralVariant,
    apply_rna_evidence_to_effects,
    make_rna_outcome,
    load_vcf,
)
from varcode.effects.effect_classes import (
    LargeDeletion,
    MultiOutcomeEffect,
    StructuralVariantEffect,
    Substitution,
)
from varcode.effect_outcomes import EffectOutcome


ensembl_grch38 = cached_release(81)
CFTR_ID = "ENST00000003084"


def _cftr():
    return ensembl_grch38.transcript_by_id(CFTR_ID)


# --------------------------------------------------------------------
# Protocol contract
# --------------------------------------------------------------------


class _StubResolver:
    """In-test resolver that returns a fixed set of outcomes for any
    (variant, transcript) pair. Useful for end-to-end tests where we
    want to confirm the apply hook fires without standing up a real
    Isovar pipeline."""

    def __init__(self, outcomes):
        self._outcomes = tuple(outcomes)

    def observed_outcomes(self, variant, transcript):
        return self._outcomes


class TestProtocol:
    def test_null_resolver_satisfies_protocol(self):
        # runtime_checkable Protocol: duck-type via observed_outcomes.
        assert isinstance(NullRNAEvidenceResolver(), RNAEvidenceResolver)

    def test_arbitrary_class_with_method_satisfies_protocol(self):
        assert isinstance(_StubResolver([]), RNAEvidenceResolver)

    def test_object_without_method_does_not_satisfy(self):
        class NotAResolver:
            def something_else(self):
                pass
        assert not isinstance(NotAResolver(), RNAEvidenceResolver)


# --------------------------------------------------------------------
# make_rna_outcome helper
# --------------------------------------------------------------------


class TestMakeRNAOutcome:
    """The factory shapes the evidence dict the way downstream code
    expects, while still letting producer-specific fields ride along
    in extra_evidence."""

    def test_minimal_call_marks_source_rna(self):
        o = make_rna_outcome(effect=Substitution)  # placeholder effect
        assert o.source == "rna"
        assert o.evidence == {}
        assert o.probability is None

    def test_well_known_fields_land_in_evidence(self):
        o = make_rna_outcome(
            effect=Substitution,
            transcript_model_id="ISOFORM_42",
            read_count=314,
            probability=0.91)
        assert o.evidence["transcript_model_id"] == "ISOFORM_42"
        assert o.evidence["read_count"] == 314
        assert o.probability == 0.91

    def test_extra_evidence_merges_alongside_well_known(self):
        o = make_rna_outcome(
            effect=Substitution,
            read_count=10,
            extra_evidence={"tpm": 4.2, "junction_id": "JUNC_1"})
        assert o.evidence["read_count"] == 10
        assert o.evidence["tpm"] == 4.2
        assert o.evidence["junction_id"] == "JUNC_1"

    def test_source_override_for_tool_specific_filtering(self):
        # Isovar / Exacto / longread_assembly producers set their own
        # source string so consumers can filter by tool.
        o = make_rna_outcome(effect=Substitution, source="isovar")
        assert o.source == "isovar"


# --------------------------------------------------------------------
# apply_rna_evidence_to_effects — no-op cases
# --------------------------------------------------------------------


class TestApplyNoOpCases:
    def test_none_resolver_returns_effects_unchanged(self):
        marker = object()
        effects = [marker]
        result = apply_rna_evidence_to_effects(effects, None)
        assert result is effects
        assert effects == [marker]

    def test_resolver_without_method_is_silent_noop(self):
        """Defensive: if a caller passes the wrong object, we don't
        crash mid-walk. Mirrors the
        ``apply_phase_resolver_to_effects`` contract."""
        class Garbage:
            pass
        marker = object()
        effects = [marker]
        result = apply_rna_evidence_to_effects(effects, Garbage())
        assert result is effects

    def test_null_resolver_leaves_outcomes_untouched(self):
        sv = StructuralVariant(
            contig="7",
            start=_cftr().start + 100,
            end=_cftr().start + 50_000,
            sv_type="DEL",
            genome=ensembl_grch38)
        effect = sv.effect_on_transcript(_cftr())
        baseline = tuple(effect.outcomes)
        apply_rna_evidence_to_effects([effect], NullRNAEvidenceResolver())
        assert tuple(effect.outcomes) == baseline


# --------------------------------------------------------------------
# apply_rna_evidence_to_effects — happy path
# --------------------------------------------------------------------


class TestApplyAttaches:
    def _sv_effect(self):
        sv = StructuralVariant(
            contig="7",
            start=_cftr().start + 100,
            end=_cftr().start + 50_000,
            sv_type="DEL",
            genome=ensembl_grch38)
        effect = sv.effect_on_transcript(_cftr())
        # Sanity: the SV annotator returns a LargeDeletion (a
        # StructuralVariantEffect, which is a MultiOutcomeEffect),
        # so it has an ``outcomes`` view that should pick up extras.
        assert isinstance(effect, LargeDeletion)
        assert isinstance(effect, MultiOutcomeEffect)
        return effect

    def test_observed_outcomes_appended_to_outcomes_view(self):
        effect = self._sv_effect()
        baseline = tuple(effect.outcomes)
        # Construct a stand-in observed outcome (the effect can be
        # whatever — for this test we recycle the LargeDeletion's own
        # most-likely candidate, since we just need a MutationEffect
        # to wrap).
        observed = make_rna_outcome(
            effect=baseline[0].effect,
            source="isovar",
            transcript_model_id="ISOFORM_A",
            read_count=42,
            probability=0.78)
        apply_rna_evidence_to_effects(
            [effect], _StubResolver([observed]))

        new_outcomes = tuple(effect.outcomes)
        # DNA-predicted outcomes preserved at the front, RNA-observed
        # appended at the end. Consumers reading by source can pick
        # whichever they want; consumers iterating in order get the
        # most-likely DNA call first (preserves existing semantics).
        assert len(new_outcomes) == len(baseline) + 1
        assert new_outcomes[:len(baseline)] == baseline
        assert new_outcomes[-1] is observed
        assert new_outcomes[-1].source == "isovar"

    def test_apply_is_idempotent_on_repeated_call_with_same_resolver(self):
        """Calling apply twice with the same resolver doesn't duplicate
        outcomes — wait, actually it DOES extend, since the apply
        function is additive by design (so a pipeline can attach Isovar
        first, Exacto after). We pin the additive contract here so
        future refactors don't silently change to dedup."""
        effect = self._sv_effect()
        baseline_len = len(tuple(effect.outcomes))
        observed = make_rna_outcome(
            effect=tuple(effect.outcomes)[0].effect, source="isovar")
        resolver = _StubResolver([observed])
        apply_rna_evidence_to_effects([effect], resolver)
        apply_rna_evidence_to_effects([effect], resolver)
        assert len(tuple(effect.outcomes)) == baseline_len + 2

    def test_mixed_sources_compose(self):
        """The intended workflow: Isovar resolver runs first, Exacto
        adds long-read disambiguation second. Consumers see both
        sources alongside the original DNA-predicted outcomes."""
        effect = self._sv_effect()
        first = make_rna_outcome(
            effect=tuple(effect.outcomes)[0].effect, source="isovar")
        second = make_rna_outcome(
            effect=tuple(effect.outcomes)[0].effect, source="exacto",
            extra_evidence={"long_read_support": True})
        apply_rna_evidence_to_effects([effect], _StubResolver([first]))
        apply_rna_evidence_to_effects([effect], _StubResolver([second]))
        sources = {o.source for o in effect.outcomes}
        assert {"varcode", "isovar", "exacto"}.issubset(sources)


# --------------------------------------------------------------------
# Non-MultiOutcomeEffect cases
# --------------------------------------------------------------------


class TestNonMultiOutcomeNotMutated:
    """Single-outcome effects (Missense, FrameShift, etc.) don't
    expose an outcomes view — replacing them with multi-outcome
    wrappers would break downstream isinstance checks. Pin that the
    apply walk leaves them untouched."""

    def test_missense_not_mutated_when_resolver_has_evidence(self):
        # Use a real point-variant effect.
        from varcode import Variant
        v = Variant(
            contig="7", start=140_753_336, ref="A", alt="T",
            genome=ensembl_grch38)  # BRAF V600E-ish locus
        effects = list(v.effects(raise_on_error=False))
        assert effects, "BRAF locus should produce at least one effect"
        # Pretend we have RNA evidence for this — apply shouldn't
        # crash and should not mutate non-MultiOutcomeEffect entries.
        observed = make_rna_outcome(
            effect=effects[0], source="isovar")
        apply_rna_evidence_to_effects(effects, _StubResolver([observed]))
        for e in effects:
            if not isinstance(e, MultiOutcomeEffect):
                assert not hasattr(e, "_extra_outcomes") or not e._extra_outcomes


# --------------------------------------------------------------------
# Wiring through Variant.effects() and VariantCollection.effects()
# --------------------------------------------------------------------


def _write_vcf(body):
    fd, path = tempfile.mkstemp(suffix=".vcf")
    with os.fdopen(fd, "w") as f:
        f.write(body)
    return path


class TestEffectsKwargWiring:
    """The user-facing API path: ``rna_resolver=`` on
    ``VariantCollection.effects()`` and ``Variant.effects()`` should
    plumb through to the apply hook."""

    def test_collection_effects_with_rna_resolver_attaches(self):
        body = (
            "##fileformat=VCFv4.2\n"
            "##reference=GRCh38\n"
            "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"type\">\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"end\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "7\t%d\t.\tA\t<DEL>\t100\tPASS\tSVTYPE=DEL;END=%d\n"
        ) % (_cftr().start + 100, _cftr().start + 50_000)
        path = _write_vcf(body)
        try:
            vc = load_vcf(
                path, genome=ensembl_grch38,
                parse_structural_variants=True)
        finally:
            os.unlink(path)
        # Capture the baseline (no resolver) and the with-resolver
        # output and compare. Each SV is annotated against multiple
        # transcripts; we just verify that *some* effect picked up
        # observed outcomes.
        baseline = vc.effects()
        baseline_sv_effects = [
            e for e in baseline
            if isinstance(e, StructuralVariantEffect)]
        assert baseline_sv_effects

        # Stub resolver returns one observed outcome per call. Use the
        # baseline's primary outcome's effect as the wrapper target.
        observed_effect = tuple(baseline_sv_effects[0].outcomes)[0].effect
        observed = make_rna_outcome(
            effect=observed_effect, source="isovar",
            transcript_model_id="ISO_1", read_count=11)
        resolver = _StubResolver([observed])

        with_evidence = vc.effects(rna_resolver=resolver)
        sv_effects = [
            e for e in with_evidence
            if isinstance(e, StructuralVariantEffect)]
        # Each SV effect should now carry the observed outcome.
        for e in sv_effects:
            sources = {o.source for o in e.outcomes}
            assert "isovar" in sources, (
                "Expected RNA-observed outcome on %r, got sources=%s"
                % (e, sources))

    def test_variant_effects_with_rna_resolver_attaches(self):
        sv = StructuralVariant(
            contig="7",
            start=_cftr().start + 100,
            end=_cftr().start + 50_000,
            sv_type="DEL",
            genome=ensembl_grch38)
        baseline = list(sv.effects(raise_on_error=False))
        sv_effects_baseline = [
            e for e in baseline if isinstance(e, StructuralVariantEffect)]
        assert sv_effects_baseline

        observed_effect = tuple(sv_effects_baseline[0].outcomes)[0].effect
        observed = make_rna_outcome(
            effect=observed_effect, source="exacto", read_count=7)
        with_evidence = list(
            sv.effects(
                raise_on_error=False,
                rna_resolver=_StubResolver([observed])))
        sv_effects = [
            e for e in with_evidence
            if isinstance(e, StructuralVariantEffect)]
        for e in sv_effects:
            sources = {o.source for o in e.outcomes}
            assert "exacto" in sources


# --------------------------------------------------------------------
# EffectOutcome ordering: DNA-predicted preserved at the front
# --------------------------------------------------------------------


class TestOrderingContract:
    """``outcomes`` returns DNA-predicted entries first, RNA-observed
    after. This preserves ``most_likely`` (which is ``candidates[0]``)
    and lets consumers that read ``outcomes[0]`` for "the primary call"
    keep working unchanged."""

    def test_dna_predicted_outcomes_remain_first(self):
        sv = StructuralVariant(
            contig="7",
            start=_cftr().start + 100,
            end=_cftr().start + 50_000,
            sv_type="DEL",
            genome=ensembl_grch38)
        effect = sv.effect_on_transcript(_cftr())
        dna_predicted = tuple(effect.outcomes)
        observed = make_rna_outcome(
            effect=dna_predicted[0].effect, source="isovar")
        apply_rna_evidence_to_effects([effect], _StubResolver([observed]))
        new = tuple(effect.outcomes)
        # First N entries are exactly the original DNA-predicted set.
        assert new[:len(dna_predicted)] == dna_predicted
        # most_likely should still be the original DNA-predicted top.
        assert effect.most_likely is effect.candidates[0]
