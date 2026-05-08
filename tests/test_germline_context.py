# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Tests for :class:`varcode.germline.GermlineContext` (#268, slice 1).

Slice 1 ships the input contract — constructors, validation,
sparseness flag, window-based lookup. No effect-side logic yet
(annotator dispatch, transcript construction with germline applied,
phase enumeration, LOH detection are subsequent slices).
"""
from __future__ import annotations

import os
import tempfile
import warnings

import pytest

from varcode import (
    GenomeBuildMismatchError,
    GermlineContext,
    SampleNotFoundError,
    Sparseness,
    Variant,
    load_vcf,
)
from varcode.variant_collection import VariantCollection


# --------------------------------------------------------------------
# Sparseness enum
# --------------------------------------------------------------------


class TestSparseness:
    """The completeness flag is the load-bearing piece for correctness:
    a downstream caller deciding whether absence-of-a-call means
    ref/ref needs to read it, so its values must be stable and
    comparable."""

    def test_distinct_values(self):
        # All four states distinguishable.
        assert len({Sparseness.COMPLETE, Sparseness.SPARSE,
                    Sparseness.HOTSPOTS_ONLY, Sparseness.EMPTY}) == 4

    def test_round_trip_via_value(self):
        # Stable string identifiers — important for serialization
        # to CSVs / JSON evidence dicts in later slices.
        for s in Sparseness:
            assert Sparseness(s.value) is s


# --------------------------------------------------------------------
# Empty context
# --------------------------------------------------------------------


class TestEmptyContext:
    def test_empty_is_falsy(self):
        ctx = GermlineContext.empty()
        assert not ctx
        assert len(ctx) == 0
        assert ctx.completeness is Sparseness.EMPTY

    def test_empty_has_no_reference(self):
        # No claim about reference build for an empty context — we
        # don't want to falsely block a cross-VCF validation just
        # because a default crept in.
        assert GermlineContext.empty().reference_name is None

    def test_empty_lookup_returns_empty_tuple(self):
        ctx = GermlineContext.empty()
        # Window lookup on an empty context is well-defined — no
        # exceptions, just empty tuple. Lets callers treat the empty
        # case the same as the no-overlap case in their loops.
        assert ctx.variants_in_window("1", 1, 1_000_000) == ()


# --------------------------------------------------------------------
# from_variants — direct construction
# --------------------------------------------------------------------


class TestFromVariants:
    def _v(self, **kwargs):
        defaults = dict(contig="1", start=100, ref="A", alt="T",
                        genome="GRCh38")
        defaults.update(kwargs)
        return Variant(**defaults)

    def test_accepts_iterable(self):
        ctx = GermlineContext.from_variants(
            [self._v()], reference_name="GRCh38")
        assert len(ctx) == 1

    def test_accepts_existing_collection(self):
        vc = VariantCollection([self._v()])
        ctx = GermlineContext.from_variants(vc, reference_name="GRCh38")
        # Reusing the same collection — no copy.
        assert ctx.variants is vc

    def test_default_completeness_is_complete(self):
        ctx = GermlineContext.from_variants(
            [self._v()], reference_name="GRCh38")
        assert ctx.completeness is Sparseness.COMPLETE

    def test_completeness_override(self):
        ctx = GermlineContext.from_variants(
            [self._v()],
            completeness=Sparseness.SPARSE,
            reference_name="GRCh38")
        assert ctx.completeness is Sparseness.SPARSE

    def test_metadata_passes_through(self):
        ctx = GermlineContext.from_variants(
            [self._v()],
            metadata={"caller": "DeepVariant", "version": "1.6"},
            reference_name="GRCh38")
        assert ctx.metadata["caller"] == "DeepVariant"

    def test_reference_inferred_from_variants_when_consistent(self):
        ctx = GermlineContext.from_variants([self._v()])
        # Variants all carry GRCh38 genome → context picks it up.
        assert ctx.reference_name == "GRCh38"


# --------------------------------------------------------------------
# Window lookup
# --------------------------------------------------------------------


class TestWindowLookup:
    """The window lookup is what slice 2's germline-application path
    will call millions of times during effect prediction — pin
    correctness on the boundary cases now while it's free to fix."""

    def _ctx(self, variants):
        return GermlineContext.from_variants(
            variants, reference_name="GRCh38")

    def _v(self, contig, start, ref="A", alt="T"):
        return Variant(contig=contig, start=start, ref=ref, alt=alt,
                       genome="GRCh38")

    def test_overlap_returns_variant(self):
        v = self._v("1", 100)
        ctx = self._ctx([v])
        assert v in ctx.variants_in_window("1", 50, 150)

    def test_inclusive_lower_bound(self):
        v = self._v("1", 100)
        ctx = self._ctx([v])
        # Window starts exactly at the variant.
        assert ctx.variants_in_window("1", 100, 200) == (v,)

    def test_inclusive_upper_bound(self):
        v = self._v("1", 100)
        ctx = self._ctx([v])
        # Window ends exactly at the variant.
        assert ctx.variants_in_window("1", 50, 100) == (v,)

    def test_outside_window_excluded(self):
        v = self._v("1", 100)
        ctx = self._ctx([v])
        assert ctx.variants_in_window("1", 200, 300) == ()
        assert ctx.variants_in_window("1", 1, 50) == ()

    def test_other_contig_excluded(self):
        ctx = self._ctx([self._v("1", 100)])
        assert ctx.variants_in_window("2", 50, 150) == ()
        assert ctx.variants_in_window("X", 50, 150) == ()

    def test_multiple_variants_in_window(self):
        vs = [self._v("1", 100), self._v("1", 110), self._v("1", 105)]
        ctx = self._ctx(vs)
        result = ctx.variants_in_window("1", 95, 115)
        assert len(result) == 3

    def test_index_is_cached(self):
        v = self._v("1", 100)
        ctx = self._ctx([v])
        # First call builds the index; cache attribute appears.
        ctx.variants_in_window("1", 50, 150)
        assert hasattr(ctx, "_index_cache")
        cached = ctx._index_cache
        # Second call returns the same index dict (no rebuild).
        ctx.variants_in_window("1", 50, 150)
        assert ctx._index_cache is cached

    def test_returned_value_is_tuple(self):
        # Tuple → callers can safely cache/share without worrying
        # about downstream mutation.
        v = self._v("1", 100)
        ctx = self._ctx([v])
        result = ctx.variants_in_window("1", 50, 150)
        assert isinstance(result, tuple)


# --------------------------------------------------------------------
# Cross-VCF validation
# --------------------------------------------------------------------


class TestValidateAgainst:
    def _v(self, genome):
        return Variant(contig="1", start=100, ref="A", alt="T",
                       genome=genome)

    def test_matching_reference_passes(self):
        germ = GermlineContext.from_variants(
            [self._v("GRCh38")], reference_name="GRCh38")
        somatic = VariantCollection([self._v("GRCh38")])
        # No exception → pass.
        germ.validate_against(somatic)

    def test_mismatched_reference_raises(self):
        germ = GermlineContext.from_variants(
            [self._v("GRCh38")], reference_name="GRCh38")
        somatic = VariantCollection([self._v("GRCh37")])
        with pytest.raises(GenomeBuildMismatchError) as exc_info:
            germ.validate_against(somatic)
        # Error carries both names so the user can see what mismatched.
        assert exc_info.value.somatic_reference == "GRCh37"
        assert exc_info.value.germline_reference == "GRCh38"

    def test_mismatch_opt_out_via_validate_reference_false(self):
        """For users who've explicitly lifted over one VCF into the
        other's build and know they're in sync. The hard error is the
        default; the opt-out is one kwarg."""
        germ = GermlineContext.from_variants(
            [self._v("GRCh38")], reference_name="GRCh38")
        somatic = VariantCollection([self._v("GRCh37")])
        # No exception even with a clear mismatch.
        germ.validate_against(somatic, validate_reference=False)

    def test_unknown_reference_skips_check(self):
        """If we couldn't determine either reference, we can't validate.
        Don't raise on uncertainty — just no-op. Slice 2's annotator
        plumbing will surface coordinate mismatches as
        ReferenceMismatchError on first contact."""
        germ = GermlineContext.from_variants([self._v("GRCh38")])
        # Manually wipe the reference (simulates a context built from
        # variants with no genome attribute).
        object.__setattr__(germ, "reference_name", None)
        somatic = VariantCollection([self._v("GRCh37")])
        # No raise — we don't know enough to claim mismatch.
        germ.validate_against(somatic)

    def test_warns_on_empty_context_with_non_empty_completeness(self):
        """A COMPLETE-flagged context with zero variants is almost
        always a wrong-file or filter-too-aggressive bug. We don't
        want this to silently pass — warn so the user investigates."""
        germ = GermlineContext.from_variants(
            [], completeness=Sparseness.COMPLETE,
            reference_name="GRCh38")
        somatic = VariantCollection([self._v("GRCh38")])
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter("always")
            germ.validate_against(somatic)
        assert any(
            "zero variants" in str(w.message)
            for w in captured), (
                "Expected zero-variants warning, got %s" % [
                    str(w.message) for w in captured])

    def test_no_warning_on_genuinely_empty_context(self):
        """An EMPTY-flagged context is the explicit fallback shape;
        the user has opted into reference-relative annotation.
        Warning here would be noise."""
        germ = GermlineContext.empty()
        somatic = VariantCollection([self._v("GRCh38")])
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter("always")
            germ.validate_against(somatic)
        assert not any(
            "zero variants" in str(w.message) for w in captured)


# --------------------------------------------------------------------
# from_germline_vcf — load path
# --------------------------------------------------------------------


def _write_vcf(body):
    fd, path = tempfile.mkstemp(suffix=".vcf")
    with os.fdopen(fd, "w") as f:
        f.write(body)
    return path


class TestFromGermlineVCF:
    def test_loads_full_call_set(self):
        body = (
            "##fileformat=VCFv4.2\n"
            "##reference=GRCh38\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t100\t.\tA\tT\t100\tPASS\t.\n"
            "1\t200\t.\tC\tG\t100\tPASS\t.\n"
        )
        path = _write_vcf(body)
        try:
            ctx = GermlineContext.from_germline_vcf(
                path, genome="GRCh38")
        finally:
            os.unlink(path)
        assert len(ctx) == 2
        assert ctx.completeness is Sparseness.COMPLETE
        assert ctx.reference_name == "GRCh38"

    def test_warns_on_zero_variant_load(self):
        """The most common 'wrong file' symptom — a VCF with no
        records. Catching this at load time saves debugging cycles
        downstream when every somatic variant comes back
        reference-relative."""
        body = (
            "##fileformat=VCFv4.2\n"
            "##reference=GRCh38\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )
        path = _write_vcf(body)
        try:
            with warnings.catch_warnings(record=True) as captured:
                warnings.simplefilter("always")
                ctx = GermlineContext.from_germline_vcf(
                    path, genome="GRCh38")
        finally:
            os.unlink(path)
        assert len(ctx) == 0
        assert any("zero variants" in str(w.message) for w in captured)


# --------------------------------------------------------------------
# from_multi_sample_vcf
# --------------------------------------------------------------------


class TestFromMultiSampleVCF:
    def _multi_sample_body(self):
        return (
            "##fileformat=VCFv4.2\n"
            "##reference=GRCh38\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"genotype\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n"
            "1\t100\t.\tA\tT\t100\tPASS\t.\tGT\t0/1\t0/1\n"
            "1\t200\t.\tC\tG\t100\tPASS\t.\tGT\t0/0\t0/1\n"
            "1\t300\t.\tG\tA\t100\tPASS\t.\tGT\t0/1\t1/1\n"
        )

    def test_extracts_calls_for_named_sample(self):
        """The NORMAL column has calls at positions 100 and 300 (GT
        non-ref) but not 200 (GT 0/0). Extracting NORMAL should give
        us those two."""
        path = _write_vcf(self._multi_sample_body())
        try:
            ctx = GermlineContext.from_multi_sample_vcf(
                path, sample="NORMAL",
                completeness=Sparseness.SPARSE,
                genome="GRCh38")
        finally:
            os.unlink(path)
        assert len(ctx) == 2
        starts = {v.start for v in ctx}
        assert starts == {100, 300}

    def test_completeness_is_required(self):
        """Multi-sample VCFs from somatic callers are sparse; pure-
        germline multi-sample VCFs are complete. Forcing the caller
        to declare prevents a class of correctness bugs."""
        path = _write_vcf(self._multi_sample_body())
        try:
            with pytest.raises(TypeError):
                # Missing required keyword-only arg: completeness=
                GermlineContext.from_multi_sample_vcf(
                    path, sample="NORMAL", genome="GRCh38")
        finally:
            os.unlink(path)

    def test_unknown_sample_raises(self):
        path = _write_vcf(self._multi_sample_body())
        try:
            with pytest.raises(SampleNotFoundError) as exc_info:
                GermlineContext.from_multi_sample_vcf(
                    path, sample="DOES_NOT_EXIST",
                    completeness=Sparseness.COMPLETE,
                    genome="GRCh38")
            # Error message lists available samples so the user sees
            # what they could have asked for.
            assert "NORMAL" in str(exc_info.value)
            assert "TUMOR" in str(exc_info.value)
        finally:
            os.unlink(path)

    def test_metadata_records_source_and_sample(self):
        path = _write_vcf(self._multi_sample_body())
        try:
            ctx = GermlineContext.from_multi_sample_vcf(
                path, sample="NORMAL",
                completeness=Sparseness.SPARSE,
                genome="GRCh38")
        finally:
            os.unlink(path)
        assert ctx.metadata["sample"] == "NORMAL"
        assert ctx.metadata["source_path"] == path


# --------------------------------------------------------------------
# Iteration / __bool__ / __len__ contracts
# --------------------------------------------------------------------


class TestDunders:
    """The dunder protocol matters for ergonomic use: ``if germline:``
    has to read intuitively, and iteration / len have to behave like
    a VariantCollection so callers can drop a context wherever a
    collection used to flow."""

    def test_bool_true_when_non_empty(self):
        ctx = GermlineContext.from_variants(
            [Variant(contig="1", start=1, ref="A", alt="T",
                     genome="GRCh38")],
            reference_name="GRCh38")
        assert bool(ctx) is True

    def test_bool_false_for_empty_complete(self):
        """A COMPLETE context with zero variants is falsy too — there's
        no germline data to apply, even though the caller intended
        completeness. The completeness flag stays for the
        rare-call-set case where downstream code wants to know
        'they meant complete but couldn't load any.'"""
        ctx = GermlineContext.from_variants(
            [], completeness=Sparseness.COMPLETE,
            reference_name="GRCh38")
        assert bool(ctx) is False

    def test_bool_false_for_empty_sparse(self):
        ctx = GermlineContext.from_variants(
            [], completeness=Sparseness.SPARSE,
            reference_name="GRCh38")
        assert bool(ctx) is False

    def test_iteration_yields_variants(self):
        v = Variant(contig="1", start=100, ref="A", alt="T",
                    genome="GRCh38")
        ctx = GermlineContext.from_variants([v], reference_name="GRCh38")
        assert list(ctx) == [v]

    def test_frozen_dataclass_immutable(self):
        """Frozen so callers can hash/share contexts safely. Mutation
        attempts should fail loudly so the immutability contract
        isn't accidentally relied on."""
        ctx = GermlineContext.empty()
        with pytest.raises(Exception):  # FrozenInstanceError or AttributeError
            ctx.completeness = Sparseness.COMPLETE
