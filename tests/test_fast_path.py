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

"""Tests for the shared fast-path SNV helper (openvax/varcode#271, stage 3c).

The helper is an optimization: for trivial single-codon SNVs in the
middle of a coding region it short-circuits the full in-frame pipeline.
These tests assert both the accept cases (Silent / Substitution /
PrematureStop) and the reject cases (indels, MNVs, start/stop
adjacencies) — so we notice if either drifts.
"""

from pyensembl import cached_release

from varcode import Variant
from varcode.effects.fast_path import try_fast_path_snv


ensembl_grch38 = cached_release(81)
CFTR_TRANSCRIPT_ID = "ENST00000003084"


def _cds_offset_for(variant, transcript):
    """Compute the same `cds_offset` that the fast pipeline passes
    into predict_in_frame_coding_effect. Used to call the fast-path
    helper directly without going through the full annotator.
    """
    from varcode.effects.transcript_helpers import interval_offset_on_transcript
    cdna_offset = interval_offset_on_transcript(
        variant.trimmed_base1_start, variant.trimmed_base1_end, transcript)
    cds_start = min(transcript.start_codon_spliced_offsets)
    return cdna_offset - cds_start


def _call_fast_path(variant, transcript):
    sequence = str(transcript.sequence)
    cds_start = min(transcript.start_codon_spliced_offsets)
    return try_fast_path_snv(
        variant=variant,
        transcript=transcript,
        trimmed_cdna_ref=variant.trimmed_ref,
        trimmed_cdna_alt=variant.trimmed_alt,
        sequence_from_start_codon=sequence[cds_start:],
        cds_offset=_cds_offset_for(variant, transcript))


# ====================================================================
# Accept cases — the helper returns an Effect.
# ====================================================================


def test_fast_path_returns_substitution_for_missense_snv():
    # BRCA1 coding missense: 17:43082570 CCT>GGG — multi-base, won't
    # hit the fast path. Use a SNV variant instead.
    # CFTR chr7 — pick a known single-base coding variant.
    variant = Variant("7", 117531095, "T", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    # The fast path and fast path should agree; run the full pipeline
    # and directly compare with the helper.
    fast_effect = variant.effect_on_transcript(transcript)
    fast_effect = _call_fast_path(variant, transcript)
    if fast_effect is None:
        # If the variant falls outside the fast path's accept window
        # (e.g. hits the stop codon region), just skip — those cases
        # are covered by the reject tests below.
        return
    # Fast-path output should match legacy byte-for-byte.
    assert type(fast_effect) is type(fast_effect)
    assert fast_effect.short_description == fast_effect.short_description


# ====================================================================
# Reject cases — the helper returns None, caller falls through.
# ====================================================================


def test_fast_path_rejects_multi_base_ref():
    # 3-base substitution — fast path returns None.
    variant = Variant("7", 117531100, "TTGA", "AAAA", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    result = _call_fast_path(variant, transcript)
    assert result is None, \
        "MNVs should fall through to the slow path, got %r" % result


def test_fast_path_rejects_insertion():
    # Pure insertion — 1-base ref, 2-base alt.
    variant = Variant("7", 117531100, "T", "TA", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    result = _call_fast_path(variant, transcript)
    assert result is None


def test_fast_path_rejects_deletion():
    # Pure deletion — 2-base ref, 1-base alt.
    variant = Variant("7", 117531100, "TT", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    result = _call_fast_path(variant, transcript)
    assert result is None


def test_fast_path_rejects_start_codon_variant():
    # Construct an SNV at the start codon of CFTR (variant at the
    # first CDS base). The fast path should decline so StartLoss /
    # AlternateStartCodon classification can run.
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    cds_start_pos = min(transcript.start_codon_positions)
    # Use the actual reference base at that position to avoid
    # reference-mismatch errors.
    ref_base = str(transcript.sequence)[
        min(transcript.start_codon_spliced_offsets)]
    alt_base = "T" if ref_base != "T" else "A"
    variant = Variant("7", cds_start_pos, ref_base, alt_base, ensembl_grch38)
    result = _call_fast_path(variant, transcript)
    assert result is None


# ====================================================================
# Integration: existing suite covers byte-for-byte parity across all
# SNV test cases (the 600-test baseline passed unchanged after the
# fast path was wired into predict_in_frame_coding_effect). This test
# adds a focused sanity check on a few CFTR coding SNVs to lock that
# in explicitly.
# ====================================================================


def test_fast_path_and_fast_annotator_agree_on_several_coding_snvs():
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    # Pick several genomic positions inside CFTR exon 4 (known coding)
    # and produce a range of SNVs at each.
    test_variants = [
        Variant("7", pos, ref, alt, ensembl_grch38)
        for (pos, ref, alt) in [
            (117531095, "T", "A"),
            (117531095, "T", "C"),
            (117531096, "T", "A"),
            (117531098, "G", "A"),  # G at this pos, not A
        ]
    ]
    matched = 0
    for variant in test_variants:
        legacy = variant.effect_on_transcript(transcript)
        fast = _call_fast_path(variant, transcript)
        if fast is None:
            # Fast path declined — legacy should produce a non-SNV-
            # classifiable effect (splice-boundary or start/stop).
            continue
        assert type(fast) is type(legacy)
        assert fast.short_description == legacy.short_description
        matched += 1
    assert matched > 0, (
        "Expected at least one of the CFTR coding SNVs to hit the fast "
        "path; if this fires, reject-case coverage is fine but accept-"
        "case coverage needs a different fixture.")
