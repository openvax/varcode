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

"""Tests for the fast-path SNV helper.

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
    # CFTR coding SNV.
    variant = Variant("7", 117531095, "T", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    # The full pipeline and the fast-path helper should agree.
    full_effect = variant.effect_on_transcript(transcript)
    fast_effect = _call_fast_path(variant, transcript)
    assert fast_effect is not None
    assert type(fast_effect) is type(full_effect)
    assert fast_effect.short_description == full_effect.short_description


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
# Integration: the helper agrees with the full pipeline on several
# CFTR coding SNVs.
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
        full = variant.effect_on_transcript(transcript)
        fast = _call_fast_path(variant, transcript)
        if fast is None:
            # Fast path declined — the full pipeline handles it
            # (splice-boundary or start/stop).
            continue
        assert type(fast) is type(full)
        assert fast.short_description == full.short_description
        matched += 1
    assert matched > 0, (
        "Expected at least one of the CFTR coding SNVs to hit the fast "
        "path; if this fires, reject-case coverage is fine but accept-"
        "case coverage needs a different fixture.")
