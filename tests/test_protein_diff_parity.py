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

"""Byte-for-byte parity harness: ``fast`` vs ``protein_diff``.

Now that ``fast`` is the default annotator (7.0.0) while ``protein_diff``
remains the substrate for the splice-outcome / germline machinery, the
two must stay reconciled — this suite is the gate that enforces it.

**What this adds over ``test_annotator_parity_adversarial.py``.** That
suite compares only the effect *class* and *short_description*. This one
additionally compares the full ``mutant_protein_sequence`` and the
frameshift/stop-loss ``shifted_sequence`` — the exact fields the
class+description comparison does NOT check, and precisely where #396
hid: ATM p.F61fs (and CFTR p.L127fs, BRCA1 p.R71fs below) matched on
class AND short_description under both annotators, yet ``protein_diff``
silently dropped the final residue of the novel C-terminus. A parity
gate that only checks class+description would not have caught it.

The corpus sweeps coding windows on both strands (CFTR +, BRCA1 -) with
SNVs / insertions / deletions, plus curated cases carrying the tricky
C-terminus shapes from #396/#397. Any disagreement must be triaged into
``EXPECTED_DIFFS`` with a reason, never silently accepted.
"""

import pytest
from pyensembl import cached_release

from varcode import Variant
from varcode.annotators.fast import FastEffectAnnotator
from varcode.annotators.protein_diff import ProteinDiffEffectAnnotator
from varcode.effects.effect_classes import FrameShift
from varcode.effects.transcript_helpers import interval_offset_on_transcript

# Selectable via `pytest -m parity`; also runs as part of the default suite.
pytestmark = pytest.mark.parity

ensembl_grch38 = cached_release(81)
CFTR_ID = "ENST00000003084"    # chr7, + strand
BRCA1_ID = "ENST00000357654"   # chr17, - strand

_FAST = FastEffectAnnotator()
_PDIFF = ProteinDiffEffectAnnotator()

# Coding windows to sweep exhaustively: (contig, transcript_id, start, end).
SWEEP_REGIONS = [
    ("7", CFTR_ID, 117531050, 117531100),    # CFTR coding exon (+ strand)
    ("17", BRCA1_ID, 43082560, 43082610),    # BRCA1 coding exon (- strand)
]

# Curated variants with known-tricky C-terminus shapes. The two
# coincidental-shared-suffix frameshifts are the #396/#397 regression:
# their novel tail ends in the reference protein's own terminal residue,
# which the pre-#396 whole-protein suffix trim silently dropped.
CURATED = [
    # (transcript_id, contig, start, ref, alt, note)
    (CFTR_ID, "7", 117531002, "GC", "G", "#396 CFTR p.L127fs (+), tail ...GYAFSLL"),
    (BRCA1_ID, "17", 43106456, "CT", "C", "#396 BRCA1 p.R71fs (-), tail ...VNLLKSY"),
]

# Documented, intentional divergences keyed by
# (transcript_id, contig, start, ref, alt) -> reason. Empty at merge
# time: any real disagreement must be triaged here, not silently
# accepted. (FrameShift vs FrameShiftTruncation is handled structurally
# by _is_known_divergence, not enumerated per-variant.)
EXPECTED_DIFFS = {}


def _is_known_divergence(fast_eff, pdiff_eff):
    """FrameShift vs FrameShiftTruncation is a documented, more-specific
    reclassification (see test_annotator_parity_adversarial.py), not a
    parity bug; the resulting protein sequences legitimately differ."""
    return ({type(fast_eff).__name__, type(pdiff_eff).__name__}
            == {"FrameShift", "FrameShiftTruncation"})


def _safe(effect, name):
    """Read an optional effect attribute without letting a lazily-computed
    property raise; missing/erroring -> None."""
    try:
        return getattr(effect, name, None)
    except Exception:
        return None


def _fingerprint(effect):
    """Full-fidelity tuple compared across annotators: class, HGVS short
    description, and the protein-level output (whole mutant protein +
    frameshift/stop-loss novel tail) that class+description alone miss."""
    return (
        type(effect).__name__,
        effect.short_description,
        _safe(effect, "mutant_protein_sequence"),
        _safe(effect, "shifted_sequence"),
    )


def _annotate(annotator, variant, transcript):
    try:
        return annotator.annotate_on_transcript(variant, transcript), None
    except Exception as exc:  # noqa: BLE001 - we compare failure symmetry
        return None, exc


def _build_corpus():
    """(label, transcript_id, variant) tuples: exhaustive coding sweeps on
    both strands plus the curated tricky cases."""
    corpus = []
    for contig, tid, start, end in SWEEP_REGIONS:
        t = ensembl_grch38.transcript_by_id(tid)
        seq = str(t.sequence)
        for pos in range(start, end):
            try:
                off = interval_offset_on_transcript(pos, pos, t)
                ref = seq[off]
            except Exception:
                continue
            for alt in ("A", "C", "G", "T"):
                if alt != ref:
                    corpus.append(("snv", tid, _mk(contig, pos, ref, alt)))
            for suffix in ("A", "AA", "AAA"):
                corpus.append(("ins", tid, _mk(contig, pos, ref, ref + suffix)))
            for del_len in (2, 3, 4):
                try:
                    off2 = interval_offset_on_transcript(pos, pos + del_len - 1, t)
                    dref = seq[off2:off2 + del_len]
                except Exception:
                    continue
                if len(dref) == del_len:
                    corpus.append(("del", tid, _mk(contig, pos, dref, dref[0])))
    for tid, contig, pos, ref, alt, _note in CURATED:
        corpus.append(("curated", tid, _mk(contig, pos, ref, alt)))
    return [entry for entry in corpus if entry[2] is not None]


def _mk(contig, pos, ref, alt):
    try:
        return Variant(contig, pos, ref, alt, ensembl_grch38)
    except Exception:
        return None


CORPUS = _build_corpus()


def test_parity_corpus_is_populated():
    """Guard against the failure mode of the old stub: a parity gate with
    an empty corpus passes vacuously. Require a substantial corpus that
    spans both strands and includes the curated #396 frameshifts."""
    assert len(CORPUS) > 500, (
        "parity corpus unexpectedly small (%d) — sweep generation broke"
        % len(CORPUS))
    tids = {tid for _label, tid, _v in CORPUS}
    assert CFTR_ID in tids and BRCA1_ID in tids, (
        "corpus must cover both strands; got %r" % tids)
    assert sum(1 for label, _t, _v in CORPUS if label == "curated") == len(CURATED)


def test_protein_diff_matches_fast_annotator_on_corpus():
    """Byte-for-byte parity (class + short_description + full mutant
    protein sequence + frameshift tail) between fast and protein_diff
    across the whole corpus."""
    failures = []
    for label, tid, variant in CORPUS:
        transcript = ensembl_grch38.transcript_by_id(tid)
        fast_eff, fast_err = _annotate(_FAST, variant, transcript)
        pdiff_eff, pdiff_err = _annotate(_PDIFF, variant, transcript)

        if fast_err or pdiff_err:
            # Only a parity gap if exactly one side raised.
            if bool(fast_err) != bool(pdiff_err):
                failures.append(
                    "%s %s %s: asymmetric error fast_err=%r pdiff_err=%r"
                    % (label, tid, variant.short_description,
                       fast_err, pdiff_err))
            continue

        key = (tid, variant.contig, variant.start, variant.ref, variant.alt)
        if key in EXPECTED_DIFFS:
            continue
        if _is_known_divergence(fast_eff, pdiff_eff):
            continue

        fast_fp = _fingerprint(fast_eff)
        pdiff_fp = _fingerprint(pdiff_eff)
        if fast_fp != pdiff_fp:
            failures.append(
                "%s %s %s:\n    fast =%r\n    pdiff=%r"
                % (label, tid, variant.short_description, fast_fp, pdiff_fp))

    assert not failures, (
        "fast vs protein_diff parity failures (%d of %d corpus entries):\n  %s"
        % (len(failures), len(CORPUS), "\n  ".join(failures[:20])))


def test_curated_frameshift_cases_exercise_full_cterminus():
    """The #396 curated cases must actually be frameshifts whose novel tail
    coincides with the reference protein's C-terminus — otherwise the pin
    would silently stop guarding the regression it was written for."""
    for tid, contig, pos, ref, alt, _note in CURATED:
        transcript = ensembl_grch38.transcript_by_id(tid)
        variant = Variant(contig, pos, ref, alt, ensembl_grch38)
        fast_eff = _FAST.annotate_on_transcript(variant, transcript)
        pdiff_eff = _PDIFF.annotate_on_transcript(variant, transcript)
        assert isinstance(fast_eff, FrameShift), (
            "%s: expected FrameShift, got %s" % (tid, type(fast_eff).__name__))
        ref_protein = transcript.protein_sequence
        mut = fast_eff.mutant_protein_sequence
        assert mut[-1] == ref_protein[-1], (
            "%s: curated case no longer shares a terminal residue with the "
            "reference protein (mut ...%s vs ref ...%s)"
            % (tid, mut[-4:], ref_protein[-4:]))
        assert fast_eff.mutant_protein_sequence == pdiff_eff.mutant_protein_sequence
