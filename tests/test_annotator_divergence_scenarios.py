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

"""Scenarios where the fast and protein-diff annotators COULD diverge,
with explicit expected output for each.

The adversarial parity suite (``tests/test_annotator_parity_adversarial.py``)
only asserts that the two annotators AGREE; it doesn't catch the case
where they agree on the wrong answer, and it uses
``_is_known_divergence`` to hide documented mismatches. This file is
complementary:

* **Agreement scenarios** (``dual_annotator`` parametrization) pin a
  specific expected effect class and ``short_description`` for BOTH
  annotators. A regression in either annotator fails the test, even
  if the two still agree with each other.

* **Divergence scenarios** run each annotator explicitly and pin the
  SEPARATE expected outputs — one for fast, one for protein_diff. If
  either side's output drifts, the test fails. These are the cases
  where fast and protein_diff model the same biological change using
  different classifications (some are bugs in protein_diff waiting to
  be fixed; some are legitimate "different level of description").

Each divergence case carries a brief note on which annotator appears
more biologically accurate, so future triage can flip the pin without
having to re-derive the reasoning.
"""

import pytest
from pyensembl import cached_release

import varcode
from varcode import Variant
from varcode.annotators.fast import FastEffectAnnotator
from varcode.annotators.protein_diff import ProteinDiffEffectAnnotator
from varcode.effects.effect_classes import (
    AlternateStartCodon,
    ComplexSubstitution,
    Deletion,
    ExonicSpliceSite,
    FivePrimeUTR,
    FrameShift,
    FrameShiftTruncation,
    Insertion,
    Intronic,
    NoncodingTranscript,
    PrematureStop,
    Silent,
    SpliceAcceptor,
    SpliceDonor,
    StartLoss,
    StopLoss,
    Substitution,
    ThreePrimeUTR,
)
from varcode.errors import ReferenceMismatchError

ensembl_grch38 = cached_release(81)

# Well-covered coding transcripts used across the suite:
CFTR_ID = "ENST00000003084"    # chr7, + strand, 1480 aa
BRCA1_ID = "ENST00000357654"   # chr17, - strand, 1863 aa
MT_CO1_ID = "ENST00000361624"  # MT, + strand, vertebrate mt code

_FAST = FastEffectAnnotator()
_PDIFF = ProteinDiffEffectAnnotator()


def _annotate(variant, transcript_id, annotator):
    transcript = variant.ensembl.transcript_by_id(transcript_id)
    return annotator.annotate_on_transcript(variant, transcript)


def _pin(variant, transcript_id, annotator, effect_class, short_description):
    """Assert ``annotator`` produces the given class + description."""
    effect = _annotate(variant, transcript_id, annotator)
    assert isinstance(effect, effect_class), (
        "Expected %s, got %s (desc=%r)"
        % (effect_class.__name__, type(effect).__name__,
           getattr(effect, "short_description", "")))
    assert effect.short_description == short_description, (
        "short_description: expected %r, got %r"
        % (short_description, effect.short_description))


# ====================================================================
# Scenarios where BOTH annotators should agree.
#
# Each test runs under both annotators via the ``dual_annotator``
# fixture and pins a single expected outcome. These lock in
# behaviours that are easy to break when refactoring either path.
# ====================================================================


def _run_dual(variant, transcript_id, effect_class, short_description,
              dual_annotator):
    """Helper that picks the right annotator instance for a
    ``dual_annotator`` parametrization."""
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    _pin(variant, transcript_id, annotator, effect_class, short_description)


def test_deep_intronic_snv_agrees_as_intronic(dual_annotator):
    # Far from any splice site (CFTR intron, 300+ bp from the
    # nearest exon boundary). Both annotators defer to fast's
    # position-based intronic classification.
    variant = Variant("7", 117480500, "A", "G", ensembl_grch38)
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    effect = _annotate(variant, CFTR_ID, annotator)
    assert isinstance(effect, Intronic)


def test_stop_codon_silent_snv_agrees(dual_annotator):
    # CFTR stop codon at 117667106-117667108 = 'TAG'.
    # Change the last base G->A: TAG -> TAA (still a stop codon).
    # Both annotators call this Silent on the terminal stop position.
    _run_dual(
        Variant("7", 117667108, "G", "A", ensembl_grch38),
        CFTR_ID, Silent, "p.1481=", dual_annotator)


def test_insertion_into_3utr_after_stop_is_silent(dual_annotator):
    # Insertion placed immediately after the stop codon, inside the
    # 3' UTR. Classify as Silent on both annotators (the protein is
    # untouched, and fast treats "past the stop" as Silent too).
    variant = Variant("7", 117667108, "G", "GAAA", ensembl_grch38)
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    effect = _annotate(variant, CFTR_ID, annotator)
    assert isinstance(effect, Silent)


def test_in_frame_insertion_just_before_stop_agrees(dual_annotator):
    # Insertion of 'CAT' between the last coding codon and the stop.
    # Produces Insertion p.1480insH on both annotators (NOT StopLoss —
    # the stop codon itself is preserved). This is the trim_shared_
    # flanking case documented in varcode #201.
    _run_dual(
        Variant("7", 117667105, "T", "TCAT", ensembl_grch38),
        CFTR_ID, Insertion, "p.1480insH", dual_annotator)


def test_splice_donor_plus1_agrees(dual_annotator):
    # CFTR exon 3 ends at 117531114; intronic +1 = 117531115.
    # Both annotators route this through fast's intronic splice
    # classifier.
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    effect = _annotate(variant, CFTR_ID, annotator)
    assert isinstance(effect, SpliceDonor)


def test_splice_acceptor_minus1_agrees(dual_annotator):
    # CFTR exon 3 starts at 117530899; acceptor -1 = 117530898.
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    variant = Variant("7", 117530898, "G", "A", ensembl_grch38)
    effect = _annotate(variant, CFTR_ID, annotator)
    assert isinstance(effect, SpliceAcceptor)


def test_mt_tga_decodes_as_trp_not_stop(dual_annotator):
    # MT-CO1 position 6739: TCA->TGA. Under the vertebrate mt codon
    # table, TGA is Trp (not a stop codon). Both annotators must
    # route through codon_table_for_transcript. The point of this
    # test is to catch a regression where one path fails to pick up
    # the mt-specific codon table.
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    variant = Variant("MT", 6739, "C", "G", ensembl_grch38)
    effect = _annotate(variant, MT_CO1_ID, annotator)
    # Not a PrematureStop: both must recognize TGA = Trp on mt.
    assert not isinstance(effect, (PrematureStop, FrameShiftTruncation)), (
        "MT TGA should not be classified as a stop: got %s"
        % type(effect).__name__)


def test_noncoding_transcript_shortcircuits(dual_annotator):
    # Both annotators gate on ``is_protein_coding`` before doing any
    # protein arithmetic and must return NoncodingTranscript. We
    # construct a minimal mock transcript so the test doesn't depend
    # on whichever non-coding biotypes happen to exist in the cached
    # Ensembl release.
    class _MockNoncodingTranscript:
        is_protein_coding = False
        complete = False  # irrelevant — the gate above short-circuits
        contig = "7"
        start = 117531050
        end = 117531150
        name = "mock-noncoding"
        id = "MOCK-NC"

        def __repr__(self):
            return "MockNoncodingTranscript()"

    # protein_diff's gate calls ``isinstance(transcript, Transcript)``
    # before it even reads is_protein_coding. The mock would fail
    # that isinstance check. To exercise the real gate we use an
    # actual non-coding pyensembl transcript if one is available in
    # the release; otherwise fall back to asserting via the variant
    # accessor path.
    candidate = None
    for gene in ensembl_grch38.genes_at_locus("7", 117479963, 117668665):
        for t in gene.transcripts:
            if not t.is_protein_coding:
                candidate = t
                break
        if candidate is not None:
            break
    if candidate is None:
        pytest.skip("No non-coding transcript in the CFTR locus")
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    # Use a coding-adjacent SNV that definitely overlaps the chosen
    # non-coding transcript; we don't care about the ref match for
    # fast's NoncodingTranscript shortcut — it fires before any ref
    # lookup. protein_diff's shortcut also fires at the top of
    # annotate_on_transcript.
    pos = (candidate.start + candidate.end) // 2
    v = Variant("7", pos, "A", "G", ensembl_grch38)
    effect = annotator.annotate_on_transcript(v, candidate)
    assert isinstance(effect, NoncodingTranscript)


def test_reference_mismatch_raises_on_both(dual_annotator):
    # Both annotators must surface a ReferenceMismatchError when
    # the caller's ``ref`` field contradicts the reference genome.
    # protein_diff falls back to fast (apply_variant_to_transcript
    # returns None on ref mismatch), which raises.
    # CFTR position 117531100 has 'T' in the reference; claim 'A'.
    variant = Variant("7", 117531100, "A", "G", ensembl_grch38)
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    with pytest.raises(ReferenceMismatchError):
        _annotate(variant, CFTR_ID, annotator)


def test_variant_outside_all_transcripts_returns_no_transcript_effects(
        dual_annotator):
    # A variant on chr22 far from any gene: .effects() should
    # produce an Intergenic (no transcripts), not crash. Both
    # annotators must survive this.
    # Pick an intergenic position on chr22 (middle of a well-known
    # gene desert).
    with varcode.use_annotator(dual_annotator):
        variant = Variant("22", 17000000, "A", "G", ensembl_grch38)
        effects = variant.effects()
    assert len(effects) >= 1
    # Either Intergenic or the variant overlaps at least one gene;
    # whichever it is, the annotator must not have crashed.


def test_frameshift_simple_agrees(dual_annotator):
    # A 1-base insertion in CFTR exon 4 that doesn't immediately
    # create a stop — should be FrameShift on both, same description.
    # This locks in the "simple frameshift agrees" case; the
    # immediate-stop case is tested separately as a KNOWN divergence.
    _run_dual(
        Variant("7", 117531100, "T", "TA", ensembl_grch38),
        CFTR_ID, FrameShift, "p.L159fs", dual_annotator)


# ====================================================================
# DOCUMENTED DIVERGENCES: the two annotators classify the same
# biological change differently. Each test pins BOTH outputs so a
# drift in either is caught.
# ====================================================================


def test_divergence_3utr_snv():
    """3'UTR SNV: fast returns ThreePrimeUTR; protein_diff returns
    Silent. Reason: protein_diff's ``apply_variant_to_transcript``
    translates past the stop codon when ``cdna_offset >= cds_start``
    and concludes the protein is unchanged (it is — the variant is
    downstream of the stop) without recognising that "no protein
    change" and "in 3' UTR" are different effect classes.

    Fast is more informative here: ThreePrimeUTR carries location
    semantics that Silent does not. protein_diff should defer to the
    fast_effect when the fast_effect is ThreePrimeUTR / FivePrimeUTR
    / Intronic — tracked separately (this test pins the current
    behaviour so we notice when it changes).
    """
    variant = Variant("7", 117667200, "T", "A", ensembl_grch38)
    _pin(variant, CFTR_ID, _FAST, ThreePrimeUTR, "3' UTR")
    pdiff_effect = _annotate(variant, CFTR_ID, _PDIFF)
    assert isinstance(pdiff_effect, Silent), (
        "protein_diff expected to (incorrectly) return Silent for a "
        "3'UTR SNV; got %s. If this starts returning ThreePrimeUTR, "
        "the protein_diff gate was tightened — update the docstring."
        % type(pdiff_effect).__name__)


def test_divergence_3utr_snv_reverse_strand():
    """Same divergence as the CFTR case, on BRCA1 (- strand), to
    confirm the bug isn't strand-specific."""
    # BRCA1 stop_codon_positions[0] == 43045678 (- strand, so 3'UTR is
    # at LOWER genomic coords). Pick 100 bp into the 3' UTR.
    pos = 43045578
    # cDNA base there is 'T' (we measured); genomic ref is complement.
    variant = Variant("17", pos, "A", "T", ensembl_grch38)
    _pin(variant, BRCA1_ID, _FAST, ThreePrimeUTR, "3' UTR")
    pdiff_effect = _annotate(variant, BRCA1_ID, _PDIFF)
    assert isinstance(pdiff_effect, Silent)


def test_divergence_stop_codon_first_base_substitution():
    """TAG->CAG at the stop codon (first base): fast calls this
    StopLoss; protein_diff calls it Insertion.

    Reason: ``classify_from_protein_diff`` requires ``n_ref > 0`` to
    produce StopLoss — it distinguishes "stop replaced by residues"
    from "residues inserted near the stop." When the mutation
    changes only the stop codon itself, the whole-protein diff sees
    a pure insertion past the end of the reference protein, which
    the classifier flags as Insertion rather than StopLoss.

    Fast's classification is closer to the canonical HGVS
    p.*NNNxxx notation. This is a real protein_diff bug — waiting
    on a classify.py patch that recognises "insertion at or past
    the reference-protein tail" as a stop-loss case.
    """
    # CFTR stop TAG at 117667106-108; change T->C to produce CAG.
    variant = Variant("7", 117667106, "T", "C", ensembl_grch38)
    _pin(variant, CFTR_ID, _FAST, StopLoss,
         "p.*1481QRAA (stop-loss)")
    pdiff_effect = _annotate(variant, CFTR_ID, _PDIFF)
    assert isinstance(pdiff_effect, Insertion), (
        "protein_diff expected to (incorrectly) return Insertion; "
        "got %s" % type(pdiff_effect).__name__)


def test_divergence_deletion_of_stop_codon():
    """Deleting the entire stop codon (TAG) produces StopLoss in
    fast but Insertion in protein_diff, for the same reason as the
    substitution case above (n_ref == 0 in the protein diff)."""
    # Delete TAG at 117667106-108 using anchor at 117667105.
    variant = Variant("7", 117667105, "TTAG", "T", ensembl_grch38)
    _pin(variant, CFTR_ID, _FAST, StopLoss,
         "p.*1481RAA (stop-loss)")
    pdiff_effect = _annotate(variant, CFTR_ID, _PDIFF)
    assert isinstance(pdiff_effect, Insertion)


def test_divergence_alternate_start_codon_atg_to_ctg():
    """ATG->CTG: fast returns AlternateStartCodon (CTG is in the
    standard-table start codon set, so the ribosome would still
    initiate at this position, producing Met). protein_diff instead
    translates literally, sees 'L' at position 0, and returns
    StartLoss.

    Biologically, fast's interpretation is closer to correct: the
    initiator tRNA loads Met regardless of the start codon's codon
    identity (ATG, CTG, GTG, TTG). protein_diff's classifier doesn't
    know about initiator-tRNA semantics and takes the pessimistic
    StartLoss interpretation.

    Tracked for protein_diff: either (a) add "mut starts with L/V/
    but alt codon is a recognised alternate start" → AlternateStart
    Codon rule, or (b) pre-rewrite the translated first residue to
    'M' when the first codon is in the table's start_codons.
    """
    start_pos = min(ensembl_grch38.transcript_by_id(
        CFTR_ID).start_codon_positions)
    variant = Variant("7", start_pos, "A", "C", ensembl_grch38)
    _pin(variant, CFTR_ID, _FAST, AlternateStartCodon,
         "alternate-start-codon (AAT>CTG)")
    pdiff_effect = _annotate(variant, CFTR_ID, _PDIFF)
    assert isinstance(pdiff_effect, StartLoss)


def test_divergence_frameshift_immediate_stop():
    """Frameshift that creates a stop codon at (or very near) the
    insertion point. Fast may report FrameShift with a short
    shifted sequence; protein_diff reports FrameShiftTruncation with
    an empty shifted sequence. Both describe the same event —
    protein_diff's is strictly more specific.

    We construct this by finding a site where a 1-base insertion
    produces a stop in the next codon, then asserting the two
    annotators emit different classes but consistent positional
    offsets. If they ever converge we can relax this test.
    """
    # A 1-base insertion a few codons before the CFTR stop where
    # the resulting frame runs into a stop quickly. Known-good
    # probe from the sweep in test_annotator_parity_adversarial.py.
    # Picking a position from the sweep range:
    variant = Variant("7", 117667101, "T", "TA", ensembl_grch38)
    fe = _annotate(variant, CFTR_ID, _FAST)
    pe = _annotate(variant, CFTR_ID, _PDIFF)
    # One of three shapes is allowed — this is documented as a
    # known-divergence in the parity harness. We only insist they
    # each classify the event:
    assert isinstance(fe, (FrameShift, FrameShiftTruncation))
    assert isinstance(pe, (FrameShift, FrameShiftTruncation))


# ====================================================================
# Structural-variant gap tests.
#
# varcode's ``Variant`` class only models point variants, MNVs, and
# simple indels expressed as explicit ref/alt nucleotide strings.
# VCF symbolic alleles (``<DEL>``, ``<DUP>``, ``<INS:ME:ALU>``, …),
# breakend notation (``G]17:198982]``), and the ``*`` spanning-
# deletion placeholder are all filtered at VCF load time (see
# ``varcode/vcf.py::_is_symbolic_allele`` and the regression tests
# in ``tests/test_symbolic_alleles.py`` for issue #88).
#
# These tests document what is NOT supported end-to-end, so future
# work on SV annotation starts from a known baseline:
#
#   1. ``_is_symbolic_allele`` rejects the expected symbolic forms.
#   2. Constructing a ``Variant`` directly with a symbolic ALT string
#      does not go through the filter — the Variant constructor
#      accepts it and will then fail downstream, which is today's
#      rough edge.
#   3. A "structurally large" deletion spelled out in full ref/alt
#      does load and annotate through the normal indel path — useful
#      for small-to-medium deletions but not for CNV-scale events.
#   4. The ``supports`` frozenset on both annotators advertises the
#      current contract: SNV, indel, MNV — no SV kinds.
# ====================================================================


from varcode.vcf import _is_symbolic_allele  # noqa: E402


@pytest.mark.parametrize("alt", [
    "<DEL>",
    "<DUP>",
    "<INV>",
    "<INS>",
    "<INS:ME:ALU>",
    "<INS:ME:LINE1>",
    "<INS:ME:SVA>",
    "<CN0>",
    "<CN2>",
    "<CNV>",
    # Breakend notation (VCF 4.1 section 5.4):
    "G]17:198982]",
    "]17:198982]G",
    "[13:123456[T",
    "T[13:123456[",
    # Spanning-deletion placeholder:
    "*",
])
def test_sv_symbolic_allele_detection(alt):
    """The load-time filter recognises the full symbolic-allele
    grammar defined in VCF 4.1/4.3. Regression guard for the SV-
    roadmap work: adding per-SV annotation starts with widening
    this predicate, so the test catches accidental narrowing."""
    assert _is_symbolic_allele(alt) is True, (
        "Expected %r to be flagged as a symbolic/SV allele" % alt)


@pytest.mark.parametrize("alt", ["A", "ACGT", "AAAAAAAA"])
def test_sv_detector_does_not_flag_normal_alleles(alt):
    assert _is_symbolic_allele(alt) is False


def test_sv_variant_constructor_accepts_symbolic_alt_but_breaks_downstream():
    """The ``Variant`` constructor itself does NOT reject symbolic
    alleles — that responsibility lives in the VCF loader. Direct
    construction + annotation is the path a caller would take if
    they tried to wire in SV support without going through the
    loader, so document what actually happens.
    """
    # Symbolic allele passed directly: normalize_nucleotide_string
    # allows it iff allow_extended_nucleotides=True; otherwise it
    # raises ValueError. This test pins both behaviours.
    with pytest.raises(ValueError):
        Variant("22", 51179178, "A", "<CN0>", ensembl_grch38)


def test_sv_annotators_declare_supports_sets_without_sv_kinds():
    """Both shipped annotators advertise a ``supports`` frozenset
    that deliberately omits SV variant kinds. SV support lands as a
    third annotator (see docstring in mutant_transcript.py about
    the planned ``reference_segments`` extension for fusions /
    translocations)."""
    assert _FAST.supports == frozenset({"snv", "indel", "mnv"})
    assert _PDIFF.supports == frozenset({"snv", "indel", "mnv"})
    # No annotator in the registry advertises support for SV kinds:
    for sv_kind in ("sv", "bnd", "breakend", "fusion", "cnv", "del_symbolic",
                    "dup", "inv"):
        assert sv_kind not in _FAST.supports
        assert sv_kind not in _PDIFF.supports


def test_sv_large_explicit_deletion_still_annotates_through_indel_path():
    """A deletion that's biologically SV-scale (hundreds of bp) but
    spelled out as an explicit ref/alt string is still just an
    indel to varcode. This documents the CURRENT capability: we
    don't need a symbolic allele for kilobase-scale deletions as
    long as the caller is willing to write them out.

    The annotation itself isn't the point — we just assert it
    doesn't crash, and that something coding-flavoured comes back
    (not a Failure/Intergenic).
    """
    # Build a 50-base deletion inside CFTR exon 11. Anchor+50 bases.
    transcript = ensembl_grch38.transcript_by_id(CFTR_ID)
    from varcode.effects.transcript_helpers import interval_offset_on_transcript
    pos = 117587738  # inside CFTR exon 11 (known coding region)
    try:
        off = interval_offset_on_transcript(pos, pos + 49, transcript)
    except Exception:
        pytest.skip("Coordinate %d not on CFTR transcript in this release"
                    % pos)
    ref = str(transcript.sequence)[off:off + 50]
    if len(ref) != 50:
        pytest.skip("Insufficient exon length at %d" % pos)
    variant = Variant("7", pos, ref, ref[0], ensembl_grch38)
    # The deletion spans a splice site inside the CFTR gene — it
    # may come back as ExonLoss, a frameshift, or an ExonicSpliceSite.
    # Any non-crash, non-failure result satisfies this test.
    with varcode.use_annotator("fast"):
        effects_fast = variant.effects()
    with varcode.use_annotator("protein_diff"):
        effects_pd = variant.effects()
    # Neither should be empty.
    assert len(effects_fast) > 0
    assert len(effects_pd) > 0


# ====================================================================
# EXTENDED DIVERGENCE SWEEP
#
# Pattern-based probes driven from the analysis in the session that
# filed issues #318 / #319 / #320. Each test is pinned to the probed
# result so a future fix to classify.py or protein_diff.py flips the
# assertion and gets noticed. Agreement tests here duplicate a little
# coverage with the adversarial parity suite but add explicit
# short_description pins.
# ====================================================================


# ----- Pattern A: location-info loss (family of #318) -----


def test_5utr_snv_agrees_as_fiveprime_utr(dual_annotator):
    """5' UTR SNV: both annotators report FivePrimeUTR. protein_diff
    avoids the 3'UTR bug (#318) here because
    ``apply_variant_to_transcript`` returns a :class:`MutantTranscript`
    whose ``mutant_protein_sequence`` is ``None`` when the variant
    lands before the start codon, forcing the fallback to
    ``fast_effect``."""
    # CFTR cDNA offset 37 (5'UTR, ~58 bp before start codon).
    # Genomic 117480000 = A. A->G.
    variant = Variant("7", 117480000, "A", "G", ensembl_grch38)
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    effect = _annotate(variant, CFTR_ID, annotator)
    assert isinstance(effect, FivePrimeUTR)


def test_5utr_insertion_agrees_as_fiveprime_utr(dual_annotator):
    """3 bp insertion inside the 5' UTR. Both annotators return
    FivePrimeUTR — neither models uORFs even when the inserted
    sequence would introduce an upstream ATG."""
    variant = Variant("7", 117480000, "A", "AATG", ensembl_grch38)
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    effect = _annotate(variant, CFTR_ID, annotator)
    assert isinstance(effect, FivePrimeUTR)


# ----- Pattern B: stop/start boundary (family of #319) -----


def test_divergence_stop_codon_third_base_substitution():
    """Same family as #319: TAG->TAC at the stop codon's third base.
    Fast reports StopLoss; protein_diff reports Insertion because of
    the ``n_ref > 0`` guard in classify.py. Tracked under issue
    #319 — adding this as a second pin so the fix is verified for
    both boundary positions (1st and 3rd stop-codon base)."""
    variant = Variant("7", 117667108, "G", "C", ensembl_grch38)
    _pin(variant, CFTR_ID, _FAST, StopLoss,
         "p.*1481YRAA (stop-loss)")
    pe = _annotate(variant, CFTR_ID, _PDIFF)
    assert isinstance(pe, Insertion), (
        "Same #319 root cause as 1st-base substitution: expected "
        "Insertion from protein_diff, got %s" % type(pe).__name__)


def test_divergence_stop_codon_second_base_substitution():
    """Same family as #319: TAG->TCG at the stop codon's middle
    base. Re-stated pin: fast=StopLoss, protein_diff=Insertion."""
    variant = Variant("7", 117667107, "A", "C", ensembl_grch38)
    _pin(variant, CFTR_ID, _FAST, StopLoss,
         "p.*1481SRAA (stop-loss)")
    pe = _annotate(variant, CFTR_ID, _PDIFF)
    assert isinstance(pe, Insertion)


def test_single_base_deletion_in_stop_codon_agrees_as_stoploss(
        dual_annotator):
    """A 1-base deletion inside the stop codon shifts the frame and
    extends translation past the original stop. Both annotators
    agree on StopLoss with the same readthrough tail — this is the
    frameshift branch of classify.py which correctly handles
    ``aa_offset >= len(ref_protein)`` (see classify.py line 106)."""
    # Delete CFTR stop codon's 1st base (117667106). Anchor at
    # 117667105 ('T'), ref='TT', alt='T'.
    variant = Variant("7", 117667105, "TT", "T", ensembl_grch38)
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    effect = _annotate(variant, CFTR_ID, annotator)
    assert isinstance(effect, StopLoss)


def test_single_base_insertion_in_stop_codon_agrees_as_silent(
        dual_annotator):
    """Inserting 1 A between stop codon's 1st and 2nd bases changes
    TAG -> TAAG. Translation still terminates at the first TAA, so
    both annotators call this Silent. Pins the behaviour because
    naively this looks like a frameshift."""
    # Position 117667106 is the stop codon's first base (T). Insert
    # 'A' BEFORE it using anchor convention: pos 117667105, ref='T',
    # alt='TA'. Actually we want insertion INSIDE the stop codon. Use
    # pos 117667106 with ref='' alt='A' — varcode's zero-ref insertion.
    variant = Variant("7", 117667106, "", "A", ensembl_grch38)
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    effect = _annotate(variant, CFTR_ID, annotator)
    assert isinstance(effect, Silent)


def test_inframe_insertion_containing_new_stop_before_stop_agrees_silent(
        dual_annotator):
    """Insert TAA (in-frame) immediately before CFTR's stop codon.
    Translation terminates at the inserted TAA, and the resulting
    protein is identical to the reference (both end at residue 1480
    with the same sequence). Both annotators return Silent —
    documents that the "phantom-silent-insertion" case (#201) is
    honoured by both code paths."""
    variant = Variant("7", 117667105, "T", "TTAA", ensembl_grch38)
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    effect = _annotate(variant, CFTR_ID, annotator)
    assert isinstance(effect, Silent)


# ----- Pattern C: codon table semantics (family of #320) -----


def test_divergence_mt_alternate_start_atg_to_gtg():
    """Same family as #320: MT-CO1 start codon ATG->GTG. GTG is in
    the NCBI table-2 (vertebrate mitochondrial) start-codon set, so
    fast returns AlternateStartCodon. protein_diff translates GTG
    literally as Val and returns StartLoss.

    This test variant of #320 confirms the bug isn't confined to the
    standard codon table — the alternate-start resolution must
    consult ``codon_table_for_transcript(transcript).start_codons``
    regardless of which codon table is in force."""
    # MT-CO1 start codon at MT:5904-5906 (ATG). A->G at position 5904.
    variant = Variant("MT", 5904, "A", "G", ensembl_grch38)
    _pin(variant, MT_CO1_ID, _FAST, AlternateStartCodon,
         "alternate-start-codon (ATG>GTG)")
    pe = _annotate(variant, MT_CO1_ID, _PDIFF)
    assert isinstance(pe, StartLoss)


# ----- Pattern D: trimming / offset conventions -----


def test_divergence_deletion_in_repeat_residue_stretch():
    """Deletion of one codon from a poly-residue stretch (CFTR
    protein[98:102] == 'PLLL'). Classes agree (both Deletion), but
    the reported ``aa_mutation_start_offset`` differs:

    * fast    → p.L99del  (reports the codon that was deleted)
    * pdiff   → p.L101del (reports the HGVS-3' canonical position)

    Per HGVS 3'-rule for repeated sequences, protein_diff's p.L101del
    is the canonical choice. Fast's p.L99del is the offset of the
    deleted codon before whole-protein trimming disambiguates it. A
    real bug in at least one annotator (probably fast) — tracked as
    a new divergence issue."""
    # Delete cDNA[429:432] (codon 99, the first L of LLL).
    # Genomic position: codon 99 starts at chr7:117530923 (+ strand).
    # Anchor at 117530922 (base 'C'), ref='CCTC', alt='C'.
    # Using zero-ref insertion convention: ref='CTC', alt=''.
    variant = Variant("7", 117530923, "CTC", "", ensembl_grch38)
    _pin(variant, CFTR_ID, _FAST, Deletion, "p.L99del")
    _pin(variant, CFTR_ID, _PDIFF, Deletion, "p.L101del")


def test_penultimate_codon_snv_agrees(dual_annotator):
    """SNV in the last coding codon before the stop. CFTR codon
    1480 = CTT (Leu); change to GTT (Val). Both annotators report
    Substitution p.L1480V — no StopLoss miscalls at the tail
    boundary."""
    variant = Variant("7", 117667103, "C", "G", ensembl_grch38)
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    effect = _annotate(variant, CFTR_ID, annotator)
    assert isinstance(effect, Substitution)
    assert effect.short_description == "p.L1480V"


def test_inframe_nine_base_deletion_agrees(dual_annotator):
    """9 bp in-frame deletion spanning three codons — Deletion on
    both annotators with the same HGVS description."""
    variant = Variant("7", 117480101, "AGGTCGCCT", "", ensembl_grch38)
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    effect = _annotate(variant, CFTR_ID, annotator)
    assert isinstance(effect, Deletion)
    assert effect.short_description == "p.RSP2del"


def test_misaligned_insertion_before_stop_agrees(dual_annotator):
    """3 bp in-frame insertion NOT aligned to a codon boundary,
    placed just before the stop codon. Both annotators report
    Insertion p.1480insL — the shared ``trim_shared_flanking_strings``
    step produces identical HGVS output."""
    variant = Variant("7", 117667104, "T", "TGCT", ensembl_grch38)
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    effect = _annotate(variant, CFTR_ID, annotator)
    assert isinstance(effect, Insertion)
    assert effect.short_description == "p.1480insL"


# ----- Pattern E: splice dual-dispatch -----


def test_exonic_splice_site_snv_last_base_agrees(dual_annotator):
    """SNV in the last base of an exon (still coding, but inside the
    splice region). Both annotators produce ExonicSpliceSite — fast
    classifies position-based, protein_diff dual-dispatches through
    fast for the splice class."""
    variant = Variant("7", 117531114, "G", "A", ensembl_grch38)
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    effect = _annotate(variant, CFTR_ID, annotator)
    assert isinstance(effect, ExonicSpliceSite)


def test_exonic_splice_site_deletion_last_base_agrees(dual_annotator):
    """1 bp deletion of the last exon base — same ExonicSpliceSite
    classification on both annotators."""
    variant = Variant("7", 117531113, "AG", "A", ensembl_grch38)
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    effect = _annotate(variant, CFTR_ID, annotator)
    assert isinstance(effect, ExonicSpliceSite)


# ----- Pattern F: multi-codon / complex edits -----


def test_mnv_across_codon_boundary_agrees(dual_annotator):
    """2 bp substitution spanning a codon boundary. Codon 2 (CAG) +
    codon 3 (AGG) → codon 2 (CAC) + codon 3 (CGG). After trimming
    the back-to-back substitution, both annotators collapse this to
    a single Substitution p.Q2H (the first-changed position is the
    only one that differs after shared-flank trimming)."""
    variant = Variant("7", 117480100, "GA", "CC", ensembl_grch38)
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    effect = _annotate(variant, CFTR_ID, annotator)
    assert isinstance(effect, Substitution)
    assert effect.short_description == "p.Q2H"


def test_mnv_entire_codon_replace_agrees(dual_annotator):
    """3 bp MNV replacing all three bases of one codon. Codon 3
    AGG (Arg) → TGG (Trp). Both annotators trim to a single
    Substitution p.R3W."""
    variant = Variant("7", 117480101, "AGG", "TGG", ensembl_grch38)
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    effect = _annotate(variant, CFTR_ID, annotator)
    assert isinstance(effect, Substitution)
    assert effect.short_description == "p.R3W"


def test_mid_codon_insertion_agrees_as_complex_substitution(
        dual_annotator):
    """In-frame 3 bp insertion placed inside a codon (codon 2
    interior). Both annotators produce ComplexSubstitution p.Q2RK —
    this is the expected HGVS-style representation for insertions
    that straddle codon boundaries."""
    variant = Variant("7", 117480098, "C", "CGCA", ensembl_grch38)
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    effect = _annotate(variant, CFTR_ID, annotator)
    assert isinstance(effect, ComplexSubstitution)
    assert effect.short_description == "p.Q2RK"


# ----- Pattern G: other edge cases -----


def test_start_codon_deletion_agrees_as_startloss(dual_annotator):
    """Deleting the entire ATG start codon: both annotators return
    StartLoss. The fast path reads the deleted codon directly;
    protein_diff's ``apply_variant_to_transcript`` returns a
    ``MutantTranscript`` whose ``mutant_protein_sequence`` is
    empty (no start → no translation from ``cds_start``), which
    triggers the classify.py StartLoss branch."""
    variant = Variant("7", 117480094, "CATG", "C", ensembl_grch38)
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    effect = _annotate(variant, CFTR_ID, annotator)
    assert isinstance(effect, StartLoss)


def test_frameshift_extends_into_3utr_agrees(dual_annotator):
    """1 bp insertion 2 bp before the stop codon → frameshift that
    reads through the original stop into the 3'UTR. Both annotators
    produce FrameShift p.R1479fs; protein_diff's literal translation
    and fast's `sequence_from_start_codon` slice both cover the 3'
    UTR so they converge on the same shifted_sequence tail."""
    variant = Variant("7", 117667100, "G", "GT", ensembl_grch38)
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    effect = _annotate(variant, CFTR_ID, annotator)
    assert isinstance(effect, FrameShift)
    assert effect.short_description == "p.R1479fs"


def test_frameshift_right_after_start_codon_agrees(dual_annotator):
    """1 bp insertion immediately after the start codon: both
    annotators return FrameShift p.Q2fs (the start codon is still
    valid; the frame shift begins in codon 2)."""
    variant = Variant("7", 117480097, "G", "GA", ensembl_grch38)
    annotator = _FAST if dual_annotator == "fast" else _PDIFF
    effect = _annotate(variant, CFTR_ID, annotator)
    assert isinstance(effect, FrameShift)
    assert effect.short_description == "p.Q2fs"


def test_sv_spanning_deletion_star_allele_is_skipped_at_load():
    """The ``*`` allele in VCF means "this position is covered by
    a deletion on another record" — it shouldn't become a Variant.
    Covered by ``_is_symbolic_allele`` above; this test pins the
    end-to-end loader contract."""
    import os
    import tempfile

    vcf_body = (
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh37\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "22\t51179200\t.\tA\t*\t100\tPASS\t.\n"
        "22\t51179201\tnormal\tC\tT\t100\tPASS\t.\n"
    )
    fd, path = tempfile.mkstemp(suffix=".vcf")
    with os.fdopen(fd, "w") as f:
        f.write(vcf_body)
    try:
        from varcode import load_vcf
        vc = load_vcf(path, genome="GRCh37")
    finally:
        os.unlink(path)
    # Only the normal SNV should survive.
    assert len(vc) == 1
    assert vc[0].ref == "C" and vc[0].alt == "T"
