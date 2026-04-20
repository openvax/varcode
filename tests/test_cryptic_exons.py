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

"""Tests for the cryptic-exon candidate enumerator (PR 11; #252).

The enumerator is deliberately coarse — motif scoring nominates many
false positives. These tests pin the *shape* of the output and the
plug-in contract for external scorers, not biological accuracy.
Real calling is downstream: SpliceAI, RNA-seq junction support,
long-read resolution.
"""

import pytest

from varcode import (
    Outcome,
    StructuralVariant,
    enumerate_cryptic_exon_candidates,
    score_acceptor,
    score_donor,
)
from varcode.effects import CrypticExonCandidate


# --------------------------------------------------------------------
# Scoring primitives
# --------------------------------------------------------------------


def test_donor_canonical_motif_scores_high():
    # MAG|GTRAGT canonical: CAG|GTAAGT
    score = score_donor("CAGGTAAGT")
    assert score >= 0.9


def test_donor_without_gt_anchor_scores_zero():
    """Anchor GT at +1/+2 is mandatory — without it, the site is
    not a canonical 5' splice site regardless of the surrounding
    bases."""
    # Anchor broken: CAG|AAAAAAA
    assert score_donor("CAGAAAAAA") == 0.0


def test_donor_low_similarity_scores_low():
    # GT anchor present but flanks are random
    score = score_donor("TTTGTCCCC")
    # Meets GT anchor so non-zero, but low match count
    assert 0.0 < score < 0.7


def test_acceptor_canonical_motif_scores_high():
    # YAG|G canonical: CAG|G
    assert score_acceptor("CAGG") >= 0.9
    assert score_acceptor("TAGG") >= 0.9


def test_acceptor_without_ag_anchor_scores_zero():
    # Anchor AG broken
    assert score_acceptor("CCCG") == 0.0


# --------------------------------------------------------------------
# Candidate enumeration
# --------------------------------------------------------------------


def test_enumerate_finds_planted_cryptic_exon():
    """Construct a sequence with a strong acceptor, then random
    middle, then a strong donor. The enumerator should report a
    candidate spanning the middle."""
    # Strong acceptor at offset 0: TAGG
    # Middle: 50 bp of filler
    # Strong donor at offset 54: CAG|GTAAGT (9 bp window)
    acceptor = "TAGG"
    middle = "AAAA" * 12   # 48 bp, so donor lands at pos 52
    donor = "CAGGTAAGT"
    seq = acceptor + middle + donor

    candidates = enumerate_cryptic_exon_candidates(
        contig="17",
        sequence=seq,
        sequence_start=1_000_000,
        min_intron_length=10,
        max_intron_length=100,
    )
    assert len(candidates) >= 1
    top = candidates[0]
    assert isinstance(top, CrypticExonCandidate)
    assert top.contig == "17"
    # Genomic coordinates offset from sequence_start
    assert top.interval_start >= 1_000_000
    assert top.interval_end > top.interval_start


def test_enumerate_empty_when_no_motifs():
    """A sequence without canonical GT / AG anchors yields nothing."""
    seq = "A" * 500
    candidates = enumerate_cryptic_exon_candidates(
        contig="1", sequence=seq, sequence_start=100)
    assert candidates == []


def test_enumerate_short_candidates_filtered_by_min_intron():
    """Exon-length filter gates tiny candidates."""
    # Acceptor immediately followed by donor — too short.
    seq = "TAGG" + "CAGGTAAGT"
    # Default min_intron_length is 30, so this tight pair filters out.
    candidates = enumerate_cryptic_exon_candidates(
        contig="1", sequence=seq, sequence_start=100)
    assert candidates == []


def test_enumerate_respects_external_score_fn():
    """An external scorer (e.g. SpliceAI wrapper) overrides the
    motif-based default. Contract: the scorer gets the window and
    kind, returns a probability. This is how the pluggable external-
    predictor integration works without varcode depending on the
    scorer's model."""
    call_log = []

    def custom_scorer(window: str, kind: str) -> float:
        call_log.append((len(window), kind))
        # Pretend every window is a strong site, to prove the scorer
        # is consulted.
        return 1.0

    seq = "A" * 100
    candidates = enumerate_cryptic_exon_candidates(
        contig="1",
        sequence=seq,
        sequence_start=1,
        score_fn=custom_scorer,
        min_intron_length=10,
        max_intron_length=200,
    )
    # The custom scorer saw donor (9-bp) and acceptor (4-bp) windows.
    assert any(length == 9 and kind == "donor"
               for length, kind in call_log)
    assert any(length == 4 and kind == "acceptor"
               for length, kind in call_log)
    # And since every site was scored 1.0, many candidates got nominated.
    assert len(candidates) > 0


def test_cryptic_exon_candidate_wraps_as_outcome():
    """A :class:`CrypticExonCandidate` can be wrapped in an
    :class:`Outcome` for the harmonized multi-outcome interface.
    This is how an SV annotator attaches a cryptic-exon outcome
    alongside its varcode-nominated `LargeDeletion` / `GeneFusion`
    outcome."""
    cand = CrypticExonCandidate(
        variant=None,
        contig="17",
        interval_start=1_000_000,
        interval_end=1_000_120,
        donor_score=0.85,
        acceptor_score=0.90)
    outcome = Outcome(
        effect=cand,
        probability=0.87,
        source="varcode-motif",
        evidence={"donor_score": 0.85, "acceptor_score": 0.90})
    assert outcome.source == "varcode-motif"
    assert outcome.probability == 0.87
    # An external rescore (SpliceAI-style) wraps the same candidate
    # with its own probability — downstream consumers filter by source.
    rescored = Outcome(
        effect=cand,
        probability=0.98,
        source="spliceai",
        evidence={"ds_dg": 0.98, "ds_al": 0.02})
    assert rescored.source == "spliceai"


def test_enumerate_from_structural_variant_uses_alt_assembly():
    """When an SV carries a long-read-derived ``alt_assembly``, the
    convenience wrapper scans that sequence directly, bypassing the
    reference."""
    # Build a fake SV with a hand-crafted assembly that plants a
    # cryptic-exon candidate.
    assembly = "TAGG" + ("A" * 50) + "CAGGTAAGT"
    sv = StructuralVariant(
        contig="1",
        start=1_000,
        end=1_001,
        sv_type="BND",
        alt="N]?:?]",
        alt_assembly=assembly,
        genome="GRCh38",
    )
    from varcode import enumerate_from_structural_variant
    candidates = enumerate_from_structural_variant(
        sv, min_intron_length=10, max_intron_length=200)
    assert len(candidates) >= 1
    assert candidates[0].contig == "1"


def test_enumerate_from_structural_variant_without_alt_or_genome():
    """No alt_assembly + no resolvable reference = empty list, not
    crash. Graceful degradation is the design contract."""
    # Use an SV whose "genome" hook doesn't expose a sequence method
    # — a pyensembl release does, but we don't guarantee that code
    # path works in tests without the full release loaded. Instead
    # verify that the function returns cleanly on a BND without
    # assembly.
    class _FakeGenome:
        pass

    sv = StructuralVariant(
        contig="99",  # nonexistent contig
        start=1_000,
        sv_type="BND",
        genome="GRCh38",
    )
    from varcode import enumerate_from_structural_variant
    # With a nonexistent contig, reference sequence fetch fails — the
    # function should return empty, not crash.
    candidates = enumerate_from_structural_variant(sv)
    assert candidates == []
