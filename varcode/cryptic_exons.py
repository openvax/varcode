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

"""Candidate cryptic-exon enumerator (PR 11; #252).

When a structural variant rearranges genomic sequence (translocation,
large inversion, deletion that brings intronic sequence into range
of a neighboring exon), the rearranged allele can contain novel
splice-site motifs that create *cryptic exons* — exons that the
reference transcript never had. varcode's job is to enumerate the
plausible candidates so downstream tools can rank them.

Approach
--------

**Varcode internal**: pure-Python consensus motif scoring. No ML
models, no external dependencies. The splice-site consensus motifs
are well characterized (Mount 1982, Shapiro & Senapathy 1987):

* 5' splice site (donor): exonic ``MAG`` + intronic ``GTRAGT`` —
  9 bp window, strict anchor on intronic ``GT``.
* 3' splice site (acceptor): intronic ``YAG`` + exonic ``G/A`` —
  4 bp core window (polypyrimidine tract scored separately in
  future refinements).

Scoring is a simple position-match ratio against the consensus
residues. This is crude — it will nominate many false-positive
candidates — but it's the right *candidate enumerator* for downstream
tools to score properly.

**External tools** (optional): callers that have SpliceAI / Pangolin
/ SQUIRLS outputs can bypass the motif scorer entirely. Pass a
custom ``score_fn`` to :func:`enumerate_candidates`, or take the
output of :func:`enumerate_candidates` and re-score each
:class:`CrypticExonCandidate` — the :class:`~varcode.Outcome`
shape from #299 carries both the effect and the external
``probability`` / ``evidence`` without needing varcode to understand
the scorer.

**RNA evidence** (optional): if the pipeline has RNA-seq, the
correct discriminator is "are there split reads that span the
candidate junction." Isovar-style tools add
``Outcome(source="isovar", evidence={"junction_reads": N})``
on top of the varcode-nominated candidate.

**Long-read assembly** (optional): a resolved rearranged allele
assembled from HiFi / ONT reads can be passed as the ``sequence``
argument directly, grounding the candidates in actual molecules
rather than inferred reference+breakpoint reconstructions.
"""

from typing import Callable, List, Optional, Tuple

from .effects.effect_classes import CrypticExonCandidate


# --------------------------------------------------------------------
# Consensus position-match models (Mount 1982, Shapiro & Senapathy
# 1987). Each position lists the nucleotides that count as a match;
# the score is matches / len(motif).
# --------------------------------------------------------------------


# 5' splice site (donor): 9-bp window = 3 exonic + 6 intronic.
# Consensus: ``MAG | GTRAGT``  (| = exon/intron boundary)
_DONOR_CONSENSUS: Tuple[frozenset, ...] = (
    frozenset({"A", "C"}),      # pos -3 exon (M = A/C)
    frozenset({"A"}),           # pos -2 exon
    frozenset({"G"}),           # pos -1 exon
    frozenset({"G"}),           # pos +1 intron (strong anchor)
    frozenset({"T"}),           # pos +2 intron (strong anchor)
    frozenset({"A", "G"}),      # pos +3 intron (R = A/G)
    frozenset({"A"}),           # pos +4 intron
    frozenset({"G"}),           # pos +5 intron
    frozenset({"T"}),           # pos +6 intron
)

# Intronic GT at positions +1/+2 is the near-invariant anchor.
_DONOR_ANCHOR_POSITIONS = (3, 4)  # 0-indexed into the 9-bp window

# 3' splice site (acceptor): 4-bp core = 3 intronic + 1 exonic.
# Consensus: ``YAG | G``  — simplified; no polypyrimidine scoring.
_ACCEPTOR_CONSENSUS: Tuple[frozenset, ...] = (
    frozenset({"C", "T"}),      # pos -3 intron (Y = C/T)
    frozenset({"A"}),           # pos -2 intron
    frozenset({"G"}),           # pos -1 intron (strong anchor)
    frozenset({"G", "A"}),      # pos +1 exon
)

# Intronic AG at positions -2/-1 is the near-invariant anchor.
_ACCEPTOR_ANCHOR_POSITIONS = (1, 2)  # 0-indexed into the 4-bp window


DONOR_WINDOW = len(_DONOR_CONSENSUS)
ACCEPTOR_WINDOW = len(_ACCEPTOR_CONSENSUS)


# --------------------------------------------------------------------
# Scoring primitives
# --------------------------------------------------------------------


def score_donor(window: str) -> float:
    """Score a 9-bp window as a candidate 5' splice site (donor).
    Returns a ratio in ``[0, 1]``. Anchor mismatches (GT at +1/+2)
    force the score to 0, since biologically a donor without GT is
    not a donor under the canonical (U2-type) spliceosome.
    """
    if len(window) != DONOR_WINDOW:
        return 0.0
    w = window.upper()
    # Anchor check: both GT positions must match.
    for pos in _DONOR_ANCHOR_POSITIONS:
        if w[pos] not in _DONOR_CONSENSUS[pos]:
            return 0.0
    matches = sum(
        1 for i, base in enumerate(w)
        if base in _DONOR_CONSENSUS[i])
    return matches / DONOR_WINDOW


def score_acceptor(window: str) -> float:
    """Score a 4-bp window as a candidate 3' splice site (acceptor).
    Returns a ratio in ``[0, 1]``. Anchor mismatches (AG at -2/-1)
    force the score to 0.
    """
    if len(window) != ACCEPTOR_WINDOW:
        return 0.0
    w = window.upper()
    for pos in _ACCEPTOR_ANCHOR_POSITIONS:
        if w[pos] not in _ACCEPTOR_CONSENSUS[pos]:
            return 0.0
    matches = sum(
        1 for i, base in enumerate(w)
        if base in _ACCEPTOR_CONSENSUS[i])
    return matches / ACCEPTOR_WINDOW


# --------------------------------------------------------------------
# Candidate enumeration
# --------------------------------------------------------------------


# Callable interface for an alternative scorer. Gets (window_sequence,
# kind) where kind is "donor" or "acceptor", returns a probability in
# [0, 1]. The default is (score_donor, score_acceptor); an integration
# can pass a SpliceAI / Pangolin wrapper.
ScoreFn = Callable[[str, str], float]


def _default_score(window: str, kind: str) -> float:
    if kind == "donor":
        return score_donor(window)
    if kind == "acceptor":
        return score_acceptor(window)
    raise ValueError(
        "Unknown splice site kind %r — expected 'donor' or 'acceptor'"
        % (kind,))


def enumerate_candidates(
        contig: str,
        sequence: str,
        sequence_start: int,
        variant=None,
        donor_threshold: float = 0.7,
        acceptor_threshold: float = 0.7,
        min_intron_length: int = 30,
        max_intron_length: int = 10_000,
        score_fn: Optional[ScoreFn] = None,
) -> List[CrypticExonCandidate]:
    """Scan ``sequence`` for candidate cryptic exons — (acceptor,
    donor) pairs that could form a new internal exon on the
    rearranged allele.

    Parameters
    ----------
    contig : str
        Chromosome / contig the sequence lives on (for reporting;
        not used for scoring).
    sequence : str
        The genomic sequence to scan — typically the region around
        an SV breakpoint, or a long-read-assembled rearranged
        allele.
    sequence_start : int
        1-based genomic start position of ``sequence[0]``. Candidate
        intervals are reported relative to this origin.
    variant : Variant or None
        Source variant for provenance on returned candidates. ``None``
        is allowed for pure-sequence scans (e.g. long-read
        assemblies not tied to a single SV).
    donor_threshold, acceptor_threshold : float
        Minimum score to accept a site. Defaults (0.7) are chosen to
        be liberal — motif-based scoring has high false-positive
        rates, so the point is to enumerate candidates for
        downstream scoring, not to be the final caller.
    min_intron_length, max_intron_length : int
        A donor following an acceptor must leave an intron in this
        length range for the candidate to be reported. Narrow
        defaults keep the output tractable; callers that want
        everything pass ``min=0, max=None``.
    score_fn : callable or None
        Override the default motif scorer. Signature
        ``(window_sequence, kind) -> probability in [0, 1]``, with
        ``kind`` one of ``"donor"``, ``"acceptor"``. Plug in a
        SpliceAI / Pangolin wrapper here.

    Returns
    -------
    list of :class:`CrypticExonCandidate`
        One entry per accepted (acceptor, donor) pair forming a
        plausible internal exon. Empty list when nothing scores
        above threshold.
    """
    scorer = score_fn if score_fn is not None else _default_score
    seq = sequence.upper()

    # First pass: collect all candidate donor and acceptor positions.
    donor_positions: List[Tuple[int, float]] = []
    acceptor_positions: List[Tuple[int, float]] = []

    # Donor: the window is 9 bp centered on the exon/intron boundary.
    # Position p means the boundary is after exon position p+2
    # (0-indexed exonic anchor).
    for i in range(len(seq) - DONOR_WINDOW + 1):
        window = seq[i:i + DONOR_WINDOW]
        score = scorer(window, "donor")
        if score >= donor_threshold:
            # Report the intron-start position (index of first intron
            # base within the local sequence) so acceptor / donor
            # distances are comparable.
            intron_start = i + 3  # after the 3 exonic bases
            donor_positions.append((intron_start, score))

    # Acceptor: 4-bp window, intronic AG at positions 1/2, exonic
    # base at position 3.
    for i in range(len(seq) - ACCEPTOR_WINDOW + 1):
        window = seq[i:i + ACCEPTOR_WINDOW]
        score = scorer(window, "acceptor")
        if score >= acceptor_threshold:
            # Report the exon-start position (index of first exonic
            # base within the local sequence).
            exon_start = i + 3
            acceptor_positions.append((exon_start, score))

    # Second pass: pair acceptors with downstream donors to produce
    # cryptic internal exon candidates.
    candidates: List[CrypticExonCandidate] = []
    max_len = max_intron_length if max_intron_length is not None else len(seq)
    for exon_start_local, acc_score in acceptor_positions:
        for intron_start_local, don_score in donor_positions:
            if intron_start_local <= exon_start_local:
                continue
            exon_length = intron_start_local - exon_start_local
            if exon_length < min_intron_length:
                continue
            if exon_length > max_len:
                continue
            candidates.append(CrypticExonCandidate(
                variant=variant,
                contig=contig,
                interval_start=sequence_start + exon_start_local,
                interval_end=sequence_start + intron_start_local,
                donor_score=don_score,
                acceptor_score=acc_score,
            ))
    # Sort by summed motif score, best first — downstream consumers
    # typically only care about the top N.
    candidates.sort(
        key=lambda c: (c.donor_score or 0) + (c.acceptor_score or 0),
        reverse=True)
    return candidates


def enumerate_from_structural_variant(
        variant,
        flank: int = 500,
        score_fn: Optional[ScoreFn] = None,
        **kwargs) -> List[CrypticExonCandidate]:
    """Convenience wrapper: enumerate candidates in the genomic
    region around an SV breakpoint.

    Reads ``flank`` bp either side of ``variant.start`` from the
    variant's genome (and the same around ``variant.mate_start`` for
    breakends). If the variant carries an ``alt_assembly`` (long-
    read resolution), scans that instead of the reference.

    Returns the combined, sorted candidate list.
    """
    candidates: List[CrypticExonCandidate] = []
    genome = getattr(variant, "genome", None) or getattr(
        variant, "ensembl", None)

    # If the SV carries a long-read-assembled allele, prefer it.
    if getattr(variant, "alt_assembly", None):
        contig = getattr(variant, "contig", "alt_assembly")
        candidates.extend(enumerate_candidates(
            contig=contig,
            sequence=variant.alt_assembly,
            sequence_start=getattr(variant, "start", 0),
            variant=variant,
            score_fn=score_fn,
            **kwargs))
        return candidates

    # Otherwise, scan reference flanking regions around breakpoints.
    breakpoints = [(variant.contig, variant.start)]
    mate_contig = getattr(variant, "mate_contig", None)
    mate_start = getattr(variant, "mate_start", None)
    if mate_contig is not None and mate_start is not None:
        breakpoints.append((mate_contig, mate_start))

    for contig, pos in breakpoints:
        if genome is None:
            continue
        try:
            ref_seq = genome.genome.sequence(
                contig, pos - flank, pos + flank)
        except Exception:
            try:
                ref_seq = genome.genome.contig_sequence(contig)[
                    max(0, pos - flank):pos + flank]
            except Exception:
                continue
        candidates.extend(enumerate_candidates(
            contig=contig,
            sequence=ref_seq,
            sequence_start=pos - flank,
            variant=variant,
            score_fn=score_fn,
            **kwargs))
    candidates.sort(
        key=lambda c: (c.donor_score or 0) + (c.acceptor_score or 0),
        reverse=True)
    return candidates
