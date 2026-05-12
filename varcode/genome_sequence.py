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

"""Reference sequence lookup against a pyensembl ``Genome`` or
:class:`varcode.Genome`.

varcode has two possible sources of reference bases at a genomic
position, and falls back between them in this order:

1. **Chromosome FASTA** (when attached) — read directly from the
   FASTA file via the ``.fasta`` attribute on :class:`varcode.Genome`.
   Covers every base on every contig the FASTA contains. Pass a
   FASTA when constructing the genome to enable this source.
2. **Transcript cDNA** — fall back to pyensembl's transcript
   sequences via ``transcript.spliced_offset()`` for any transcript
   covering the position. Reverse-complements for ``-`` strand
   transcripts so the result is always on the ``+`` strand. Covers
   every exonic base of every transcript varcode loaded (which is
   always available — pyensembl ships transcript FASTAs by default).
3. *Nothing covers it.* Return ``""`` — the caller decides whether
   empty is an error or just "no shift possible / no candidates /
   etc."

Internal API: feature code calls :func:`reference_base` or
:func:`reference_range`. User-facing API: :class:`varcode.Genome`'s
methods (or its construction-time ``fasta=`` kwarg) and the same
module-level functions for ad-hoc use.

Tracked in openvax/varcode#372. Upstream pyensembl support for the
chromosome FASTA itself is tracked in openvax/pyensembl#337.
"""

import warnings
from typing import Any, Optional

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANTGCAN")


# --------------------------------------------------------------------
# Public lookup helpers
# --------------------------------------------------------------------


def reference_base(genome: Any, contig: str, position: int) -> str:
    """Return the ``+`` strand reference base at ``(contig, position)``.

    1-based inclusive coordinates matching pyensembl. Returns ``""``
    when no source covers this position.

    Works on any genome shape — :class:`varcode.Genome` (chromosome
    FASTA available via ``.fasta``), bare ``pyensembl.Genome``
    (transcript cDNA only), or anything duck-typed for
    ``transcripts_at_locus``.

    See module docstring for the full fallback order.
    """
    fasta = getattr(genome, "fasta", None)
    if fasta is not None:
        from_fasta = _fasta_range(fasta, contig, position, position)
        if from_fasta:
            return from_fasta
    return _base_from_genome_transcripts(genome, contig, position)


def reference_range(genome: Any, contig: str, start: int, end: int) -> str:
    """Return the ``+`` strand reference sequence over ``[start, end]``.

    1-based inclusive coordinates. Returns ``""`` if **any** position
    in the range is not covered by the chosen source — the lookup is
    all-or-nothing, so callers see a partial range only when they ask
    for one explicitly. (For best-effort partial answers, iterate
    :func:`reference_base`.)

    Transcript-cDNA path first tries a single-transcript range lookup
    (cheap: one locus query, one slice). Falls back to per-position
    scanning only when no single transcript spans the whole range.

    See module docstring for the full fallback order.
    """
    if end < start:
        raise ValueError(
            "reference_range: end=%d < start=%d on %r"
            % (end, start, contig))

    fasta = getattr(genome, "fasta", None)
    if fasta is not None:
        from_fasta = _fasta_range(fasta, contig, start, end)
        if from_fasta and len(from_fasta) == end - start + 1:
            return from_fasta

    span = _range_from_single_transcript(genome, contig, start, end)
    if span is not None:
        return span

    bases = []
    for pos in range(start, end + 1):
        base = _base_from_genome_transcripts(genome, contig, pos)
        if not base:
            return ""
        bases.append(base)
    return "".join(bases)


# --------------------------------------------------------------------
# Source readers (also consumed by varcode.Genome construction)
# --------------------------------------------------------------------


def _fasta_range(fasta: Any, contig: str, start: int, end: int) -> str:
    """Read a range from an attached chromosome FASTA.

    Returns ``""`` for missing contigs, out-of-range slices, or any
    pyfaidx/adapter error. Bases come back uppercase (soft-masking
    dropped — see module docstring)."""
    try:
        slicer = fasta[contig]
    except Exception:
        # Defensive: pyfaidx raises KeyError for missing contigs but
        # custom adapters may raise other types. Treat any failure
        # to resolve the contig as "not covered."
        return ""
    try:
        span = slicer[start - 1:end]   # 1-based inclusive -> 0-based half-open
    except Exception:
        return ""
    raw = getattr(span, "seq", None)
    if raw is None:
        raw = span if isinstance(span, str) else str(span)
    return raw.upper()


def _base_from_genome_transcripts(
        genome: Any, contig: str, position: int) -> str:
    """Walk the genome's transcripts at ``position`` and return the
    first ``+`` strand base any of them yields.

    Returns ``""`` when no transcript covers ``position``. When
    multiple transcripts cover it and disagree on the base, warns
    once per ``(contig, position)`` and returns the first answer."""
    try:
        transcripts = genome.transcripts_at_locus(contig, position, position)
    except Exception:
        return ""

    answer = ""
    for t in transcripts:
        base = _read_transcript_base(t, position)
        if not base:
            continue
        if answer and base != answer:
            _warn_transcript_disagreement(contig, position, answer, base)
            break
        if not answer:
            answer = base
    return answer


def _plus_strand_slice(transcript: Any, start: int, end: int) -> str:
    """Return the ``+`` strand bases over ``[start, end]`` from one
    transcript's spliced cDNA, or ``""`` if the requested range isn't
    fully on this transcript's exons. 1-based inclusive.

    Strand handling lives here so callers don't have to reason about
    it. For ``+`` strand transcripts the slice maps directly. For
    ``-`` strand transcripts the cDNA is in transcript orientation
    (the reverse-complement of the ``+`` strand), so the slice is
    complemented and then reversed to come back in ``+`` strand
    5'→3' order.

    Returns ``""`` rather than raising for any non-coverage condition
    (intronic, past the end of the transcript, range that spans an
    exon boundary, etc.) — callers loop over multiple transcripts and
    need a uniform "this one didn't work" signal.
    """
    try:
        start_offset = transcript.spliced_offset(start)
        end_offset = transcript.spliced_offset(end)
    except (ValueError, KeyError):
        return ""
    seq = transcript.sequence
    if seq is None:
        return ""
    lo, hi = sorted((start_offset, end_offset))
    if hi >= len(seq):
        return ""
    slice_seq = seq[lo:hi + 1].upper()
    if len(slice_seq) != end - start + 1:
        return ""
    if transcript.strand == "-":
        slice_seq = slice_seq.translate(_COMPLEMENT)[::-1]
    return slice_seq


def _read_transcript_base(transcript: Any, position: int) -> str:
    """Single-position specialization of :func:`_plus_strand_slice`."""
    return _plus_strand_slice(transcript, position, position)


def _range_from_single_transcript(
        genome: Any, contig: str, start: int, end: int) -> Optional[str]:
    """Try to satisfy a range read from one transcript that fully
    covers ``[start, end]``. Returns the ``+`` strand sequence on
    success, or ``None`` if no single transcript spans the range
    (caller falls back to per-position).
    """
    try:
        transcripts = genome.transcripts_at_locus(contig, start, end)
    except Exception:
        return None
    for t in transcripts:
        span = _plus_strand_slice(t, start, end)
        if span:
            return span
    return None


# --------------------------------------------------------------------
# Internal: transcript-disagreement warn-once
# --------------------------------------------------------------------


_DISAGREEMENT_WARNED = set()


def _warn_transcript_disagreement(contig, position, first_base, other_base):
    """Warn once per ``(contig, position)`` when transcripts overlapping
    a single position disagree on the ``+`` strand base.

    Range scans cross many positions; without dedup, an
    annotation-skew region would emit a separate warning for every
    base. The dedup is process-wide because the underlying transcript
    annotations are too — re-warning across calls would still produce
    duplicate signals for the same biological discrepancy.
    """
    key = (contig, position)
    if key in _DISAGREEMENT_WARNED:
        return
    _DISAGREEMENT_WARNED.add(key)
    warnings.warn(
        "reference_base: transcripts overlapping %s:%d disagree on "
        "the + strand base (%r vs %r). Returning the first answer; "
        "the discrepancy may indicate annotation skew or an "
        "alt-locus position. This warning fires at most once per "
        "(contig, position)."
        % (contig, position, first_base, other_base),
        stacklevel=4)
