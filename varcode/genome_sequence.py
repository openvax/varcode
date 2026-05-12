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

"""Tiered reference sequence lookup against a pyensembl ``Genome`` or
:class:`varcode.Genome`.

The genome object is the reference. varcode's :class:`Genome` wrapper
adds an optional chromosome FASTA (the ``.fasta`` attribute);
``pyensembl.Genome`` has only transcript and protein FASTAs. Both
work here — the dispatch reads ``getattr(genome, 'fasta', None)``
and falls through tiers:

1. **Tier 1** — if the genome has a ``.fasta`` attached, read from it.
2. **Tier 2** — fall back to pyensembl transcript cDNA via
   ``transcript.spliced_offset()`` for any transcript covering the
   position. Reverse-complements for ``-`` strand transcripts so the
   result is always on the ``+`` strand.
3. **Tier 3** — return ``""``. The caller decides whether empty is an
   error or just "no shift possible / no candidates / etc."

Internal API: feature code calls :func:`reference_base` or
:func:`reference_range`. User-facing API: :class:`varcode.Genome`'s
methods (or its construction-time ``fasta=`` kwarg) and the same
module-level functions for ad-hoc use.

Tracked in openvax/varcode#372. Upstream pyensembl support for the
chromosome FASTA itself is tracked in openvax/pyensembl#337.
"""

import warnings

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANTGCAN")


# --------------------------------------------------------------------
# Public lookup helpers
# --------------------------------------------------------------------


def reference_base(genome, contig, position):
    """Return the ``+`` strand reference base at ``(contig, position)``.

    1-based inclusive coordinates matching pyensembl. Returns ``""``
    when no tier of the lookup covers this position.

    Works on any genome shape — :class:`varcode.Genome` (Tier 1
    available via ``.fasta``), bare ``pyensembl.Genome`` (Tier 2
    only), or anything duck-typed for ``transcripts_at_locus``.

    See module docstring for the tiered fallback rules.
    """
    fasta = getattr(genome, "fasta", None)
    if fasta is not None:
        tier1 = _tier1_range(fasta, contig, position, position)
        if tier1:
            return tier1
    return _tier2_base(genome, contig, position)


def reference_range(genome, contig, start, end):
    """Return the ``+`` strand reference sequence over ``[start, end]``.

    1-based inclusive coordinates. Returns ``""`` if **any** position
    in the range is not covered by the current tier — the lookup is
    all-or-nothing, so callers see a partial range only when they ask
    for one explicitly. (For best-effort partial answers, iterate
    :func:`reference_base`.)

    Tier 2 first tries a single-transcript range lookup (cheap: one
    locus query, one slice). Falls back to per-position scanning only
    when no single transcript spans the whole range.

    See module docstring for the tiered fallback rules.
    """
    if end < start:
        raise ValueError(
            "reference_range: end=%d < start=%d on %r"
            % (end, start, contig))

    fasta = getattr(genome, "fasta", None)
    if fasta is not None:
        tier1 = _tier1_range(fasta, contig, start, end)
        if tier1 and len(tier1) == end - start + 1:
            return tier1

    span = _tier2_range_via_single_transcript(genome, contig, start, end)
    if span is not None:
        return span

    bases = []
    for pos in range(start, end + 1):
        base = _tier2_base(genome, contig, pos)
        if not base:
            return ""
        bases.append(base)
    return "".join(bases)


# --------------------------------------------------------------------
# Tier implementations (also consumed by varcode.Genome construction)
# --------------------------------------------------------------------


def _tier1_range(fasta, contig, start, end):
    """Single-range FASTA lookup. Returns ``""`` for missing contigs,
    out-of-range slices, or any pyfaidx/adapter error. Bases come back
    uppercase (soft-masking dropped — see module docstring)."""
    try:
        slicer = fasta[contig]
    except (KeyError, ValueError, Exception):
        return ""
    try:
        span = slicer[start - 1:end]   # 1-based inclusive -> 0-based half-open
    except (IndexError, ValueError, Exception):
        return ""
    raw = getattr(span, "seq", None)
    if raw is None:
        raw = span if isinstance(span, str) else str(span)
    return raw.upper()


def _tier2_base(genome, contig, position):
    """Pyensembl-transcript fallback for a single ``+`` strand base.

    Returns ``""`` when no transcript covers ``position``. When
    multiple transcripts cover it and disagree on the base, logs a
    warning and returns the first answer."""
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
            warnings.warn(
                "reference_base: transcripts overlapping %s:%d "
                "disagree on the + strand base (%r vs %r). Returning "
                "the first answer; the discrepancy may indicate "
                "annotation skew or an alt-locus position."
                % (contig, position, answer, base),
                stacklevel=4)
            break
        if not answer:
            answer = base
    return answer


def _plus_strand_slice(transcript, start, end):
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


def _read_transcript_base(transcript, position):
    """Single-position specialization of :func:`_plus_strand_slice`."""
    return _plus_strand_slice(transcript, position, position)


def _tier2_range_via_single_transcript(genome, contig, start, end):
    """Try to satisfy a Tier 2 range from one transcript that fully
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
