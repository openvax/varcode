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

"""Reference sequence lookup against a pyensembl :class:`Genome`.

varcode's reference is the pyensembl ``Genome`` object the user passes
through ``genome=``. By default that genome ships transcript and
protein FASTAs but not the chromosome FASTA — so any feature that
needs raw bases at arbitrary genomic positions (intronic, intergenic,
flanking) has to go through this module's tiered lookup:

1. **Tier 1** — if the genome has an attached chromosome FASTA via
   :func:`attach_genome_fasta`, read from it.
2. **Tier 2** — fall back to pyensembl transcript cDNA via
   ``transcript.spliced_offset()`` for any transcript covering the
   position. Reverse-complements for ``-`` strand transcripts so the
   result is always on the ``+`` strand.
3. **Tier 3** — return ``""``. The caller decides whether empty is an
   error or just "no shift possible / no candidates / etc."

Internal API: feature code calls :func:`reference_base` or
:func:`reference_range`. User-facing API: :func:`attach_genome_fasta`.

Tracked in openvax/varcode#372.

Notes
-----

* Returned sequence is always uppercase. Soft-masked (lowercase) repeat
  annotations from the FASTA are normalized away — varcode treats all
  bases uniformly. Callers that need the soft-masking signal should
  read the FASTA directly.
* :func:`attach_genome_fasta` mutates the genome object in place.
  Because pyensembl typically returns cached singletons
  (``cached_release(N)`` is shared across the process), the attachment
  is visible to every caller holding that genome reference. This is
  usually what you want; if you need a FASTA-less view, construct a
  separate :class:`pyensembl.Genome` instance for it.
"""

import warnings

_VARCODE_FASTA_ATTR = "_varcode_fasta"
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANTGCAN")


# --------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------


def attach_genome_fasta(genome, fasta_or_path, *, verify=True):
    """Attach a chromosome FASTA to a pyensembl ``Genome`` so varcode
    can read raw reference bases at arbitrary positions.

    The default ``pyensembl install`` ships transcript and protein
    FASTAs only. Features that need intronic, intergenic, or flanking
    bases (left-alignment, cryptic-exon scoring, sequence-aware
    splice prediction) require this attachment.

    Mutates ``genome`` in place. Returns it for use as a pseudo-
    constructor:

    >>> g = varcode.attach_genome_fasta(EnsemblRelease(81), fasta)  # doctest: +SKIP

    Parameters
    ----------
    genome : pyensembl.Genome
        The genome to extend. **Mutated in place** — see module docstring
        re: pyensembl's cached singletons.
    fasta_or_path : str, path-like, FASTA-shaped object, or None
        One of:

        * A path string (``"/path/to/GRCh38.fa"``) — opened with
          pyfaidx internally. pyfaidx becomes a runtime dependency
          *only* on this path.
        * A ``pyfaidx.Fasta`` instance — used as-is.
        * Any object supporting ``fa[contig][start:end]`` and
          returning either a string or an object with ``.seq``.
          Validated by a probe call (see notes).
        * ``None`` — detach. Removes any previously attached FASTA.
    verify : bool, optional
        When True (default), spot-check the attached FASTA against an
        exonic position varcode already knows (via Tier 2) and warn if
        the bases disagree. Catches mislabeled FASTAs (GRCh37 attached
        to a GRCh38 release). Set False to skip if you have an
        intentionally-divergent setup.

    Returns
    -------
    pyensembl.Genome
        The same ``genome`` object, for chaining.

    Raises
    ------
    TypeError
        If ``fasta_or_path`` is not a path, ``None``, or a FASTA-shaped
        object that responds to a probe lookup.

    Notes
    -----

    The "FASTA-shaped object" check is a probe, not a class check:
    we attempt ``fasta[contig][0:1]`` against a contig the genome
    reports. This rejects strings, lists, and bare dicts (which all
    have ``__getitem__`` but don't honor the slice protocol pyfaidx
    uses) and accepts pyfaidx, the in-tree ``_DictBackedFasta`` test
    helper, and any custom mmap/S3-backed storage with the same shape.
    """
    if fasta_or_path is None:
        # Detach.
        if hasattr(genome, _VARCODE_FASTA_ATTR):
            delattr(genome, _VARCODE_FASTA_ATTR)
        return genome

    fasta = _resolve_fasta(fasta_or_path, genome=genome)
    setattr(genome, _VARCODE_FASTA_ATTR, fasta)

    if verify:
        _verify_fasta_against_transcripts(genome, fasta)

    return genome


def reference_base(genome, contig, position):
    """Return the ``+`` strand reference base at ``(contig, position)``.

    1-based inclusive coordinates matching pyensembl. Returns ``""``
    when no tier of the lookup covers this position.

    See module docstring for the tiered fallback rules.
    """
    fasta = getattr(genome, _VARCODE_FASTA_ATTR, None)
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
    locus query, one slice into ``transcript.sequence``). Falls back
    to per-position scanning only when no single transcript spans the
    whole range.

    See module docstring for the tiered fallback rules.
    """
    if end < start:
        raise ValueError(
            "reference_range: end=%d < start=%d on %r"
            % (end, start, contig))

    fasta = getattr(genome, _VARCODE_FASTA_ATTR, None)
    if fasta is not None:
        tier1 = _tier1_range(fasta, contig, start, end)
        if tier1 and len(tier1) == end - start + 1:
            return tier1
        # Tier 1 returned partial or empty — fall through to Tier 2.

    # Tier 2 fast path: a single transcript that fully covers the range.
    span = _tier2_range_via_single_transcript(genome, contig, start, end)
    if span is not None:
        return span

    # Tier 2 slow path: walk each position. If any is uncovered,
    # all-or-nothing returns empty.
    bases = []
    for pos in range(start, end + 1):
        base = _tier2_base(genome, contig, pos)
        if not base:
            return ""
        bases.append(base)
    return "".join(bases)


# --------------------------------------------------------------------
# Internal helpers
# --------------------------------------------------------------------


def _resolve_fasta(arg, genome):
    """Validate and normalize a FASTA argument.

    Path strings open via pyfaidx; objects are probed against a known
    contig from the genome (we don't trust ``hasattr`` alone — strings
    and dicts have ``__getitem__`` but don't honor the slice protocol).
    """
    import os
    if isinstance(arg, (str, bytes, os.PathLike)):
        try:
            import pyfaidx
        except ImportError as e:
            raise ImportError(
                "attach_genome_fasta needs pyfaidx to open a FASTA "
                "path. Install with `pip install pyfaidx`, or pass a "
                "pre-opened FASTA object that supports "
                "fa[contig][start:end]."
            ) from e
        return pyfaidx.Fasta(os.fspath(arg))

    # Probe duck-typed object: it must support fa[contig][slice].
    probe_contig = _probe_contig_for_validation(genome)
    if probe_contig is None:
        # We can't validate — genome reports no contigs. Accept the
        # object on faith; if it's wrong, the first real lookup fails.
        return arg
    try:
        slicer = arg[probe_contig]
        _ = slicer[0:1]
    except Exception as e:
        raise TypeError(
            "attach_genome_fasta: object of type %s does not honor "
            "the pyfaidx-style protocol fa[contig][start:end] "
            "(probe against contig %r failed: %s). Pass a path, a "
            "pyfaidx.Fasta, or an object with the same indexing shape."
            % (type(arg).__name__, probe_contig, e))
    return arg


def _probe_contig_for_validation(genome):
    """Pick a contig name from the genome for probing a candidate FASTA.
    Returns ``None`` if the genome can't supply one."""
    for accessor in ("contigs", "all_contigs"):
        getter = getattr(genome, accessor, None)
        if callable(getter):
            try:
                contigs = list(getter())
            except Exception:
                contigs = []
            if contigs:
                return contigs[0]
    return None


def _verify_fasta_against_transcripts(genome, fasta):
    """Spot-check that the attached FASTA matches what pyensembl says
    about a few exonic positions. Warns (does not raise) on mismatch
    so users can override with ``verify=False`` if they know what
    they're doing.

    Picks up to three transcripts and reads their first CDS-start base
    from both Tier 1 (FASTA) and Tier 2 (transcript cDNA). Mismatches
    almost always mean the FASTA is for a different build than the
    pyensembl release.
    """
    try:
        # Reuse the genome's transcript iteration; bail quietly if
        # pyensembl can't enumerate (no DB built yet, etc.).
        transcripts = list(genome.transcripts())[:50]
    except Exception:
        return

    checked = 0
    mismatches = []
    for t in transcripts:
        if checked >= 3:
            break
        if not t.exons:
            continue
        try:
            exon = t.exons[0]
            position = exon.start + 1   # 1-based; safely inside the exon
            tier2 = _tier2_base(genome, t.contig, position)
            tier1 = _tier1_range(fasta, t.contig, position, position)
        except Exception:
            continue
        if not tier1 or not tier2:
            continue
        checked += 1
        if tier1 != tier2:
            mismatches.append((t.contig, position, tier1, tier2))

    if mismatches:
        details = "; ".join(
            "%s:%d FASTA=%r vs transcript=%r" % m for m in mismatches)
        warnings.warn(
            "attach_genome_fasta: FASTA contents disagree with "
            "pyensembl transcript sequence at %d/%d probed positions "
            "(%s). The attached FASTA may be for a different build "
            "than the pyensembl release (e.g. GRCh37 vs GRCh38). "
            "Pass verify=False to suppress this check."
            % (len(mismatches), checked, details),
            stacklevel=3)


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
    # pyfaidx returns a Sequence with .seq; bare-string adapters return
    # a string directly.
    raw = getattr(span, "seq", None)
    if raw is None:
        raw = span if isinstance(span, str) else str(span)
    return raw.upper()


def _tier2_base(genome, contig, position):
    """Pyensembl-transcript fallback for a single + strand base.

    Returns ``""`` when no transcript covers ``position``. When
    multiple transcripts cover it and disagree on the base, logs a
    debug-level warning and returns the first one's answer (callers
    that care about disagreement can intercept the warning).
    """
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


def _read_transcript_base(transcript, position):
    """Return the + strand base at ``position`` from one transcript's
    spliced cDNA, or ``""`` if the position isn't on any exon of this
    transcript."""
    try:
        offset = transcript.spliced_offset(position)
    except (ValueError, KeyError):
        return ""
    seq = transcript.sequence
    if seq is None or offset >= len(seq):
        return ""
    base = seq[offset].upper()
    if transcript.strand == "-":
        base = base.translate(_COMPLEMENT)
    return base


def _tier2_range_via_single_transcript(genome, contig, start, end):
    """Try to satisfy a Tier 2 range from one transcript that fully
    covers ``[start, end]``. Returns the + strand sequence on success,
    or ``None`` if no single transcript spans the range (caller falls
    back to per-position).
    """
    try:
        transcripts = genome.transcripts_at_locus(contig, start, end)
    except Exception:
        return None

    for t in transcripts:
        try:
            start_offset = t.spliced_offset(start)
            end_offset = t.spliced_offset(end)
        except (ValueError, KeyError):
            continue
        seq = t.sequence
        if seq is None:
            continue
        lo, hi = sorted((start_offset, end_offset))
        if hi >= len(seq):
            continue
        slice_seq = seq[lo:hi + 1].upper()
        if len(slice_seq) != end - start + 1:
            # Shouldn't happen for a well-formed transcript, but guard
            # against pyensembl returning a truncated sequence.
            continue
        if t.strand == "-":
            slice_seq = slice_seq.translate(_COMPLEMENT)[::-1]
        return slice_seq
    return None
