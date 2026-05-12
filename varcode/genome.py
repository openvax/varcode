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

"""Declarative wrapper for a pyensembl ``Genome`` + optional chromosome FASTA.

varcode's reference is the pyensembl ``Genome`` object. By default it
ships transcript and protein FASTAs but not the chromosome FASTA, so
features needing raw genomic bases (intronic, intergenic, flanking)
fall back to transcript-cDNA coverage only. :class:`Genome` is the
construction-time composition point: pass a release identifier plus
an optional ``fasta``, and downstream lookups see both tiers.

Acts as a drop-in for ``pyensembl.Genome`` — every pyensembl attribute
is accessible via :py:meth:`__getattr__` delegation, so existing code
that calls ``genome.transcripts_at_locus(...)``,
``genome.reference_name``, etc. works unchanged.

This module is a stopgap until ``pyensembl.Genome`` natively supports
the chromosome FASTA (tracked in openvax/pyensembl#337). Once that
lands, this wrapper either becomes a one-line delegate or evaporates.
Designed so the migration is straightforward — the API surface here
(``Genome.sequence()``, ``Genome.fasta``, the construction-time
composition pattern) is intentionally identical to the proposed
pyensembl-native API.

Tracked in openvax/varcode#372.
"""

import warnings


class Genome:
    """Pyensembl ``Genome`` + optional chromosome FASTA.

    Parameters
    ----------
    ensembl_release : int, str, pyensembl.Genome, or Genome
        Release identifier — passes through :func:`varcode.reference.infer_genome`,
        so the same shapes accepted by ``load_vcf(genome=...)`` work here.
        Passing another :class:`Genome` rewraps idempotently (inheriting
        ``fasta`` unless overridden).
    fasta : str, path-like, FASTA-shaped object, or None, optional
        Optional chromosome FASTA:

        * Path string — opened with pyfaidx (only required on this path).
        * ``pyfaidx.Fasta`` — used as-is.
        * Any object supporting ``fa[contig][start:end]`` and returning
          a string or an object with ``.seq``.
        * ``None`` — no chromosome-level access; features fall back to
          transcript cDNA (Tier 2 only).
    verify : bool, optional
        When True (default) and ``fasta`` is provided, spot-check a few
        exonic positions against pyensembl's transcript cDNA to catch
        mislabeled FASTAs (e.g. GRCh37 attached to a GRCh38 release).
        Iterates ``self.transcripts()`` which on a fresh process may
        trigger lazy DB construction; pass ``verify=False`` to defer
        that cost.

    Examples
    --------

    >>> import varcode
    >>> g = varcode.Genome(81, fasta="/path/to/GRCh38.fa")  # doctest: +SKIP
    >>> vc = varcode.load_vcf("tumor.vcf", genome=g)        # doctest: +SKIP
    >>> g.sequence("7", 117_480_000, 117_480_050)           # doctest: +SKIP
    'AC...'

    Notes
    -----

    Returned sequence is always uppercase. Soft-masked (lowercase)
    repeat annotations from the FASTA are dropped — varcode treats all
    bases uniformly. Callers that need the soft-masking signal should
    read the FASTA directly.
    """

    def __init__(self, ensembl_release=None, *, fasta=None, verify=True):
        # Late imports break a circular dep: reference.py imports
        # pyensembl; we don't want this module imported during varcode
        # package construction in unusual orders.
        from .reference import infer_genome

        fasta_was_provided = fasta is not None
        if isinstance(ensembl_release, Genome):
            # Idempotent rewrap. Inherits ._ensembl always; inherits
            # .fasta unless the caller is providing a fresh one.
            self._ensembl = ensembl_release._ensembl
            if fasta_was_provided:
                self.fasta = _resolve_fasta(fasta, self._ensembl)
            else:
                self.fasta = ensembl_release.fasta
        elif ensembl_release is None:
            # Tests / placeholder use — no genome to delegate to.
            # Most consumers will fail downstream when they try to
            # look up transcripts; intentional, no silent surprise.
            self._ensembl = None
            self.fasta = (_resolve_fasta(fasta, None)
                          if fasta_was_provided else None)
        else:
            self._ensembl, _ = infer_genome(ensembl_release)
            self.fasta = (_resolve_fasta(fasta, self._ensembl)
                          if fasta_was_provided else None)

        # Verify only when the caller freshly provided a FASTA.
        # Inheriting a verified FASTA on rewrap doesn't need a re-check.
        if (verify and fasta_was_provided
                and self.fasta is not None
                and self._ensembl is not None):
            _verify_fasta_against_transcripts(self._ensembl, self.fasta)

        # Per-instance warn-once for missing-reference signals from
        # consumers (cryptic_exons, etc.). Bound to this Genome's
        # lifetime — no module-level cache to invalidate.
        self._missing_reference_warned = False

    def __getattr__(self, name):
        """Delegate everything else to the wrapped pyensembl Genome.

        Called only when normal attribute lookup fails, so explicit
        attributes (``_ensembl``, ``fasta``, ``_missing_reference_warned``)
        still win. Underscore-prefixed names are not delegated to keep
        Python internals (pickling, copying) from accidentally pulling
        private state from the wrapped object.
        """
        if name.startswith("_"):
            raise AttributeError(name)
        ensembl = self.__dict__.get("_ensembl")
        if ensembl is None:
            raise AttributeError(
                "varcode.Genome has no wrapped pyensembl Genome; "
                "cannot resolve attribute %r" % name)
        return getattr(ensembl, name)

    def __repr__(self):
        ref = getattr(self._ensembl, "reference_name", "?")
        fasta_repr = "no FASTA" if self.fasta is None else "FASTA attached"
        return "varcode.Genome(reference_name=%r, %s)" % (ref, fasta_repr)

    # -- chromosome FASTA API (mirrors openvax/pyensembl#337) ----------

    def sequence(self, contig, start, end):
        """Chromosome FASTA sequence on the ``+`` strand.

        1-based inclusive coordinates. Returns ``""`` when no FASTA is
        attached or the contig / range isn't covered. Mirrors the
        ``pyensembl.Genome.sequence()`` shape proposed upstream so the
        eventual migration is mechanical.
        """
        if self.fasta is None:
            return ""
        from .genome_sequence import _tier1_range
        return _tier1_range(self.fasta, contig, start, end)

    # -- tiered lookup convenience -------------------------------------

    def reference_base(self, contig, position):
        """Tiered ``+`` strand base lookup (FASTA → transcript cDNA → ``""``).

        Convenience method delegating to
        :func:`varcode.genome_sequence.reference_base`.
        """
        from .genome_sequence import reference_base
        return reference_base(self, contig, position)

    def reference_range(self, contig, start, end):
        """Tiered ``+`` strand range lookup. All-or-nothing — returns
        ``""`` if any position in the range is uncovered by the chosen
        tier."""
        from .genome_sequence import reference_range
        return reference_range(self, contig, start, end)


# --------------------------------------------------------------------
# Internal helpers (shared with varcode.genome_sequence). Centralized
# here because Genome construction is the only public path that builds
# them up.
# --------------------------------------------------------------------


def _resolve_fasta(arg, ensembl_genome):
    """Validate and normalize a FASTA argument.

    Path strings open via pyfaidx; objects are probed against a known
    contig from the genome (``hasattr(arg, '__getitem__')`` alone is
    too permissive — strings, lists, etc. all have it but don't honor
    the slice protocol).
    """
    import os
    if isinstance(arg, (str, bytes, os.PathLike)):
        try:
            import pyfaidx
        except ImportError as e:
            raise ImportError(
                "varcode.Genome(fasta=<path>) needs pyfaidx. Install "
                "with `pip install pyfaidx`, or pass a pre-opened FASTA "
                "object that supports fa[contig][start:end]."
            ) from e
        return pyfaidx.Fasta(os.fspath(arg))

    probe_contig = _probe_contig_for_validation(ensembl_genome)
    if probe_contig is None:
        # Can't validate (no genome, or genome reports no contigs).
        # Accept on faith; the first real lookup will fail loudly if
        # the object is malformed.
        return arg
    try:
        slicer = arg[probe_contig]
        _ = slicer[0:1]
    except Exception as e:
        raise TypeError(
            "varcode.Genome: object of type %s does not honor the "
            "pyfaidx-style protocol fa[contig][start:end] (probe "
            "against contig %r failed: %s). Pass a path, a pyfaidx.Fasta, "
            "or an object with the same indexing shape."
            % (type(arg).__name__, probe_contig, e))
    return arg


def _probe_contig_for_validation(genome):
    """Pick a contig name from the genome for probing a candidate FASTA.
    Returns ``None`` if no genome was attached or it can't enumerate."""
    if genome is None:
        return None
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


def _verify_fasta_against_transcripts(ensembl_genome, fasta):
    """Spot-check that the attached FASTA matches pyensembl's transcript
    cDNA at a few exonic positions. Warns (does not raise) on mismatch
    so users with intentionally-divergent setups can override with
    ``verify=False``.
    """
    from .genome_sequence import _read_transcript_base, _tier1_range

    try:
        transcripts = list(ensembl_genome.transcripts())[:50]
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
            position = exon.start + 1
            tier2 = _read_transcript_base(t, position)
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
            "varcode.Genome: FASTA contents disagree with pyensembl "
            "transcript sequence at %d/%d probed positions (%s). "
            "The attached FASTA may be for a different build than "
            "the pyensembl release (e.g. GRCh37 vs GRCh38). Pass "
            "verify=False to suppress this check."
            % (len(mismatches), checked, details),
            stacklevel=3)
