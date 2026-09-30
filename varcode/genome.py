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

Varcode inherits optional native reference DNA from the wrapped
pyensembl ``Genome``. An explicit ``fasta=`` overrides that source.
Without chromosome DNA, tiered lookups fall back to transcript cDNA.
Construction never opens or downloads native reference DNA.

Acts as a drop-in for ``pyensembl.Genome`` — every pyensembl attribute
is accessible via :py:meth:`__getattr__` delegation, so existing code
that calls ``genome.transcripts_at_locus(...)``,
``genome.reference_name``, etc. works unchanged.

Tracked in openvax/varcode#372 and #488.
"""

import os
import warnings
from typing import Any, Optional

# Module-level imports of the lookup helpers. Functions in this file
# call into them on hot paths; per-call imports would add measurable
# overhead in tight loops (e.g. left-alignment shifts in #369).
from .genome_sequence import (
    _fasta_range,
    _read_transcript_base,
    reference_base as _reference_base_lookup,
    reference_range as _reference_range_lookup,
)
from .reference import infer_genome


class Genome:
    """Pyensembl ``Genome`` + optional chromosome FASTA.

    Two source layers under one object:

    * ``self.fasta`` — optional chromosome FASTA. When attached,
      :meth:`sequence` reads from it directly; :meth:`reference_base`
      and :meth:`reference_range` prefer it before falling back to
      transcript cDNA.
    * Wrapped pyensembl ``Genome`` — transcript annotations + cDNA.
      Always present; provides the fallback for exonic positions
      when no FASTA is attached.

    The two methods have intentionally different fall-through semantics:

    * :meth:`sequence` — chromosome FASTA *only*. Delegates to native
      PyEnsembl DNA when configured, including its missing-DNA and
      invalid-interval errors. Returns ``""`` for genomes without DNA
      configured. Explicit FASTA overrides retain the legacy permissive
      lookup behavior.
    * :meth:`reference_base` / :meth:`reference_range` — tiered.
      FASTA first when attached; otherwise transcript cDNA. Use these
      when you want "whatever varcode can tell you" about a position.

    Parameters
    ----------
    ensembl_release :
        Release identifier — passes through :func:`varcode.reference.infer_genome`,
        so the same shapes accepted by ``load_vcf(genome=...)`` work here
        (int, reference-name string, ``pyensembl.Genome``).
        Passing another :class:`Genome` rewraps idempotently (inheriting
        ``fasta`` unless overridden).
    fasta :
        Optional chromosome FASTA:

        * Path string — opened with pyfaidx (only required on this path).
        * ``pyfaidx.Fasta`` — used as-is.
        * Any object supporting ``fa[contig][start:end]`` and returning
          a string or an object with ``.seq``.
        * ``None`` — inherit native PyEnsembl DNA (or a rewrapped
          Genome's override). Without DNA, features fall back to cDNA.

        Native and pre-opened readers are borrowed; their owner controls
        their lifetime. :meth:`close` closes only a reader this wrapper
        opened from an explicit path. Rewrapping borrows that reader, so
        keep its owning wrapper open while using it.
    verify :
        When True (default) and ``fasta`` is provided, spot-check a few
        exonic positions against pyensembl's transcript cDNA to catch
        mislabeled FASTAs (e.g. GRCh37 attached to a GRCh38 release).
        Iterates ``self.transcripts()`` which on a fresh process may
        trigger lazy DB construction; pass ``verify=False`` to defer
        that cost. Inherited native DNA is not opened or spot-checked.

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

    Equality and hashing fall back to ``object`` identity — two
    ``Genome`` instances wrapping the same pyensembl release compare
    unequal. This is intentional: a wrapper with an attached FASTA
    and one without are different objects in any sense that matters
    for varcode features, even if they share the same ``reference_name``.
    """

    def __init__(
            self,
            ensembl_release: Any = None,
            *,
            fasta: Optional[Any] = None,
            verify: bool = True):
        fasta_was_provided = fasta is not None
        self._fasta_override = None
        self._owns_fasta = False
        if isinstance(ensembl_release, Genome):
            self._ensembl = ensembl_release._ensembl
            # Copy only the override, without opening a lazy native reader.
            self._fasta_override = ensembl_release._fasta_override
        else:
            if ensembl_release is None:
                raise ValueError(
                    "varcode.Genome requires an ensembl_release argument "
                    "(int release number, reference-name string, or a "
                    "pyensembl.Genome / varcode.Genome instance).")
            self._ensembl, _ = infer_genome(ensembl_release)

        if fasta_was_provided:
            self._fasta_override = _resolve_fasta(fasta, self._ensembl)
            self._owns_fasta = isinstance(fasta, (str, bytes, os.PathLike))

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

    @property
    def fasta(self):
        """Explicit FASTA override, otherwise the lazy native DNA reader.

        Native access may index a local FASTA, but never downloads DNA.
        Returns None when native DNA is unconfigured or uninstalled.
        Assigning None removes the override and resumes native access.
        """
        if self._fasta_override is not None:
            return self._fasta_override
        return getattr(self._ensembl, "fasta", None)

    @fasta.setter
    def fasta(self, reader):
        if reader is self._fasta_override:
            return
        self.close()
        self._fasta_override = reader

    def close(self):
        """Close only a FASTA opened by this wrapper from an explicit path.

        Native genomes and borrowed readers remain open. Repeated calls
        are harmless. After closing an owned reader, assign a new reader
        (or None to resume native access) before making further lookups.
        """
        if self._owns_fasta:
            self._fasta_override.close()
            self._owns_fasta = False

    def __getattr__(self, name):
        """Delegate everything else to the wrapped pyensembl Genome.

        Called only when normal attribute lookup fails, so explicit
        attributes (``_ensembl``, ``fasta``, ``_missing_reference_warned``)
        still win. Underscore-prefixed names are not delegated — this
        prevents infinite recursion when Python's copy / pickle
        machinery probes for ``__deepcopy__`` / ``__reduce__`` /
        ``__getstate__`` etc., and keeps the wrapper's private state
        from accidentally pulling private state from the wrapped
        object.
        """
        if name.startswith("_"):
            raise AttributeError(name)
        ensembl = self.__dict__.get("_ensembl")
        if ensembl is None:
            raise AttributeError(
                "varcode.Genome has no wrapped pyensembl Genome; "
                "cannot resolve attribute %r" % name)
        return getattr(ensembl, name)

    def __repr__(self) -> str:
        ref = getattr(self._ensembl, "reference_name", "?")
        if self._fasta_override is not None:
            fasta_repr = "FASTA attached"
        elif getattr(self._ensembl, "requires_genome_fasta", False):
            fasta_repr = "native FASTA configured"
        else:
            fasta_repr = "no FASTA"
        return "varcode.Genome(reference_name=%r, %s)" % (ref, fasta_repr)

    def __dir__(self):
        """Include attributes of the wrapped pyensembl Genome so
        ``dir(genome)`` and IDE auto-completion surface the full
        delegated API, not just the wrapper's own attributes."""
        own = set(super().__dir__())
        ensembl = self.__dict__.get("_ensembl")
        if ensembl is not None:
            own.update(dir(ensembl))
        return sorted(own)

    # -- chromosome FASTA API ------------------------------------------

    def sequence(self, contig: str, start: int, end: int) -> str:
        """Chromosome FASTA sequence on the ``+`` strand.

        Uses 1-based inclusive coordinates and returns uppercase DNA.
        Configured native DNA delegates to PyEnsembl: missing DNA raises
        ``MissingGenomeFastaError``; invalid coordinates or absent contigs
        raise ``ValueError``. Install native DNA explicitly with
        ``download_genome_fasta()`` before querying remote sources.

        Without configured DNA, returns ``""``. Explicit overrides use
        the legacy permissive reader (unreadable intervals return ``""``).
        Does not fall back to cDNA; use :meth:`reference_base` or
        :meth:`reference_range` for tiered lookup.
        """
        if (self._fasta_override is None
                and getattr(self._ensembl, "requires_genome_fasta", False)):
            return self._ensembl.sequence(contig, start, end)
        fasta = self.fasta
        if fasta is None:
            return ""
        return _fasta_range(fasta, contig, start, end)

    # -- tiered lookup convenience -------------------------------------

    def reference_base(self, contig: str, position: int) -> str:
        """Tiered ``+`` strand base lookup (FASTA → transcript cDNA → ``""``).

        Convenience method delegating to
        :func:`varcode.genome_sequence.reference_base`.
        """
        return _reference_base_lookup(self, contig, position)

    def reference_range(self, contig: str, start: int, end: int) -> str:
        """Tiered ``+`` strand range lookup. All-or-nothing — returns
        ``""`` if any position in the range is uncovered by the chosen
        source."""
        return _reference_range_lookup(self, contig, start, end)


# --------------------------------------------------------------------
# Internal helpers
# --------------------------------------------------------------------


def _resolve_fasta(arg: Any, ensembl_genome: Any) -> Any:
    """Validate and normalize a FASTA argument.

    Path strings open via pyfaidx; objects are probed against a known
    contig from the genome (``hasattr(arg, '__getitem__')`` alone is
    too permissive — strings, lists, etc. all have it but don't honor
    the slice protocol).
    """
    if isinstance(arg, (str, bytes, os.PathLike)):
        try:
            import pyfaidx
        except ImportError as e:
            raise ImportError(
                "varcode.Genome(fasta=<path>) needs pyfaidx. Install "
                "with `pip install pyfaidx`, or pass a pre-opened FASTA "
                "object that supports fa[contig][start:end]."
            ) from e
        path_str = os.fspath(arg)
        try:
            return pyfaidx.Fasta(path_str)
        except (FileNotFoundError, OSError) as e:
            # Frame the path error as a varcode problem so users see
            # "the FASTA you passed to Genome" rather than a bare
            # pyfaidx traceback that they have to translate back.
            raise FileNotFoundError(
                "varcode.Genome(fasta=%r): could not open FASTA at this "
                "path (%s: %s). Check the path and that the file is "
                "readable; pyfaidx needs the FASTA and a matching .fai "
                "index alongside it (the index is created automatically "
                "on first open if writable)."
                % (path_str, type(e).__name__, e)) from e

    probe_contig = _probe_contig_for_validation(ensembl_genome)
    if probe_contig is None:
        # Can't validate (genome reports no contigs). Accept on faith;
        # the first real lookup will fail loudly if the object is malformed.
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


def _probe_contig_for_validation(genome: Any) -> Optional[str]:
    """Pick a contig name from the genome for probing a candidate FASTA.
    Returns ``None`` if the genome can't enumerate."""
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


def _verify_fasta_against_transcripts(ensembl_genome: Any, fasta: Any) -> None:
    """Spot-check that the attached FASTA matches pyensembl's transcript
    cDNA at a few exonic positions. Warns (does not raise) on mismatch
    so users with intentionally-divergent setups can override with
    ``verify=False``.
    """
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
            transcript_base = _read_transcript_base(t, position)
            fasta_base = _fasta_range(fasta, t.contig, position, position)
        except Exception:
            continue
        if not fasta_base or not transcript_base:
            continue
        checked += 1
        if fasta_base != transcript_base:
            mismatches.append((t.contig, position, fasta_base, transcript_base))

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
