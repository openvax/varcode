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

from serializable import Serializable

from .contigs import _RCRS_LENGTH, ContigResolver, _is_grch37, _is_t2t_chm13

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


class Genome(Serializable):
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

        Its contig names may be spelled differently from the
        annotation's: a FASTA with UCSC names (``chr22``, ``chrM``)
        serves an Ensembl annotation (``22``, ``MT``). A name the FASTA
        lists is read as is; otherwise its one spelling of that contig,
        matched as :class:`~varcode.ContigResolver` matches them. A
        differently spelled mitochondrion is read only when it is the
        annotation's: under GRCh37, a ``chrM`` of the rCRS length
        (16,569 bp; hg19's own chrM is not the rCRS); under T2T-CHM13,
        never, since its FASTAs hold one of two mitochondria of that
        length. Readers that can't list their names are read by the
        annotation's name.

        Native and pre-opened readers are borrowed; their owner controls
        their lifetime. :meth:`close` closes only a reader this wrapper
        opened from an explicit path. Rewrapping borrows that reader, so
        keep its owning wrapper open while using it.
    verify :
        When True (default) and ``fasta`` is provided, spot-check a few
        exonic positions against pyensembl's transcript cDNA to catch
        mislabeled FASTAs (e.g. GRCh37 attached to a GRCh38 release).
        A FASTA that lists its names but names none of the annotation's
        contigs, even respelled, raises ``ValueError``.
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
            _check_fasta_names_annotated_contigs(self)
            _verify_fasta_against_transcripts(self)

        # Per-instance warn-once for missing-reference signals from
        # consumers (cryptic_exons, etc.). Bound to this Genome's
        # lifetime — no module-level cache to invalidate.
        self._missing_reference_warned = False

    # A DNA wrapper keeps its existing identity semantics, independently of
    # Serializable's value-based default implementations.
    __eq__ = object.__eq__
    __hash__ = object.__hash__

    def to_dict(self):
        """Archive the exact annotation and an explicit FASTA path, never an open reader.

        Borrowed in-memory readers have no portable installation to reopen and
        raise instead of silently losing chromosome DNA. ORF/proteome snapshots
        retain their sequence and reference metadata independently of this API.
        """
        fasta = self._fasta_override
        path = None if fasta is None else getattr(fasta, "filename", None)
        if fasta is not None and not isinstance(path, (str, bytes, os.PathLike)):
            raise TypeError("Cannot archive a Genome with an in-memory FASTA reader; "
                            "use a filesystem FASTA or a native PyEnsembl DNA source")
        return dict(ensembl_release=self._ensembl,
                    fasta=None if path is None else os.path.abspath(os.fsdecode(path)),
                    verify=False)

    __str__ = object.__str__

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

    def __setstate__(self, state):
        """Restore wrappers saved before native DNA inheritance as borrowers."""
        self.__dict__.update(state)
        if "fasta" in state:
            self._fasta_override = self.__dict__.pop("fasta")
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
        return _fasta_range(fasta, self._fasta_contig(contig), start, end)

    def _fasta_contig(self, contig):
        """The attached FASTA's own name for annotation ``contig`` (#721).

        Each name is resolved once per reader; see the ``fasta``
        parameter for the rules. Absent or ambiguous names are read as
        given, which finds nothing.
        """
        fasta = self.fasta
        spellings = self.__dict__.get("_fasta_spellings")
        if spellings is None or spellings[0] is not fasta:
            names = _fasta_names(fasta)
            resolver = None if names is None else ContigResolver(
                names, reference_name=_contig_reference_name(self._ensembl))
            spellings = self._fasta_spellings = (fasta, resolver, {})
        _, resolver, found = spellings
        name = found.get(contig)
        if name is None:
            name = found[contig] = _fasta_spelling(resolver, fasta, contig)
        return name

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

    # Probe a name the reader lists, if it lists any: its spelling of the
    # annotation's names (chr1 for 1) need not match.
    names = _fasta_names(arg)
    probe_contig = names[0] if names else _probe_contig_for_validation(ensembl_genome)
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


def _fasta_names(fasta: Any) -> Optional[list]:
    """The contig names a FASTA reader lists (``pyfaidx.Fasta`` and
    mappings do), or None for readers that can't list them."""
    keys = getattr(fasta, "keys", None)
    if not callable(keys):
        return None
    try:
        return [name for name in keys() if isinstance(name, str)]
    except Exception:
        return None


def _contig_reference_name(genome: Any) -> Optional[str]:
    """The assembly whose contig spelling rules apply, as ``load_vcf`` reads it."""
    name = getattr(genome, "assembly", None) or getattr(genome, "reference_name", None)
    return name if isinstance(name, str) and name.strip() else None


def _fasta_spelling(resolver: Optional[ContigResolver], fasta: Any, contig: Any) -> Any:
    """``contig`` as the FASTA spells it, or ``contig`` itself when the
    FASTA lists it, can't list its names, or has no single match."""
    if resolver is None:
        return contig
    try:
        resolution = resolver.resolve(contig)
    except (TypeError, ValueError):
        return contig
    if resolution.status != "alias":
        return contig
    if resolution.normalized_name == "MT" and not _fasta_mitochondrion_is_annotations(
            resolver.reference_name, contig, fasta, resolution.resolved_name):
        return contig
    return resolution.resolved_name


def _fasta_mitochondrion_is_annotations(
        reference_name: Optional[str], contig: Any, fasta: Any, fasta_name: str) -> bool:
    """Whether the FASTA's ``fasta_name`` (``chrM``, say) holds the
    annotation's mitochondrion ``contig`` (``MT``).

    hg19's chrM (NC_001807, 16,571 bp) is not the rCRS that Ensembl's
    GRCh37 names MT, and its positions are offset from it (#707). Under
    GRCh37, a FASTA's chrM therefore serves MT only at the rCRS length,
    as in GENCODE's GRCh37 (#710), and hg19's chrM, which variants keep
    as ``chrM``, reads no other spelling. T2T-CHM13 FASTAs hold either
    CHM13's own mitochondrion or the rCRS, both 16,569 bp, with positions
    576 to 577 bp apart (#722), so a respelling never chooses between
    them. Elsewhere (hg38, mm10, canFam3) a FASTA's chrM is the
    annotation's MT, as a VCF's chrM becomes MT.
    """
    if _is_t2t_chm13(reference_name):
        return False
    if not _is_grch37(reference_name):
        return True
    if str(contig).upper() != "MT":
        return False
    try:
        return len(fasta[fasta_name]) == _RCRS_LENGTH
    except Exception:
        return False


def _check_fasta_names_annotated_contigs(genome: "Genome") -> None:
    """Raise when an attached FASTA that lists its names names none of
    the annotation's contigs, even respelled: every DNA lookup would
    find nothing, silently (#721)."""
    names = _fasta_names(genome.fasta)
    try:
        contigs = list(genome._ensembl.contigs())
    except Exception:
        return
    if not names or not contigs:
        return
    listed = set(names)
    if any(genome._fasta_contig(contig) in listed for contig in contigs):
        return
    raise ValueError(
        "varcode.Genome: the FASTA names none of the %d contigs in the %s "
        "annotation, even allowing for spellings such as chr22 for 22. "
        "It has %s; the annotation has %s. Attach DNA of this assembly "
        "whose contig names are the annotation's or their UCSC spellings."
        % (len(contigs), _contig_reference_name(genome._ensembl) or "attached",
           _some(names), _some(contigs)))


def _some(names: list, shown: int = 5) -> str:
    listed = ", ".join(repr(name) for name in names[:shown])
    return listed + (", …" if len(names) > shown else "")


def _verify_fasta_against_transcripts(genome: "Genome") -> None:
    """Spot-check that the attached FASTA matches pyensembl's transcript
    cDNA at a few exonic positions. Warns (does not raise) on mismatch
    so users with intentionally-divergent setups can override with
    ``verify=False``.
    """
    try:
        transcripts = list(genome._ensembl.transcripts())[:50]
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
            fasta_base = _fasta_range(
                genome.fasta, genome._fasta_contig(t.contig), position, position)
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
