"""Reference-aware contig spelling conversion and vocabulary resolution.

These functions do not establish sequence/coordinate equivalence, perform
liftover, or open reference data. Callers supply the reference vocabulary.
"""

from collections.abc import Iterable
from dataclasses import dataclass
import re
from typing import Literal

from pyensembl.locus import normalize_chromosome


ContigResolutionStatus = Literal[
    "exact", "alias", "ambiguous", "missing", "reference_unavailable"]


def _check_reference_name(reference_name: str | None) -> None:
    if reference_name is not None:
        if not isinstance(reference_name, str):
            raise TypeError("reference_name must be a string or None")
        if not reference_name.strip():
            raise ValueError("reference_name cannot be empty")


def _convert_ucsc_contig_name(contig: str, reference_name: str | None) -> str:
    """Convert spelling only; also used by Variant with normalization disabled."""
    match = re.fullmatch(
        r"chr(?:Un|[0-9XYM]+)_([A-Za-z]+[0-9]+)v([0-9]+)(?:_(?:random|alt|fix))?",
        contig)
    if match:
        return "%s.%s" % (match[1].upper(), match[2])
    # UCSC hg19.chromAlias.txt assigns these unversioned names to version 1.
    # Never infer that version from a target vocabulary or a newer assembly.
    if reference_name is not None and reference_name.strip().lower() in {
            "grch37", "hg19", "b37", "ncbi37"}:
        match = re.fullmatch(r"chr(?:Un|[0-9XYM]+)_(gl000[0-9]{3})(?:_random)?", contig)
        if match and 191 <= int(match[1][-3:]) <= 249:
            return match[1].upper() + ".1"
    if contig.startswith("chr") and "_" not in contig:
        contig = contig[3:]
    return "MT" if contig == "M" else contig


def normalize_contig_name(
        name: str | int, *, reference_name: str | None = None,
        convert_ucsc: bool = True) -> str:
    """Return a normalized spelling without querying a reference dataset.

    Parameters
    ----------
    name : str or int
        Chromosome/contig name. Uses PyEnsembl's integer and case rules.
    reference_name : str, optional
        Assembly spelling context. GRCh37/hg19 (also b37/NCBI37) enables
        unversioned UCSC GL aliases. Other names do not imply that rule.
        For an annotation genome, prefer its ``assembly`` attribute when
        available; dated annotations may use an accession as reference_name.
    convert_ucsc : bool
        Convert chr prefixes, M/MT and accession-bearing UCSC names (default
        True). False performs only PyEnsembl's case/integer normalization.

    Returns
    -------
    str
        Normalized name, with explicit accession versions retained. This is
        not a claim that the name exists in any reference. Use a resolver to
        select the actual spelling in an annotation, BAM or FASTA vocabulary.
    """
    _check_reference_name(reference_name)
    if convert_ucsc and isinstance(name, str):
        canonical_chr_match = re.fullmatch(r"chr([0-9]+|[xym]|mt)", name, flags=re.IGNORECASE)
        if canonical_chr_match:
            name = "chr" + canonical_chr_match[1].upper()
    contig = normalize_chromosome(name)
    return _convert_ucsc_contig_name(contig, reference_name) if convert_ucsc else contig


def _resolution_key(name: str | int, reference_name: str | None) -> str:
    """Build a spelling key without applying annotation's zero-contig sentinel.

    ``0`` is a valid external reference name (including SAM SN). Resolving
    that name does not make it a valid Variant/annotation chromosome.
    """
    if (isinstance(name, str) and name == "0") or (type(name) is int and name == 0):
        return "0"
    return normalize_contig_name(name, reference_name=reference_name)


@dataclass(frozen=True)
class ContigResolution:
    """A name lookup with original spelling and explicit resolution state.

    ``resolved_name`` is an actual target-vocabulary name for ``exact`` or
    ``alias`` results, otherwise None. ``candidates`` contains the selected
    name or all ambiguous matches in sorted order. ``normalized_name`` is
    only a comparison key and need not occur in the target vocabulary.
    """

    original_name: str | int
    normalized_name: str
    resolved_name: str | None
    reference_name: str | None
    status: ContigResolutionStatus
    candidates: tuple[str, ...] = ()

    @property
    def resolved(self) -> bool:
        return self.resolved_name is not None

    @property
    def diagnostic(self) -> str | None:
        """Failure detail for caller-side aggregation; no warning is emitted."""
        if self.resolved:
            return None
        context = " for reference %r" % self.reference_name if self.reference_name is not None else ""
        if self.status == "reference_unavailable":
            return "Cannot resolve contig %r%s: reference names unavailable" % (self.original_name, context)
        if self.status == "ambiguous":
            return "Ambiguous contig %r%s: %s" % (
                self.original_name, context, ", ".join(repr(name) for name in self.candidates))
        return "Contig %r (normalized %r) is absent from the supplied names%s" % (
            self.original_name, self.normalized_name, context)


class ContigResolver:
    """Index one reference vocabulary for repeated exact-first name lookups.

    Parameters
    ----------
    available_contigs : iterable of str or None
        Actual names, e.g. ``bam.references``, ``genome.contigs()`` or FASTA
        keys. The iterable is consumed once and copied; duplicates are ignored.
        None means unavailable reference names. An empty iterable is an empty
        reference vocabulary, so lookups are missing rather than unavailable.
    reference_name : str, optional
        Assembly spelling context, as in :func:`normalize_contig_name`.

    Notes
    -----
    An exact supplied spelling wins even when other aliases are present.
    Otherwise exactly one normalized match is required. Ambiguity is returned
    explicitly, never resolved by iteration order. There is no global vocabulary cache,
    warning emission, reference download, or automatic assembly inference.
    Canonical chr/mitochondrial spellings accept case variants. Other names
    retain the public normalization rules; no viral/decoy synonym catalogue
    is inferred. The valid external name ``0`` is accepted even though
    annotation/Variant chromosome normalization rejects it.
    """

    def __init__(self, available_contigs: Iterable[str] | None, *, reference_name: str | None = None):
        _check_reference_name(reference_name)
        if isinstance(available_contigs, (str, bytes)):
            raise TypeError("available_contigs must be an iterable of names, not a string")
        self._reference_name = reference_name
        self._available_contigs: frozenset[str] | None = None
        self._aliases: dict[str, tuple[str, ...]] = {}
        if available_contigs is not None:
            aliases: dict[str, set[str]] = {}
            names: set[str] = set()
            for name in available_contigs:
                if not isinstance(name, str):
                    raise TypeError("available_contigs must contain only string names")
                key = _resolution_key(name, reference_name)
                names.add(name)
                aliases.setdefault(key, set()).add(name)
            self._available_contigs = frozenset(names)
            self._aliases = {key: tuple(sorted(values)) for key, values in aliases.items()}

    @property
    def reference_name(self) -> str | None:
        return self._reference_name

    @property
    def available_contigs(self) -> frozenset[str] | None:
        return self._available_contigs

    def resolve(self, name: str | int) -> ContigResolution:
        """Resolve a name, preserving both the supplied and target spellings."""
        key = _resolution_key(name, self.reference_name)
        candidates: tuple[str, ...] = ()
        status: ContigResolutionStatus
        if self.available_contigs is None:
            status = "reference_unavailable"
        elif str(name) in self.available_contigs:
            status, candidates = "exact", (str(name),)
        else:
            candidates = self._aliases.get(key, ())
            status = "alias" if len(candidates) == 1 else "ambiguous" if candidates else "missing"
        resolved = candidates[0] if status in ("exact", "alias") else None
        return ContigResolution(name, key, resolved, self.reference_name, status, candidates)


def resolve_contig_name(
        name: str | int, available_contigs: Iterable[str] | None, *,
        reference_name: str | None = None) -> ContigResolution:
    """Resolve one name; reuse :class:`ContigResolver` for repeated lookups."""
    return ContigResolver(available_contigs, reference_name=reference_name).resolve(name)
