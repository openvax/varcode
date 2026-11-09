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


_GRCH37_NAMES = frozenset({"grch37", "hg19", "b37", "ncbi37"})


def _is_grch37(reference_name: str | None) -> bool:
    return reference_name is not None and reference_name.strip().lower() in _GRCH37_NAMES


def _is_hg19_mitochondrion(
        contig: str | int, reference_name: str | None, chrm_is_rcrs: bool = False) -> bool:
    """True for ``chrM`` or ``M`` (any case) under GRCh37/hg19, unless the
    input declares its chrM to be the rCRS (``chrm_is_rcrs``, #710).

    hg19's chrM is NC_001807, not the rCRS (NC_012920) that Ensembl's GRCh37
    names MT, and its positions are offset from it: hg19 chrM:12995 is rCRS
    MT:12994. A variant there must not be renamed onto MT (#707). hg38's chrM
    is the rCRS, and ``chrMT`` is a chr-prefixed rCRS ``MT``.
    """
    return (not chrm_is_rcrs and _is_grch37(reference_name)
            and str(contig).upper() in ("CHRM", "M"))


# The rCRS (NC_012920) is 16,569 bp; hg19's chrM (NC_001807) is 16,571 bp.
_RCRS_LENGTH = 16569


_T2T_NAMES = frozenset({"t2t-chm13v2.0", "gca_009914755.4"})


def _is_t2t(reference_name: str | None) -> bool:
    return reference_name is not None and reference_name.strip().lower() in _T2T_NAMES


def _may_alias_mitochondrion(
        annotation_contig: str, sequence_contig: str, reference_name: str | None,
        sequence_length: int | None) -> bool:
    """Whether a sequence source's ``sequence_contig`` may be read as the
    annotation's ``annotation_contig`` when only the mitochondrion's spelling
    differs (``chrM`` or ``M`` against ``MT``).

    Such a pair can name two different molecules:

    * GRCh37: hg19's ``chrM`` is NC_001807 (16,571 bp), offset from the rCRS
      that Ensembl names ``MT`` (#707). A source's ``chrM`` is read as ``MT``
      only at the rCRS length, 16,569 bp, as in GENCODE's GRCh37 (#710).
    * T2T-CHM13: a FASTA's ``chrM`` is CHM13's own mitochondrion (CP068254.1)
      in ``hs1`` and ``chm13v2.0.fa``, but the rCRS (NC_012920.1) in the
      ``_rCRS`` analysis sets. Both are 16,569 bp, numbered 576 to 577 bp
      apart, so the spellings are never paired (#721).

    Any other pair of spellings, and any other reference, is unaffected.
    """
    annotation_family, sequence_family = (_mitochondrion_spelling(annotation_contig),
                                          _mitochondrion_spelling(sequence_contig))
    if None in (annotation_family, sequence_family) or annotation_family == sequence_family:
        return True
    if _is_t2t(reference_name):
        return False
    if _is_grch37(reference_name):
        return annotation_family == "MT" and sequence_length == _RCRS_LENGTH
    return True


def _mitochondrion_spelling(contig: str) -> str | None:
    """``"chrM"`` for ``chrM``/``M``, ``"MT"`` for ``MT``/``chrMT``, else None."""
    name = str(contig).upper()
    return "chrM" if name in ("CHRM", "M") else "MT" if name in ("MT", "CHRMT") else None


# Appended to errors about a GRCh37 chrM, which Ensembl's annotation lacks.
_HG19_MITOCHONDRION_NOTE = (
    "hg19's chrM (NC_001807) is not the rCRS that Ensembl names MT, and its "
    "positions are offset from it; give rCRS positions on contig MT (if this "
    "chrM is already the rCRS, as in GENCODE's GRCh37, pass chrm_is_rcrs=True; "
    "load_vcf does when the VCF header gives chrM the rCRS length, 16569)")


def _annotation_contig_name(
        contig: str | int, reference_name: str | None, convert_ucsc: bool = True,
        chrm_is_rcrs: bool = False) -> str:
    """The name of ``contig`` in Ensembl annotation of ``reference_name``.

    As :func:`normalize_contig_name`, except that hg19's mitochondrion has no
    Ensembl name, so UCSC conversion leaves it ``chrM`` rather than renaming
    it onto the rCRS ``MT`` (#707), unless ``chrm_is_rcrs`` says this chrM is
    the rCRS (#710). In a genome ORFdata pins, conversion uses ORFdata's
    contig table: T2T-CHM13v2.0's ``chr1``, ``CP068277.2`` and
    ``NC_060925.1`` are ``1``, and its ``chrM`` is ``MT``, CHM13's own
    mitochondrion, or, with ``chrm_is_rcrs``, ``NC_012920.1``, the rCRS
    (#723, #722). The public function, and so :class:`ContigResolver`,
    only matches spellings and still pairs them.
    """
    if convert_ucsc:
        # A genome ORFdata pins, such as T2T-CHM13v2.0, names each contig's
        # UCSC, GenBank and RefSeq spellings (#723).
        from .pinned_references import pinned_contig_name
        pinned = pinned_contig_name(contig, reference_name, chrm_is_rcrs)
        if pinned is not None:
            return pinned
    if convert_ucsc and _is_hg19_mitochondrion(contig, reference_name, chrm_is_rcrs):
        return "chrM"
    return normalize_contig_name(contig, reference_name=reference_name, convert_ucsc=convert_ucsc)


def _convert_ucsc_contig_name(contig: str, reference_name: str | None) -> str:
    """Convert spelling only; also used by Variant with normalization disabled."""
    match = re.fullmatch(
        r"chr(?:Un|[0-9XYM]+)_([A-Za-z]+[0-9]+)v([0-9]+)(?:_(?:random|alt|fix))?",
        contig)
    if match:
        return "%s.%s" % (match[1].upper(), match[2])
    # UCSC hg19.chromAlias.txt assigns these unversioned names to version 1.
    # Never infer that version from a target vocabulary or a newer assembly.
    if _is_grch37(reference_name):
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
