"""Pinned rapid-release genomes (such as T2T-CHM13v2.0), delegated to PyEnsembl.

PyEnsembl 2.36+ natively catalogues checksum-pinned rapid-release assemblies
(``pyensembl.pinned_assemblies()``, ``pyensembl.PinnedGenome``,
``pyensembl.contig_aliases()``), including human ``T2T-CHM13v2.0``
(``GCA_009914755.4`` / ``hs1``) and its UCSC/GenBank/RefSeq contig names and
mitochondrion disambiguation (``CP068254.1`` native ``MT`` vs ``NC_012920.1``
rCRS in ``_rCRS`` analysis sets).

This module bridges PyEnsembl's native pinned-assembly API into Varcode's
reference/contig resolution, while retaining an offline fallback to an existing
ORFdata genome cache when a user has not yet run
``pyensembl install --reference-name T2T-CHM13v2.0``.
"""

from pyensembl import (
    PinnedGenome,
    RCRS_ACCESSION as RCRS,
    find_pinned_reference_name,
    is_rcrs_reference,
    mitochondrion_note as _pyensembl_mitochondrion_note,
    normalize_pinned_contig,
    pinned_contig_lengths as _pyensembl_contig_lengths,
    pinned_contig_spellings as _pyensembl_contig_spellings,
    resolve_pinned_assembly,
)
from pyensembl.pinned_genomes import lookup_pinned_alias


class PinnedGenomeNotInstalledError(FileNotFoundError, ValueError):
    """A catalogued pinned genome is not installed in PyEnsembl (or ORFdata)."""


def find_pinned_name(text):
    """The longest catalogued pinned reference name in ``text``, or None.

    Names are matched ignoring case, dashes and spaces, and only as whole
    words: ``hs1`` is found in ``hs1.fa`` and ``/refs/hs1/``, never in
    ``ths1`` or ``hs10``. The longest name wins, so ``chm13v2.0_maskedY_rCRS.fa``
    names the rCRS analysis set rather than ``chm13v2.0``.
    """
    return find_pinned_reference_name(text)


def pinned_name(name):
    """PyEnsembl's alias row for ``name`` (ignoring case/punctuation), as a dict, or None."""
    row = lookup_pinned_alias(name)
    if row is None:
        return None
    return dict(
        name=row.name,
        genome=row.assembly,
        mitochondrion=row.mitochondrion,
        note=row.note,
    )


def check_supported(name):
    """Raise ValueError for an unsupported version of a pinned assembly, such as T2T-CHM13v1.0."""
    row = lookup_pinned_alias(name)
    if row is not None and row.assembly is None:
        raise ValueError("%s is not a supported reference: %s" % (name, row.note))


def pinned_assembly(name):
    """The canonical assembly name of a pinned genome for ``name``, e.g.
    ``T2T-CHM13v2.0`` for ``hs1`` or ``chm13v2.0``; None for other names and
    for other versions such as T2T-CHM13v1.0 (#730)."""
    if not isinstance(name, str):
        return None
    row = lookup_pinned_alias(name)
    if row is None or row.assembly is None:
        return None
    return row.reference_name


def names_rcrs(name):
    """Whether files called ``name`` carry the rCRS as their mitochondrion."""
    return is_rcrs_reference(name)


def is_pinned_assembly(reference_name):
    """Whether ``reference_name`` names a catalogued pinned assembly or alias (#730)."""
    return pinned_assembly(reference_name) is not None


def _not_installed_error_type():
    try:
        import orfdata
    except ImportError:
        return PinnedGenomeNotInstalledError

    class _CompatNotInstalledError(orfdata.NotInstalledError, PinnedGenomeNotInstalledError):
        pass

    return _CompatNotInstalledError


def load_pinned_genome(name, release=None, source=None):
    """Open a catalogued pinned genome offline via PyEnsembl (with ORFdata cache fallback).

    ``release`` must match the pinned annotation's version, e.g. ``2022_07``.
    When not yet installed in PyEnsembl, falls back to an existing local
    ORFdata installation if present. If neither cache has the genome, raises an
    error (subclassing ``ValueError`` and ``FileNotFoundError``) naming
    ``pyensembl install --reference-name <reference_name>``.
    """
    check_supported(name)
    row = lookup_pinned_alias(name)
    if row is None or row.assembly is None:
        raise ValueError("Unknown pinned reference: %r" % (name,))
    pin = resolve_pinned_assembly(row.name)
    if release is not None and str(release) != str(pin.annotation_version):
        raise ValueError(
            "%s has one pinned annotation, %s %s (%s); release %s isn't available"
            % (pin.reference_name, pin.annotation_name, pin.annotation_version, pin.name, release)
        )

    genome = PinnedGenome.cached(pin.name)
    if genome.installed():
        return genome

    # Backwards-compatibility fallback: reuse an existing local ORFdata install
    # when PyEnsembl's cache doesn't have the pinned genome yet.
    try:
        import orfdata
    except ImportError:
        orfdata = None
    if orfdata is not None:
        try:
            return orfdata.load_genome(row.name)
        except orfdata.NotInstalledError:
            pass

    found = "" if source in (None, row.name) else " (from %r)" % (source,)
    if row.name == pin.reference_name:
        lead = "The annotation of %s%s is not installed." % (row.name, found)
    else:
        lead = "%s%s is %s, whose annotation is not installed." % (
            row.name, found, pin.reference_name
        )
    install_cmd = genome.install_string()
    err_cls = _not_installed_error_type()
    raise err_cls(
        "%s Install it once with '%s'." % (lead, install_cmd)
    )


def pinned_contig_name(contig, reference_name, chrm_is_rcrs=False):
    """The annotation's name for ``contig`` in a pinned assembly, or None."""
    if not is_pinned_assembly(reference_name):
        return None
    return normalize_pinned_contig(contig, reference_name, chrm_is_rcrs=bool(chrm_is_rcrs))


def pinned_spellings(reference_name, chrm_is_rcrs=False, mitochondrion=True):
    """Each annotation contig of a pinned assembly to all its spellings."""
    return _pyensembl_contig_spellings(
        reference_name,
        chrm_is_rcrs=bool(chrm_is_rcrs),
        include_mitochondrion=bool(mitochondrion),
    )


def pinned_contig_lengths(reference_name, chrm_is_rcrs=False):
    """The annotation's name for each contig to its length."""
    return _pyensembl_contig_lengths(reference_name, chrm_is_rcrs=bool(chrm_is_rcrs))


def mitochondrion_note(contig, reference_name):
    """Why a pinned assembly's mitochondrion has no annotation, or None."""
    if not is_pinned_assembly(reference_name) or str(contig) not in ("MT", "chrM", "M", RCRS):
        return None
    canonical = pinned_assembly(reference_name)
    if str(contig) == RCRS:
        return _pyensembl_mitochondrion_note(contig, canonical)
    return (
        "%s's own mitochondrion has no annotation: Ensembl's annotation of %s has no "
        "mitochondrial genes, and its position 1 is what the rCRS numbers 577, so it can't be "
        "annotated as the rCRS. Only the chm13v2.0_maskedY_rCRS analysis set's chrM is the rCRS, "
        "annotated as %s; load_vcf reads that from the file name, and chrm_is_rcrs=True says it "
        "for a Variant whose chrM is the rCRS" % (canonical, canonical, RCRS)
    )
