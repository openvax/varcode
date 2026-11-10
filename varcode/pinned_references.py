"""Genomes ORFdata pins, such as T2T-CHM13v2.0, found by the names their files go by.

PyEnsembl has no T2T-CHM13 release. ORFdata pins Ensembl's rapid-release
annotation of T2T-CHM13v2.0 as ``t2t-chm13``, and publishes the names its
files go by (``orfdata.reference_names()``: ``hs1``, ``chm13v2.0_maskedY_rCRS``,
``GCA_009914755.4``, ...) and how each file names its contigs
(``orfdata.contig_aliases()``). Varcode reads both from ORFdata, so the two
packages accept the same names; nothing here lists them again.

A T2T file's ``chrM`` is one of two molecules, both 16,569 bp: CHM13's own
mitochondrion (CP068254.1), which Ensembl names ``MT``, or, in the
``_rCRS`` analysis set, the rCRS (NC_012920.1), which begins 576 bases
earlier in its numbering. Only the file's name tells them apart. ORFdata
annotates the rCRS as its own contig, ``NC_012920.1``, with Ensembl GRCh38's
MT genes, and CHM13's own ``MT`` has no annotation.
"""

import re
from functools import lru_cache

# The revised Cambridge Reference Sequence: GRCh38's MT.
RCRS = "NC_012920.1"


def _normalize(name):
    return name.strip().lower().replace("-", "_").replace(" ", "_")


@lru_cache(maxsize=1)
def _names():
    """ORFdata's names, longest first: (whole-word pattern, name, genome, mitochondrion, note)."""
    from orfdata import reference_names

    rows = []
    for row in reference_names().itertuples(index=False):
        key = _normalize(row.name)
        pattern = re.compile(r"(?<![a-z0-9])%s(?![a-z0-9])" % re.escape(key))
        rows.append((pattern, row.name, row.genome, row.mitochondrion, row.note))
    return tuple(sorted(rows, key=lambda row: -len(row[1])))


def find_pinned_name(text):
    """The longest of ORFdata's reference names in ``text``, or None.

    Names are matched ignoring case, dashes and spaces, and only as whole
    words: ``hs1`` is found in ``hs1.fa`` and ``/refs/hs1/``, never in
    ``ths1`` or ``hs10``. The longest name wins, so ``chm13v2.0_maskedY_rCRS.fa``
    names the rCRS analysis set rather than ``chm13v2.0``.
    """
    normalized = _normalize(text)
    for pattern, name, *_ in _names():
        if pattern.search(normalized):
            return name
    return None


def pinned_name(name):
    """ORFdata's row for exactly ``name`` (ignoring case), as a dict, or None."""
    key = _normalize(name)
    for pattern, found, genome, mitochondrion, note in _names():
        if _normalize(found) == key:
            return dict(name=found, genome=genome, mitochondrion=mitochondrion, note=note)
    return None


def check_supported(name):
    """Raise ValueError for a name of another version of a pinned assembly, such as T2T-CHM13v1.0."""
    row = pinned_name(name)
    if row is not None and row["genome"] is None:
        raise ValueError("%s is not a supported reference: %s" % (name, row["note"]))


def pinned_assembly(name):
    """The assembly of the genome ORFdata pins for ``name``, e.g.
    ``T2T-CHM13v2.0`` for ``hs1``; None for other names and for other
    versions such as T2T-CHM13v1.0."""
    row = pinned_name(name)
    if row is None or row["genome"] is None:
        return None
    return next((reference for reference, pin in _pins().items() if pin.name == row["genome"]), None)


def names_rcrs(name):
    """Whether files called ``name`` carry the rCRS as their mitochondrion."""
    row = pinned_name(name)
    return row is not None and row["mitochondrion"] == RCRS


@lru_cache(maxsize=1)
def _pins():
    from orfdata import pinned_genomes

    return {genome.reference_name: genome for genome in pinned_genomes()}


def is_pinned_assembly(reference_name):
    """Whether ``reference_name`` is the assembly of a genome ORFdata pins, e.g. ``T2T-CHM13v2.0``."""
    return isinstance(reference_name, str) and reference_name in _pins()


def load_pinned_genome(name, release=None, source=None):
    """Open the genome ORFdata pins for ``name``, offline.

    ``release`` must be the pinned annotation's version, e.g. ``2022_07``.
    A genome that isn't installed raises ORFdata's ``NotInstalledError``
    (both a ``ValueError`` and a ``FileNotFoundError``), which names the
    command that installs it; ``source``, the text ``name`` was found in, is
    added to its message.
    """
    import orfdata

    check_supported(name)
    row = pinned_name(name)
    pin = next(g for g in orfdata.pinned_genomes() if g.name == row["genome"])
    if release is not None and str(release) != str(pin.annotation_version):
        raise ValueError(
            "%s has one pinned annotation, %s %s (ORFdata's %s); release %s isn't available"
            % (pin.reference_name, pin.annotation_name, pin.annotation_version, pin.name, release))
    try:
        return orfdata.load_genome(row["name"])
    except orfdata.NotInstalledError as error:
        found = "" if source in (None, row["name"]) else " (from %r)" % (source,)
        if row["name"] == pin.reference_name:
            lead = "The annotation of %s%s comes from ORFdata." % (row["name"], found)
        else:
            lead = "%s%s is %s, whose annotation comes from ORFdata." % (
                row["name"], found, pin.reference_name)
        raise type(error)("%s %s" % (lead, error)) from error


@lru_cache(maxsize=None)
def _contigs(reference_name, chrm_is_rcrs):
    """Every spelling of each contig (UCSC, GenBank, RefSeq and the
    annotation's own) to the annotation's name, and each name's length."""
    from orfdata import contig_aliases

    pin = _pins()[reference_name]
    frame = contig_aliases(pin.name)
    rows = [row for row in frame.to_dict("records")]
    if chrm_is_rcrs:
        # This file's chrM is the rCRS, not the assembly's own mitochondrion.
        own = [row for row in rows if row["ucsc"] in ("chrM", "chrMT")]
        rows = [row for row in rows if row not in own]
        contig = (pin.mitochondria.get(RCRS) or {}).get("contig") or dict(
            name=RCRS, ucsc="chrM", genbank="J01415.2", refseq=RCRS, length=16569)
        rows.append(contig)
    names, lengths = {}, {}
    for row in rows:
        for alias in (row["name"], row["ucsc"], row["genbank"], row["refseq"]):
            if alias:
                names[alias] = row["name"]
        if row["ucsc"] == "chrM":
            names["M"] = row["name"]
        lengths[row["name"]] = int(row["length"])
    return names, lengths


def pinned_contig_name(contig, reference_name, chrm_is_rcrs=False):
    """The annotation's name for ``contig`` in a pinned assembly, or None.

    ``chr22``, ``CP068256.2`` and ``NC_060946.1`` are all T2T-CHM13v2.0's
    ``22``. ``chrM`` is ``MT``, CHM13's own mitochondrion, unless
    ``chrm_is_rcrs`` says it is the rCRS, which is ``NC_012920.1``. None
    when ``reference_name`` isn't a pinned assembly or the name is unknown.
    """
    if not is_pinned_assembly(reference_name):
        return None
    names = _contigs(reference_name, bool(chrm_is_rcrs))[0]
    contig = str(contig)
    if contig in names:
        return names[contig]
    # The canonical chr spellings in any case: chrx, CHR1, chrm.
    match = re.fullmatch(r"chr([0-9]+|[xym]|mt)", contig, flags=re.IGNORECASE)
    if match:
        return names.get("chr" + match[1].upper())
    return names.get(contig.upper())


def pinned_spellings(reference_name, chrm_is_rcrs=False, mitochondrion=True):
    """Each annotation contig of a pinned assembly to all its spellings.

    ``mitochondrion=False`` leaves out the mitochondrion, for a file whose
    name doesn't say which one its ``chrM`` is.
    """
    names, lengths = _contigs(reference_name, bool(chrm_is_rcrs))
    mitochondria = {names.get("chrM")} if not mitochondrion else set()
    spellings = {}
    for alias, name in names.items():
        if name not in mitochondria:
            spellings.setdefault(name, set()).add(alias)
    return spellings


def pinned_contig_lengths(reference_name, chrm_is_rcrs=False):
    """The annotation's name for each contig to its length, from ORFdata's table."""
    return _contigs(reference_name, bool(chrm_is_rcrs))[1]


def mitochondrion_note(contig, reference_name):
    """Why a pinned assembly's mitochondrion has no annotation, or None."""
    if not is_pinned_assembly(reference_name) or str(contig) not in ("MT", "chrM", "M", RCRS):
        return None
    pin = _pins()[reference_name]
    if str(contig) == RCRS:
        return ("this %s annotation has no genes on the rCRS; ORFdata's %s annotates them with "
                "Ensembl GRCh38's MT genes" % (reference_name, pin.name))
    return ("%s's own mitochondrion has no annotation: Ensembl's annotation of %s has no "
            "mitochondrial genes, and its position 1 is what the rCRS numbers 577, so it can't be "
            "annotated as the rCRS. Only the chm13v2.0_maskedY_rCRS analysis set's chrM is the rCRS, "
            "annotated as %s; load_vcf reads that from the file name, and chrm_is_rcrs=True says it "
            "for a Variant whose chrM is the rCRS" % (reference_name, reference_name, RCRS))
