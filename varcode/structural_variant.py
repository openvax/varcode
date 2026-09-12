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

"""Structural-variant types (#252 / #257 / #259 / #264 / #305).

varcode's core :class:`~varcode.Variant` models point variants,
MNVs, and simple indels as explicit ref/alt nucleotide strings.
Structural variants — deletions, duplications, inversions, large
insertions, translocations / breakends — don't fit that shape; they
were previously filtered at VCF load time by :func:`_is_symbolic_allele`
(#88, #264).

This module introduces :class:`StructuralVariant`, a :class:`Variant`
subclass that preserves the SV-specific fields (type, end position,
breakend mate, confidence intervals) instead of dropping them. The
class is deliberately minimal — it describes *what the VCF said*, not
*what the biological consequence is*. Effect annotation happens in
:class:`~varcode.annotators.structural_variant.StructuralVariantAnnotator`
(PR 10).

Design notes for future integrations
------------------------------------

The ``StructuralVariant`` shape is intentionally open to these
downstream integrations without forcing them into the core class:

* **Full-genome / personalized-genome support** — a caller with a
  patient-specific reference FASTA can pass a different ``genome=``
  to the SV, just like for point variants. varcode's pyensembl
  genome hook is agnostic to whether the reference is GRCh38 or a
  custom assembly.

* **Long-read-assembly resolution** — callers that have resolved the
  rearranged allele (e.g. via Shasta, hifiasm, or a targeted local
  assembly) can attach the assembled sequence as
  :attr:`alt_assembly`. The SV annotator can prefer that over
  inferring the rearranged allele from breakpoint coordinates alone.

* **Short-read RNA evidence** — evidence from RNA-seq (Isovar-style
  junction reads, exon-skipping counts) attaches through the
  :class:`~varcode.effect_candidates.EffectCandidate` ``evidence`` dict rather than
  fields on the variant itself. This keeps variants immutable and
  lets multiple pieces of evidence accumulate per outcome.

* **External SV callers' extra fields** — CIPOS / CIEND / HOMLEN /
  SVMETHOD etc. go into :attr:`info`, an open-ended dict. The named
  fields (:attr:`sv_type`, :attr:`end`, etc.) are the only ones the
  annotator requires; the rest ride along for downstream consumers.
"""

from typing import Any, Mapping, NamedTuple, Optional, Tuple

from .variant import Variant


# Canonical SV-type codes. Matches VCF 4.3 §5.4 plus common extensions.
SV_TYPES = frozenset({
    "DEL",     # deletion
    "DUP",     # duplication
    "INV",     # inversion
    "INS",     # insertion (often mobile element; see INS:ME:* below)
    "CNV",     # copy number variant (unspecified direction)
    "BND",     # breakend (half of a translocation / complex rearrangement)
})


class Breakend(NamedTuple):
    """One end of a novel adjacency.

    ``position`` is the last base the rearrangement keeps on that side,
    and ``keeps`` says which side that is: ``"left"`` keeps the bases up
    to and including ``position``, ``"right"`` those from ``position``
    on. ``keeps`` is ``None`` when the record didn't say — a breakend
    ALT varcode couldn't read.
    """

    contig: str
    position: int
    keeps: Optional[str] = None


def typed_event_from_breakends(a, b):
    """The DEL / DUP / INV span that two paired breakend records
    describe, or ``None`` when they don't describe one.

    Callers such as esvee and GRIDSS write deletions, duplications and
    inversions as a pair of breakend records labeled ``SVTYPE=DEL`` and
    so on. Given both halves, with contig names already normalized,
    return ``(sv_type, start, end)`` with the coordinates the equivalent
    ``<DEL>`` / ``<DUP>`` / ``<INV>`` record would have, where ``start``
    is the base before the event. Returns ``None`` unless both halves
    agree on a same-contig DEL, DUP or INV label, point at each other,
    and keep sides that fit that label.
    """
    from .sv_allele_parser import breakend_sides
    labels = {
        str((variant.info or {}).get("svtype") or "").upper()
        for variant in (a, b)}
    if len(labels) != 1:
        return None
    (label,) = labels
    if label not in ("DEL", "DUP", "INV"):
        return None
    if a.contig != b.contig or a.start == b.start:
        return None
    if (a.mate_contig, a.mate_start) != (b.contig, b.start):
        return None
    if (b.mate_contig, b.mate_start) != (a.contig, a.start):
        return None
    sides = {}
    for variant in (a, b):
        variant_sides = breakend_sides(variant.symbolic_alt)
        if variant_sides is None:
            return None
        sides[variant.start] = variant_sides[0]
    low, high = sorted((a.start, b.start))
    low_keeps_left = sides[low] == "left"
    high_keeps_left = sides[high] == "left"
    if label == "DEL" and low_keeps_left and not high_keeps_left:
        # Adjacent breakpoints delete nothing — that's an insertion
        # point, not a deletion.
        if high - low > 1:
            return ("DEL", low, high - 1)
    elif label == "DUP" and not low_keeps_left and high_keeps_left:
        return ("DUP", low - 1, high)
    elif label == "INV" and low_keeps_left == high_keeps_left:
        if low_keeps_left:
            return ("INV", low, high)
        return ("INV", low - 1, high - 1)
    return None


class StructuralVariant(Variant):
    """A structural variant — deletion, duplication, inversion,
    insertion, CNV, or breakend — too large or too complex to
    represent as a simple ref/alt nucleotide pair.

    Subclasses :class:`Variant` so ``isinstance(v, Variant)`` still
    works; downstream code that handles variant kinds generically
    (effect collections, serialization) sees a :class:`Variant` and
    the shared contract still applies. The SV-specific fields
    (:attr:`sv_type`, :attr:`end`, breakend mate fields) live here
    and are consulted by SV-aware code.

    The SV position model:

    * :attr:`start` — 1-based start of the affected region
      (matches VCF POS).
    * :attr:`end` — 1-based *inclusive* end. For a DEL/DUP/INV/CNV
      this is the SV endpoint on the same contig. For an INS it
      equals start (insertions are zero-width in reference coords).
      For a BND, ``end == start`` and the other breakpoint lives in
      :attr:`mate_contig` / :attr:`mate_start`.

    Parameters
    ----------
    contig : str
        Chromosome of the (first) breakpoint.
    start : int
        1-based start position.
    sv_type : str
        One of :data:`SV_TYPES`.
    end : int, optional
        1-based inclusive end position. Defaults to ``start`` for
        zero-width SVs (INS, BND).
    alt : str, optional
        Original ALT field from the VCF — ``<DEL>``, ``<INS:ME:ALU>``,
        ``G]17:198982]``, etc. Kept so round-trip to VCF is possible.
        Defaults to ``"<{sv_type}>"``.
    ref : str, optional
        Original REF base (usually one nucleotide, the anchor).
        Defaults to ``"N"``.
    mate_contig : str, optional
        For BND: the mate breakpoint's chromosome. Normalized the same
        way as ``contig`` (e.g. "chr4" -> "4" when converting UCSC names).
    mate_start : int, optional
        For BND: the mate breakpoint's position.
    mate_orientation : str, optional
        For BND: one of ``"[["``, ``"[]"``, ``"][``, ``"]]"``
        encoding the VCF 4.1 breakend strand + direction shorthand
        (first bracket = preceding; second = following). See VCF
        §5.4 for the full grammar.
    ci_start : (int, int), optional
        Confidence interval around ``start`` (VCF CIPOS).
    ci_end : (int, int), optional
        Confidence interval around ``end`` (VCF CIEND).
    alt_assembly : str, optional
        Caller-supplied assembled sequence of the rearranged allele.
        Hook for long-read / targeted-assembly pipelines. The SV
        annotator can prefer this over inferring from breakpoints.
    info : Mapping[str, Any], optional
        Open-ended bag for extra VCF INFO fields the core class
        doesn't model (HOMLEN, SVMETHOD, MATEID, etc.). Kept as
        a Mapping so callers can pass whatever shape their caller
        produces.
    genome, ensembl, normalize_contig_names, convert_ucsc_contig_names
        Same meaning as :class:`Variant`.
    """

    __slots__ = (
        "sv_type",
        "mate_contig",
        "mate_start",
        "mate_orientation",
        "ci_start",
        "ci_end",
        "alt_assembly",
        "info",
        "_sv_alt",
        "_junctions",
        "_locus_transcripts",
    )

    def __init__(
            self,
            contig: str,
            start: int,
            sv_type: str,
            end: Optional[int] = None,
            alt: Optional[str] = None,
            ref: str = "N",
            mate_contig: Optional[str] = None,
            mate_start: Optional[int] = None,
            mate_orientation: Optional[str] = None,
            ci_start: Optional[Tuple[int, int]] = None,
            ci_end: Optional[Tuple[int, int]] = None,
            alt_assembly: Optional[str] = None,
            info: Optional[Mapping[str, Any]] = None,
            genome=None,
            ensembl=None,
            normalize_contig_names: bool = True,
            convert_ucsc_contig_names=None):
        if sv_type not in SV_TYPES:
            raise ValueError(
                "Unknown sv_type %r (expected one of %s)"
                % (sv_type, sorted(SV_TYPES)))
        if end is None:
            end = start

        # Initialize the base Variant with a placeholder ref/alt so the
        # nucleotide-normalization path doesn't reject <DEL>-style
        # symbolic alleles. The original symbolic ALT is preserved in
        # ``_sv_alt``; :attr:`alt` returns it via a property override.
        Variant.__init__(
            self,
            contig=contig,
            start=start,
            ref=ref if ref else "N",
            alt="A" if ref == "N" else "A",  # placeholder, overridden below
            genome=genome,
            ensembl=ensembl,
            allow_extended_nucleotides=True,
            normalize_contig_names=normalize_contig_names,
            convert_ucsc_contig_names=convert_ucsc_contig_names)

        # Override end position — base Variant ignores end for SNVs
        # but we need it as an explicit SV endpoint.
        self.end = int(end)

        # SV-specific fields.
        self.sv_type = sv_type
        self.mate_contig = (
            self._normalize_contig_name(mate_contig)
            if mate_contig is not None else None)
        self.mate_start = int(mate_start) if mate_start is not None else None
        self.mate_orientation = mate_orientation
        self.ci_start = tuple(ci_start) if ci_start is not None else None
        self.ci_end = tuple(ci_end) if ci_end is not None else None
        self.alt_assembly = alt_assembly
        self.info = dict(info) if info is not None else {}

        # Preserve the original symbolic ALT string so round-tripping
        # and downstream consumers see what the VCF said.
        self._sv_alt = alt if alt is not None else "<%s>" % sv_type

        # Caches, like Variant's overlapping-gene / transcript caches:
        # the junction list derived from the fields above, and
        # protein-coding transcripts looked up per locus by the SV
        # annotator (one breakend's mate is queried once per variant
        # rather than once per transcript at the other end).
        self._junctions = None
        self._locus_transcripts = {}

    @property
    def is_structural(self) -> bool:
        return True

    @property
    def symbolic_alt(self) -> str:
        """The original symbolic / breakend ALT string from the VCF."""
        return self._sv_alt

    @property
    def junctions(self) -> Tuple[Tuple[Breakend, Breakend], ...]:
        """The novel adjacencies this variant creates, each a pair of
        :class:`Breakend` ends.

        One junction for a breakend record (its own position joined to
        its mate), and one for a deletion or tandem duplication. A
        symbolic inversion has two, since it joins both of its ends; an
        inversion built from one breakend pair has only the junction
        that pair observed. Empty for insertions, CNVs and single
        breakends, which have no second end to join.

        Positions follow the VCF symbolic-allele convention that
        ``start`` is the base before the event, so a deletion joins
        ``start`` (keeping the left side) to ``end + 1`` (keeping the
        right side), and a tandem duplication joins ``end`` to
        ``start + 1``.
        """
        if self._junctions is None:
            self._junctions = self._compute_junctions()
        return self._junctions

    def _compute_junctions(self):
        # Local import: sv_allele_parser imports this module.
        from .sv_allele_parser import breakend_sides
        sides = breakend_sides(self._sv_alt)
        contig, start, end = self.contig, self.start, self.end
        if self.sv_type == "BND":
            if self.mate_contig is None or self.mate_start is None:
                return ()
            this_side, mate_side = sides if sides else (None, None)
            return ((
                Breakend(contig, start, this_side),
                Breakend(self.mate_contig, self.mate_start, mate_side)),)
        if self.sv_type == "DEL":
            return ((
                Breakend(contig, start, "left"),
                Breakend(contig, end + 1, "right")),)
        if self.sv_type == "DUP":
            return ((
                Breakend(contig, start + 1, "right"),
                Breakend(contig, end, "left")),)
        if self.sv_type == "INV":
            keeping_left = (
                Breakend(contig, start, "left"),
                Breakend(contig, end, "left"))
            keeping_right = (
                Breakend(contig, start + 1, "right"),
                Breakend(contig, end + 1, "right"))
            if sides == ("left", "left"):
                return (keeping_left,)
            if sides == ("right", "right"):
                return (keeping_right,)
            return (keeping_left, keeping_right)
        return ()

    @property
    def breakpoints(self) -> Tuple[Tuple[str, int], ...]:
        """Distinct ``(contig, position)`` breakpoints of this variant's
        junctions, for consumers that read sequence around them
        (cryptic-exon and splice-window scans). Falls back to the
        variant's own start when it has no junction."""
        loci = []
        for junction in self.junctions:
            for breakend in junction:
                locus = (breakend.contig, breakend.position)
                if locus not in loci:
                    loci.append(locus)
        return tuple(loci) or ((self.contig, self.start),)

    @property
    def length(self) -> Optional[int]:
        """Length of the SV span in reference coordinates, when
        well-defined. ``None`` for breakends (the span depends on
        the mate, which may be on another contig)."""
        if self.sv_type == "BND":
            return None
        if self.sv_type == "INS":
            # Pure insertion: the length is the inserted sequence's
            # length, not an interval on the reference. Callers that
            # need the ref-coord span see 0 (INS is zero-width on ref).
            return 0
        return self.end - self.start + 1

    @property
    def short_description(self) -> str:
        if self.sv_type == "BND":
            if self.mate_contig is None:
                # Single breakend, or a symbolic BND with no mate.
                return "BND(%s:%d)" % (self.contig, self.start)
            mate = "%s:%s" % (self.mate_contig, self.mate_start)
            return "BND(%s:%d -> %s)" % (self.contig, self.start, mate)
        return "%s(%s:%d-%d)" % (
            self.sv_type, self.contig, self.start, self.end)

    def __str__(self) -> str:
        return (
            "StructuralVariant(contig=%r, start=%d, end=%d, "
            "sv_type=%r, alt=%r, reference_name=%r)"
        ) % (
            self.contig, self.start, self.end,
            self.sv_type, self._sv_alt, self.reference_name)

    def __repr__(self) -> str:
        return str(self)

    def _sv_identity(self):
        """What the record said, beyond the base variant's locus. The
        base :class:`Variant` identity compares ref/alt, which are
        placeholders on an SV, so two different SV records at one
        position would otherwise compare equal and
        ``load_vcf(distinct=True)`` would drop one."""
        return (
            self.sv_type,
            self.end,
            self._sv_alt,
            self.mate_contig,
            self.mate_start,
            self.mate_orientation,
            self.alt_assembly,
            self.ci_start,
            self.ci_end)

    def __eq__(self, other) -> bool:
        if self is other:
            return True
        if not isinstance(other, StructuralVariant):
            return False
        return (
            Variant.__eq__(self, other)
            and self._sv_identity() == other._sv_identity())

    def __hash__(self) -> int:
        return Variant.__hash__(self)
