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
  :class:`~varcode.outcomes.Outcome` ``evidence`` dict rather than
  fields on the variant itself. This keeps variants immutable and
  lets multiple pieces of evidence accumulate per outcome.

* **External SV callers' extra fields** — CIPOS / CIEND / HOMLEN /
  SVMETHOD etc. go into :attr:`info`, an open-ended dict. The named
  fields (:attr:`sv_type`, :attr:`end`, etc.) are the only ones the
  annotator requires; the rest ride along for downstream consumers.
"""

from typing import Any, Mapping, Optional, Tuple

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
        For BND: the mate breakpoint's chromosome.
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
        self.mate_contig = mate_contig
        self.mate_start = int(mate_start) if mate_start is not None else None
        self.mate_orientation = mate_orientation
        self.ci_start = tuple(ci_start) if ci_start is not None else None
        self.ci_end = tuple(ci_end) if ci_end is not None else None
        self.alt_assembly = alt_assembly
        self.info = dict(info) if info is not None else {}

        # Preserve the original symbolic ALT string so round-tripping
        # and downstream consumers see what the VCF said.
        self._sv_alt = alt if alt is not None else "<%s>" % sv_type

    @property
    def is_structural(self) -> bool:
        return True

    @property
    def symbolic_alt(self) -> str:
        """The original symbolic / breakend ALT string from the VCF."""
        return self._sv_alt

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
