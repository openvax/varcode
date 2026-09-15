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

"""Coordinate-aware genomic layouts for realized transcript products.

The existing :class:`~varcode.mutant_transcript.TranscriptEdit` model uses
reference-cDNA offsets.  That is a good compact representation for ordinary
point variants, but it cannot express an exon duplicated twice, an inverted
segment, or a cross-contig adjacency without losing where each base came
from.  This module keeps the rearranged molecule as a short list of lazy
reference intervals.  Sequence is fetched only when a consumer materializes
the layout, so a transcript spanning a very large intron does not allocate
one Python object per genomic base.

Coordinates are one-based and inclusive.  Segment order is molecule order;
``strand`` says how the corresponding reference interval is traversed.
"""

from dataclasses import dataclass, field
from typing import Callable, Optional, Tuple

from .nucleotides import reverse_complement


class SequenceUnavailable(ValueError):
    """Raised when a lazy reference segment is materialized without FASTA."""


@dataclass(frozen=True)
class LayoutSegment:
    """One contiguous piece of a rearranged genomic molecule.

    Reference-derived segments carry ``contig``, ``start`` and ``end``.
    Inserted sequence has all three set to ``None``.  When ``sequence`` is
    present it is already in molecule order; otherwise it is fetched lazily
    and reverse-complemented when ``strand == "-"``.
    """

    contig: Optional[str]
    start: Optional[int]
    end: Optional[int]
    strand: str = "+"
    sequence: Optional[str] = None
    origin_kind: str = "reference"
    source_variant: Optional[object] = field(
        default=None, compare=False, repr=False)

    def __post_init__(self):
        has_coordinates = self.contig is not None
        if has_coordinates != (self.start is not None and self.end is not None):
            raise ValueError(
                "contig, start and end must be all set or all None")
        if self.strand not in ("+", "-"):
            raise ValueError("strand must be '+' or '-', got %r" % self.strand)
        if has_coordinates:
            if self.start < 1 or self.end < self.start:
                raise ValueError(
                    "invalid one-based interval %r:%r-%r" % (
                        self.contig, self.start, self.end))
            if self.sequence is not None and len(self.sequence) != self.length:
                raise ValueError(
                    "segment sequence has length %d, expected %d" % (
                        len(self.sequence), self.length))
        elif not self.sequence:
            raise ValueError("inserted segments require non-empty sequence")

    @property
    def length(self):
        if self.start is None:
            return len(self.sequence)
        return self.end - self.start + 1

    def genomic_position(self, offset):
        """Reference position of a base offset, or ``None`` if inserted."""
        if offset < 0 or offset >= self.length:
            raise IndexError(offset)
        if self.start is None:
            return None
        if self.strand == "+":
            return self.start + offset
        return self.end - offset

    def slice(self, start, end):
        """Return the half-open molecule-order slice ``[start:end]``."""
        if start < 0 or end < start or end > self.length:
            raise IndexError((start, end))
        if start == end:
            return None
        sequence = self.sequence[start:end] if self.sequence is not None else None
        if self.start is None:
            return LayoutSegment(
                contig=None,
                start=None,
                end=None,
                strand=self.strand,
                sequence=sequence,
                origin_kind=self.origin_kind,
                source_variant=self.source_variant)
        first = self.genomic_position(start)
        last = self.genomic_position(end - 1)
        return LayoutSegment(
            contig=self.contig,
            start=min(first, last),
            end=max(first, last),
            strand=self.strand,
            sequence=sequence,
            origin_kind=self.origin_kind,
            source_variant=self.source_variant)

    def reversed(self):
        """Return this segment traversed in the opposite direction."""
        sequence = (
            reverse_complement(self.sequence)
            if self.sequence is not None else None)
        return LayoutSegment(
            contig=self.contig,
            start=self.start,
            end=self.end,
            strand="-" if self.strand == "+" else "+",
            sequence=sequence,
            origin_kind=self.origin_kind,
            source_variant=self.source_variant)

    def materialize(self, sequence_provider=None):
        """Return bases in molecule order."""
        if self.sequence is not None:
            return self.sequence
        if sequence_provider is None:
            raise SequenceUnavailable(
                "No genomic sequence provider for %s:%d-%d" % (
                    self.contig, self.start, self.end))
        sequence = sequence_provider(
            self.contig, self.start, self.end).upper()
        if len(sequence) != self.length:
            raise ValueError(
                "Sequence provider returned %d bases for %s:%d-%d; "
                "expected %d" % (
                    len(sequence), self.contig, self.start, self.end,
                    self.length))
        return reverse_complement(sequence) if self.strand == "-" else sequence


def _coalesce(segments):
    """Merge adjacent compatible lazy reference segments."""
    result = []
    for segment in segments:
        if segment is None or segment.length == 0:
            continue
        if result:
            previous = result[-1]
            same_metadata = (
                previous.contig == segment.contig
                and previous.strand == segment.strand
                and previous.origin_kind == segment.origin_kind
                and previous.source_variant is segment.source_variant
                and previous.sequence is None
                and segment.sequence is None)
            contiguous = False
            if same_metadata:
                contiguous = (
                    previous.end + 1 == segment.start
                    if segment.strand == "+"
                    else segment.end + 1 == previous.start)
            if same_metadata and contiguous:
                result[-1] = LayoutSegment(
                    contig=previous.contig,
                    start=min(previous.start, segment.start),
                    end=max(previous.end, segment.end),
                    strand=previous.strand,
                    origin_kind=previous.origin_kind,
                    source_variant=previous.source_variant)
                continue
        result.append(segment)
    return tuple(result)


@dataclass(frozen=True)
class GenomicLayout:
    """A rearranged molecule represented by ordered interval segments."""

    segments: Tuple[LayoutSegment, ...]
    sequence_provider: Optional[Callable] = field(
        default=None, compare=False, repr=False)

    def __post_init__(self):
        object.__setattr__(self, "segments", _coalesce(tuple(self.segments)))

    @classmethod
    def from_interval(
            cls, contig, start, end, strand="+", sequence_provider=None):
        """Create a lazy reference layout for one genomic interval."""
        return cls(
            segments=(LayoutSegment(contig, start, end, strand=strand),),
            sequence_provider=sequence_provider)

    @classmethod
    def from_transcript(
            cls, transcript, flank=50, sequence_provider=None):
        """Create a transcript-oriented layout around ``transcript``."""
        start = max(1, transcript.start - flank)
        end = transcript.end + flank
        strand = "-" if transcript.on_backward_strand else "+"
        return cls.from_interval(
            transcript.contig, start, end, strand, sequence_provider)

    @property
    def length(self):
        return sum(segment.length for segment in self.segments)

    def materialize(self):
        """Render the complete molecule in layout order."""
        return "".join(
            segment.materialize(self.sequence_provider)
            for segment in self.segments)

    def origins(self):
        """Return ``(contig, position, origin_kind)`` for every base.

        This intentionally expands only on request.  Core layout editing and
        exon projection operate on intervals and remain proportional to the
        number of segments rather than the genomic span.
        """
        return tuple(
            (segment.contig, segment.genomic_position(i), segment.origin_kind)
            for segment in self.segments
            for i in range(segment.length))

    def reversed(self):
        """Reverse-complement the molecule while preserving base origins."""
        return GenomicLayout(
            tuple(segment.reversed() for segment in reversed(self.segments)),
            self.sequence_provider)

    def _locate(self, contig, position, allowed_origin_kinds=None):
        """Locate the first matching reference base in molecule order."""
        offset = 0
        for index, segment in enumerate(self.segments):
            allowed = (
                allowed_origin_kinds is None
                or segment.origin_kind in allowed_origin_kinds)
            if (allowed and segment.contig == contig
                    and segment.start <= position <= segment.end):
                within = (
                    position - segment.start
                    if segment.strand == "+"
                    else segment.end - position)
                return offset + within, index, within, segment
            offset += segment.length
        raise KeyError((contig, position))

    def _split_at(self, offset):
        """Split into segment tuples before and after a molecule offset."""
        if offset < 0 or offset > self.length:
            raise IndexError(offset)
        left = []
        right = []
        consumed = 0
        for segment in self.segments:
            segment_end = consumed + segment.length
            if segment_end <= offset:
                left.append(segment)
            elif consumed >= offset:
                right.append(segment)
            else:
                cut = offset - consumed
                left.append(segment.slice(0, cut))
                right.append(segment.slice(cut, segment.length))
            consumed = segment_end
        return _coalesce(left), _coalesce(right)

    def slice(self, start, end):
        """Return the half-open molecule-order interval ``[start:end]``."""
        if start < 0 or end < start or end > self.length:
            raise IndexError((start, end))
        _, from_start = self._split_at(start)
        middle = GenomicLayout(from_start, self.sequence_provider)
        selected, _ = middle._split_at(end - start)
        return GenomicLayout(selected, self.sequence_provider)

    def replace(self, start, end, replacement=()):
        """Replace a half-open molecule interval with segments."""
        before, remainder = self._split_at(start)
        after_layout = GenomicLayout(remainder, self.sequence_provider)
        _, after = after_layout._split_at(end - start)
        return GenomicLayout(
            before + tuple(replacement) + after,
            self.sequence_provider)

    def _span_offsets(self, contig, start, end):
        allowed = frozenset(("reference", "alternate"))
        start_location = self._locate(contig, start, allowed)
        end_location = self._locate(contig, end, allowed)
        start_offset = start_location[0]
        end_offset = end_location[0]
        low = min(start_offset, end_offset)
        high = max(start_offset, end_offset) + 1
        strand = start_location[3].strand
        if end_location[3].strand != strand:
            raise ValueError(
                "Variant span crosses differently oriented segments")
        return low, high, strand

    def apply_point_variant(self, variant, validate_reference=True):
        """Apply one normalized nucleotide allele to its first occurrence."""
        contig = variant.contig
        ref = variant.trimmed_ref.upper()
        alt = variant.trimmed_alt.upper()
        start = variant.trimmed_base1_start
        if not ref:
            anchor, _, _, segment = self._locate(
                contig, start, frozenset(("reference", "alternate")))
            insert_at = anchor + 1 if segment.strand == "+" else anchor
            sequence = alt if segment.strand == "+" else reverse_complement(alt)
            inserted = LayoutSegment(
                None, None, None, sequence=sequence,
                origin_kind="inserted", source_variant=variant)
            return self.replace(insert_at, insert_at, (inserted,))

        end = variant.trimmed_base1_end
        low, high, strand = self._span_offsets(contig, start, end)
        existing = self.slice(low, high).materialize()
        expected = ref if strand == "+" else reverse_complement(ref)
        if validate_reference and existing.upper() != expected:
            raise ValueError(
                "Reference allele %r does not match realized layout %r "
                "at %s:%d-%d" % (expected, existing, contig, start, end))
        replacement = _point_replacement_segments(
            contig=contig,
            start=start,
            ref=ref,
            alt=alt,
            strand=strand,
            source_variant=variant)
        return self.replace(low, high, replacement)

    def apply_structural_variant(self, variant):
        """Apply a DEL, DUP, INV or sequence-resolved INS.

        ``affected_start`` / ``affected_end`` are authoritative.  They keep
        VCF padding anchors out of the altered span for typed breakend pairs.
        BND is handled by :func:`join_breakends`, since it needs two source
        layouts.
        """
        kind = variant.sv_type
        if kind == "BND":
            raise ValueError("Use join_breakends for BND variants")
        if kind == "CNV":
            raise ValueError("CNV direction is unresolved")
        if kind == "INS":
            sequence = variant.alt_assembly
            if not sequence:
                raise SequenceUnavailable(
                    "INS requires alt_assembly to realize sequence")
            anchor, _, _, segment = self._locate(
                variant.contig, variant.start,
                frozenset(("reference", "alternate")))
            insert_at = anchor + 1 if segment.strand == "+" else anchor
            if segment.strand == "-":
                sequence = reverse_complement(sequence)
            inserted = LayoutSegment(
                None, None, None, sequence=sequence,
                origin_kind="inserted", source_variant=variant)
            return self.replace(insert_at, insert_at, (inserted,))

        start = variant.affected_start
        end = variant.affected_end
        low, high, strand = self._span_offsets(variant.contig, start, end)
        affected = self.slice(low, high)
        if kind == "DEL":
            return self.replace(low, high)
        if kind == "INV":
            inverted = tuple(
                LayoutSegment(
                    contig=segment.contig,
                    start=segment.start,
                    end=segment.end,
                    strand="-" if segment.strand == "+" else "+",
                    sequence=(
                        reverse_complement(segment.sequence)
                        if segment.sequence is not None else None),
                    origin_kind="inverted",
                    source_variant=variant)
                for segment in reversed(affected.segments))
            return self.replace(low, high, inverted)
        if kind == "DUP":
            duplicated = tuple(
                LayoutSegment(
                    contig=segment.contig,
                    start=segment.start,
                    end=segment.end,
                    strand=segment.strand,
                    sequence=segment.sequence,
                    origin_kind="duplicated",
                    source_variant=variant)
                for segment in affected.segments)
            insert_at = high if strand == "+" else low
            return self.replace(insert_at, insert_at, duplicated)
        raise ValueError("Unsupported structural variant type %r" % kind)

    def apply_variants(self, variants, validate_reference=True):
        """Apply variants in descending genomic order.

        Point variants at a higher coordinate are applied before an
        overlapping SV.  Thus a later deletion removes them and a later
        duplication copies the realized allele, matching one-molecule
        haplotype semantics.
        """
        def coordinate(variant):
            if getattr(variant, "is_structural", False):
                return variant.affected_start
            return variant.trimmed_base1_start

        layout = self
        for variant in sorted(variants, key=coordinate, reverse=True):
            if getattr(variant, "is_structural", False):
                layout = layout.apply_structural_variant(variant)
            else:
                layout = layout.apply_point_variant(
                    variant, validate_reference=validate_reference)
        return layout

    def retain_genomic_side(self, contig, position, keeps):
        """Keep the genomic left or right side of a breakend coordinate."""
        if keeps not in ("left", "right"):
            raise ValueError("breakend keeps must be 'left' or 'right'")
        kept = []
        for segment in self.segments:
            if segment.contig != contig:
                continue
            if keeps == "left":
                overlap_start = segment.start
                overlap_end = min(segment.end, position)
            else:
                overlap_start = max(segment.start, position)
                overlap_end = segment.end
            if overlap_start > overlap_end:
                continue
            if segment.strand == "+":
                first = overlap_start - segment.start
                last = overlap_end - segment.start + 1
            else:
                first = segment.end - overlap_end
                last = segment.end - overlap_start + 1
            kept.append(segment.slice(first, last))
        return GenomicLayout(tuple(kept), self.sequence_provider)

    def endpoint_origin(self, first):
        """Return ``(contig, position)`` for the first or last base."""
        segment = self.segments[0] if first else self.segments[-1]
        offset = 0 if first else segment.length - 1
        return segment.contig, segment.genomic_position(offset)


def _point_replacement_segments(
        contig, start, ref, alt, strand, source_variant):
    """Build small origin-aware segments for a point-variant ALT."""
    aligned = min(len(ref), len(alt))
    bases = []
    for index, base in enumerate(alt):
        origin = start + index if index < aligned else None
        bases.append((base, origin))
    if strand == "-":
        bases = [
            (reverse_complement(base), origin)
            for base, origin in reversed(bases)]

    segments = []
    inserted = []
    for base, origin in bases:
        if origin is None:
            inserted.append(base)
            continue
        if inserted:
            segments.append(LayoutSegment(
                None, None, None, sequence="".join(inserted),
                origin_kind="inserted", source_variant=source_variant))
            inserted = []
        segments.append(LayoutSegment(
            contig, origin, origin, strand=strand, sequence=base,
            origin_kind="alternate", source_variant=source_variant))
    if inserted:
        segments.append(LayoutSegment(
            None, None, None, sequence="".join(inserted),
            origin_kind="inserted", source_variant=source_variant))
    return tuple(segments)


def join_breakends(layout_a, breakend_a, layout_b, breakend_b):
    """Join two retained breakend sides into one oriented layout.

    Each retained side is reversed when necessary so ``breakend_a`` is the
    last base of the left piece and ``breakend_b`` the first base of the
    right piece.  This makes the novel adjacency explicit and works for
    same-contig rearrangements as well as cross-contig gene fusions.
    """
    left = layout_a.retain_genomic_side(
        breakend_a.contig, breakend_a.position, breakend_a.keeps)
    right = layout_b.retain_genomic_side(
        breakend_b.contig, breakend_b.position, breakend_b.keeps)
    if not left.segments or not right.segments:
        raise ValueError("Breakend retained an empty layout")
    if left.endpoint_origin(first=False) != (
            breakend_a.contig, breakend_a.position):
        left = left.reversed()
    if right.endpoint_origin(first=True) != (
            breakend_b.contig, breakend_b.position):
        right = right.reversed()
    if left.endpoint_origin(first=False) != (
            breakend_a.contig, breakend_a.position):
        raise ValueError("First breakend is not at the retained-side boundary")
    if right.endpoint_origin(first=True) != (
            breakend_b.contig, breakend_b.position):
        raise ValueError("Second breakend is not at the retained-side boundary")
    provider = layout_a.sequence_provider or layout_b.sequence_provider
    return GenomicLayout(left.segments + right.segments, provider)


def project_intervals(layout, intervals):
    """Project reference intervals onto a layout without expanding bases.

    Returns one :class:`GenomicLayout` per maximal matching run, in molecule
    order.  Inverted segments do not match by default because their
    ``origin_kind`` is not sense-strand transcript sequence.
    """
    intervals_by_contig = {}
    for contig, start, end in intervals:
        intervals_by_contig.setdefault(contig, []).append((start, end))
    pieces = []
    for segment in layout.segments:
        if segment.origin_kind == "inverted" or segment.contig is None:
            continue
        for start, end in intervals_by_contig.get(segment.contig, ()):
            overlap_start = max(start, segment.start)
            overlap_end = min(end, segment.end)
            if overlap_start > overlap_end:
                continue
            if segment.strand == "+":
                first = overlap_start - segment.start
                last = overlap_end - segment.start + 1
            else:
                first = segment.end - overlap_end
                last = segment.end - overlap_start + 1
            pieces.append(GenomicLayout(
                (segment.slice(first, last),), layout.sequence_provider))
    return tuple(pieces)
