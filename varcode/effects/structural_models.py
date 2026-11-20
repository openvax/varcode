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

"""Structural event normalization and ``MutantTranscript`` builders.

Maps local span events (``DEL``, ``DUP``, ``INV``), junctions, translocations,
and caller-supplied assemblies onto :class:`~varcode.MutantTranscript` models
populated with :class:`~varcode.ReferenceSegment` entries.
"""

from dataclasses import dataclass, replace

from ..circular_origin import TranscriptFrame, feature_overlaps, transcript_exons
from ..mutant_transcript import AssembledAllele, MutantTranscript, ReferenceSegment
from ..structural_variant import affected_span, junction_span_details
from ..sv_allele_parser import breakend_local_side
from .effect_helpers import exon_length
from .selenocysteine import segment_selenocysteine
from .transcript_helpers import has_complete_coding_sequence, poly_a_completion


def exon_cdna_ranges(transcript):
    """Yield ``(exon, cdna_start, cdna_end)`` in transcript order.

    Parameters
    ----------
    transcript : Transcript-like
        Spliced transcript.

    Yields
    ------
    tuple of (Exon, int, int)
        Exon object and its half-open ``[cdna_start, cdna_end)`` offsets into
        ``transcript.sequence``.
    """
    offset = 0
    for exon in transcript_exons(transcript):
        length = exon_length(exon)
        yield exon, offset, offset + length
        offset += length


def merge_adjacent_ranges(ranges):
    """Merge abutting sorted half-open ``[(start, end), ...]`` intervals.

    Parameters
    ----------
    ranges : iterable of tuple of (int, int)
        Intervals sorted by ``start``.

    Returns
    -------
    list of tuple of (int, int)
        Merged non-empty intervals with adjacent ``end == next_start`` joined.
    """
    merged = []
    for s, e in ranges:
        if s >= e:
            continue
        if merged and merged[-1][1] == s:
            merged[-1] = (merged[-1][0], e)
        else:
            merged.append((s, e))
    return merged


def feature_contains_breakend(feature, breakend):
    """Return whether ``feature`` spans ``breakend`` (including circular origins)."""
    return feature_overlaps(
        feature, breakend.contig, breakend.position, breakend.position)


def retains_five_prime_end(transcript, side):
    """Return whether keeping ``side`` (``"left"``/``"right"``) retains ``transcript``'s 5' end."""
    return (side == "left") == (transcript.strand == "+")


def cached_variant_lookup(variant, key, compute):
    """Return ``compute()``, memoized on ``variant._annotation_cache`` under ``key``."""
    cache = getattr(variant, "_annotation_cache", None)
    if cache is None:
        return compute()
    if key not in cache:
        cache[key] = compute()
    return cache[key]


@dataclass(frozen=True)
class LocalEvent:
    """Reference span changed by a local DEL/DUP/INV/INS/CNV model and its reported variant.

    A span record is its own event. For the transcripts of a gene that contains
    both ends of a breakend's junction, the event is the DEL, DUP, or INV the
    junction's kept sides describe (#551).
    """

    variant: object
    sv_type: str
    start: int
    end: int
    affected_start: int
    affected_end: int
    junctions: tuple

    @classmethod
    def of(cls, variant):
        start = variant.start
        end = getattr(variant, "end", start)
        return cls(
            variant,
            variant.sv_type,
            start,
            end,
            getattr(variant, "affected_start", start),
            getattr(variant, "affected_end", end),
            tuple(getattr(variant, "junctions", ())),
        )

    @property
    def contig(self):
        return self.variant.contig

    @property
    def alt_assembly(self):
        return getattr(self.variant, "alt_assembly", None)

    @property
    def breakpoints(self):
        """Distinct ``(contig, position)`` junction ends."""
        loci = []
        for junction in self.junctions:
            for breakend in junction:
                if (breakend.contig, breakend.position) not in loci:
                    loci.append((breakend.contig, breakend.position))
        return tuple(loci) or ((self.contig, self.start),)

    def junction_inserted_sequence(self, junction):
        return self.variant.junction_inserted_sequence(junction)


def local_event(variant, transcript):
    """Return the :class:`LocalEvent` ``variant`` makes on ``transcript``, or ``None``."""
    if variant.sv_type != "BND":
        return LocalEvent.of(variant)
    gene = transcript.gene
    frame = TranscriptFrame(transcript)
    for junction in variant.junctions:
        if not all(feature_contains_breakend(gene, end) for end in junction):
            continue
        details = junction_span_details(tuple(
            end._replace(position=frame.position(end.position))
            for end in junction))
        if details is not None:
            return LocalEvent(variant, *details, junctions=(junction,))
    return None


def event_evidence(event):
    """Return ``{"sv_type": ...}`` plus ``"junction_event"`` when typed from a breakend."""
    evidence = {"sv_type": event.variant.sv_type}
    if event.sv_type != event.variant.sv_type:
        evidence["junction_event"] = event.sv_type
    return evidence


def exonic_bases_before(transcript, position):
    """Return ``(cdna_offset, position_is_exonic)`` for ``position`` on ``transcript``."""
    frame = TranscriptFrame(transcript)
    return frame_bases_before(frame, frame.position(position))


def frame_bases_before(frame, position):
    """Return ``(cdna_offset, position_is_exonic)`` for ``position`` in ``frame``."""
    offset = 0
    for start, end in frame.exon_intervals:
        length = end - start + 1
        if frame.backward:
            if start > position:
                offset += length
            elif end >= position:
                return offset + (end - position), True
            else:
                return offset, False
        else:
            if end < position:
                offset += length
            elif start <= position:
                return offset + (position - start), True
            else:
                return offset, False
    return offset, False


def cdna_cut(transcript, position, keep_five_prime):
    """Return the cDNA slice boundary where a junction at ``position`` cuts ``transcript``.

    The 5' side keeps ``sequence[:cut]`` and the 3' side keeps ``sequence[cut:]``.
    An exonic breakpoint keeps that base on whichever side retains it; an intronic
    breakpoint snaps to the adjacent exon boundary.
    """
    offset, exonic = exonic_bases_before(transcript, position)
    return offset + 1 if keep_five_prime and exonic else offset


def deleted_cdna_ranges(variant, transcript):
    """Return sorted half-open cDNA spans removed by ``variant`` on ``transcript``."""
    frame = TranscriptFrame(transcript)
    ranges = []
    for low, high in frame.images(*affected_span(variant)):
        first, last = (high, low) if transcript.on_backward_strand else (low, high)
        start, _ = frame_bases_before(frame, first)
        end, exonic = frame_bases_before(frame, last)
        ranges.append((start, end + int(exonic)))
    return tuple(sorted(ranges))


def contained_image(variant, transcript):
    """Return the single frame image of ``variant`` contained within ``transcript``, or ``None``."""
    frame = TranscriptFrame(transcript)
    inside = [
        (low, high) for low, high in frame.images(*affected_span(variant))
        if frame.start <= low and high <= frame.end]
    return inside[0] if len(inside) == 1 else None


def cdna_ranges_within_sv(variant, transcript):
    """Return merged half-open cDNA body ranges for a local DUP/INV, or ``[]`` if nonlocal."""
    if any(not feature_contains_breakend(transcript, end)
           for junction in variant.junctions for end in junction):
        return []
    image = contained_image(variant, transcript)
    if image is None:
        return []
    sv_start, sv_end = image
    inside = []
    reverse = transcript.on_backward_strand
    frame_exons = TranscriptFrame(transcript).exon_intervals
    for (ex_s, ex_e), (_, c_s, c_e) in zip(
            frame_exons, exon_cdna_ranges(transcript)):
        if ex_e < sv_start or ex_s > sv_end:
            continue
        clipped_s = max(ex_s, sv_start)
        clipped_e = min(ex_e, sv_end)
        if reverse:
            cdna_lo = c_s + (ex_e - clipped_e)
            cdna_hi = c_s + (ex_e - clipped_s + 1)
        else:
            cdna_lo = c_s + (clipped_s - ex_s)
            cdna_hi = c_s + (clipped_e - ex_s + 1)
        inside.append((cdna_lo, cdna_hi))
    return merge_adjacent_ranges(sorted(inside))


def junction_insertion(inserted_sequence, transcript_positions):
    """Return ``(insertion_segments, evidence)`` for oriented junction bases."""
    evidence = {}
    insertion = ()
    if inserted_sequence is None:
        evidence["sequence_status"] = "unresolved_junction_insertion"
    elif inserted_sequence:
        exonic = [any(e.start <= pos <= e.end for e in t.exons)
                  for t, pos in transcript_positions]
        evidence["junction_inserted_sequence"] = inserted_sequence
        if all(exonic):
            insertion = (ReferenceSegment(
                source=AssembledAllele(inserted_sequence), start=0,
                end=len(inserted_sequence), label="junction_insertion"),)
            evidence["junction_insertion_status"] = "retained"
        elif any(exonic):
            evidence["sequence_status"] = "unresolved_insertion_retention"
        else:
            evidence["junction_insertion_status"] = "excluded_by_reference_splicing"
    return insertion, evidence


def local_junction_insertion(variant, transcript):
    """Orient a local DEL/DUP junction in ``transcript`` direction and resolve inserted bases."""
    junction, = variant.junctions
    if variant.junction_inserted_sequence(junction) == "":
        return (), {}
    if not retains_five_prime_end(transcript, junction[0].keeps):
        junction = junction[::-1]
    inserted = variant.junction_inserted_sequence(junction)
    if inserted and not all(feature_contains_breakend(transcript, end) for end in junction):
        return (), {"junction_inserted_sequence": inserted,
                    "sequence_status": "unresolved_insertion_retention"}
    return junction_insertion(
        inserted, ((transcript, end.position) for end in junction))


def build_alt_assembly_mutant_transcript(variant, transcript):
    """Wrap ``variant.alt_assembly`` as a :class:`MutantTranscript` when populated."""
    assembly = getattr(variant, "alt_assembly", None)
    if not assembly:
        return None
    return MutantTranscript.from_sequence(
        assembly,
        reference_transcript=transcript,
        annotator_name="structural_variant",
        evidence={"source": "alt_assembly"},
        label="alt_assembly",
    )


def build_deletion_mutant_transcript(variant, transcript):
    """Build a :class:`MutantTranscript` for a local ``DEL`` structural variant."""
    ranges = deleted_cdna_ranges(variant, transcript)
    insertion, evidence = local_junction_insertion(variant, transcript)
    full_len = len(str(transcript.sequence))
    segments = []
    kept_from = 0
    for index, (start, end) in enumerate(ranges):
        if start > kept_from:
            segments.append(
                ReferenceSegment(transcript, kept_from, start, label="del_kept"))
        if index == len(ranges) - 1:
            segments.extend(insertion)
        kept_from = max(kept_from, end)
    if kept_from < full_len:
        segments.append(
            ReferenceSegment(transcript, kept_from, full_len, label="del_kept"))
    segments = tuple(segments)
    joined = None if "sequence_status" in evidence else "".join(
        s.source.sequence[s.start:s.end] for s in segments)
    return MutantTranscript(
        reference_transcript=transcript,
        reference_segments=segments,
        cdna_sequence=joined,
        evidence=evidence,
        annotator_name="structural_variant",
    )


def build_duplication_mutant_transcript(variant, transcript):
    """Build a tandem-duplication :class:`MutantTranscript` for a local ``DUP`` SV."""
    inside = cdna_ranges_within_sv(variant, transcript)
    if not inside:
        return None
    full_cdna = str(transcript.sequence)
    full_len = len(full_cdna)
    last_inside = inside[-1][1]
    body_segments = tuple(
        ReferenceSegment(
            source=transcript, start=s, end=e, strand="+", label="dup_body_copy")
        for s, e in inside)
    insertion, evidence = local_junction_insertion(variant, transcript)
    segments = (
        ReferenceSegment(
            source=transcript, start=0, end=last_inside,
            strand="+", label="pre_dup_including_body"),
    ) + insertion + body_segments + (
        ReferenceSegment(
            source=transcript, start=last_inside, end=full_len,
            strand="+", label="post_dup"),
    )
    joined = None if "sequence_status" in evidence else "".join(
        s.source.sequence[s.start:s.end] for s in segments)
    return MutantTranscript(
        reference_transcript=transcript,
        reference_segments=segments,
        cdna_sequence=joined,
        evidence=evidence,
        annotator_name="structural_variant",
    )


def build_inversion_mutant_transcript(variant, transcript):
    """Build a segment-mapped :class:`MutantTranscript` for a local ``INV`` SV."""
    inside = cdna_ranges_within_sv(variant, transcript)
    if not inside:
        return None
    first_inside_start = inside[0][0]
    last_inside_end = inside[-1][1]
    full_len = len(str(transcript.sequence))
    segments = (
        ReferenceSegment(
            source=transcript, start=0, end=first_inside_start,
            strand="+", label="pre_inv"),
    ) + tuple(
        ReferenceSegment(
            source=transcript, start=s, end=e,
            strand="-", label="inv_body")
        for s, e in inside
    ) + (
        ReferenceSegment(
            source=transcript, start=last_inside_end, end=full_len,
            strand="+", label="post_inv"),
    )
    return MutantTranscript(
        reference_transcript=transcript,
        reference_segments=segments,
        cdna_sequence=None,
        annotator_name="structural_variant",
    )


def translate_fused_cdna(
        fused_cdna, five_prime_transcript, five_prime_len,
        selenocysteine=(), three_prime_transcript=None):
    """Translate ``fused_cdna`` from the 5' partner's retained CDS start (#336, #701)."""
    from .codon_tables import (
        codon_table_for_transcript,
        initiator_as_met,
        translate_sequence,
    )
    if not has_complete_coding_sequence(five_prime_transcript):
        return None
    try:
        cds_start = min(five_prime_transcript.start_codon_spliced_offsets)
    except (AttributeError, ValueError):
        return None
    if cds_start >= five_prime_len:
        return None
    codon_table = codon_table_for_transcript(five_prime_transcript)
    coding = fused_cdna[cds_start:]
    if three_prime_transcript is not None:
        coding += poly_a_completion(three_prime_transcript, coding)
    truncated = coding[:(len(coding) // 3) * 3]
    try:
        protein = translate_sequence(
            truncated, codon_table=codon_table, to_stop=True,
            selenocysteine={pos - cds_start for pos in selenocysteine})
    except ValueError:
        return None
    return initiator_as_met(protein, coding[:3], codon_table)


def build_fusion_mutant_transcript(
        reference_transcript,
        five_prime_transcript, five_prime_position,
        three_prime_transcript, three_prime_position,
        inserted_sequence="",
        three_prime_start=None):
    """Build a :class:`MutantTranscript` for a :class:`~varcode.effects.GeneFusion` (#336)."""
    five_prime_cdna = five_prime_transcript.sequence
    three_prime_cdna = three_prime_transcript.sequence
    if not five_prime_cdna or not three_prime_cdna:
        return None
    five_prime_cdna = str(five_prime_cdna)
    three_prime_cdna = str(three_prime_cdna)
    five_prime_end = cdna_cut(five_prime_transcript, five_prime_position, True)
    if three_prime_start is None:
        three_prime_start = cdna_cut(three_prime_transcript, three_prime_position, False)
    segments = (
        ReferenceSegment(
            source=five_prime_transcript, start=0, end=five_prime_end,
            strand="+", label="5p_partner"),
        ReferenceSegment(
            source=three_prime_transcript,
            start=three_prime_start, end=len(three_prime_cdna),
            strand="+", label="3p_partner"),
    )
    insertion, evidence = junction_insertion(
        inserted_sequence,
        ((five_prime_transcript, five_prime_position),
         (three_prime_transcript, three_prime_position)),
    )
    segments = segments[:1] + insertion + segments[1:]
    if "sequence_status" in evidence:
        return MutantTranscript(
            reference_transcript=reference_transcript,
            reference_segments=segments,
            evidence=evidence,
            annotator_name="structural_variant",
        )
    fused_cdna = "".join(s.source.sequence[s.start:s.end] for s in segments)
    model = MutantTranscript(
        reference_transcript=reference_transcript,
        reference_segments=segments,
        cdna_sequence=fused_cdna,
        evidence=evidence,
        annotator_name="structural_variant",
    )
    decoded, uncertain, _ = segment_selenocysteine(model) or (set(), set(), set())
    return replace(
        model,
        mutant_protein_sequence=translate_fused_cdna(
            fused_cdna, five_prime_transcript, five_prime_end, decoded | uncertain,
            three_prime_transcript=three_prime_transcript),
    )


def build_translocation_mutant_transcript(variant, transcript):
    """Build a retained-reference-fragment :class:`MutantTranscript` for a BND."""
    side = breakend_local_side(variant.symbolic_alt)
    if (side is None or str(variant.contig) != str(transcript.contig)
            or not TranscriptFrame(transcript).overlaps(variant.start, variant.start)):
        return None
    keep_five_prime = retains_five_prime_end(transcript, side)
    cut = cdna_cut(transcript, variant.start, keep_five_prime)
    full_len = len(str(transcript.sequence))
    start, end = (0, cut) if keep_five_prime else (cut, full_len)
    segments = (
        ReferenceSegment(
            source=transcript, start=start, end=end, strand="+",
            label="translocation_5p" if keep_five_prime else "translocation_3p"),
    )
    return MutantTranscript(
        reference_transcript=transcript,
        reference_segments=segments,
        cdna_sequence=None,
        evidence={"sequence_status": "retained_reference_fragment"},
        annotator_name="structural_variant",
    )
