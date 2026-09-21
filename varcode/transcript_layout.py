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

"""Projection and translation of coordinate-aware genomic layouts."""

from dataclasses import dataclass

from .effect_hypotheses import RealizedTranscriptProduct
from .cryptic_exons import (
    ACCEPTOR_WINDOW,
    DONOR_WINDOW,
    score_acceptor,
    score_donor,
)
from .genomic_layout import SequenceUnavailable
from .effects.classify import classify_from_protein_diff
from .effects.codon_tables import (
    codon_table_for_transcript,
    translate_sequence,
)
from .effects.selenocysteine import layout_selenocysteine
from .effects.effect_classes import (
    AlternateStartCodon,
    FivePrimeUTR,
    Intronic,
    Silent,
    StartLoss,
    ThreePrimeUTR,
)
from .mutant_transcript import MutantTranscript


@dataclass(frozen=True)
class ExonRun:
    """One occurrence of a sense-strand reference exon in a layout."""

    key: tuple
    exon_number: int
    occurrence: int
    layout_start: int
    layout_end: int
    reference_start: int
    reference_end: int
    first_genomic_position: int
    last_genomic_position: int

    @property
    def length(self):
        return self.layout_end - self.layout_start


@dataclass(frozen=True)
class _ExonFragment:
    exon_number: int
    layout_start: int
    layout_end: int
    reference_start: int
    reference_end: int
    first_genomic_position: int
    last_genomic_position: int
    copy_group: str


def _ordered_exons(transcript):
    return tuple(sorted(
        transcript.exons,
        key=lambda exon: exon.start,
        reverse=transcript.on_backward_strand))


def _gap_is_inserted(layout, start, end):
    if start == end:
        return True
    gap = layout.slice(start, end)
    return bool(gap.segments) and all(
        segment.contig is None for segment in gap.segments)


def build_exon_runs(transcript, layout):
    """Project sense-strand exon occurrences onto a genomic layout.

    Point edits and internal deletions remain part of one exon run.  A
    duplicated exon is a distinct occurrence, while inverted (antisense)
    sequence is intentionally excluded.
    """
    exons = _ordered_exons(transcript)
    transcript_strand = "-" if transcript.on_backward_strand else "+"
    fragments = []
    layout_offset = 0
    for segment in layout.segments:
        if (segment.contig == transcript.contig
                and segment.strand == transcript_strand
                and segment.origin_kind != "inverted"):
            for exon_number, exon in enumerate(exons, 1):
                overlap_start = max(segment.start, exon.start)
                overlap_end = min(segment.end, exon.end)
                if overlap_start > overlap_end:
                    continue
                if segment.strand == "+":
                    within_start = overlap_start - segment.start
                    within_end = overlap_end - segment.start + 1
                    first_position = overlap_start
                    last_position = overlap_end
                else:
                    within_start = segment.end - overlap_end
                    within_end = segment.end - overlap_start + 1
                    first_position = overlap_end
                    last_position = overlap_start
                copy_group = (
                    "duplicated"
                    if segment.origin_kind == "duplicated" else "original")
                fragments.append(_ExonFragment(
                    exon_number=exon_number,
                    layout_start=layout_offset + within_start,
                    layout_end=layout_offset + within_end,
                    reference_start=exon.start,
                    reference_end=exon.end,
                    first_genomic_position=first_position,
                    last_genomic_position=last_position,
                    copy_group=copy_group))
        layout_offset += segment.length

    fragments.sort(key=lambda fragment: fragment.layout_start)
    merged = []
    for fragment in fragments:
        if merged:
            previous = merged[-1]
            if (previous.exon_number == fragment.exon_number
                    and previous.copy_group == fragment.copy_group
                    and _gap_is_inserted(
                        layout, previous.layout_end, fragment.layout_start)):
                merged[-1] = _ExonFragment(
                    exon_number=previous.exon_number,
                    layout_start=previous.layout_start,
                    layout_end=fragment.layout_end,
                    reference_start=previous.reference_start,
                    reference_end=previous.reference_end,
                    first_genomic_position=previous.first_genomic_position,
                    last_genomic_position=fragment.last_genomic_position,
                    copy_group=previous.copy_group)
                continue
        merged.append(fragment)

    occurrences = {}
    runs = []
    for fragment in merged:
        occurrence = occurrences.get(fragment.exon_number, 0) + 1
        occurrences[fragment.exon_number] = occurrence
        runs.append(ExonRun(
            key=(fragment.exon_number, occurrence),
            exon_number=fragment.exon_number,
            occurrence=occurrence,
            layout_start=fragment.layout_start,
            layout_end=fragment.layout_end,
            reference_start=fragment.reference_start,
            reference_end=fragment.reference_end,
            first_genomic_position=fragment.first_genomic_position,
            last_genomic_position=fragment.last_genomic_position))
    return tuple(runs)


def _is_canonical_junction(transcript, left, right):
    if left.occurrence != 1 or right.occurrence != 1:
        return False
    if right.exon_number != left.exon_number + 1:
        return False
    left_exon = _ordered_exons(transcript)[left.exon_number - 1]
    right_exon = _ordered_exons(transcript)[right.exon_number - 1]
    if transcript.on_backward_strand:
        return (
            left.last_genomic_position == left_exon.start
            and right.first_genomic_position == right_exon.end)
    return (
        left.last_genomic_position == left_exon.end
        and right.first_genomic_position == right_exon.start)


def junction_signature(transcript, runs):
    """Return non-reference exon joins in transcript order."""
    signature = []
    for left, right in zip(runs, runs[1:]):
        if _is_canonical_junction(transcript, left, right):
            continue
        signature.append((
            transcript.contig,
            left.last_genomic_position,
            transcript.contig,
            right.first_genomic_position,
        ))
    return tuple(signature)


def _start_codon_first_position(transcript):
    positions = tuple(transcript.start_codon_positions)
    if not positions:
        return None
    return max(positions) if transcript.on_backward_strand else min(positions)


def _start_offset(layout, selected_runs, transcript):
    target = _start_codon_first_position(transcript)
    if target is None:
        return None
    cdna_offset = 0
    for run in selected_runs:
        run_layout = layout.slice(run.layout_start, run.layout_end)
        for index, (_, contig, position, kind) in enumerate(
                run_layout.base(i) for i in range(run_layout.length)):
            if (contig == transcript.contig and position == target
                    and kind != "inverted"):
                return cdna_offset + index
        cdna_offset += run.length
    return None


def _translate_from_start(transcript, cdna, start, pieces):
    """Return protein, start presence and codon at the mapped CDS start."""
    if start is None:
        return "", False, None
    coding = cdna[start:]
    start_codon = coding[:3].upper()
    codon_table = codon_table_for_transcript(transcript)
    if start_codon not in codon_table.start_codons:
        return "", False, start_codon
    coding = coding[:len(coding) // 3 * 3]
    protein = translate_sequence(
        coding,
        codon_table=codon_table,
        to_stop=True,
        selenocysteine={pos - start for pos in
                        layout_selenocysteine(transcript, pieces)})
    # NCBI's recognized initiation codons encode Met at the CDS start,
    # even when their internal translation is another residue (e.g. CTG).
    return "M" + protein[1:], True, start_codon


def realize_exon_path(transcript, layout, kept_run_keys=None):
    """Splice selected exon runs and translate from the mapped CDS start."""
    runs = build_exon_runs(transcript, layout)
    if kept_run_keys is None:
        selected = runs
    else:
        kept = set(kept_run_keys)
        selected = tuple(run for run in runs if run.key in kept)
    run_layouts = tuple(
        layout.slice(run.layout_start, run.layout_end) for run in selected)
    cdna = "".join(run_layout.materialize() for run_layout in run_layouts)
    genomic_positions = tuple(
        position
        for run_layout in run_layouts
        for _, position, _ in run_layout.origins())
    start = _start_offset(layout, selected, transcript)
    protein, start_present, start_codon = _translate_from_start(
        transcript, cdna, start, run_layouts)
    return RealizedTranscriptProduct(
        transcript=transcript,
        cdna_sequence=cdna,
        protein_sequence=protein,
        junction_signature=junction_signature(transcript, selected),
        start_codon_present=start_present,
        evidence={
            "start_codon": start_codon,
            "exon_runs": tuple(run.key for run in selected),
            "genomic_positions": genomic_positions,
        })


def _best_cryptic_boundary(layout, canonical_offset, side, scan_flank=50):
    """Return the best non-canonical interbase offset and motif score."""
    scan_start = max(0, canonical_offset - scan_flank)
    scan_end = min(layout.length, canonical_offset + scan_flank)
    sequence = layout.slice(scan_start, scan_end).materialize()
    scorer = score_donor if side == "donor" else score_acceptor
    window = DONOR_WINDOW if side == "donor" else ACCEPTOR_WINDOW
    canonical_local = canonical_offset - scan_start
    best = None
    for index in range(len(sequence) - window + 1):
        score = scorer(sequence[index:index + window])
        if score <= 0:
            continue
        boundary = index + 3
        if boundary == canonical_local:
            continue
        candidate = (scan_start + boundary, score)
        if best is None or candidate[1] > best[1]:
            best = candidate
    return best


def _coalesce_intervals(intervals):
    result = []
    for start, end in sorted(intervals):
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(result[-1][1], end))
        else:
            result.append((start, end))
    return tuple(result)


def _junctions_for_bounds(transcript, layout, selected, bounds):
    signature = []
    for left, right in zip(selected, selected[1:]):
        left_end = bounds[left.key][1]
        right_start = bounds[right.key][0]
        if left_end >= right_start:
            continue
        left_contig, left_position, _ = layout.origin(left_end - 1)
        right_contig, right_position, _ = layout.origin(right_start)
        canonical = (
            left_end == left.layout_end
            and right_start == right.layout_start
            and _is_canonical_junction(transcript, left, right))
        if not canonical:
            signature.append((
                left_contig, left_position, right_contig, right_position))
    return tuple(signature)


def realize_splice_plan(
        transcript, layout, plan, statuses, scan_flank=50):
    """Realize exon-skip, retention and cryptic choices on one layout.

    Returns ``(product, unresolved_mechanisms)``. Sequence-free tier 0 can
    still realize canonical and exon-skip paths; sequence-dependent choices
    remain explicit instead of being guessed.
    """
    runs = build_exon_runs(transcript, layout)
    run_by_key = {run.key: run for run in runs}
    selected = [run_by_key[key] for key in plan.kept_runs if key in run_by_key]
    bounds = {run.key: [run.layout_start, run.layout_end] for run in selected}
    status_by_key = {status.key: status for status in statuses}
    unresolved = []
    resolved_choices = []

    for key, option in plan.choices:
        if option.mechanism in ("normal", "exon_skip"):
            continue
        status = status_by_key.get(key)
        run = run_by_key.get(status.run_key) if status is not None else None
        if run is None or run.key not in bounds:
            unresolved.append(option.mechanism)
            continue
        index = selected.index(run)
        if option.mechanism == "intron_retention":
            if key.side == "donor" and index + 1 < len(selected):
                bounds[run.key][1] = selected[index + 1].layout_start
            elif key.side == "acceptor" and index > 0:
                bounds[run.key][0] = selected[index - 1].layout_end
            else:
                unresolved.append(option.mechanism)
            continue
        if option.mechanism == "cryptic":
            canonical = (
                run.layout_end if key.side == "donor" else run.layout_start)
            try:
                cryptic = _best_cryptic_boundary(
                    layout, canonical, key.side, scan_flank)
            except (SequenceUnavailable, KeyError, ValueError):
                cryptic = None
            if cryptic is None:
                unresolved.append(option.mechanism)
                continue
            boundary, score = cryptic
            if key.side == "donor":
                bounds[run.key][1] = boundary
            else:
                bounds[run.key][0] = boundary
            _, position, _ = layout.origin(
                boundary - 1 if key.side == "donor" else boundary)
            resolved_choices.append((key, option.mechanism, position, score))
            continue
        unresolved.append(option.mechanism)

    if unresolved:
        return realize_exon_path(transcript, layout, plan.kept_runs), tuple(
            unresolved)

    intervals = _coalesce_intervals(tuple(
        tuple(bounds[run.key]) for run in selected))
    pieces = tuple(layout.slice(start, end) for start, end in intervals)
    try:
        cdna = "".join(piece.materialize() for piece in pieces)
    except SequenceUnavailable:
        sequence_dependent = tuple(
            option.mechanism
            for _, option in plan.choices
            if option.mechanism in ("intron_retention", "cryptic"))
        return realize_exon_path(
            transcript, layout, plan.kept_runs), sequence_dependent
    genomic_positions = tuple(
        position
        for piece in pieces
        for _, position, _ in piece.origins())

    target = _start_codon_first_position(transcript)
    start_offset = None
    if target is not None:
        consumed = 0
        for piece in pieces:
            origins = piece.origins()
            for index, (contig, position, kind) in enumerate(origins):
                if (contig == transcript.contig and position == target
                        and kind != "inverted"):
                    start_offset = consumed + index
                    break
            if start_offset is not None:
                break
            consumed += piece.length
    protein, start_present, start_codon = _translate_from_start(
        transcript, cdna, start_offset, pieces)
    evidence = {
        "start_codon": start_codon,
        "exon_runs": tuple(run.key for run in selected),
        "genomic_positions": genomic_positions,
        "splice_choices": tuple(
            (key, option.mechanism) for key, option in plan.choices),
        "resolved_splice_choices": tuple(resolved_choices),
    }
    return RealizedTranscriptProduct(
        transcript=transcript,
        cdna_sequence=cdna,
        protein_sequence=protein,
        junction_signature=_junctions_for_bounds(
            transcript, layout, selected, bounds),
        start_codon_present=start_present,
        evidence=evidence), ()


def _variant_interval(variant):
    if getattr(variant, "is_structural", False):
        return variant.affected_start, variant.affected_end
    return variant.trimmed_base1_start, variant.trimmed_base1_end


def _unchanged_product_effect(variant, transcript, mutant):
    """Classify by genomic location when the translated protein is unchanged."""
    start, end = _variant_interval(variant)
    exons = _ordered_exons(transcript)
    overlapping_positions = []
    for exon in exons:
        overlap_start = max(start, exon.start)
        overlap_end = min(end, exon.end)
        if overlap_start <= overlap_end:
            overlapping_positions.extend((overlap_start, overlap_end))
    if not overlapping_positions:
        return Intronic(variant, transcript, nearest_exon=None,
                        distance_to_exon=0)

    offsets = []
    for position in overlapping_positions:
        try:
            offsets.append(transcript.spliced_offset(position))
        except (KeyError, ValueError):
            pass
    if not offsets:
        return Intronic(variant, transcript, nearest_exon=None,
                        distance_to_exon=0)
    first_offset = min(offsets)
    cds_start = min(transcript.start_codon_spliced_offsets)
    stop_offsets = tuple(transcript.stop_codon_spliced_offsets)
    cds_end = max(stop_offsets) + 1 if stop_offsets else (
        cds_start + 3 * len(str(transcript.protein_sequence)))
    # On the reverse strand, inserted bases precede their genomic anchor.
    # Anchoring on the retained first CDS base leaves the insertion in the UTR.
    if (first_offset < cds_start
            or (getattr(variant, "is_insertion", False)
                and transcript.on_backward_strand and first_offset == cds_start)):
        return FivePrimeUTR(variant, transcript)
    if first_offset >= cds_end:
        return ThreePrimeUTR(variant, transcript)
    aa_pos = (first_offset - cds_start) // 3
    reference_protein = str(transcript.protein_sequence)
    aa_ref = (
        reference_protein[aa_pos]
        if 0 <= aa_pos < len(reference_protein) else "")
    effect = Silent(variant, transcript, aa_pos=aa_pos, aa_ref=aa_ref)
    mutant_origins = mutant.evidence.get("genomic_positions", ())
    if mutant_origins:
        effect.excluded_from_mrna = not any(
            start <= position <= end for position in mutant_origins
            if position is not None)
    return effect


def classify_products(variant, transcript, baseline, mutant):
    """Classify a mutant product against the patient's baseline product."""
    if baseline.start_codon_present and not mutant.start_codon_present:
        return StartLoss(variant, transcript)
    if baseline.protein_sequence == mutant.protein_sequence:
        ref_codon = baseline.evidence.get("start_codon")
        alt_codon = mutant.evidence.get("start_codon")
        if (baseline.start_codon_present and mutant.start_codon_present
                and ref_codon and alt_codon and ref_codon != alt_codon):
            return AlternateStartCodon(
                variant, transcript, ref_codon=ref_codon, alt_codon=alt_codon)
        return _unchanged_product_effect(variant, transcript, mutant)
    mutant_transcript = MutantTranscript(
        reference_transcript=transcript,
        cdna_sequence=mutant.cdna_sequence,
        mutant_protein_sequence=mutant.protein_sequence,
        annotator_name="transcript_model")
    return classify_from_protein_diff(
        variant=variant,
        transcript=transcript,
        ref_protein=baseline.protein_sequence or "",
        mut_protein=mutant.protein_sequence or "",
        length_delta=(
            len(mutant.cdna_sequence or "")
            - len(baseline.cdna_sequence or "")),
        mutant_transcript=mutant_transcript)
