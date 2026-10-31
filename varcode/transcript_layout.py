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
from .effects.classify import _attach_transcript_models, classify_from_protein_diff
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
    Unresolved,
)


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
    if transcript.on_backward_strand:
        return (
            left.last_genomic_position == left.reference_start
            and right.first_genomic_position == right.reference_end)
    return (
        left.last_genomic_position == left.reference_end
        and right.first_genomic_position == right.reference_start)


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
    # Recognized and source-declared initiators encode Met at this start,
    # even when their internal translation is another residue (e.g. CTG/GTG).
    return "M" + protein[1:], True, start_codon


def _assemble_product(transcript, layout, intervals, junctions, evidence):
    """Assemble and translate one selected path without losing base origins."""
    pieces = tuple(layout.slice(start, end) for start, end in intervals)
    cdna = "".join(piece.materialize() for piece in pieces)
    origins = tuple(origin for piece in pieces for origin in piece.origins())
    target = _start_codon_first_position(transcript)
    start = next((i for i, (contig, position, kind) in enumerate(origins)
                  if contig == transcript.contig and position == target
                  and kind != "inverted"), None) if target is not None else None
    protein, start_present, start_codon = _translate_from_start(
        transcript, cdna, start, pieces)
    return RealizedTranscriptProduct(
        transcript=transcript,
        cdna_sequence=cdna,
        protein_sequence=protein,
        junction_signature=junctions,
        start_codon_present=start_present,
        evidence={
            **evidence,
            "start_codon": start_codon,
            "start_codon_offset": start,
            "genomic_positions": tuple(position for _, position, _ in origins),
            # None alone cannot distinguish a known insertion from unmapped
            # external RNA. Discovery consumes this explicit distinction.
            "inserted_base_offsets": tuple(i for i, (_, _, kind) in enumerate(origins)
                                            if kind == "inserted"),
        })


def realize_exon_path(transcript, layout, kept_run_keys=None):
    """Splice selected exon runs and translate through the shared assembler."""
    runs = build_exon_runs(transcript, layout)
    kept = set(kept_run_keys) if kept_run_keys is not None else None
    selected = runs if kept_run_keys is None else tuple(
        run for run in runs if run.key in kept)
    return _assemble_product(
        transcript, layout,
        tuple((run.layout_start, run.layout_end) for run in selected),
        junction_signature(transcript, selected),
        {"exon_runs": tuple(run.key for run in selected)})


def _best_cryptic_boundary(layout, canonical_offset, side, valid_interval, scan_flank=50):
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
        if not valid_interval[0] <= candidate[0] <= valid_interval[1]:
            continue
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
            if (key.side == "donor" and index + 1 == len(selected)
                    or key.side == "acceptor" and index == 0):
                unresolved.append(option.mechanism)
                continue
            valid_interval = (
                (bounds[run.key][0] + 1, bounds[selected[index + 1].key][0] - 1)
                if key.side == "donor" else
                (bounds[selected[index - 1].key][1] + 1, bounds[run.key][1] - 1))
            canonical = (
                run.layout_end if key.side == "donor" else run.layout_start)
            try:
                cryptic = _best_cryptic_boundary(
                    layout, canonical, key.side, valid_interval, scan_flank)
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

    evidence = {
        "exon_runs": tuple(run.key for run in selected),
        "splice_choices": tuple(
            (key, option.mechanism) for key, option in plan.choices),
        "resolved_splice_choices": tuple(resolved_choices),
    }

    def unavailable(mechanisms):
        mechanisms = tuple(mechanisms)
        return RealizedTranscriptProduct(
            transcript, None, None,
            evidence={**evidence, "sequence_status": "unavailable",
                      "unresolved_mechanisms": mechanisms}), mechanisms

    if unresolved:
        return unavailable(unresolved)
    intervals = _coalesce_intervals(tuple(
        tuple(bounds[run.key]) for run in selected))
    try:
        product = _assemble_product(
            transcript, layout, intervals,
            _junctions_for_bounds(transcript, layout, selected, bounds), evidence)
    except SequenceUnavailable:
        sequence_dependent = tuple(
            option.mechanism for _, option in plan.choices
            if option.mechanism in ("intron_retention", "cryptic"))
        if not sequence_dependent:
            raise
        return unavailable(sequence_dependent)
    return product, ()


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
    last_offset = max(offsets)
    cds_start = min(transcript.start_codon_spliced_offsets)
    stop_offsets = tuple(transcript.stop_codon_spliced_offsets)
    cds_end = max(stop_offsets) + 1 if stop_offsets else (
        cds_start + 3 * len(str(transcript.protein_sequence)))
    mutant_origins = mutant.evidence.get("genomic_positions", ())
    mapped_offsets = [i for i, position in enumerate(mutant_origins)
                      if position is not None and start <= position <= end]
    protein = str(transcript.protein_sequence)
    if mapped_offsets and mutant.evidence.get("start_codon_offset") is not None:
        first_offset = min(mapped_offsets)
        last_offset = max(mapped_offsets)
        cds_start = mutant.evidence["start_codon_offset"]
        protein = mutant.protein_sequence
        cds_end = cds_start + 3 * (len(protein) + 1)
    # On the reverse strand, inserted bases precede their genomic anchor.
    # Anchoring on the retained first CDS base leaves the insertion in the UTR.
    if (first_offset < cds_start
            or (getattr(variant, "is_insertion", False)
                and transcript.on_backward_strand and first_offset == cds_start)):
        return FivePrimeUTR(variant, transcript)
    if first_offset >= cds_end:
        return ThreePrimeUTR(variant, transcript)
    aa_pos = (first_offset - cds_start) // 3
    aa_end = (min(last_offset, cds_end - 1) - cds_start) // 3 + 1
    aa_ref = protein[aa_pos:aa_end] if 0 <= aa_pos < len(protein) else ""
    effect = Silent(variant, transcript, aa_pos=aa_pos, aa_ref=aa_ref)
    if mutant_origins:
        effect.excluded_from_mrna = not any(
            start <= position <= end for position in mutant_origins
            if position is not None)
    return effect


def _unknown_codons(product):
    """Mapped ambiguous codons; equal X strings alone cannot prove equality."""
    start = product.evidence.get("start_codon_offset")
    if start is None or product.cdna_sequence is None:
        return ()
    origins = product.evidence.get("genomic_positions", ())
    return tuple((origins[start + 3 * i:start + 3 * i + 3],
                  product.cdna_sequence[start + 3 * i:start + 3 * i + 3])
                 for i, residue in enumerate(product.protein_sequence or "")
                 if residue == "X")


def classify_products(variant, transcript, baseline, mutant):
    """Classify a mutant product against the patient's baseline product."""
    baseline_model = baseline.to_mutant_transcript()
    mutant_model = mutant.to_mutant_transcript()

    def finish(effect):
        return _attach_transcript_models(
            effect, mutant_model, baseline_model,
            baseline_model.mutant_protein_sequence)

    if (baseline.protein_sequence is None or mutant.protein_sequence is None
            or baseline.start_codon_present is False):
        unresolved = tuple(mutant.evidence.get("unresolved_mechanisms", ())) + tuple(
            baseline.evidence.get("unresolved_mechanisms", ()))
        return finish(Unresolved(
            variant, transcript, mechanism=unresolved[0] if unresolved else "transcript_comparison",
            reason=mutant.evidence.get("unresolved_reason") or
            "A translated baseline and mutant product are required for comparison"))
    if baseline.start_codon_present and not mutant.start_codon_present:
        return finish(StartLoss(variant, transcript, mutant_transcript=mutant_model,
                         ref_codon=baseline.evidence.get("start_codon"),
                         ref_codon_offset=baseline.evidence.get("start_codon_offset")))
    if "X" in baseline.protein_sequence or "X" in mutant.protein_sequence:
        before, after = baseline.protein_sequence, mutant.protein_sequence
        offset = next((i for i, (a, b) in enumerate(zip(before, after)) if a != b),
                      min(len(before), len(after)))
        ambiguous = (before[offset:offset + 1] == "X" or after[offset:offset + 1] == "X"
                     or before == after and _unknown_codons(baseline) != _unknown_codons(mutant))
        if ambiguous:
            return finish(Unresolved(
                variant, transcript, mechanism="ambiguous_coding_sequence",
                reason="Ambiguous codons do not establish the first protein change"))
    if baseline.protein_sequence == mutant.protein_sequence:
        ref_codon = baseline.evidence.get("start_codon")
        alt_codon = mutant.evidence.get("start_codon")
        if (baseline.start_codon_present and mutant.start_codon_present
                and ref_codon and alt_codon and ref_codon != alt_codon):
            return finish(AlternateStartCodon(
                variant, transcript, ref_codon=ref_codon, alt_codon=alt_codon))
        return finish(_unchanged_product_effect(variant, transcript, mutant))
    def coding_sequence(product):
        origins = product.evidence.get("genomic_positions", ())
        target = _start_codon_first_position(transcript)
        start = origins.index(target) if target in origins else 0
        return (product.cdna_sequence or "")[start:], origins[start:]

    reference_coding, reference_origins = coding_sequence(baseline)
    mutant_coding, mutant_origins = coding_sequence(mutant)
    reference_offsets = {p: i for i, p in enumerate(reference_origins) if p is not None}
    shifted_codons = [min(i, reference_offsets[p]) // 3
                      for i, p in enumerate(mutant_origins)
                      if p in reference_offsets and (i - reference_offsets[p]) % 3]
    first_changed = next((i for i, (a, b) in enumerate(zip(reference_coding, mutant_coding))
                          if a != b), min(len(reference_coding), len(mutant_coding)))
    frame_change = min(shifted_codons) if shifted_codons else first_changed // 3
    return classify_from_protein_diff(
        variant=variant,
        transcript=transcript,
        ref_protein=baseline.protein_sequence or "",
        mut_protein=mutant.protein_sequence or "",
        length_delta=len(mutant_coding) - len(reference_coding),
        frame_change_aa_offset=frame_change,
        mutant_transcript=mutant_model,
        baseline_transcript=baseline_model,
        coding_start=mutant.evidence.get("start_codon_offset"))
