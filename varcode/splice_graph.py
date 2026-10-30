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

"""Splice-site state and graph-valid mechanism enumeration on layouts."""

from dataclasses import dataclass

from .effect_hypotheses import (
    SpliceAxis,
    SpliceOption,
    SpliceSiteKey,
    enumerate_splice_plans,
)
from .genomic_layout import SequenceUnavailable
from .errors import HypothesisLimitError
from .nucleotides import reverse_complement
from .transcript_layout import build_exon_runs


@dataclass(frozen=True)
class SpliceSiteStatus:
    """State of one donor or acceptor on a realized genomic layout."""

    key: SpliceSiteKey
    run_key: tuple
    strength: str
    reason: str
    context: tuple = ()


def _base(layout, offset):
    """Return base metadata, tolerating absent tier-0 sequence."""
    try:
        base, contig, position, kind = layout.base(offset)
    except SequenceUnavailable:
        base = None
        contig, position, kind = layout.origin(offset)
    return base, contig, position, kind


def _at_expected_origin(layout, offset, contig, position):
    if offset < 0 or offset >= layout.length:
        return False
    actual_contig, actual_position, kind = layout.origin(offset)
    return (
        actual_contig == contig
        and actual_position == position
        and kind not in ("inserted", "inverted"))


def _motif_changed(transcript, layout, offsets, canonical):
    """Compare edited bases with the annotated site, including GC/AT-AC.

    Annotation establishes an unedited site's existence even without genomic
    FASTA. Canonical bases are only a fallback for explicitly edited positions
    whose reference bases are unavailable, never a test of an unedited site.
    """
    for offset, fallback in zip(offsets, canonical):
        contig, position, kind = layout.origin(offset)
        if kind in ("reference", "duplicated"):
            continue
        base = _base(layout, offset)[0]
        reference = fallback
        # VCF REF is evidence even when no chromosome FASTA is attached.
        segment = layout.slice(offset, offset + 1).segments[0]
        source = segment.source_variant
        if source is not None and not getattr(source, "is_structural", False):
            source_offset = position - source.trimmed_base1_start
            reference = source.trimmed_ref[source_offset:source_offset + 1] or fallback
            if transcript.on_backward_strand:
                reference = reverse_complement(reference)
        if layout.sequence_provider is not None:
            try:
                reference = layout.sequence_provider(contig, position, position)
            except SequenceUnavailable:
                pass
            else:
                if transcript.on_backward_strand:
                    reference = reverse_complement(reference)
        if base is not None and base.upper() != reference.upper():
            return True
    return False


def _site_context(layout, run, side):
    """Comparable sequence/origin context for coupling unchanged sites."""
    offsets = (range(run.layout_end - 3, run.layout_end + 6)
               if side == "donor" else
               range(run.layout_start - 3, run.layout_start + 1))
    return tuple(
        _base(layout, offset)[:3] if 0 <= offset < layout.length else None
        for offset in offsets)


def donor_status(transcript, layout, run):
    """Return ``(strength, reason)`` for an exon occurrence's donor."""
    direction = -1 if transcript.on_backward_strand else 1
    expected = run.reference_start if direction == -1 else run.reference_end
    last_offset = run.layout_end - 1
    contig, position, _ = layout.origin(last_offset)
    if contig != transcript.contig or position != expected:
        return "lost", "exon end removed"
    canonical_offsets = (last_offset + 1, last_offset + 2)
    if not all(
            _at_expected_origin(
                layout, offset, transcript.contig,
                expected + direction * distance)
            for offset, distance in zip(canonical_offsets, (1, 2))):
        return "lost", "donor +1/+2 removed"
    if _motif_changed(transcript, layout, canonical_offsets, "GT"):
        return "lost", "donor motif changed"
    for distance in range(3, 7):
        offset = last_offset + distance
        if not _at_expected_origin(
                layout, offset, transcript.contig,
                expected + direction * distance):
            return "weak_intronic", "donor +%d changed" % distance
        if layout.origin(offset)[2] not in ("reference", "duplicated"):
            return "weak_intronic", "donor +%d changed" % distance
    for distance in range(3):
        offset = last_offset - distance
        if not _at_expected_origin(
                layout, offset, transcript.contig,
                expected - direction * distance):
            return "weak_exonic", "exonic donor region changed"
        if layout.origin(offset)[2] not in ("reference", "duplicated"):
            reason = (
                "last exon base changed" if distance == 0
                else "exonic donor region changed")
            return "weak_exonic", reason
    return "intact", ""


def acceptor_status(transcript, layout, run):
    """Return ``(strength, reason)`` for an exon occurrence's acceptor."""
    direction = -1 if transcript.on_backward_strand else 1
    expected = run.reference_end if direction == -1 else run.reference_start
    first_offset = run.layout_start
    contig, position, _ = layout.origin(first_offset)
    if contig != transcript.contig or position != expected:
        return "lost", "exon start removed"
    canonical_offsets = (first_offset - 2, first_offset - 1)
    if not all(
            _at_expected_origin(
                layout, offset, transcript.contig,
                expected - direction * distance)
            for offset, distance in zip(canonical_offsets, (2, 1))):
        return "lost", "acceptor -2/-1 removed"
    if _motif_changed(transcript, layout, canonical_offsets, "AG"):
        return "lost", "acceptor motif changed"
    offset = first_offset - 3
    if (not _at_expected_origin(
            layout, offset, transcript.contig, expected - direction * 3)
            or layout.origin(offset)[2] not in ("reference", "duplicated")):
        return "weak_intronic", "acceptor -3 changed"
    if layout.origin(first_offset)[2] not in ("reference", "duplicated"):
        return "weak_exonic", "first exon base changed"
    return "intact", ""


def disrupted_splice_sites(transcript, layout, runs=None):
    """Return every non-intact donor/acceptor in the realized layout."""
    if runs is None:
        runs = build_exon_runs(transcript, layout)
    n_exons = len(transcript.exons)
    statuses = []
    for run in runs:
        if run.exon_number < n_exons:
            strength, reason = donor_status(transcript, layout, run)
            if strength != "intact":
                statuses.append(SpliceSiteStatus(
                    SpliceSiteKey("donor", run.exon_number, run.occurrence),
                    run.key, strength, reason, _site_context(layout, run, "donor")))
        if run.exon_number > 1:
            strength, reason = acceptor_status(transcript, layout, run)
            if strength != "intact":
                statuses.append(SpliceSiteStatus(
                    SpliceSiteKey(
                        "acceptor", run.exon_number, run.occurrence),
                    run.key, strength, reason, _site_context(layout, run, "acceptor")))
    return tuple(sorted(statuses, key=lambda status: status.key))


_MECHANISM_ORDER = {
    ("donor", "lost"): (
        "exon_skip", "intron_retention", "cryptic", "normal"),
    ("acceptor", "lost"): (
        "exon_skip", "intron_retention", "cryptic", "normal"),
    ("donor", "weak_intronic"): (
        "normal", "exon_skip", "intron_retention", "cryptic"),
    ("acceptor", "weak_intronic"): (
        "normal", "exon_skip", "intron_retention", "cryptic"),
    ("donor", "weak_exonic"): (
        "normal", "exon_skip", "cryptic", "intron_retention"),
    ("acceptor", "weak_exonic"): (
        "normal", "exon_skip", "cryptic", "intron_retention"),
}


def splice_axes(transcript, statuses, run_keys=None):
    """Convert disrupted sites into ordinal, uncalibrated graph choices."""
    axes = []
    n_exons = len(transcript.exons)
    run_keys = tuple(run_keys or ())
    for status in statuses:
        mechanisms = _MECHANISM_ORDER[(status.key.side, status.strength)]
        options = []
        for rank, mechanism in enumerate(mechanisms):
            if (mechanism == "exon_skip"
                    and status.key.exon_number in (1, n_exons)):
                continue
            kwargs = {
                "mechanism": mechanism,
                "ordinal_rank": rank,
                "evidence": {
                    "strength": status.strength,
                    "reason": status.reason,
                },
            }
            if mechanism == "exon_skip":
                kwargs["skipped_runs"] = (status.run_key,)
            else:
                required = [status.run_key]
                if mechanism == "intron_retention" and status.run_key in run_keys:
                    run_index = run_keys.index(status.run_key)
                    neighbor_index = (
                        run_index + 1
                        if status.key.side == "donor" else run_index - 1)
                    if 0 <= neighbor_index < len(run_keys):
                        required.append(run_keys[neighbor_index])
                kwargs["required_runs"] = tuple(required)
            options.append(SpliceOption(**kwargs))
        axes.append(SpliceAxis(
            key=status.key,
            anchor_runs=(status.run_key,),
            options=tuple(options)))
    return tuple(axes)


def enumerate_layout_splice_plans(transcript, layout, max_plans=64):
    """Detect disrupted sites and enumerate only graph-valid plans."""
    runs = build_exon_runs(transcript, layout)
    statuses = disrupted_splice_sites(transcript, layout, runs)
    run_keys = tuple(run.key for run in runs)
    axes = splice_axes(transcript, statuses, run_keys)
    if not axes:
        return enumerate_splice_plans((), run_keys, 1)
    return enumerate_splice_plans(
        axes, run_keys, max_plans=max_plans)


def enumerate_paired_splice_plans(
        transcript, baseline_statuses, mutant_statuses,
        baseline_run_keys, mutant_run_keys, shared_keys, max_plans=64):
    """Preserve uncertain baseline choices and couple unchanged active sites.

    Each side is graph-valid independently. Shared sites use the same mechanism
    when present in both plans. If a rewrite removes the site on one side, the
    other side retains all compatible choices instead of an arbitrary default.
    The limit applies to complete comparisons, never a truncated candidate set.
    """
    baseline_plans = enumerate_splice_plans(
        splice_axes(transcript, baseline_statuses, baseline_run_keys),
        baseline_run_keys, max_plans=max_plans)
    mutant_plans = enumerate_splice_plans(
        splice_axes(transcript, mutant_statuses, mutant_run_keys),
        mutant_run_keys, max_plans=max_plans)
    pairs = []
    for mutant in mutant_plans:
        mutant_choices = dict(mutant.choices)
        for baseline in baseline_plans:
            baseline_choices = dict(baseline.choices)
            active_shared = shared_keys & baseline_choices.keys() & mutant_choices.keys()
            if any(baseline_choices[key].mechanism != mutant_choices[key].mechanism
                   for key in active_shared):
                continue
            pairs.append((baseline, mutant))
            if len(pairs) > max_plans:
                raise HypothesisLimitError(max_plans)
    return tuple(pairs)
