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
from .transcript_layout import build_exon_runs


@dataclass(frozen=True)
class SpliceSiteStatus:
    """State of one donor or acceptor on a realized genomic layout."""

    key: SpliceSiteKey
    run_key: tuple
    strength: str
    reason: str


def _base(layout, offset):
    """Return base metadata, tolerating absent tier-0 sequence."""
    try:
        base, contig, position, kind = layout.base(offset)
    except SequenceUnavailable:
        base, contig, position, kind = None, None, None, None
        # Coordinates and origin kind are still available without sequence.
        consumed = 0
        for segment in layout.segments:
            if offset < consumed + segment.length:
                within = offset - consumed
                contig = segment.contig
                position = segment.genomic_position(within)
                kind = segment.origin_kind
                break
            consumed += segment.length
    return base, contig, position, kind


def _at_expected_origin(layout, offset, contig, position):
    if offset < 0 or offset >= layout.length:
        return False
    _, actual_contig, actual_position, kind = _base(layout, offset)
    return (
        actual_contig == contig
        and actual_position == position
        and kind not in ("inserted", "inverted"))


def _motif(layout, offsets, expected):
    bases = [_base(layout, offset)[0] for offset in offsets]
    # Tier 0 knows canonical reference motifs from the annotation. Substitute
    # that expected reference base only where sequence is unavailable; any
    # explicit alternate base still participates in the comparison.
    realized = "".join(
        expected[index] if base is None else base.upper()
        for index, base in enumerate(bases))
    return realized == expected


def donor_status(transcript, layout, run):
    """Return ``(strength, reason)`` for an exon occurrence's donor."""
    direction = -1 if transcript.on_backward_strand else 1
    expected = run.reference_start if direction == -1 else run.reference_end
    last_offset = run.layout_end - 1
    _, contig, position, _ = _base(layout, last_offset)
    if contig != transcript.contig or position != expected:
        return "lost", "exon end removed"
    canonical_offsets = (last_offset + 1, last_offset + 2)
    if not all(
            _at_expected_origin(
                layout, offset, transcript.contig,
                expected + direction * distance)
            for offset, distance in zip(canonical_offsets, (1, 2))):
        return "lost", "donor +1/+2 removed"
    if not _motif(layout, canonical_offsets, "GT"):
        return "lost", "donor GT changed"
    for distance in range(3, 7):
        offset = last_offset + distance
        if not _at_expected_origin(
                layout, offset, transcript.contig,
                expected + direction * distance):
            return "weak_intronic", "donor +%d changed" % distance
        if _base(layout, offset)[3] not in ("reference", "duplicated"):
            return "weak_intronic", "donor +%d changed" % distance
    for distance in range(3):
        offset = last_offset - distance
        if not _at_expected_origin(
                layout, offset, transcript.contig,
                expected - direction * distance):
            return "weak_exonic", "exonic donor region changed"
        if _base(layout, offset)[3] not in ("reference", "duplicated"):
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
    _, contig, position, _ = _base(layout, first_offset)
    if contig != transcript.contig or position != expected:
        return "lost", "exon start removed"
    canonical_offsets = (first_offset - 2, first_offset - 1)
    if not all(
            _at_expected_origin(
                layout, offset, transcript.contig,
                expected - direction * distance)
            for offset, distance in zip(canonical_offsets, (2, 1))):
        return "lost", "acceptor -2/-1 removed"
    if not _motif(layout, canonical_offsets, "AG"):
        return "lost", "acceptor AG changed"
    offset = first_offset - 3
    if (not _at_expected_origin(
            layout, offset, transcript.contig, expected - direction * 3)
            or _base(layout, offset)[3] not in ("reference", "duplicated")):
        return "weak_intronic", "acceptor -3 changed"
    if _base(layout, first_offset)[3] not in ("reference", "duplicated"):
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
                    run.key, strength, reason))
        if run.exon_number > 1:
            strength, reason = acceptor_status(transcript, layout, run)
            if strength != "intact":
                statuses.append(SpliceSiteStatus(
                    SpliceSiteKey(
                        "acceptor", run.exon_number, run.occurrence),
                    run.key, strength, reason))
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
