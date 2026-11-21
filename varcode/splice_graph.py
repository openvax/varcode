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

from .circular_origin import exon_number_before_origin, origin_partners
from .effect_hypotheses import (
    SpliceAxis,
    SpliceOption,
    SpliceSiteKey,
    enumerate_splice_plans,
)
from .genomic_layout import SequenceUnavailable
from .errors import HypothesisLimitError
from .nucleotides import reverse_complement
from .transcript_layout import (
    DEFAULT_CRYPTIC_SPLICE_POLICY,
    build_exon_runs,
    enumerate_cryptic_boundaries,
)


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
            source_offset = (
                layout.genomic(contig, position) - source.trimmed_base1_start)
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
    # The exons on either side of a circular origin meet without an intron.
    exon_before_origin = exon_number_before_origin(transcript)
    statuses = []
    for run in runs:
        if run.exon_number < n_exons and run.exon_number != exon_before_origin:
            strength, reason = donor_status(transcript, layout, run)
            if strength != "intact":
                statuses.append(SpliceSiteStatus(
                    SpliceSiteKey("donor", run.exon_number, run.occurrence),
                    run.key, strength, reason, _site_context(layout, run, "donor")))
        if run.exon_number > 1 and run.exon_number - 1 != exon_before_origin:
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


def cryptic_candidates_for_status(
        layout, runs, status, policy=DEFAULT_CRYPTIC_SPLICE_POLICY):
    """Enumerate distinct graph-valid cryptic boundaries for one disrupted site.

    Returns ``tuple[(boundary_offset, genomic_position, motif_score), ...]`` or
    ``()`` when flanking genomic sequence is unavailable or no boundary meets
    ``policy``.
    """
    run_by_key = {run.key: run for run in runs}
    run = run_by_key.get(status.run_key)
    if run is None:
        return ()
    index = runs.index(run)
    if status.key.side == "donor":
        if index + 1 >= len(runs):
            return ()
        valid_interval = (run.layout_start + 1, runs[index + 1].layout_start - 1)
        canonical = run.layout_end
    else:
        if index == 0:
            return ()
        valid_interval = (runs[index - 1].layout_end + 1, run.layout_end - 1)
        canonical = run.layout_start
    if valid_interval[0] > valid_interval[1]:
        return ()
    try:
        candidates = enumerate_cryptic_boundaries(
            layout, canonical, status.key.side, valid_interval, policy=policy
        )
    except (SequenceUnavailable, KeyError, ValueError):
        return ()
    resolved = []
    for boundary, score in candidates:
        contig, position, _ = layout.origin(
            boundary - 1 if status.key.side == "donor" else boundary
        )
        genomic_position = layout.genomic(contig, position)
        resolved.append((boundary, genomic_position, score))
    return tuple(resolved)


def splice_axes(
        transcript, statuses, run_keys=None, *,
        layout=None, runs=None, cryptic_policy=None):
    """Convert disrupted sites into ordinal, uncalibrated graph choices."""
    axes = []
    n_exons = len(transcript.exons)
    run_keys = tuple(run_keys or ())
    policy = cryptic_policy if cryptic_policy is not None else DEFAULT_CRYPTIC_SPLICE_POLICY
    if layout is not None and runs is None:
        runs = build_exon_runs(transcript, layout)
    for status in statuses:
        mechanisms = _MECHANISM_ORDER[(status.key.side, status.strength)]
        # Exons meeting at a circular origin are one exon on the circle and
        # are skipped together.
        skipped_exons = origin_partners(transcript, status.key.exon_number)
        options = []
        option_rank = 0
        for mechanism in mechanisms:
            if (mechanism == "exon_skip"
                    and (1 in skipped_exons or n_exons in skipped_exons)):
                continue
            base_evidence = {
                "strength": status.strength,
                "reason": status.reason,
            }
            if mechanism == "cryptic" and layout is not None and policy.max_sites_per_junction > 1:
                cryptic_sites = cryptic_candidates_for_status(
                    layout, runs, status, policy=policy
                )
                if cryptic_sites:
                    for site_rank, (boundary, genomic_pos, score) in enumerate(cryptic_sites):
                        options.append(SpliceOption(
                            mechanism="cryptic",
                            ordinal_rank=option_rank,
                            required_runs=(status.run_key,),
                            evidence={
                                **base_evidence,
                                "cryptic_boundary_offset": boundary,
                                "cryptic_genomic_position": genomic_pos,
                                "cryptic_motif_score": score,
                                "cryptic_rank": site_rank,
                                "discovery_status": "heuristic_dna_motif",
                                "cryptic_search_policy": policy.as_evidence(),
                            },
                        ))
                        option_rank += 1
                    continue
            kwargs = {
                "mechanism": mechanism,
                "ordinal_rank": option_rank,
                "evidence": (
                    {**base_evidence, "cryptic_search_policy": policy.as_evidence()}
                    if mechanism == "cryptic"
                    else base_evidence
                ),
            }
            if mechanism == "exon_skip":
                occurrence = status.run_key[1]
                kwargs["skipped_runs"] = tuple(
                    key for key in ((number, occurrence) for number in skipped_exons)
                    if key == status.run_key or key in run_keys)
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
            option_rank += 1
        axes.append(SpliceAxis(
            key=status.key,
            anchor_runs=(status.run_key,),
            options=tuple(options)))
    return tuple(axes)


def enumerate_layout_splice_plans(
        transcript, layout, max_plans=64, *, cryptic_policy=None):
    """Detect disrupted sites and enumerate only graph-valid plans."""
    runs = build_exon_runs(transcript, layout)
    statuses = disrupted_splice_sites(transcript, layout, runs)
    run_keys = tuple(run.key for run in runs)
    axes = splice_axes(
        transcript, statuses, run_keys,
        layout=layout, runs=runs, cryptic_policy=cryptic_policy,
    )
    if not axes:
        return enumerate_splice_plans((), run_keys, 1)
    return enumerate_splice_plans(
        axes, run_keys, max_plans=max_plans)


def _shared_choice_matches(baseline_option, mutant_option):
    if baseline_option.mechanism != mutant_option.mechanism:
        return False
    if baseline_option.mechanism == "cryptic":
        return (
            baseline_option.evidence.get("cryptic_genomic_position")
            == mutant_option.evidence.get("cryptic_genomic_position")
        )
    return True


def enumerate_paired_splice_plans(
        transcript, baseline_statuses, mutant_statuses,
        baseline_run_keys, mutant_run_keys, shared_keys, max_plans=64, *,
        baseline_layout=None, mutant_layout=None,
        baseline_runs=None, mutant_runs=None, cryptic_policy=None):
    """Preserve uncertain baseline choices and couple unchanged active sites.

    Each side is graph-valid independently. Shared sites use the same mechanism
    (and cryptic boundary genomic position when enumerated) when present in
    both plans. If a rewrite removes the site on one side, the other side
    retains all compatible choices instead of an arbitrary default. The limit
    applies to complete comparisons, never a truncated candidate set.
    """
    baseline_plans = enumerate_splice_plans(
        splice_axes(
            transcript, baseline_statuses, baseline_run_keys,
            layout=baseline_layout, runs=baseline_runs, cryptic_policy=cryptic_policy,
        ),
        baseline_run_keys, max_plans=max_plans)
    mutant_plans = enumerate_splice_plans(
        splice_axes(
            transcript, mutant_statuses, mutant_run_keys,
            layout=mutant_layout, runs=mutant_runs, cryptic_policy=cryptic_policy,
        ),
        mutant_run_keys, max_plans=max_plans)
    pairs = []
    for mutant in mutant_plans:
        mutant_choices = dict(mutant.choices)
        for baseline in baseline_plans:
            baseline_choices = dict(baseline.choices)
            active_shared = shared_keys & baseline_choices.keys() & mutant_choices.keys()
            if any(not _shared_choice_matches(baseline_choices[key], mutant_choices[key])
                   for key in active_shared):
                continue
            pairs.append((baseline, mutant))
            if len(pairs) > max_plans:
                raise HypothesisLimitError(max_plans)
    return tuple(pairs)
