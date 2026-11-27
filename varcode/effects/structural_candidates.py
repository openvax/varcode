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

"""Breakpoint splice-window and cryptic-exon candidate attachment for SVs."""

from ..circular_origin import TranscriptFrame
from ..structural_variant import affected_span
from .structural_models import cached_variant_lookup, event_evidence


def breakpoint_splice_window(transcript, breakpoint_pos):
    """Return ``(splice_cls, nearest_exon, distance)`` when ``breakpoint_pos`` hits an intronic splice window (#341).

    Splice windows (intronic, matching the point-variant classifier):

    * <= 2 bp on the donor side -> :class:`SpliceDonor`
    * <= 2 bp on the acceptor side -> :class:`SpliceAcceptor`
    * 3-6 bp on the donor side -> :class:`IntronicSpliceSite`
    * 3 bp on the acceptor side -> :class:`IntronicSpliceSite`
    """
    from .effect_classes import IntronicSpliceSite, SpliceAcceptor, SpliceDonor

    best_exon = None
    best_distance = None
    best_before = True
    for exon in transcript.exons:
        if exon.start <= breakpoint_pos <= exon.end:
            return None
        dist_start = exon.start - breakpoint_pos
        dist_end = breakpoint_pos - exon.end
        distance = dist_start if dist_start > 0 else dist_end
        if best_distance is None or distance < best_distance:
            best_distance = distance
            best_exon = exon
            if exon.strand == "+":
                best_before = breakpoint_pos < exon.start
            else:
                best_before = breakpoint_pos > exon.end
    if best_exon is None or best_distance is None:
        return None
    if best_distance <= 2:
        cls = SpliceAcceptor if best_before else SpliceDonor
        return (cls, best_exon, best_distance)
    if not best_before and best_distance <= 6:
        return (IntronicSpliceSite, best_exon, best_distance)
    if best_before and best_distance <= 3:
        return (IntronicSpliceSite, best_exon, best_distance)
    return None


def cryptic_exon_candidate(cryptic, *, base_evidence=None):
    """Wrap a :class:`~varcode.effects.CrypticExonCandidate` as a ``varcode_motif`` candidate."""
    from ..effect_candidates import EffectCandidate

    donor = cryptic.donor_score or 0.0
    acceptor = cryptic.acceptor_score or 0.0
    evidence = {
        **(dict(base_evidence) if base_evidence else {}),
        "donor_score": cryptic.donor_score,
        "acceptor_score": cryptic.acceptor_score,
        "motif_score": (donor + acceptor) / 2.0,
        "interval_start": cryptic.interval_start,
        "interval_end": cryptic.interval_end,
    }
    return EffectCandidate(
        effect=cryptic,
        source="varcode_motif",
        evidence=evidence,
    )


def enumerate_and_attach_splice_outcomes(event, transcript, effect):
    """Synthesize and attach ``varcode_splice`` candidates when SV breakpoints hit splice windows (#341)."""
    from ..effect_candidates import EffectCandidate
    from ..splice_outcomes import enumerate_splice_outcomes
    from .effect_classes import NormalSplicing, StructuralVariantEffect

    if not isinstance(effect, StructuralVariantEffect):
        return
    affected_start, affected_end = affected_span(event)
    if (affected_start, affected_end) != (event.start, event.end):
        bp_positions = {position for _, position in event.breakpoints}
    else:
        bp_positions = {event.start, event.end}
    sv_evidence = event_evidence(event)
    attached = []
    seen_exons = set()
    genomic = TranscriptFrame(transcript).genomic
    for bp in bp_positions:
        window = breakpoint_splice_window(transcript, genomic(bp))
        if window is None:
            continue
        splice_cls, exon, distance = window
        if exon.exon_id in seen_exons:
            continue
        seen_exons.add(exon.exon_id)
        try:
            synthetic = splice_cls(
                variant=event.variant,
                transcript=transcript,
                nearest_exon=exon,
                distance_to_exon=distance,
            )
            splice_set = enumerate_splice_outcomes(synthetic)
        except (AttributeError, KeyError, ValueError):
            continue
        for candidate in splice_set.candidates:
            if isinstance(candidate.effect, NormalSplicing):
                continue
            attached.append(EffectCandidate(
                effect=candidate.effect,
                source="varcode_splice",
                evidence={**dict(candidate.evidence), **sv_evidence},
            ))
    if attached:
        effect.add_candidates(attached)


def enumerate_and_attach_cryptics(variant, effect):
    """Enumerate cryptic-exon candidates around ``variant`` and attach them to ``effect`` (#337)."""
    from ..cryptic_exons import enumerate_from_structural_variant
    from .effect_classes import StructuralVariantEffect

    if not isinstance(effect, StructuralVariantEffect):
        return

    def enumerate_candidates():
        try:
            return enumerate_from_structural_variant(variant)
        except (AttributeError, KeyError, ValueError, OSError):
            return []

    candidates = cached_variant_lookup(
        variant, "cryptic_candidates", enumerate_candidates)
    if candidates:
        base_evidence = dict(effect.candidate_evidence or {})
        sv_type = getattr(variant, "sv_type", None)
        if sv_type is not None:
            base_evidence["sv_type"] = sv_type
        effect.add_candidates(
            cryptic_exon_candidate(c, base_evidence=base_evidence)
            for c in candidates
        )
