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

"""Fusion partner discovery and breakend mate classification helpers."""

from dataclasses import replace
import warnings

from ..circular_origin import transcript_exons, transcripts_at_locus
from .effect_classes import GeneFusion
from .effect_helpers import exon_length
from .structural_models import (
    build_fusion_mutant_transcript,
    cached_variant_lookup,
    feature_contains_breakend,
    retains_five_prime_end,
)

# How far past a junction end a sense-oriented coding transcript may start and
# still be a 3' fusion partner: LINX's limit for fusions other than known pairs
# and high-impact promiscuous partners (#550).
MAX_UPSTREAM_PARTNER_DISTANCE = 10_000


def warn_on_unknown_breakend_orientation(variant):
    """Warn when a breakend with a mate lacks readable kept-side orientation."""
    warnings.warn(
        "Structural prediction can't read the breakend orientation "
        "of %s from ALT %r; treating the annotated transcript as the "
        "5' fusion partner. Use a VCF breakend ALT such as N[17:100[ "
        "or pass StructuralVariant.alt_assembly." % (
            variant.short_description, variant.symbolic_alt),
        stacklevel=3,
    )


def junction_ends(variant, transcript):
    """Return ``(near, far)`` pairs where ``transcript`` spans ``near`` and its gene omits ``far``."""
    gene = transcript.gene
    return [
        (near, far)
        for junction in getattr(variant, "junctions", ())
        for near, far in (junction, junction[::-1])
        if feature_contains_breakend(transcript, near)
        and not feature_contains_breakend(gene, far)
    ]


def unannotated_partner_contig(variant, transcript):
    """Return the first unannotated mate contig leaving ``transcript``'s gene, or ``None`` (#709)."""
    annotated = variant.genome.contigs()
    for _, far in junction_ends(variant, transcript):
        if str(far.contig) not in annotated:
            return str(far.contig)
    return None


def coding_transcripts_at(variant, contig, position):
    """Return protein-coding transcripts overlapping ``contig:position`` (cached per variant)."""
    def lookup():
        try:
            transcripts = transcripts_at_locus(
                variant.genome, str(contig), int(position), int(position))
        except ValueError:
            return []
        return [t for t in transcripts if t.is_protein_coding]

    return cached_variant_lookup(
        variant, ("coding_transcripts", str(contig), int(position)), lookup)


def fusion_partners(variant, transcript, near, far, transcript_is_five_prime):
    """Yield compatible sense-to-sense protein-coding fusion partners at ``far``."""
    for candidate in coding_transcripts_at(variant, far.contig, far.position):
        if candidate.gene_id == transcript.gene_id:
            continue
        if feature_contains_breakend(candidate.gene, near):
            continue
        if far.keeps is None or (
                retains_five_prime_end(candidate, far.keeps) != transcript_is_five_prime):
            yield candidate


def fusion_across_junction(variant, transcript, assembly):
    """Build primary :class:`GeneFusion` candidates across compatible junctions and isoforms."""
    def keeps_five_prime(breakend):
        return breakend.keeps is None or retains_five_prime_end(transcript, breakend.keeps)

    fusions = []
    seen = set()
    for near, far in sorted(
            junction_ends(variant, transcript),
            key=lambda pair: not keeps_five_prime(pair[0])):
        five_prime = keeps_five_prime(near)
        for partner in fusion_partners(variant, transcript, near, far, five_prime):
            key = (near, far, partner.id)
            if key in seen:
                continue
            seen.add(key)
            if five_prime:
                five, five_position = transcript, near.position
                three, three_position = partner, far.position
            else:
                five, five_position = partner, far.position
                three, three_position = transcript, near.position
            junction = (near, far) if five_prime else (far, near)
            inserted = (variant.junction_inserted_sequence(junction)
                        if hasattr(variant, "junction_inserted_sequence") else None)
            fusions.append(GeneFusion(
                variant=variant,
                transcript=transcript,
                partner_transcript=partner,
                five_prime_transcript=five,
                three_prime_transcript=three,
                mutant_transcript=assembly or build_fusion_mutant_transcript(
                    transcript, five, five_position, three, three_position,
                    inserted_sequence=inserted),
            ))
    if not fusions:
        return None
    primary = fusions[0]
    primary._attach_primary_effects(fusions[1:])
    return primary


def upstream_partners(variant, far):
    """Return ``(transcript, distance)`` 3' partners starting within :data:`MAX_UPSTREAM_PARTNER_DISTANCE`."""
    def lookup():
        position = int(far.position)
        forward = far.keeps == "right"
        low, high = (
            (position + 1, position + MAX_UPSTREAM_PARTNER_DISTANCE)
            if forward else
            (max(1, position - MAX_UPSTREAM_PARTNER_DISTANCE), position - 1)
        )
        try:
            transcripts = transcripts_at_locus(
                variant.genome, str(far.contig), low, high)
        except ValueError:
            return []
        partners = []
        for transcript in transcripts:
            first = transcript_exons(transcript)[0]
            distance = first.start - position if forward else position - first.end
            if (transcript.is_protein_coding and len(transcript.exons) > 1
                    and transcript.strand == ("+" if forward else "-")
                    and 0 < distance <= MAX_UPSTREAM_PARTNER_DISTANCE):
                partners.append((transcript, distance))
        return partners

    return cached_variant_lookup(
        variant,
        ("upstream_partners", str(far.contig), int(far.position), far.keeps),
        lookup,
    )


def upstream_fusions(variant, transcript, assembly):
    """Return :class:`GeneFusion` candidates with 3' partners starting just past ``far`` (#550)."""
    fusions = []
    for near, far in junction_ends(variant, transcript):
        if (near.keeps is None or far.keeps is None
                or not retains_five_prime_end(transcript, near.keeps)):
            continue
        inserted = (variant.junction_inserted_sequence((near, far))
                    if hasattr(variant, "junction_inserted_sequence") else None)
        for partner, distance in upstream_partners(variant, far):
            if (partner.gene_id == transcript.gene_id
                    or feature_contains_breakend(partner.gene, near)):
                continue
            model = assembly or build_fusion_mutant_transcript(
                transcript, transcript, near.position, partner, far.position,
                inserted_sequence=inserted,
                three_prime_start=exon_length(transcript_exons(partner)[0]),
            )
            if model is not None and not assembly:
                model = replace(model, evidence={
                    **dict(model.evidence or {}),
                    "three_prime_breakend": "upstream",
                    "three_prime_upstream_distance": distance,
                })
            fusions.append(GeneFusion(
                variant=variant,
                transcript=transcript,
                partner_transcript=partner,
                five_prime_transcript=transcript,
                three_prime_transcript=partner,
                mutant_transcript=model,
            ))
    return fusions
