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

"""Structural-variant classification shared by ``fast`` and ``transcript_model``.

Mutant-transcript construction lives in :mod:`varcode.effects.structural_models`,
partner discovery in :mod:`varcode.effects.structural_fusions`, and splice/cryptic
candidate attachment in :mod:`varcode.effects.structural_candidates`.
"""

from dataclasses import replace

from ..circular_origin import TranscriptFrame
from ..mutant_transcript import AssembledAllele, TranscriptEdit, retained_cds_start
from ..structural_variant import affected_span
from ..sv_allele_parser import breakend_sides
from .classify import classify_from_protein_diff
from .codon_tables import codon_table_for_transcript, translate_sequence
from .effect_classes import (
    ExonLoss,
    FivePrimeUTR,
    Intergenic,
    Intronic,
    NoncodingTranscript,
    StartLoss,
    StructuralVariantEffect,
    ThreePrimeUTR,
    TranslocationToIntergenic,
    TranslocationToUnannotatedContig,
    Unresolved,
)
from . import structural_candidates
from .selenocysteine import segment_selenocysteine
from .sequence_change import structural_sequence_changes
from .structural_candidates import (
    breakpoint_splice_window,
    enumerate_and_attach_cryptics,
    enumerate_and_attach_splice_outcomes,
)
from .structural_fusions import (
    MAX_UPSTREAM_PARTNER_DISTANCE,
    coding_transcripts_at,
    fusion_across_junction,
    fusion_partners,
    junction_ends,
    unannotated_partner_contig,
    upstream_fusions,
    upstream_partners,
    warn_on_unknown_breakend_orientation,
)
from .structural_models import (
    LocalEvent,
    build_alt_assembly_mutant_transcript,
    build_deletion_mutant_transcript,
    build_duplication_mutant_transcript,
    build_fusion_mutant_transcript,
    build_inversion_mutant_transcript,
    build_translocation_mutant_transcript,
    cdna_cut,
    cdna_ranges_within_sv,
    contained_image,
    deleted_cdna_ranges,
    event_evidence,
    feature_contains_breakend,
    local_event,
    retains_five_prime_end,
    translate_fused_cdna,
)
from .transcript_helpers import (
    has_complete_coding_sequence,
    poly_a_completion,
    three_prime_utr_offset,
)

# Existing assembled-SV pickles use this legacy module path.
_AssembledAllele = AssembledAllele


def deletes_whole_transcript(variant, transcript):
    """Return whether ``variant`` is an unassembled DEL/CN0 removing all of ``transcript``."""
    if (variant.sv_type not in ("DEL", "CN0") or variant.alt_assembly
            or str(variant.contig) != str(transcript.contig)):
        return False
    frame = TranscriptFrame(transcript)
    covered = sum(
        max(0, min(high, frame.end) - max(low, frame.start) + 1)
        for low, high in frame.images(*affected_span(variant)))
    return covered >= frame.end - frame.start + 1


def overlapping_exons(variant, transcript):
    """Return exons of ``transcript`` (in transcript order) overlapping ``variant``."""
    if str(variant.contig) != str(transcript.contig):
        return []
    frame = TranscriptFrame(transcript)
    images = frame.images(*affected_span(variant))
    return [
        exon for exon, (start, end) in zip(frame.exons, frame.exon_intervals)
        if any(start <= high and low <= end for low, high in images)
    ]


def intronic_or_intergenic(event, transcript):
    """Classify an SV overlapping zero exons as :class:`Intronic` or :class:`Intergenic`."""
    image = (contained_image(event, transcript)
             if str(event.contig) == str(transcript.contig) else None)
    if image is not None:
        var_start = image[0]
        frame = TranscriptFrame(transcript)
        distance, nearest_exon = min(
            (min(abs(start - var_start), abs(end - var_start)), index)
            for index, (start, end) in enumerate(frame.exon_intervals))
        nearest_exon = frame.exons[nearest_exon]
        return Intronic(
            variant=event.variant,
            transcript=transcript,
            nearest_exon=nearest_exon,
            distance_to_exon=distance,
        )
    return Intergenic(variant=event.variant)


def classify_structural_transcript(event, transcript, model, affected):
    """Translate a mapped local DEL/DUP model and classify its protein consequence."""
    variant = event.variant

    def unresolved(reason):
        return Unresolved(variant, transcript, "structural_transcript", reason), model

    if deletes_whole_transcript(event, transcript):
        return ExonLoss(variant, transcript, exons=tuple(transcript.exons)), model
    if model is None or model.cdna_sequence is None:
        return unresolved("No complete spliced sequence is available")
    if not has_complete_coding_sequence(transcript):
        if event.sv_type in ("DEL", "CN0") and not event.alt_assembly:
            return ExonLoss(variant, transcript, exons=tuple(affected)), model
        return unresolved("Reference transcript has no complete annotated CDS")
    if event.alt_assembly or model.reference_segments is None or model.edits:
        return unresolved("The assembled sequence has no mapped annotated CDS")

    cds_start = retained_cds_start(model, transcript)
    if cds_start is None:
        return StartLoss(variant, transcript, mutant_transcript=model), model
    reference_start = min(transcript.start_codon_spliced_offsets)
    sequence = model.cdna_sequence
    mapped_sec = segment_selenocysteine(model)
    if mapped_sec is None:
        return unresolved("Transcript segments do not map the spliced sequence")
    decoded, uncertain, _ = mapped_sec
    if uncertain:
        effect, _ = unresolved("Selenocysteine decoding depends on the retained SECIS")
        effect.modifies_coding_sequence, effect.modifies_protein_sequence = (
            structural_sequence_changes(StructuralVariantEffect(
                variant, transcript, mutant_transcript=model)))
        return effect, model
    coding = sequence[cds_start:]
    if set(coding.upper()) - set("ACGT"):
        return unresolved("The coding sequence contains ambiguous bases")
    table = codon_table_for_transcript(transcript)
    coding += poly_a_completion(transcript, coding)
    protein = translate_sequence(
        coding[:len(coding) // 3 * 3], codon_table=table, to_stop=True,
        selenocysteine={pos - cds_start for pos in decoded})
    protein = "M" + protein[1:] if protein else protein
    model = replace(
        model,
        mutant_protein_sequence=protein,
        evidence={**dict(model.evidence or {}), "cds_start": cds_start},
    )

    reference = str(transcript.sequence)
    inserted = "".join(
        s.source.sequence[s.start:s.end]
        for s in model.reference_segments
        if s.label == "junction_insertion"
    )
    if event.sv_type in ("DEL", "CN0"):
        ranges = deleted_cdna_ranges(event, transcript)
        edits = tuple(
            TranscriptEdit(start, end, inserted if index == len(ranges) - 1 else "", variant)
            for index, (start, end) in enumerate(ranges))
    elif event.sv_type == "DUP":
        inside = cdna_ranges_within_sv(event, transcript)
        start = inside[-1][1]
        edits = (TranscriptEdit(
            start, start, inserted + "".join(reference[s:e] for s, e in inside), variant),)
    else:
        return unresolved("No mapped coding edit is available")

    stop = three_prime_utr_offset(transcript)
    reference_protein = str(transcript.protein_sequence)
    reference_protein = "M" + reference_protein[1:]
    after_stop = all(edit.cdna_start >= stop for edit in edits)
    if all(edit.cdna_end <= reference_start for edit in edits):
        effect = FivePrimeUTR(variant, transcript)
    elif after_stop and protein == reference_protein:
        effect = ThreePrimeUTR(variant, transcript)
    else:
        diff_model = replace(model, reference_segments=None, edits=edits)
        effect = classify_from_protein_diff(
            variant, transcript, reference_protein, protein,
            length_delta=sum(
                edit.length_delta for edit in edits if edit.cdna_start < stop),
            mutant_transcript=diff_model,
        )
        effect._original_protein_sequence = reference_protein
        if after_stop:
            effect._modifies_coding_sequence = False
    return effect, model


def local_consequence(event, transcript, affected, model):
    """Wrap a classified local SV transcript model with splice-ambiguity evidence."""
    if (event.sv_type == "INV" and model is not None and not event.alt_assembly
            and any(event.junction_inserted_sequence(junction) != ""
                    for junction in event.junctions)):
        model = replace(
            model, cdna_sequence=None, mutant_protein_sequence=None,
            evidence={**dict(model.evidence or {}),
                      "sequence_status": "unresolved_junction_insertion"})
    effect, model = classify_structural_transcript(event, transcript, model, affected)
    frame = TranscriptFrame(transcript)
    images = frame.images(*affected_span(event))
    ambiguous = any(
        start < low <= end or start <= high < end
        for start, end in frame.exon_intervals
        for low, high in images
        if start <= high and low <= end)
    evidence = {
        **event_evidence(event),
        "splice_ambiguous": ambiguous,
        "assumption": "reference_splicing",
    }
    if event.sv_type == "DUP":
        evidence["structure_assumption"] = "tandem_duplication"
    if model is not None:
        model = replace(model, evidence={**dict(model.evidence or {}), **evidence})
    effect.mutant_transcript = model
    effect.affected_exons = tuple(affected)
    result = StructuralVariantEffect(
        event.variant, transcript, primary_effects=(effect,), mutant_transcript=model)
    result.affected_exons = tuple(affected)
    result.candidate_evidence = evidence
    return result


def annotate_deletion(event, transcript, assembly=None):
    """Classify the surviving spliced transcript after a ``DEL`` event."""
    affected = overlapping_exons(event, transcript)
    if not affected and not any(
            event.junction_inserted_sequence(junction) != ""
            and any(feature_contains_breakend(exon, end)
                    for end in junction for exon in transcript.exons)
            for junction in event.junctions):
        return intronic_or_intergenic(event, transcript)
    return local_consequence(
        event, transcript, affected,
        assembly or build_deletion_mutant_transcript(event, transcript))


def annotate_duplication(event, transcript, assembly=None):
    """Classify a local ``DUP`` event on ``transcript``."""
    affected = overlapping_exons(event, transcript)
    if not affected:
        return intronic_or_intergenic(event, transcript)
    return local_consequence(
        event, transcript, affected,
        assembly or build_duplication_mutant_transcript(event, transcript))


def annotate_inversion(event, transcript, assembly=None):
    """Classify a local ``INV`` event on ``transcript``."""
    affected = overlapping_exons(event, transcript)
    if not affected:
        return intronic_or_intergenic(event, transcript)
    return local_consequence(
        event, transcript, affected,
        assembly or build_inversion_mutant_transcript(event, transcript))


def annotate_insertion(event, transcript, assembly=None):
    """Classify a local ``INS`` or ``CNV`` event on ``transcript``."""
    affected = overlapping_exons(event, transcript)
    if not affected:
        return intronic_or_intergenic(event, transcript)
    return local_consequence(event, transcript, affected, assembly)


def annotate_breakend(variant, transcript, assembly=None):
    """Classify a breakend whose junction leaves ``transcript``'s gene."""
    if variant.junctions:
        if assembly is None and breakend_sides(variant.symbolic_alt) is None:
            warn_on_unknown_breakend_orientation(variant)
        fusion = fusion_across_junction(variant, transcript, assembly)
        if fusion is not None:
            return fusion
    model = assembly or build_translocation_mutant_transcript(variant, transcript)
    partner_contig = unannotated_partner_contig(variant, transcript)
    if partner_contig is not None:
        return TranslocationToUnannotatedContig(
            variant=variant,
            transcript=transcript,
            partner_contig=partner_contig,
            mutant_transcript=model,
        )
    return TranslocationToIntergenic(
        variant=variant, transcript=transcript, mutant_transcript=model)


SPAN_EFFECT_BUILDERS = {
    "DEL": annotate_deletion,
    "DUP": annotate_duplication,
    "INV": annotate_inversion,
    "CNV": annotate_insertion,
    "INS": annotate_insertion,
}


def predict_structural_variant_effect(variant, transcript):
    """Classify structural ``variant`` on ``transcript``."""
    from pyensembl import Transcript
    from ..orf import ORFTranscriptView

    if not isinstance(transcript, (Transcript, ORFTranscriptView)):
        raise TypeError(
            "Expected a reference or ORF transcript, got %s" % type(transcript))

    sv_type = getattr(variant, "sv_type", None)
    if (not getattr(variant, "is_structural", False)
            or (sv_type != "BND" and sv_type not in SPAN_EFFECT_BUILDERS)):
        return NotImplemented

    if not transcript.is_protein_coding:
        return NoncodingTranscript(variant, transcript)

    assembly = build_alt_assembly_mutant_transcript(variant, transcript)
    event = local_event(variant, transcript)
    if event is None:
        effect = annotate_breakend(variant, transcript, assembly)
    else:
        span_effect = SPAN_EFFECT_BUILDERS[event.sv_type](event, transcript, assembly)
        effect = fusion_across_junction(variant, transcript, assembly)
        if effect is None:
            effect = span_effect
        elif isinstance(span_effect, StructuralVariantEffect):
            effect._attach_primary_effects((span_effect,))
    if isinstance(effect, StructuralVariantEffect):
        effect._attach_primary_effects(
            upstream_fusions(variant, transcript, assembly))

    structural_candidates.enumerate_and_attach_cryptics(variant, effect)
    structural_candidates.enumerate_and_attach_splice_outcomes(
        event or LocalEvent.of(variant), transcript, effect)
    if (type(effect) is StructuralVariantEffect
            and len(effect.candidates) == 1
            and not effect.candidate_evidence.get("splice_ambiguous")
            and not isinstance(effect.most_likely_effect, Unresolved)):
        return effect.most_likely_effect
    return effect


__all__ = [
    "LocalEvent",
    "MAX_UPSTREAM_PARTNER_DISTANCE",
    "SPAN_EFFECT_BUILDERS",
    "annotate_breakend",
    "annotate_deletion",
    "annotate_duplication",
    "annotate_insertion",
    "annotate_inversion",
    "breakpoint_splice_window",
    "build_alt_assembly_mutant_transcript",
    "build_deletion_mutant_transcript",
    "build_duplication_mutant_transcript",
    "build_fusion_mutant_transcript",
    "build_inversion_mutant_transcript",
    "build_translocation_mutant_transcript",
    "cdna_cut",
    "cdna_ranges_within_sv",
    "classify_structural_transcript",
    "coding_transcripts_at",
    "contained_image",
    "deleted_cdna_ranges",
    "deletes_whole_transcript",
    "enumerate_and_attach_cryptics",
    "enumerate_and_attach_splice_outcomes",
    "event_evidence",
    "feature_contains_breakend",
    "fusion_across_junction",
    "fusion_partners",
    "intronic_or_intergenic",
    "junction_ends",
    "local_consequence",
    "local_event",
    "overlapping_exons",
    "predict_structural_variant_effect",
    "retains_five_prime_end",
    "translate_fused_cdna",
    "unannotated_partner_contig",
    "upstream_fusions",
    "upstream_partners",
    "warn_on_unknown_breakend_orientation",
]
