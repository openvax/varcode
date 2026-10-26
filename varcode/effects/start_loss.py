"""Retain nucleotide context without inferring a new translation initiation site."""

from dataclasses import replace

from ..mutant_transcript import (
    MutantTranscript, _AssembledAllele, _mutant_cds_start, _resolve_variant_edit,
)


def start_loss_details(variant, transcript, model=None):
    """Map the reference start into an existing product or one validated point edit."""
    try:
        offsets = transcript.start_codon_spliced_offsets
        reference_start = min(offsets) if offsets else None
    except (AttributeError, ValueError):
        reference_start = None
    try:
        reference_sequence = transcript.sequence
    except (AttributeError, ValueError):
        reference_sequence = None
    reference_codon = None
    if reference_start is not None and reference_sequence is not None:
        reference_codon = str(reference_sequence)[reference_start:reference_start + 3].upper()
        if len(reference_codon) != 3:
            reference_codon = None
    result = dict(ref_codon=reference_codon, alt_codon=None,
                  ref_codon_offset=reference_start, alt_codon_offset=None,
                  codon_change_status="edit_unresolved", mutant_transcript=None)
    if ((model is None or not model.edits) and variant is not None
            and not getattr(variant, "is_structural", False)):
        if reference_sequence is not None:
            resolved = _resolve_variant_edit(variant, transcript, str(reference_sequence))
            if resolved is not None:
                edit, _ = resolved
                edited = (str(reference_sequence)[:edit.cdna_start] + edit.alt_bases
                          + str(reference_sequence)[edit.cdna_end:])
                if model is None:
                    model = MutantTranscript(
                        reference_transcript=transcript, edits=(edit,),
                        cdna_sequence=edited, annotator_name="varcode")
                elif model.cdna_sequence == edited:
                    # Only describe a realized product by this single edit when
                    # it accounts for the whole sequence (no other cis/splice edits).
                    model = replace(model, edits=(edit,), reference_segments=None)
    if model is None:
        return result
    if (not model.edits and model.reference_segments is None
            and model.cdna_sequence is not None
            and model.cdna_sequence != reference_sequence):
        model = replace(model, reference_segments=MutantTranscript.from_sequence(
            model.cdna_sequence).reference_segments)
    # Literal translation from the old frame is not a known mutant protein.
    model = replace(model, mutant_protein_sequence=None)
    result["mutant_transcript"] = model
    if model.cdna_sequence is None:
        result["codon_change_status"] = "sequence_unavailable"
        return result
    evidence = model.evidence or {}
    if "start_codon_offset" in evidence:
        # The realized path has already mapped (or lost) the original start.
        start = evidence["start_codon_offset"]
    elif reference_start is None:
        return result
    elif model.reference_segments is not None:
        if any(isinstance(s.source, _AssembledAllele) for s in model.reference_segments):
            return result  # External sequence alone does not map the original start.
        starts, offset = [], 0
        for segment in model.reference_segments:
            if (segment.source == transcript and segment.strand == "+"
                    and segment.start <= reference_start < segment.end):
                starts.append(offset + reference_start - segment.start)
            offset += segment.length
        if len(starts) > 1:
            result["codon_change_status"] = "start_ambiguous"
            return result
        start = starts[0] if starts else None
    elif model.edits:
        removed = any(e.cdna_start <= reference_start < e.cdna_end and e.length_delta
                      for e in model.edits)
        start = None if removed else _mutant_cds_start(transcript, model.edits)
    else:
        return result  # An unmapped assembly cannot acquire a guessed start.
    if start is None:
        result["codon_change_status"] = "start_removed"
        return result
    codon = model.cdna_sequence[start:start + 3].upper()
    if len(codon) != 3:
        result["codon_change_status"] = "sequence_unavailable"
        return result
    result.update(alt_codon=codon, alt_codon_offset=start, codon_change_status="mapped")
    return result
