"""Sequence-change flags for existing SV predictions, not a new SV model."""

from .codon_tables import codon_table_for_transcript, translate_sequence


def modification_status(effect, attribute):
    """Any changed candidate wins; otherwise unknown wins over unchanged.

    Walk by identity because SVs include themselves among their candidates.
    Read their local status directly, avoiding recursive property calls.
    """
    from .effect_classes import MultiOutcomeEffect, StructuralVariantEffect

    pending, seen = [effect], set()
    unknown = False
    evaluated = False
    while pending:
        current = pending.pop()
        if id(current) in seen:
            continue
        seen.add(id(current))
        candidates = (current.candidates
                      if isinstance(current, MultiOutcomeEffect) else ())
        if isinstance(current, StructuralVariantEffect):
            if any(candidate.effect is current for candidate in candidates):
                evaluated = True
                coding, protein = structural_sequence_changes(current)
                status = protein if attribute == "modifies_protein_sequence" else coding
            else:
                # An explicit primary_effects tuple replaces the local result.
                status = False if candidates else None
        else:
            evaluated = True
            status = getattr(current, attribute)
        if status is True:
            return True
        unknown |= status is None
        pending.extend(candidate.effect for candidate in candidates)
    return None if unknown or not evaluated else False


def _retained_start(model, transcript):
    """Map all three annotated start bases, including split segment boundaries."""
    positions = [set(), set(), set()]
    start = min(transcript.start_codon_spliced_offsets)
    offset = 0
    for segment in model.reference_segments or ():
        if segment.source == transcript and segment.strand == "+":
            for i in range(3):
                if segment.start <= start + i < segment.end:
                    positions[i].add(offset + start + i - segment.start)
        offset += segment.length
    return next((p for p in sorted(positions[0])
                 if p + 1 in positions[1] and p + 2 in positions[2]), None)


def _coding_sequence(model, initiator):
    """Read a known ORF, or the retained reference ORF; never search for one.

    Return None for unmapped assemblies, partial fragments, missing sequence,
    overlaid edits whose coordinates haven't been mapped, or an incomplete ORF.
    """
    if model.cdna_sequence is None:
        return None
    evidence = model.evidence or {}
    start, end = evidence.get("cds_start"), evidence.get("cds_end")
    if start is None or end is None:
        if not initiator.complete or model.edits:
            return None
        start = _retained_start(model, initiator)
        end = None
    if start is None:
        return None
    sequence = model.cdna_sequence[start:end].upper()
    table = codon_table_for_transcript(initiator)
    if sequence[:3] not in table.start_codons:
        return None
    for i in range(3, len(sequence) - 2, 3):
        if sequence[i:i + 3] in table.stop_codons:
            coding = sequence[:i + 3]
            return coding if not set(coding) - set("ACGT") else None
    return None


def structural_sequence_changes(effect):
    """Return local (coding, protein) changes as True / False / None.

    Compare against the transcript being annotated, including on a fusion's
    3' side. The 5' partner supplies the initiation site and translation table.
    These are sequence predictions, not evidence of expression or translation.
    """
    from .effect_classes import LargeDeletion

    transcript = effect.transcript
    model = effect.mutant_transcript
    coding_status = protein_status = None
    if model is not None:
        protein = model.mutant_protein_sequence
        reference_protein = transcript.protein_sequence
        if protein is not None and reference_protein is not None:
            # An empty protein is a known absence, not missing data.
            if "X" not in protein and "X" not in reference_protein:
                protein_status = protein != reference_protein
        initiator = getattr(effect, "five_prime_transcript", transcript)
        coding = _coding_sequence(model, initiator)
        if coding is not None:
            reference_coding = transcript.coding_sequence
            if reference_coding is not None:
                coding_status = coding != reference_coding.upper()
            if protein_status is None and reference_protein is not None:
                # Initiator methionine also applies to alternative start codons.
                protein = "M" + translate_sequence(
                    coding[3:], codon_table=codon_table_for_transcript(initiator))
                if "X" not in reference_protein:
                    protein_status = protein != reference_protein

    # A pure deletion can establish CDS/start loss without a translated model.
    # Do not apply DNA-span inference to a supplied assembly or an RNA outcome.
    if (isinstance(effect, LargeDeletion)
            and not getattr(effect.variant, "alt_assembly", None)
            and (model is None or model.annotator_name == "structural_variant")):
        from .structural import _affected_span
        start, end = _affected_span(effect.variant)
        try:
            coding_ranges = transcript.coding_sequence_position_ranges
        except ValueError:
            # A protein-coding biotype need not have annotated CDS features.
            return coding_status, protein_status
        coding_overlap = any(start <= hi and lo <= end
                             for lo, hi in coding_ranges)
        if coding_status is None:
            coding_status = coding_overlap
        if protein_status is None:
            if not coding_overlap:
                protein_status = False
            elif (transcript.contains_start_codon
                  and any(start <= pos <= end for pos in transcript.start_codon_positions)):
                protein_status = True

    if protein_status is True and coding_status is None:
        coding_status = True
    return coding_status, protein_status
