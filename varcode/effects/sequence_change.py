"""Sequence-change flags for existing SV predictions, not a new SV model."""

import operator

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


def _residue(codon, table, initiator):
    # As for a complete ORF, any start codon at the initiator reads as Met.
    if initiator and codon in table.start_codons:
        return "M"
    return table.forward_table.get(codon, "*")


def _partial_sequence_changes(model, transcript, table):
    """Compare only the mapped, in-frame codons of a partial observation.

    ``cds_start`` / ``cds_end`` place the producer's reading frame in the
    observed cDNA; transcript segments give each base's reference offset.
    A codon is compared only when its three bases map contiguously onto one
    reference CDS codon. A difference there proves a local change. Equal
    codons cannot show that unobserved sequence is unchanged, and missing
    coverage is not a deletion or truncation, so this never returns False.
    """
    evidence = model.evidence or {}
    try:
        # Any integer type (e.g. numpy); not missing or fractional bounds.
        start = operator.index(evidence.get("cds_start"))
        end = operator.index(evidence.get("cds_end"))
    except TypeError:
        return None, None
    sequence = model.cdna_sequence
    segments = model.reference_segments or ()
    if (sequence is None or model.edits or not transcript.complete
            or not 0 <= start < end <= len(sequence)
            or sum(s.length for s in segments) != len(sequence)):
        # Unmapped edits, or segments that don't render this sequence,
        # leave observed bases without reference coordinates.
        return None, None
    reference_offsets = {}
    offset = 0
    for segment in segments:
        if segment.source == transcript and segment.strand == "+":
            for pos in range(max(start, offset), min(end, offset + segment.length)):
                reference_offsets[pos] = segment.start + pos - offset
        offset += segment.length
    reference_start = min(transcript.start_codon_spliced_offsets)
    reference_coding = transcript.coding_sequence.upper()
    coding_status = None
    for pos in range(start, end - 2, 3):
        codon = sequence[pos:pos + 3].upper()
        ref = reference_offsets.get(pos)
        k = None if ref is None else ref - reference_start
        if (k is not None and k % 3 == 0 and 0 <= k < len(reference_coding)
                and reference_offsets.get(pos + 1) == ref + 1
                and reference_offsets.get(pos + 2) == ref + 2):
            reference_codon = reference_coding[k:k + 3]
            if codon != reference_codon and not set(codon + reference_codon) - set("ACGT"):
                coding_status = True
                # Only an observed ORF that begins here can initiate here.
                observed = _residue(codon, table, k == 0 and pos == start)
                if observed != _residue(reference_codon, table, k == 0):
                    return True, True
        if codon in table.stop_codons:
            break  # Later bases are not translated in this frame.
    return coding_status, None


def structural_sequence_changes(effect):
    """Return local (coding, protein) changes as True / False / None.

    Compare against the transcript being annotated, including on a fusion's
    3' side. The 5' partner supplies the initiation site and translation table.
    These are sequence predictions, not evidence of expression or translation.
    """
    from .effect_classes import LargeDeletion

    transcript = effect.transcript
    model = effect.mutant_transcript
    initiator = getattr(effect, "five_prime_transcript", transcript)
    coding_status = protein_status = None
    if model is not None and (model.evidence or {}).get(
            "protein_completeness", "start_to_stop") != "start_to_stop":
        # A fragment is not a full-length protein; a shorter one isn't truncated.
        coding_status, protein_status = _partial_sequence_changes(
            model, transcript, codon_table_for_transcript(initiator))
    elif model is not None:
        protein = model.mutant_protein_sequence
        reference_protein = transcript.protein_sequence
        if protein is not None and reference_protein is not None:
            # An empty protein is a known absence, not missing data.
            if "X" not in protein and "X" not in reference_protein:
                protein_status = protein != reference_protein
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
