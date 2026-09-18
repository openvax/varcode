"""Sequence-change flags for existing SV predictions, not a new SV model."""

import operator

from .codon_tables import codon_table_for_transcript, translate_sequence
from .selenocysteine import reference_selenocysteine, segment_selenocysteine


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


def _orf_bounds(model, initiator):
    """A known ORF, or the retained reference start; never search for one.

    Return ``(start, end)`` with ``end`` None when unbounded, or None for
    unmapped assemblies, partial fragments, missing sequence, or overlaid
    edits whose coordinates haven't been mapped.
    """
    if model.cdna_sequence is None:
        return None
    evidence = model.evidence or {}
    start, end = evidence.get("cds_start"), evidence.get("cds_end")
    if start is None or end is None:
        if not initiator.complete or model.edits:
            return None
        start, end = _retained_start(model, initiator), None
    return None if start is None else (start, end)


def _read_orf(sequence, start, end, table, selenocysteine):
    """The ORF from ``start`` through its stop codon, or None if incomplete.

    TGA codons at the offsets in ``selenocysteine`` are read as Sec. A
    producer's ``end`` assumed they terminate, so reading through lifts it.
    """
    sequence = sequence[start:].upper()
    limit = len(sequence) if end is None else end - start
    if sequence[:3] not in table.start_codons:
        return None
    for i in range(3, len(sequence) - 2, 3):
        if start + i in selenocysteine:
            limit = len(sequence)
        elif i + 3 > limit:
            return None
        elif sequence[i:i + 3] in table.stop_codons:
            coding = sequence[:i + 3]
            return coding if not set(coding) - set("ACGT") else None
    return None


def _initiated(protein):
    # Any start codon initiates as Met; Ensembl writes CTG/TTG initiators literally.
    return "M" + protein[1:] if protein else protein


def _orf_changes(model, transcript, reference_protein, table, start, end, selenocysteine):
    """Compare the ORF read with one Sec decoding against the reference."""
    coding = _read_orf(model.cdna_sequence, start, end, table, selenocysteine)
    if coding is None:
        return None, None
    reference_coding = transcript.coding_sequence
    coding_status = None if reference_coding is None else coding != reference_coding.upper()
    if reference_protein is None or "X" in reference_protein:
        return coding_status, None
    residues = translate_sequence(coding, codon_table=table, to_stop=False,
                                  selenocysteine={pos - start for pos in selenocysteine})
    # Initiator methionine also applies to alternative start codons.
    protein = "M" + residues[1:-1]
    return coding_status, protein != reference_protein


def _complete_sequence_changes(model, transcript, initiator):
    """Compare a start-to-stop prediction with the whole reference.

    CDS bases are compared with every Sec codon read through; how Sec is
    decoded affects only the protein. Where that is unknown, evaluate both
    readings and keep the protein flag only when they agree.
    """
    supplied = None
    protein = _initiated(model.mutant_protein_sequence)
    reference_protein = _initiated(transcript.protein_sequence)
    if (protein is not None and reference_protein is not None
            and "X" not in protein and "X" not in reference_protein
            # Ending exactly at a Sec residue is a decoding choice, not a change.
            and not (reference_protein.startswith(protein)
                     and reference_protein[len(protein):len(protein) + 1] == "U")):
        # An empty protein is a known absence, not missing data.
        supplied = protein != reference_protein
    bounds = _orf_bounds(model, initiator)
    if bounds is None:
        return None, supplied
    start, end = bounds
    mapped = segment_selenocysteine(model)
    if mapped is not None and any(s.source == transcript for s in model.reference_segments or ()):
        decoded, uncertain, terminating = mapped
    else:
        # Without coordinates, a TGA at a reference Sec index of an ORF read
        # from this transcript's own start may be Sec.
        decoded, uncertain, terminating = set(), set(), set()
        if initiator == transcript and reference_protein:
            sequence = model.cdna_sequence
            uncertain = {pos for pos in (start + 3 * i for i, aa in enumerate(reference_protein)
                                         if aa == "U")
                         if sequence[pos:pos + 3].upper() == "TGA"}
    table = codon_table_for_transcript(initiator)
    coding_status, _ = _orf_changes(model, transcript, reference_protein, table, start, end,
                                    decoded | uncertain | terminating)
    proteins = {_orf_changes(model, transcript, reference_protein, table, start, end,
                             selenocysteine)[1]
                for selenocysteine in {frozenset(decoded), frozenset(decoded | uncertain)}}
    protein_status = proteins.pop() if len(proteins) == 1 else None
    if uncertain:
        # A supplied protein assumed one decoding of these codons.
        return coding_status, protein_status
    return coding_status, supplied if supplied is not None else protein_status


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
    reference_sec = reference_selenocysteine(transcript)
    decoded = segment_selenocysteine(model)[0]
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
                reference = ("U" if ref in reference_sec
                             else _residue(reference_codon, table, k == 0))
                if observed != reference:
                    return True, True
        if codon in table.stop_codons and pos not in decoded:
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
        coding_status, protein_status = _complete_sequence_changes(
            model, transcript, initiator)

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
                # Deleting a selenoprotein's 3' UTR may remove its SECIS element.
                stops = reference_selenocysteine(transcript) and transcript.stop_codon_positions
                if not stops or (end < min(stops) if transcript.strand == "+"
                                 else start > max(stops)):
                    protein_status = False
            elif (transcript.contains_start_codon
                  and any(start <= pos <= end for pos in transcript.start_codon_positions)):
                protein_status = True

    if protein_status is True and coding_status is None:
        coding_status = True
    return coding_status, protein_status
