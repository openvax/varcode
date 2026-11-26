"""Sequence-change flags for existing SV predictions, not a new SV model."""

import operator

from ..effect_traversal import get_effect_candidates
from ..mutant_transcript import retained_cds_start
from ..structural_variant import affected_span

from .codon_tables import codon_table_for_transcript, translate_sequence
from .selenocysteine import reference_selenocysteine, segment_selenocysteine
from .transcript_helpers import (
    has_complete_coding_sequence,
    poly_a_completion,
    polyadenylation_completes_stop,
)


def modification_status(effect, attribute):
    """Any changed candidate wins; otherwise unknown wins over unchanged.

    Walk by identity because SVs include themselves among their candidates.
    Read their local status directly, avoiding recursive property calls.
    """
    from .effect_classes import StructuralVariantEffect

    pending, seen = [effect], set()
    unknown = False
    evaluated = False
    while pending:
        current = pending.pop()
        if id(current) in seen:
            continue
        seen.add(id(current))
        candidates = get_effect_candidates(current)
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
        if not has_complete_coding_sequence(initiator) or model.edits:
            return None
        start, end = retained_cds_start(model, initiator), None
    return None if start is None else (start, end)


def _terminal_poly_a_completion(model, coding):
    """Poly(A) bases finishing ``coding`` when ``model`` reaches a poly(A)-stopped 3' end."""
    segments = model.reference_segments or ()
    if not segments:
        return ""
    last = segments[-1]
    source_seq = getattr(last.source, "sequence", None)
    if (last.strand == "+"
            and source_seq is not None
            and last.end == len(str(source_seq))
            and polyadenylation_completes_stop(last.source)):
        return poly_a_completion(last.source, coding)
    return ""


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


def _proper_reference_suffix(sequence, reference):
    """Could missing 5' coverage alone explain this shorter sequence?"""
    return (bool(sequence) and reference is not None
            and sequence != reference and reference.endswith(sequence))


def _orf_changes(model, transcript, reference_protein, table, start, end, selenocysteine,
                 partial_start=False, initiator=None):
    """Compare the ORF read with one Sec decoding against the reference."""
    tail = "" if end is not None else _terminal_poly_a_completion(model, model.cdna_sequence[start:])
    coding = _read_orf(model.cdna_sequence + tail, start, end, table, selenocysteine)
    if coding is None:
        return None, None
    reference_coding = transcript.coding_sequence
    if reference_coding is not None:
        reference_coding = reference_coding + poly_a_completion(transcript, reference_coding)
    coding_status = None if reference_coding is None else coding != reference_coding.upper()
    if partial_start and _proper_reference_suffix(
            coding, reference_coding.upper() if reference_coding is not None else None):
        coding_status = None
    if reference_protein is None or "X" in reference_protein:
        return coding_status, None
    residues = translate_sequence(coding, codon_table=table, to_stop=False,
                                  selenocysteine={pos - start for pos in selenocysteine})
    # Initiator methionine also applies to alternative start codons.
    protein = "M" + residues[1:-1]
    if partial_start and _proper_reference_suffix(protein, reference_protein):
        return _partial_sequence_changes(model, transcript, table, initiator=initiator)
    return coding_status, protein != reference_protein


def _complete_sequence_changes(model, transcript, initiator):
    """Compare a start-to-stop prediction with the whole reference.

    CDS bases are compared with every Sec codon read through; how Sec is
    decoded affects only the protein. Where that is unknown, evaluate both
    readings and keep the protein flag only when they agree.

    ORF boundaries do not establish transcript completeness. An unanchored
    observed reference suffix is compared as partial when its producer marks
    the transcript's completeness unknown.
    """
    supplied = None
    protein = _initiated(model.mutant_protein_sequence)
    bounds = _orf_bounds(model, initiator)
    if bounds is None and protein is None:
        return None, None
    reference_protein = _initiated(transcript.protein_sequence)
    partial_start = (model.evidence or {}).get(
        "sequence_status") == "observed_model_completeness_unknown" and not (
            bounds is not None and not model.edits and has_complete_coding_sequence(initiator)
            and sum(s.length for s in model.reference_segments or ()) == len(model.cdna_sequence)
            and bounds[0] == retained_cds_start(model, initiator))
    table = codon_table_for_transcript(initiator)
    if partial_start and _proper_reference_suffix(protein, reference_protein):
        # An observed internal Met-to-stop ORF may only be missing 5' coverage.
        # Keep its ORF boundaries; use mapped codons to establish any change.
        return _partial_sequence_changes(model, transcript, table, initiator=initiator)
    if (protein is not None and reference_protein is not None
            and "X" not in protein and "X" not in reference_protein
            # Ending exactly at a Sec residue is a decoding choice, not a change.
            and not (reference_protein.startswith(protein)
                     and reference_protein[len(protein):len(protein) + 1] == "U")):
        # An empty protein is a known absence, not missing data.
        supplied = protein != reference_protein
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
    coding_status, _ = _orf_changes(model, transcript, reference_protein, table, start, end,
                                    decoded | uncertain | terminating, partial_start,
                                    initiator=initiator)
    proteins = {_orf_changes(model, transcript, reference_protein, table, start, end,
                             selenocysteine, partial_start, initiator=initiator)[1]
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


def _in_frame_anchor(segments, target, start, end, translated_positions):
    """Return ``ref - pos`` for the first mapped in-frame CDS codon of ``target``."""
    if target is None or not has_complete_coding_sequence(target):
        return None
    reference_offsets = {}
    offset = 0
    for segment in segments:
        if segment.source == target and segment.strand == "+":
            for pos in range(max(start, offset), min(end, offset + segment.length)):
                reference_offsets[pos] = segment.start + pos - offset
        offset += segment.length
    if not reference_offsets:
        return None
    reference_start = min(target.start_codon_spliced_offsets)
    reference_coding_len = (
        len(target.coding_sequence)
        + len(poly_a_completion(target, target.coding_sequence))
    )
    for pos in translated_positions:
        ref = reference_offsets.get(pos)
        if (ref is not None and (ref - reference_start) % 3 == 0
                and reference_start <= ref < reference_start + reference_coding_len
                and reference_offsets.get(pos + 1) == ref + 1
                and reference_offsets.get(pos + 2) == ref + 2):
            return ref - pos
    return None


def _unanchored_fusion_changes(
    sequence,
    start,
    translated_positions,
    reference_coding,
    reference_start,
    reference_sec,
    table,
    transcript_table,
    decoded,
    uncertain,
):
    """Compare a 5'-partner-anchored partial fusion ORF against a 3' partner.

    Used when the reading frame is anchored on the 5' partner ``initiator``
    and the 3' partner ``transcript`` has no in-frame mapped codon (for
    example, an out-of-frame 3'-partner fusion or a 5'-partner-only read).
    Stops before the first ambiguous codon because it could hide a stop codon.
    """
    unambiguous = []
    obs_codons = []
    for pos in translated_positions:
        codon = sequence[pos:pos + 3].upper()
        if set(codon) - set("ACGT"):
            break
        unambiguous.append(pos)
        obs_codons.append(codon)
    n_obs = len(obs_codons)
    if n_obs == 0:
        return None, None
    n_ref = len(reference_coding) // 3
    ref_codons = [reference_coding[3 * j:3 * j + 3] for j in range(n_ref)]
    coding_matches = n_obs <= n_ref and any(
        all(
            bool(set(ref_codons[j + i]) - set("ACGT")) or obs_codons[i] == ref_codons[j + i]
            for i in range(n_obs)
        )
        for j in range(n_ref - n_obs + 1)
    )
    coding_status = None if coding_matches else True

    obs_residues = []
    for pos, codon in zip(unambiguous, obs_codons):
        if pos in uncertain:
            obs_residues.append({"U", "*"})
        elif pos in decoded:
            obs_residues.append({"U"})
        elif pos == start and codon in table.start_codons:
            obs_residues.append({"M", _residue(codon, table, False)})
        else:
            obs_residues.append({_residue(codon, table, False)})
    ref_residues = [
        None if set(ref_codons[j]) - set("ACGT")
        else ("U" if reference_start + 3 * j in reference_sec
              else _residue(ref_codons[j], transcript_table, j == 0))
        for j in range(n_ref)
    ]
    protein_matches = n_obs <= n_ref and any(
        all(
            ref_residues[j + i] is None or ref_residues[j + i] in obs_residues[i]
            for i in range(n_obs)
        )
        for j in range(n_ref - n_obs + 1)
    )
    protein_status = None if protein_matches else True
    return coding_status, protein_status


def _partial_sequence_changes(model, transcript, table, initiator=None):
    """Compare observed codons using mapped, in-frame reference anchors.

    ``cds_start`` / ``cds_end`` place the producer's reading frame in the
    observed cDNA; transcript segments give each base's reference offset.
    A mapped codon anchors comparison across an observed junction or indel:
    subsequent/preceding codons can be compared with the reference in that
    frame, even when their bases come from another source.

    Policy for upstream (``k < 0``) codons and 3'-partner fusions:

    - On the same transcript (``initiator`` is ``None`` or ``transcript``),
      codons upstream of the annotated start (``k < 0``) are ignored: an
      upstream start codon in a partial observation (e.g. ``partial_end``)
      does not by itself establish translation initiation there rather than at
      the annotated start.
    - On a fusion's 3' partner (``initiator != transcript``), an in-frame
      mapped CDS codon on ``initiator`` anchors the fusion's reading frame:
      unambiguous upstream codons reading continuously into ``transcript``'s
      CDS establish an N-terminal fusion change, and when ``transcript`` has
      no in-frame codon (e.g. an out-of-frame 3'-partner fusion), the
      observed codons and peptide are checked against all codon-aligned
      windows of ``transcript``'s reference CDS/protein.

    Equal codons cannot show that unobserved sequence is unchanged, so this
    never returns False or infers loss from missing ends.
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
    if (sequence is None or model.edits or not has_complete_coding_sequence(transcript)
            or not 0 <= start < end <= len(sequence)
            or sum(s.length for s in segments) != len(sequence)):
        # Unmapped edits, or segments that don't render this sequence,
        # leave observed bases without reference coordinates.
        return None, None
    reference_start = min(transcript.start_codon_spliced_offsets)
    reference_coding = (
        transcript.coding_sequence
        + poly_a_completion(transcript, transcript.coding_sequence)
    ).upper()
    reference_sec = reference_selenocysteine(transcript)
    decoded, uncertain, _ = segment_selenocysteine(model)
    translated_positions = []
    for pos in range(start, end - 2, 3):
        codon = sequence[pos:pos + 3].upper()
        translated_positions.append(pos)
        if codon in table.stop_codons and pos not in decoded:
            break
    # The first in-frame codon of this transcript anchors its continuation.
    # Later offsets may jump across a deletion/duplication; those jumps are
    # precisely why comparing each mapped codon only with itself is insufficient.
    anchor = _in_frame_anchor(segments, transcript, start, end, translated_positions)
    initiator_anchor = (
        _in_frame_anchor(segments, initiator, start, end, translated_positions)
        if initiator is not None and initiator != transcript
        else None
    )
    transcript_table = codon_table_for_transcript(transcript)
    if anchor is None:
        if initiator_anchor is None:
            return None, None
        return _unanchored_fusion_changes(
            sequence,
            start,
            translated_positions,
            reference_coding,
            reference_start,
            reference_sec,
            table,
            transcript_table,
            decoded,
            uncertain,
        )
    coding_status = None
    upstream_extension = False
    ambiguous_upstream = False
    for pos in translated_positions:
        codon = sequence[pos:pos + 3].upper()
        # A junction or frameshift can depart from the reference without
        # changing any of the individually mapped codons. The retained frame
        # still establishes a coordinate for its observed continuation.
        ref = pos + anchor
        k = ref - reference_start
        if k < 0:
            if initiator_anchor is not None:
                if set(codon) - set("ACGT"):
                    ambiguous_upstream = True
                else:
                    upstream_extension = True
        elif 0 <= k < len(reference_coding):
            reference_codon = reference_coding[k:k + 3]
            if codon != reference_codon and not set(codon + reference_codon) - set("ACGT"):
                coding_status = True
                # Only an observed ORF that begins here can initiate here.
                observed = ({"U", "*"} if pos in uncertain else
                            {"U"} if pos in decoded else
                            {_residue(codon, table, k == 0 and pos == start)})
                reference = ("U" if ref in reference_sec
                             else _residue(reference_codon, transcript_table, k == 0))
                if reference not in observed:
                    return True, True
    if upstream_extension and not ambiguous_upstream:
        return True, True
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
            model, transcript, codon_table_for_transcript(initiator),
            initiator=initiator)
    elif model is not None:
        coding_status, protein_status = _complete_sequence_changes(
            model, transcript, initiator)

    # A pure deletion can establish CDS/start loss without a translated model.
    # Do not apply DNA-span inference to a supplied assembly or an RNA outcome.
    if (isinstance(effect, LargeDeletion)
            and not getattr(effect.variant, "alt_assembly", None)
            and (model is None or model.annotator_name == "structural_variant")):
        start, end = affected_span(effect.variant)
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
