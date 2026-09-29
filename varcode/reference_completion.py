"""Explicit reference-completion hypotheses for mapped RNA observations."""

from copy import deepcopy
from dataclasses import replace
from itertools import product

from .effect_candidates import EffectCandidate
from .effects.codon_tables import codon_table_for_transcript, translate_sequence
from .effects.effect_classes import GeneFusion, TranslocationToIntergenic
from .effects.selenocysteine import segment_selenocysteine
from .effects.sequence_change import _retained_start
from .mutant_transcript import ReferenceSegment
from .nucleotides import reverse_complement


_SOURCE = "varcode_reference_completion"
_ENDS = ("five_prime", "three_prime")


def _genomic_pieces(transcript, start, end):
    """Inclusive genomic pieces for a cDNA interval, in transcript order."""
    offset = 0
    for lo, hi in sorted(transcript.exon_intervals, reverse=transcript.strand == "-"):
        length = hi - lo + 1
        a, b = max(start, offset) - offset, min(end, offset + length) - offset
        if a < b:
            yield ((lo + a, lo + b - 1) if transcript.strand == "+"
                   else (hi - a, hi - b + 1))
        offset += length


def _remap(segments, source, target):
    """A compatible isoform keeps one cDNA shift across all mapped donor runs."""
    if any(getattr(source, key, None) != getattr(target, key, None)
           for key in ("gene_id", "contig", "strand", "genome")):
        return None
    sequence = getattr(target, "sequence", None)
    if not sequence:
        return None
    shift = None
    mapped = []
    for segment in segments:
        if segment.source != source:
            mapped.append(segment)
            continue
        if segment.strand != "+":
            return None
        if target == source:
            start = segment.start
        else:
            offsets = []
            try:
                for first, last in _genomic_pieces(source, segment.start, segment.end):
                    a, b = target.spliced_offset(first), target.spliced_offset(last)
                    if b - a != abs(last - first):
                        return None
                    if offsets and a != offsets[-1][1] + 1:
                        return None
                    offsets.append((a, b))
            except (AttributeError, ValueError):
                return None
            if not offsets or offsets[-1][1] - offsets[0][0] + 1 != segment.length:
                return None
            start = offsets[0][0]
        delta = start - segment.start
        if shift is not None and shift != delta:
            return None
        shift = delta
        end = start + segment.length
        if (start < 0 or end > len(sequence)
                or sequence[start:end].upper() != source.sequence[segment.start:segment.end].upper()
                or set(sequence[start:end].upper()) - set("ACGT")):
            return None
        mapped.append(replace(segment, source=target, start=start, end=end))
    return tuple(mapped)


def _end_options(model, source, choices, end, enabled):
    segments = model.reference_segments
    terminal = segments[0 if end == "five_prime" else -1]
    if not enabled or source is None or terminal.source != source or terminal.strand != "+":
        return [(source, segments, None)]
    options, seen = [], set()
    for transcript in (source,) if choices is None else choices:
        identifier = getattr(transcript, "id", None)
        if identifier is None:
            raise ValueError("Completion transcripts must have an id")
        if identifier in seen:
            continue
        seen.add(identifier)
        mapped = _remap(segments, source, transcript)
        if mapped is None:
            continue
        edge = mapped[0 if end == "five_prime" else -1]
        start, stop = ((0, edge.start) if end == "five_prime"
                       else (edge.end, len(transcript.sequence)))
        extension = (ReferenceSegment(transcript, start, stop,
                                      label="assumed_reference_" + end)
                     if start < stop else None)
        options.append((transcript, mapped, extension))
    return options


def _translate(model, transcript):
    """Use the retained annotated start, never an ORF search over added bases."""
    evidence = dict(model.evidence)
    evidence.update(protein_status="not_determined", protein_completeness="unknown")
    try:
        starts = sorted(getattr(transcript, "start_codon_spliced_offsets", ()))
    except ValueError:
        # PyEnsembl raises when start-codon annotation is absent or partial.
        starts = []
    if len(starts) != 3 or starts != list(range(starts[0], starts[0] + 3)):
        return replace(model, evidence=evidence)
    start = _retained_start(model, transcript)
    sequence = model.cdna_sequence
    table = codon_table_for_transcript(transcript)
    if start is None or sequence[start:start + 3] not in table.start_codons:
        return replace(model, evidence=evidence)
    decoded, uncertain, _ = segment_selenocysteine(model)
    sec = decoded | uncertain
    stop = next((i for i in range(start + 3, len(sequence) - 2, 3)
                 if sequence[i:i + 3] in table.stop_codons and i not in sec), None)
    end = len(sequence) if stop is None else stop + 3
    coding = sequence[start:end - (end - start) % 3]
    if set(coding) - set("ACGT"):
        evidence["protein_status"] = "unresolved_ambiguous_coding_sequence"
        return replace(model, evidence=evidence)
    protein = translate_sequence(coding, codon_table=table,
                                 selenocysteine={p - start for p in sec})
    protein = "M" + protein[1:]
    observed_start = evidence["observed_model_evidence"].get("cds_start")
    relation = "unknown"
    if isinstance(observed_start, int) and not isinstance(observed_start, bool):
        delta = observed_start + evidence["reference_completion"]["observed_span"][0] - start
        relation = "in_frame" if delta % 3 == 0 else "out_of_frame"
    evidence.update(
        cds_start=start, cds_end=end, initiation_status="assumed_reference_start",
        protein_status="predicted_from_reference_completion",
        protein_completeness="partial_end" if stop is None else "start_to_stop",
        protein_translation_table=table.ncbi_id,
        observed_orf_frame=relation,
        selenocysteine_decoded_offsets=sorted(decoded),
        selenocysteine_uncertain_offsets=sorted(uncertain))
    return replace(model, mutant_protein_sequence=protein, evidence=evidence)


def reference_completion_hypotheses(
        candidate, *, five_prime_transcripts=None, three_prime_transcripts=None,
        ends=("five_prime", "three_prime")):
    """Complete mapped RNA ends as separate, explicitly assumed hypotheses.

    Parameters
    ----------
    candidate : EffectCandidate
        Observed GeneFusion or TranslocationToIntergenic with materialized
        cDNA and reference segments. The candidate is never modified.
    five_prime_transcripts, three_prime_transcripts : iterable or None
        Reference isoforms to consider at each end. None uses that end's
        annotated transcript. Alternatives must agree with the observed mapped
        donor runs in genomic coordinates, sequence and splice path. An empty
        iterable yields no hypotheses when that end can be completed.
    ends : iterable of str
        Ends to complete: ``five_prime``, ``three_prime``, or both (default).

    Returns
    -------
    tuple of EffectCandidate
        New candidates with source ``varcode_reference_completion``. At least
        one reference base is added per hypothesis. Unmapped ends are left
        alone. Explicit ``rna_five_prime_complete`` / ``rna_three_prime_complete``
        evidence set to True prevents extension of that end. ORF completeness
        and read ends do not establish RNA ends.

        Evidence records half-open observed/assumed spans in completed cDNA,
        selected transcripts, and nested original evidence/protein. Read counts
        are not promoted to support for the inferred sequence. Translation uses
        the retained annotated 5' start and existing selenocysteine rules; it
        does not establish biological initiation, translation or expression.
        Return these alongside observations, never as RNA-observed outcomes.

    Raises
    ------
    ValueError
        For unsupported effects, malformed layouts, invalid end evidence, or
        an already reference-completed input.
    """
    if isinstance(ends, str):
        raise ValueError("ends must contain five_prime and/or three_prime")
    ends = set(ends)
    if ends - set(_ENDS):
        raise ValueError("ends must contain five_prime and/or three_prime")
    effect = candidate.effect
    if not isinstance(effect, (GeneFusion, TranslocationToIntergenic)):
        raise ValueError("Reference completion requires an observed fusion outcome")
    model = effect.mutant_transcript
    if (candidate.source == _SOURCE or (model and (model.evidence or {}).get(
            "sequence_status") == "reference_completion_hypothesis")):
        raise ValueError("Cannot complete an already reference-completed hypothesis")
    if (model is None or not model.cdna_sequence or not model.reference_segments
            or model.edits):
        raise ValueError("Reference completion requires a materialized segment layout without edits")
    rendered = []
    for segment in model.reference_segments:
        sequence = getattr(segment.source, "sequence", None)
        if sequence is None or segment.length <= 0 or segment.end > len(sequence):
            raise ValueError("Invalid observed reference segment")
        piece = sequence[segment.start:segment.end]
        rendered.append(piece if segment.strand == "+" else reverse_complement(piece))
    if "".join(rendered).upper() != model.cdna_sequence.upper():
        raise ValueError("Reference segments do not render the observed cDNA")
    enabled = {}
    for end in _ENDS:
        flags = [(evidence or {}).get("rna_" + end + "_complete")
                 for evidence in (candidate.evidence, model.evidence)]
        if any(flag is not None and not isinstance(flag, bool) for flag in flags):
            raise ValueError("RNA end completeness must be bool or None")
        enabled[end] = end in ends and True not in flags
    first = getattr(effect, "five_prime_transcript", effect.transcript)
    last = getattr(effect, "three_prime_transcript", None)
    left = _end_options(model, first, five_prime_transcripts, "five_prime", enabled["five_prime"])
    right = _end_options(model, last, three_prime_transcripts, "three_prime", enabled["three_prime"])
    hypotheses = []
    for (five, left_segments, prefix), (three, right_segments, suffix) in product(left, right):
        if prefix is None and suffix is None:
            continue
        # Both remappings apply to different genes, preserving all opaque bases.
        observed = tuple(r if s.source == last else l for s, l, r in zip(
            model.reference_segments, left_segments, right_segments))
        segments = ((prefix,) if prefix else ()) + observed + ((suffix,) if suffix else ())
        offset = prefix.length if prefix else 0
        observed_end = offset + len(model.cdna_sequence)
        assumed = []
        for segment, start, end in ((prefix, 0, offset),
                                    (suffix, observed_end, observed_end + (suffix.length if suffix else 0))):
            if segment:
                assumed.append(dict(start=start, end=end, transcript_id=segment.source.id,
                                    reference_start=segment.start, reference_end=segment.end,
                                    label=segment.label))
        evidence = dict(
            sequence_status="reference_completion_hypothesis",
            observed_source=candidate.source,
            observed_evidence=deepcopy(dict(candidate.evidence)),
            observed_model_evidence=deepcopy(dict(model.evidence or {})),
            observed_protein_sequence=model.mutant_protein_sequence,
            reference_completion=dict(
                observed_span=[offset, observed_end], assumed_spans=assumed,
                five_prime_transcript_id=five.id,
                three_prime_transcript_id=three.id if three else None))
        sequence = ((prefix.source.sequence[prefix.start:prefix.end] if prefix else "")
                    + model.cdna_sequence
                    + (suffix.source.sequence[suffix.start:suffix.end] if suffix else ""))
        completed = _translate(replace(
            model, reference_transcript=five, reference_segments=segments,
            cdna_sequence=sequence.upper(), mutant_protein_sequence=None,
            annotator_name=_SOURCE, evidence=evidence), five)
        if isinstance(effect, GeneFusion):
            result = GeneFusion(effect.variant, five, three, mutant_transcript=completed)
        else:
            result = TranslocationToIntergenic(effect.variant, five, mutant_transcript=completed)
        result.candidate_evidence = completed.evidence
        hypotheses.append(EffectCandidate(result, source=_SOURCE, evidence=completed.evidence))
    return tuple(hypotheses)
