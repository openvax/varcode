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


"""Intron-retention and cryptic splice-site candidate builders."""

from .effects.effect_classes import (
    ExonicSpliceSite,
    IntronRetention,
    IntronicSpliceSite,
    SpliceAcceptor,
    SpliceDonor,
)
from .effects.effect_helpers import exon_length
from .effects.transcript_helpers import has_complete_coding_sequence
from .mutant_transcript import (
    MutantTranscript,
    TranscriptEdit,
    resolve_variant_edit,
    translate_from_cds,
)
from .splice_exon_builders import (
    aa_fields_from_protein_diff,
    affected_exon,
    exon_start_offset_in_transcript,
    intronic_splice_side_is_acceptor,
    wrap_splice_candidate,
)


def _adjacent_intron_coords(transcript, exon, side):
    """Return ``(intron_start, intron_end)`` (1-based inclusive) of the
    intron on ``side`` of ``exon`` within ``transcript``, or ``None``
    if the exon has no intron on that side (e.g. first or last exon).

    ``side`` is ``"donor"`` (transcript-3' of exon → intron after it
    in transcript order) or ``"acceptor"`` (transcript-5' of exon).
    Handles forward + reverse strand transcripts.
    """
    exons = list(transcript.exons)
    try:
        idx = next(
            i for i, ex in enumerate(exons) if ex.exon_id == exon.exon_id)
    except StopIteration:
        return None
    reverse = transcript.on_backward_strand
    if side == "donor":
        adjacent_idx = idx + 1
    elif side == "acceptor":
        adjacent_idx = idx - 1
    else:
        return None
    if adjacent_idx < 0 or adjacent_idx >= len(exons):
        return None
    adjacent = exons[adjacent_idx]
    # Intron genomic coords are between the two exons on the reference.
    if reverse:
        # Transcript order: higher genomic → lower. Intron on donor
        # side sits below exon; acceptor side above.
        if side == "donor":
            intron_start = adjacent.end + 1
            intron_end = exon.start - 1
        else:
            intron_start = exon.end + 1
            intron_end = adjacent.start - 1
    else:
        if side == "donor":
            intron_start = exon.end + 1
            intron_end = adjacent.start - 1
        else:
            intron_start = adjacent.end + 1
            intron_end = exon.start - 1
    if intron_end < intron_start:
        return None
    return (intron_start, intron_end)


def _build_intron_retention_mutant_transcript(
        variant, transcript, exon, side, genomic_sequence):
    """Construct a :class:`MutantTranscript` for an intron-retention
    outcome using caller-supplied genomic sequence (#296).

    Returns ``None`` if the intron can't be resolved (no adjacent
    exon, provider raised, or transcript is incomplete).
    """
    if not has_complete_coding_sequence(transcript):
        return None
    coords = _adjacent_intron_coords(transcript, exon, side)
    if coords is None:
        return None
    intron_start, intron_end = coords
    try:
        intron_seq = genomic_sequence(
            transcript.contig, intron_start, intron_end)
    except Exception:
        return None
    if not intron_seq:
        return None
    intron_seq = intron_seq.upper()
    intron_seq = _apply_point_variant_to_genomic_sequence(
        intron_seq, intron_start, variant)
    if transcript.on_backward_strand:
        # Caller returns forward-strand genomic bases; for a reverse
        # strand transcript we need the reverse complement to get
        # transcript-order cDNA.
        from .nucleotides import reverse_complement
        intron_seq = reverse_complement(intron_seq)
    # Insertion point in the cDNA: immediately after the exon in
    # transcript order when the donor side is retained; immediately
    # before it when the acceptor side is retained.
    try:
        exon_start_offset = exon_start_offset_in_transcript(transcript, exon)
    except ValueError:
        return None
    length = exon_length(exon)
    if side == "donor":
        insert_at = exon_start_offset + length
    else:
        insert_at = exon_start_offset
    full_cdna = str(transcript.sequence)
    retention_edit = TranscriptEdit(
        cdna_start=insert_at,
        cdna_end=insert_at,
        alt_bases=intron_seq,
        source_variant=variant,
    )
    originating_edit = _point_variant_transcript_edit(
        variant, transcript, full_cdna)
    if originating_edit is False:
        return None
    edits = [retention_edit]
    if originating_edit is not None:
        edits.append(originating_edit)
    composed = _compose_transcript_edits(full_cdna, edits)
    if composed is None:
        return None
    mutant_cdna, edits = composed
    # Translate from the canonical start codon; premature stop inside
    # the retained intron terminates the ORF.
    mut_protein = translate_from_cds(mutant_cdna, transcript, edits)
    return MutantTranscript(
        reference_transcript=transcript,
        edits=edits,
        cdna_sequence=mutant_cdna,
        mutant_protein_sequence=mut_protein,
        annotator_name="splice_outcomes",
    )


def _splice_side_for_effect(splice_effect):
    """Return ``"donor"`` / ``"acceptor"`` for the disrupted side, or
    ``None`` if the splice class doesn't distinguish (falls through
    to a single-side default in the caller)."""
    if isinstance(splice_effect, SpliceDonor):
        return "donor"
    if isinstance(splice_effect, SpliceAcceptor):
        return "acceptor"
    if isinstance(splice_effect, IntronicSpliceSite):
        # The existing intronic_splice_side_is_acceptor helper keys
        # the side decision off the variant's position relative to
        # the exon.
        return "acceptor" if (
            intronic_splice_side_is_acceptor(splice_effect)) else "donor"
    if isinstance(splice_effect, ExonicSpliceSite):
        # ExonicSpliceSite sits at the end of an exon — default donor.
        return "donor"
    return None


def build_intron_retention_candidate(splice_effect, genomic_sequence=None):
    """Intron retention: the spliceosome fails and the intron is
    transcribed through to the mature mRNA.

    With a ``genomic_sequence`` provider, materializes the retained-
    intron :class:`MutantTranscript` and fills in the protein vocab
    (#296). Without a provider, the :class:`IntronRetention` is
    returned in its "predicted" state — ``retained_intron_start`` /
    ``retained_intron_end`` / ``side`` are still populated (from the
    adjacent-exon geometry) so consumers know *which* intron, but
    the protein fields are ``None``.
    """
    transcript = getattr(splice_effect, "transcript", None)
    exon = affected_exon(splice_effect)
    side = _splice_side_for_effect(splice_effect)

    # Compute the retained-intron coords from the adjacent-exon
    # geometry, independent of provider availability. When the
    # geometry isn't resolvable (e.g. first/last exon, transcript
    # missing the adjacent exon), pass None — the IntronRetention
    # instance carries None rather than a sentinel zero so consumers
    # can tell "unknown" from "actually coord 0".
    coords = None
    if transcript is not None and exon is not None and side is not None:
        coords = _adjacent_intron_coords(transcript, exon, side)
    intron_start = coords[0] if coords is not None else None
    intron_end = coords[1] if coords is not None else None

    common = dict(
        variant=splice_effect.variant,
        transcript=transcript,
        splice_signal=splice_effect,
        retained_intron_start=intron_start,
        retained_intron_end=intron_end,
        side=side,
    )

    if (genomic_sequence is not None and transcript is not None
            and exon is not None and side is not None):
        mt = _build_intron_retention_mutant_transcript(
            splice_effect.variant, transcript, exon, side, genomic_sequence)
        if mt is not None and mt.mutant_protein_sequence is not None:
            length_delta = len(mt.cdna_sequence) - len(
                str(transcript.sequence))
            aa_fields = aa_fields_from_protein_diff(
                splice_effect.variant, transcript, mt, length_delta)
            return wrap_splice_candidate(
                IntronRetention(mutant_transcript=mt, **common, **aa_fields),
            )

    # Predicted state: no protein math.
    return wrap_splice_candidate(
        IntronRetention(**common),
    )


_CRYPTIC_SCAN_FLANK = 50
"""Default number of bases to scan on each side of the canonical
splice boundary when looking for cryptic donor/acceptor sites.
Cryptic activations within ~50 bp are the usual clinical observations;
consumers wanting a wider scan pass ``cryptic_scan_flank`` through
:func:`enumerate_splice_outcomes`."""


def _apply_point_variant_to_genomic_sequence(
        sequence, sequence_start, variant, *, keep_origins=False):
    """Apply a nucleotide variant to forward-strand genomic sequence.

    ``sequence_start`` is the 1-based coordinate of ``sequence[0]``.
    Structural variants and variants outside the window are left alone.
    When the supplied sequence does not match the variant's reference
    allele (as can happen with synthetic test providers), it is also left
    alone rather than inventing a haplotype.

    If ``keep_origins`` is true, return ``(sequence, origins)`` where each
    origin is the corresponding 1-based reference position, or ``None``
    for an inserted base. This lets cryptic-site offsets be translated
    back to genomic coordinates after an indel changes sequence length.
    """
    origins = list(range(sequence_start, sequence_start + len(sequence)))
    if getattr(variant, "is_structural", False):
        return (sequence, origins) if keep_origins else sequence

    ref = getattr(variant, "trimmed_ref", "")
    alt = getattr(variant, "trimmed_alt", "")
    start = getattr(variant, "trimmed_base1_start", None)
    if start is None:
        return (sequence, origins) if keep_origins else sequence

    if ref:
        offset = start - sequence_start
        if offset < 0 or offset + len(ref) > len(sequence):
            return (sequence, origins) if keep_origins else sequence
        if sequence[offset:offset + len(ref)].upper() != ref.upper():
            return (sequence, origins) if keep_origins else sequence
        aligned = min(len(ref), len(alt))
        alt_origins = [start + i for i in range(aligned)]
        alt_origins.extend([None] * (len(alt) - aligned))
        sequence = sequence[:offset] + alt + sequence[offset + len(ref):]
        origins[offset:offset + len(ref)] = alt_origins
    elif alt:
        # Normalized insertions occur after ``start``.
        offset = start - sequence_start + 1
        if offset < 0 or offset > len(sequence):
            return (sequence, origins) if keep_origins else sequence
        sequence = sequence[:offset] + alt + sequence[offset:]
        origins[offset:offset] = [None] * len(alt)

    return (sequence, origins) if keep_origins else sequence


def _compose_transcript_edits(full_cdna, edits):
    """Apply non-conflicting reference-cDNA edits and return their result."""
    edits = sorted(edits, key=lambda edit: (edit.cdna_start, edit.cdna_end))
    for left, right in zip(edits, edits[1:]):
        if left.cdna_end > right.cdna_start:
            return None
    mutant_cdna = full_cdna
    for edit in reversed(edits):
        mutant_cdna = (
            mutant_cdna[:edit.cdna_start]
            + edit.alt_bases
            + mutant_cdna[edit.cdna_end:]
        )
    return mutant_cdna, tuple(edits)


def _point_variant_transcript_edit(variant, transcript, full_cdna):
    """Return the originating exonic edit.

    ``None`` means the allele is wholly intronic. ``False`` means it
    overlaps an exon but cannot be represented safely as one cDNA edit.
    """
    if getattr(variant, "is_structural", False):
        return None
    resolved = resolve_variant_edit(variant, transcript, full_cdna)
    if resolved is not None:
        return resolved[0]
    start = getattr(variant, "trimmed_base1_start", None)
    end = getattr(variant, "trimmed_base1_end", start)
    if start is not None and any(
            start <= exon.end and end >= exon.start
            for exon in transcript.exons):
        return False
    return None


def _without_edit_removed_by(edit, remove_start, remove_end):
    """Drop an originating edit swallowed by a splice-boundary truncation.

    Return ``False`` for a partial overlap, which cannot be composed safely
    in reference-cDNA coordinates.
    """
    if edit is None:
        return None
    if edit.cdna_start == edit.cdna_end:
        if remove_start < edit.cdna_start < remove_end:
            return None
        return edit
    if edit.cdna_end <= remove_start or edit.cdna_start >= remove_end:
        return edit
    if remove_start <= edit.cdna_start and edit.cdna_end <= remove_end:
        return None
    return False


def _canonical_boundary_offset(origins, exon, side, reverse):
    """Locate a reference splice boundary in a mutated scan layout.

    The returned offset is an interbase (Python-slice) offset in
    transcript order: donor exonic bases end there and acceptor exonic
    bases begin there.
    """
    if side == "donor":
        exonic_pos = exon.start if reverse else exon.end
        intronic_pos = exonic_pos - 1 if reverse else exonic_pos + 1
        if exonic_pos in origins:
            return origins.index(exonic_pos) + 1
        if intronic_pos in origins:
            return origins.index(intronic_pos)
    else:
        exonic_pos = exon.end if reverse else exon.start
        intronic_pos = exonic_pos + 1 if reverse else exonic_pos - 1
        if exonic_pos in origins:
            return origins.index(exonic_pos)
        if intronic_pos in origins:
            return origins.index(intronic_pos) + 1
    return None


def _cryptic_boundary_genomic_position(origins, offset, side):
    """Map an interbase cryptic boundary to its adjacent exonic base."""
    index = offset - 1 if side == "donor" else offset
    if index < 0 or index >= len(origins):
        return None
    return origins[index]


def _best_cryptic_site(sequence, canonical_offset, kind):
    """Scan ``sequence`` for the highest-scoring cryptic splice site
    of ``kind`` (``"donor"`` or ``"acceptor"``), ignoring the
    canonical position itself (#296).

    ``sequence`` is in transcript order (5'→3'). ``canonical_offset``
    and the returned candidate offset are interbase (Python-slice)
    offsets: donor exonic bases end at the offset and acceptor exonic
    bases begin there. Returns ``(offset, score)`` or ``None`` when no
    above-threshold site is found.
    """
    from .cryptic_exons import (
        ACCEPTOR_WINDOW,
        DONOR_WINDOW,
        score_acceptor,
        score_donor,
    )
    scorer = score_donor if kind == "donor" else score_acceptor
    window = DONOR_WINDOW if kind == "donor" else ACCEPTOR_WINDOW
    # Donor window is 3 exonic + 6 intronic, so the interbase junction
    # follows the third exonic base at window_start + 3.
    #
    # Acceptor window is 3 intronic + 1 exonic, so the interbase
    # junction before the first exonic base is also window_start + 3.
    junction_offset_in_window = 3
    best = None
    for i in range(len(sequence) - window + 1):
        sub = sequence[i:i + window]
        score = scorer(sub)
        if score <= 0:
            continue
        junction = i + junction_offset_in_window
        # Skip the canonical site itself — we're looking for alternatives.
        if junction == canonical_offset:
            continue
        if best is None or score > best[1]:
            best = (junction, score)
    return best


def _build_cryptic_site_mutant_transcript(
        variant, transcript, exon, side, genomic_sequence,
        scan_flank=_CRYPTIC_SCAN_FLANK):
    """Construct a :class:`MutantTranscript` where the canonical
    splice boundary on ``side`` (``"donor"`` or ``"acceptor"``) is
    replaced by the best-scoring cryptic site within ``scan_flank``
    bp on either side (#296).

    Returns ``(mt, score, cryptic_genomic_pos)`` or ``None`` when no
    cryptic site is found or the transcript is incomplete.
    """
    if not has_complete_coding_sequence(transcript):
        return None
    reverse = transcript.on_backward_strand
    # Determine the scan region around the canonical boundary.
    # Donor boundary = exon's transcript-3' end in genomic coords:
    #   forward: exon.end; reverse: exon.start.
    # Acceptor boundary = exon's transcript-5' start:
    #   forward: exon.start; reverse: exon.end.
    if side == "donor":
        boundary_pos = exon.start if reverse else exon.end
        # On forward strand, scan extends into the intron (higher coords)
        # by scan_flank and into the exon (lower coords) by scan_flank.
        # On reverse strand, intron is at LOWER coords.
        if reverse:
            scan_start = boundary_pos - scan_flank
            scan_end = boundary_pos + scan_flank
        else:
            scan_start = boundary_pos - scan_flank
            scan_end = boundary_pos + scan_flank
    else:  # acceptor
        boundary_pos = exon.end if reverse else exon.start
        scan_start = boundary_pos - scan_flank
        scan_end = boundary_pos + scan_flank
    # Fetch forward-strand sequence and realize the originating allele
    # before scoring. Otherwise a destroyed canonical motif remains in
    # the scan and can be selected as its own cryptic replacement (#422).
    try:
        forward_seq = genomic_sequence(
            transcript.contig, scan_start, scan_end).upper()
    except Exception:
        return None
    if not forward_seq:
        return None
    forward_seq, origins = _apply_point_variant_to_genomic_sequence(
        forward_seq, scan_start, variant, keep_origins=True)
    # For reverse-strand transcript, transcript-order sequence is the
    # reverse complement of the forward-strand region.
    from .nucleotides import reverse_complement
    if reverse:
        scan_seq = reverse_complement(forward_seq)
        origins.reverse()
    else:
        scan_seq = forward_seq
    canonical_offset = _canonical_boundary_offset(
        origins, exon, side, reverse)
    if canonical_offset is None:
        return None
    result = _best_cryptic_site(scan_seq, canonical_offset, kind=side)
    if result is None:
        return None
    cryptic_offset_in_seq, score = result
    cryptic_genomic_pos = _cryptic_boundary_genomic_position(
        origins, cryptic_offset_in_seq, side)
    # A boundary inside newly inserted sequence has no single genomic
    # base coordinate. The current effect model cannot represent it.
    if cryptic_genomic_pos is None:
        return None
    # Build the mutant cDNA by replacing the canonical exon boundary
    # with the cryptic one.
    #
    # Delta in exon length (transcript-order bases):
    #   donor:   cryptic - canonical (positive = extends exon)
    #   acceptor: canonical - cryptic (positive = extends exon at 5' end)
    transcript_direction = -1 if reverse else 1
    if side == "donor":
        exon_length_delta = (
            cryptic_genomic_pos - boundary_pos) * transcript_direction
    else:
        exon_length_delta = (
            boundary_pos - cryptic_genomic_pos) * transcript_direction
    try:
        exon_start_in_tx = exon_start_offset_in_transcript(transcript, exon)
    except ValueError:
        return None
    exon_len = exon_length(exon)
    full_cdna = str(transcript.sequence)
    if side == "donor":
        # Donor change shifts the exon's 3' end. Junction in transcript
        # cDNA is at exon_start_in_tx + exon_len.
        junction_cdna = exon_start_in_tx + exon_len
    else:
        # Acceptor change shifts the exon's 5' start.
        junction_cdna = exon_start_in_tx
    originating_edit = _point_variant_transcript_edit(
        variant, transcript, full_cdna)
    if originating_edit is False:
        return None
    if exon_length_delta >= 0:
        # Exon is extended → insert intron bases from the scan window.
        if side == "donor":
            # Added bases sit immediately after the canonical junction.
            # In transcript order they come from scan_seq between
            # canonical_offset and cryptic_offset_in_seq.
            added = scan_seq[canonical_offset:cryptic_offset_in_seq]
            insert_at = junction_cdna
        else:
            # Acceptor extension: intron bases go before the canonical
            # junction.
            added = scan_seq[cryptic_offset_in_seq:canonical_offset]
            insert_at = junction_cdna
        splice_edit = TranscriptEdit(
            cdna_start=insert_at,
            cdna_end=insert_at,
            alt_bases=added,
            source_variant=variant,
        )
    else:
        # Exon is truncated → remove bases from the exon's end (donor)
        # or start (acceptor).
        removed = -exon_length_delta
        if side == "donor":
            remove_start = junction_cdna - removed
            remove_end = junction_cdna
        else:
            remove_start = junction_cdna
            remove_end = junction_cdna + removed
        if originating_edit is not None:
            originating_edit = _without_edit_removed_by(
                originating_edit, remove_start, remove_end)
            if originating_edit is False:
                return None
        splice_edit = TranscriptEdit(
            cdna_start=remove_start,
            cdna_end=remove_end,
            alt_bases="",
            source_variant=variant,
        )
    edits = [splice_edit]
    if originating_edit is not None:
        edits.append(originating_edit)
    composed = _compose_transcript_edits(full_cdna, edits)
    if composed is None:
        return None
    mutant_cdna, edits = composed
    mut_protein = translate_from_cds(mutant_cdna, transcript, edits)
    mt = MutantTranscript(
        reference_transcript=transcript,
        edits=edits,
        cdna_sequence=mutant_cdna,
        mutant_protein_sequence=mut_protein,
        annotator_name="splice_outcomes",
    )
    return (mt, score, cryptic_genomic_pos)


def build_cryptic_splice_candidate(
        splice_effect, mechanism_cls, genomic_sequence=None,
        scan_flank=_CRYPTIC_SCAN_FLANK):
    """Cryptic donor/acceptor candidate. When ``genomic_sequence`` is
    provided (#296), scan the canonical splice boundary's flanking
    region on BOTH exon and intron sides for the best-scoring weak
    motif and materialize a real :class:`MutantTranscript` that uses
    the cryptic site in place of the disrupted canonical one.

    Without a provider (or when no motif scores above threshold),
    falls back to the predicted state — :class:`CrypticDonor` /
    :class:`CrypticAcceptor` with protein fields ``None`` and the
    cryptic-site fields (``cryptic_genomic_position``,
    ``motif_score``, ``exon_length_delta``) also ``None``.

    ``mechanism_cls`` is :class:`CrypticDonor` or
    :class:`CrypticAcceptor`; ``direction`` is derived from the
    class.
    """
    direction = mechanism_cls.direction
    transcript = getattr(splice_effect, "transcript", None)
    exon = affected_exon(splice_effect)

    common = dict(
        variant=splice_effect.variant,
        transcript=transcript,
        splice_signal=splice_effect,
        affected_exon=exon,
    )

    if (genomic_sequence is not None and transcript is not None
            and exon is not None):
        result = _build_cryptic_site_mutant_transcript(
            splice_effect.variant, transcript, exon, direction,
            genomic_sequence, scan_flank=scan_flank)
        if result is not None:
            mt, score, cryptic_pos = result
            length_delta = len(mt.cdna_sequence) - len(
                str(transcript.sequence))
            if mt.mutant_protein_sequence is not None:
                aa_fields = aa_fields_from_protein_diff(
                    splice_effect.variant, transcript, mt, length_delta)
                return wrap_splice_candidate(
                    mechanism_cls(
                        cryptic_genomic_position=cryptic_pos,
                        motif_score=score,
                        exon_length_delta=length_delta,
                        mutant_transcript=mt,
                        **common,
                        **aa_fields),
                )
            # MT built but protein didn't translate — partial info.
            return wrap_splice_candidate(
                mechanism_cls(
                    cryptic_genomic_position=cryptic_pos,
                    motif_score=score,
                    exon_length_delta=length_delta,
                    mutant_transcript=mt,
                    **common),
            )

    # Predicted state.
    return wrap_splice_candidate(
        mechanism_cls(**common),
    )




__all__ = [
    "build_cryptic_splice_candidate",
    "build_intron_retention_candidate",
]
