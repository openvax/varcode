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

"""Internal structural-effect prediction shared by the default and transcript model.

The public default calls these helpers directly; this module does not select or
instantiate annotators. Structural classification, sequence construction and
candidate enumeration live here.

Mutant transcripts retain the historical ``structural_variant`` builder
provenance; their containing effect collection records the selected annotator.
"""

from .effect_classes import (
    GeneFusion,
    Intergenic,
    Intronic,
    Inversion,
    LargeDeletion,
    LargeDuplication,
    NoncodingTranscript,
    StructuralVariantEffect,
    TranslocationToIntergenic,
)
from .effect_helpers import exon_length
from ..mutant_transcript import MutantTranscript, ReferenceSegment
# Existing assembled-SV pickles use this private module path.
from ..mutant_transcript import _AssembledAllele  # noqa: F401
from ..sv_allele_parser import breakend_sides, _breakend_local_side


# --------------------------------------------------------------------
# MutantTranscript builders (#335).
#
# Map each SV classification to a :class:`MutantTranscript` populated
# with :class:`ReferenceSegment` entries that describe the
# rearranged allele. ``cdna_sequence`` is filled in where it's
# derivable from pyensembl-cached transcript cDNA alone (DEL, and
# DUP / INV where only whole exons are rearranged); otherwise left
# None and resolved by downstream helpers (#336, #338).
# --------------------------------------------------------------------


def _exon_cdna_ranges(transcript):
    """Yield ``(exon, cdna_start, cdna_end)`` for each exon of
    ``transcript`` in transcript order. ``cdna_start`` / ``cdna_end``
    are offsets into ``transcript.sequence``.
    """
    offset = 0
    for exon in transcript.exons:
        length = exon_length(exon)
        yield exon, offset, offset + length
        offset += length


def _merge_adjacent_ranges(ranges):
    """Merge ``[(s, e), ...]`` pairs that abut (``e_i == s_{i+1}``)
    so the resulting segment list doesn't contain redundant splits.
    Assumes the input is already sorted by start.
    """
    merged = []
    for s, e in ranges:
        if s >= e:
            continue
        if merged and merged[-1][1] == s:
            merged[-1] = (merged[-1][0], e)
        else:
            merged.append((s, e))
    return merged


def _affected_span(variant):
    """Inclusive reference bases changed by a structural variant."""
    start = getattr(variant, "affected_start", variant.start)
    end = getattr(variant, "affected_end", getattr(
        variant, "end", variant.start))
    return min(start, end), max(start, end)


def _cdna_ranges_kept_after_deletion(variant, transcript):
    """Compute the cDNA ranges of ``transcript`` that survive after
    the genomic deletion described by ``variant``.

    Returns a merged list of ``(cdna_start, cdna_end)`` tuples in
    transcript order. Handles both forward and reverse strand
    transcripts and the case where the deletion cuts through the
    middle of an exon.
    """
    del_start, del_end = _affected_span(variant)
    kept = []
    reverse = transcript.on_backward_strand
    for exon, c_s, c_e in _exon_cdna_ranges(transcript):
        ex_s, ex_e = exon.start, exon.end
        # Exon entirely outside the deletion.
        if ex_e < del_start or ex_s > del_end:
            kept.append((c_s, c_e))
            continue
        # Exon fully covered by the deletion.
        if del_start <= ex_s and del_end >= ex_e:
            continue
        # Partial overlap: preserve the genomic bases outside the
        # deletion and map them to cDNA. The direction of the map
        # depends on strand — on the forward strand cDNA position
        # 0 of the exon corresponds to ``ex_s``; on the reverse
        # strand it corresponds to ``ex_e``.
        if reverse:
            if del_end < ex_e:
                kept.append((c_s, c_s + (ex_e - del_end)))
            if del_start > ex_s:
                kept.append((c_e - (del_start - ex_s), c_e))
        else:
            if del_start > ex_s:
                kept.append((c_s, c_s + (del_start - ex_s)))
            if del_end < ex_e:
                kept.append((c_s + (del_end + 1 - ex_s), c_e))
    return _merge_adjacent_ranges(sorted(kept))


def _cdna_ranges_within_sv(variant, transcript):
    """Body ranges for a local DUP/INV model, or [] if not local.

    Clipping a cross-boundary event to this transcript invents a local
    rearrangement (#405). Every junction end must be in this transcript,
    not merely in its gene. Fusion and supplied-assembly paths are separate.
    """
    if any(not _contains(transcript, end)
           for junction in variant.junctions for end in junction):
        return []
    sv_start, sv_end = _affected_span(variant)
    inside = []
    reverse = transcript.on_backward_strand
    for exon, c_s, c_e in _exon_cdna_ranges(transcript):
        ex_s, ex_e = exon.start, exon.end
        if ex_e < sv_start or ex_s > sv_end:
            continue
        # Clip exon to [sv_start, sv_end].
        clipped_s = max(ex_s, sv_start)
        clipped_e = min(ex_e, sv_end)
        if reverse:
            # cDNA offsets within exon grow as genomic position falls.
            cdna_lo = c_s + (ex_e - clipped_e)
            cdna_hi = c_s + (ex_e - clipped_s + 1)
        else:
            cdna_lo = c_s + (clipped_s - ex_s)
            cdna_hi = c_s + (clipped_e - ex_s + 1)
        inside.append((cdna_lo, cdna_hi))
    return _merge_adjacent_ranges(sorted(inside))


def _breakpoint_splice_window(transcript, breakpoint_pos):
    """Return ``(splice_cls, nearest_exon, distance)`` if
    ``breakpoint_pos`` falls within a canonical splice window on
    ``transcript``; otherwise ``None`` (#341).

    Splice windows (intronic, matching the point-variant classifier):

    * ≤ 2 bp from exon boundary on the donor side (transcript-3' of
      exon) → :class:`SpliceDonor`
    * ≤ 2 bp from exon boundary on the acceptor side
      (transcript-5' of exon) → :class:`SpliceAcceptor`
    * 3-6 bp on the donor side → :class:`IntronicSpliceSite`
    * 3 bp on the acceptor side → :class:`IntronicSpliceSite`

    Returns distances in genomic-coord bases so downstream code can
    construct the synthesized effect with matching ``distance_to_exon``.
    Exonic-splice-site synthesis is deliberately out of scope here —
    SV breakpoints rarely land strictly inside exons, and the
    exonic-splice case has extra alternate-effect machinery that
    doesn't apply to an SV's already-disruptive primary outcome.
    """
    from .effect_classes import (
        IntronicSpliceSite,
        SpliceAcceptor,
        SpliceDonor,
    )
    best_exon = None
    best_distance = None
    best_before = True
    for exon in transcript.exons:
        # Strictly intronic scope only — skip the exon itself.
        if exon.start <= breakpoint_pos <= exon.end:
            return None
        dist_start = exon.start - breakpoint_pos  # >0 if breakpoint before exon
        dist_end = breakpoint_pos - exon.end      # >0 if breakpoint after exon
        distance = dist_start if dist_start > 0 else dist_end
        if best_distance is None or distance < best_distance:
            best_distance = distance
            best_exon = exon
            # "before_exon" in transcript direction:
            #   forward strand: breakpoint_pos < exon.start
            #   reverse strand: breakpoint_pos > exon.end
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


def _retains_five_prime_end(transcript, side):
    """Whether keeping ``side`` of a breakpoint inside ``transcript``
    keeps the transcript's 5' end: the left side on the forward
    strand, the right side on the reverse strand."""
    return (side == "left") == (transcript.strand == "+")


def _contains(feature, breakend):
    """Whether ``feature`` (a transcript or gene) spans ``breakend``."""
    return (str(feature.contig) == str(breakend.contig)
            and feature.start <= breakend.position <= feature.end)


def _cached(variant, key, compute):
    """``compute()``, memoized in the variant's annotation cache under
    ``key``. Results that depend only on the variant (mate-locus
    transcripts, cryptic-exon candidates) are then computed once per
    variant rather than once per transcript it's annotated on."""
    cache = getattr(variant, "_annotation_cache", None)
    if cache is None:
        return compute()
    if key not in cache:
        cache[key] = compute()
    return cache[key]


def _exonic_bases_before(transcript, position):
    """``(offset, position_is_exonic)`` for ``position`` in
    ``transcript``: how many cDNA bases lie transcript-5' of it, and
    whether the base at ``position`` is itself exonic.

    "Transcript-5' of" means a lower genomic coordinate on the forward
    strand and a higher one on the reverse strand.
    """
    offset = 0
    reverse = transcript.on_backward_strand
    for exon in transcript.exons:
        length = exon_length(exon)
        if reverse:
            if exon.start > position:
                offset += length
            elif exon.end >= position:
                return offset + (exon.end - position), True
            else:
                return offset, False
        else:
            if exon.end < position:
                offset += length
            elif exon.start <= position:
                return offset + (position - exon.start), True
            else:
                return offset, False
    return offset, False


def _cdna_cut(transcript, position, keep_five_prime):
    """The cDNA offset where a junction at ``position`` cuts
    ``transcript``: the 5' side keeps ``sequence[:cut]`` and the 3' side
    ``sequence[cut:]``.

    A breakend names the last base its side keeps, so an exonic
    breakpoint keeps that base on whichever side retains it. An
    intronic breakpoint snaps to the exon boundary: the 5' side ends
    with the last exon transcript-5' of the breakpoint, and the 3' side
    starts with the first exon transcript-3' of it.
    """
    offset, exonic = _exonic_bases_before(transcript, position)
    return offset + 1 if keep_five_prime and exonic else offset


def _warn_on_unknown_breakend_orientation(variant):
    """Warn when a breakend with a mate has no well-formed breakend
    ALT, so the kept sides are unknown and the annotator falls back
    to treating the annotated transcript as the 5' fusion partner."""
    import warnings
    warnings.warn(
        "Structural prediction can't read the breakend orientation "
        "of %s from ALT %r; treating the annotated transcript as the "
        "5' fusion partner. Use a VCF breakend ALT such as N[17:100[ "
        "or pass StructuralVariant.alt_assembly." % (
            variant.short_description, variant.symbolic_alt),
        stacklevel=3)


def _build_alt_assembly_mutant_transcript(variant, transcript):
    """If ``variant.alt_assembly`` is populated, build a
    single-:class:`ReferenceSegment` :class:`MutantTranscript` wrapping
    the assembled allele (#338). Returns ``None`` otherwise.

    The assembly is preferred over breakpoint-reconstructed cDNA
    because it reflects the molecule actually observed rather than
    the inference from reference + breakpoints. Downstream consumers
    read ``mutant_transcript.cdna_sequence`` uniformly.
    """
    assembly = getattr(variant, "alt_assembly", None)
    if not assembly:
        return None
    return MutantTranscript.from_sequence(
        assembly,
        reference_transcript=transcript,
        annotator_name="structural_variant",
        evidence={"source": "alt_assembly"}, label="alt_assembly")


def _build_deletion_mutant_transcript(variant, transcript):
    """Build a :class:`MutantTranscript` for a ``<DEL>`` SV.

    The cDNA of the surviving transcript is the concatenation of
    cDNA ranges outside the deleted span. cDNA is derivable entirely
    from pyensembl-cached transcript sequence — no genomic FASTA
    needed. ``mutant_protein_sequence`` is left to downstream
    consumers (translation requires knowing whether the CDS start
    survives, frame preservation across the junction, etc.).
    """
    kept = _cdna_ranges_kept_after_deletion(variant, transcript)
    if not kept:
        # Whole transcript deleted — express as a single zero-length
        # segment so the MutantTranscript still has a well-defined
        # shape (consumers see reference_segments=() and cdna="").
        return MutantTranscript(
            reference_transcript=transcript,
            reference_segments=(),
            cdna_sequence="",
            annotator_name="structural_variant")
    segments = tuple(
        ReferenceSegment(
            source=transcript, start=s, end=e, strand="+", label="del_kept")
        for s, e in kept)
    cdna = str(transcript.sequence)
    joined = "".join(cdna[s:e] for s, e in kept)
    return MutantTranscript(
        reference_transcript=transcript,
        reference_segments=segments,
        cdna_sequence=joined,
        annotator_name="structural_variant")


def _build_duplication_mutant_transcript(variant, transcript):
    """Build a :class:`MutantTranscript` for a ``<DUP>`` SV.

    Tandem-duplication interpretation: the duplicated body is
    inserted once between the surviving pre- and post-duplication
    transcript segments, yielding a transcript cDNA that carries an
    additional copy of the exonic ranges inside the SV span.
    """
    inside = _cdna_ranges_within_sv(variant, transcript)
    if not inside:
        return None
    full_cdna = str(transcript.sequence)
    full_len = len(full_cdna)
    last_inside = inside[-1][1]
    body_segments = tuple(
        ReferenceSegment(
            source=transcript, start=s, end=e,
            strand="+", label="dup_body_copy")
        for s, e in inside)
    segments = (
        ReferenceSegment(
            source=transcript, start=0, end=last_inside,
            strand="+", label="pre_dup_including_body"),
    ) + body_segments + (
        ReferenceSegment(
            source=transcript, start=last_inside, end=full_len,
            strand="+", label="post_dup"),
    )
    dup_body = "".join(full_cdna[s:e] for s, e in inside)
    joined = full_cdna[:last_inside] + dup_body + full_cdna[last_inside:]
    return MutantTranscript(
        reference_transcript=transcript,
        reference_segments=segments,
        cdna_sequence=joined,
        annotator_name="structural_variant")


def _build_inversion_mutant_transcript(variant, transcript):
    """Build a :class:`MutantTranscript` for an ``<INV>`` SV.

    The inverted exonic body is represented as one or more
    ``strand='-'`` :class:`ReferenceSegment` entries; consumers that
    want the assembled cDNA apply reverse-complement per segment.
    ``cdna_sequence`` is intentionally left None here — correctly
    assembling an inversion across exon boundaries requires
    care and is deferred to a follow-up (the shape lands here so
    downstream tools can see the inverted segment layout).
    """
    inside = _cdna_ranges_within_sv(variant, transcript)
    if not inside:
        return None
    first_inside_start = inside[0][0]
    last_inside_end = inside[-1][1]
    full_len = len(str(transcript.sequence))
    segments = (
        ReferenceSegment(
            source=transcript, start=0, end=first_inside_start,
            strand="+", label="pre_inv"),) + tuple(
        ReferenceSegment(
            source=transcript, start=s, end=e,
            strand="-", label="inv_body")
        for s, e in inside) + (
        ReferenceSegment(
            source=transcript, start=last_inside_end, end=full_len,
            strand="+", label="post_inv"),)
    return MutantTranscript(
        reference_transcript=transcript,
        reference_segments=segments,
        cdna_sequence=None,
        annotator_name="structural_variant")


def _translate_fused_cdna(fused_cdna, five_prime_transcript, five_prime_len):
    """Translate ``fused_cdna`` from the 5' partner's CDS start when
    that start codon lies in the retained 5' portion (#336).

    Returns the translated protein (stopping at the first stop
    codon) or ``None`` when the CDS start is past the breakpoint or
    translation otherwise fails.
    """
    from .codon_tables import (
        codon_table_for_transcript,
        translate_sequence,
    )
    if not five_prime_transcript.complete:
        return None
    try:
        cds_start = min(five_prime_transcript.start_codon_spliced_offsets)
    except (AttributeError, ValueError):
        return None
    if cds_start >= five_prime_len:
        # Start codon is past the breakpoint — the fusion loses it.
        return None
    codon_table = codon_table_for_transcript(five_prime_transcript)
    coding = fused_cdna[cds_start:]
    truncated = coding[:(len(coding) // 3) * 3]
    try:
        return translate_sequence(
            truncated, codon_table=codon_table, to_stop=True)
    except ValueError:
        return None


def _build_fusion_mutant_transcript(
        reference_transcript,
        five_prime_transcript, five_prime_position,
        three_prime_transcript, three_prime_position):
    """Build a :class:`MutantTranscript` for a :class:`GeneFusion`
    (#336).

    The fused cDNA is the 5' partner's cDNA up to its breakpoint
    followed by the 3' partner's from its breakpoint on, each keeping
    the base at its breakpoint. When the 5' partner's start codon lies
    in the retained portion, the protein is translated through the
    junction.

    ``reference_transcript`` is the transcript being annotated, which
    may be either partner, so an effect and its mutant transcript
    describe the same transcript; each partner's contribution is in the
    segments. The annotator assigns the 5' / 3' roles from breakend
    orientation and strand before calling this.
    """
    five_prime_cdna = str(five_prime_transcript.sequence)
    three_prime_cdna = str(three_prime_transcript.sequence)
    five_prime_end = _cdna_cut(
        five_prime_transcript, five_prime_position, True)
    three_prime_start = _cdna_cut(
        three_prime_transcript, three_prime_position, False)
    fused_cdna = (
        five_prime_cdna[:five_prime_end]
        + three_prime_cdna[three_prime_start:])
    segments = (
        ReferenceSegment(
            source=five_prime_transcript, start=0, end=five_prime_end,
            strand="+", label="5p_partner"),
        ReferenceSegment(
            source=three_prime_transcript,
            start=three_prime_start, end=len(three_prime_cdna),
            strand="+", label="3p_partner"),
    )
    return MutantTranscript(
        reference_transcript=reference_transcript,
        reference_segments=segments,
        cdna_sequence=fused_cdna,
        mutant_protein_sequence=_translate_fused_cdna(
            fused_cdna, five_prime_transcript, five_prime_end),
        annotator_name="structural_variant")


def _build_translocation_mutant_transcript(variant, transcript):
    """Describe only the retained local reference fragment of a BND.

    The other side and full mutant cDNA/protein remain unknown. Without
    local orientation or a breakpoint on this transcript, no fragment can
    be established. Supplied assemblies are handled by the caller.
    """
    side = _breakend_local_side(variant.symbolic_alt)
    if (side is None or str(variant.contig) != str(transcript.contig)
            or not transcript.start <= variant.start <= transcript.end):
        return None
    keep_five_prime = _retains_five_prime_end(transcript, side)
    cut = _cdna_cut(transcript, variant.start, keep_five_prime)
    full_len = len(str(transcript.sequence))
    start, end = (0, cut) if keep_five_prime else (cut, full_len)
    segments = (
        ReferenceSegment(
            source=transcript, start=start, end=end, strand="+",
            label="translocation_5p" if keep_five_prime else "translocation_3p"),)
    return MutantTranscript(
        reference_transcript=transcript,
        reference_segments=segments,
        cdna_sequence=None,
        evidence={"sequence_status": "retained_reference_fragment"},
        annotator_name="structural_variant")


def predict_structural_variant_effect(variant, transcript):
    """Classify ``variant`` on ``transcript``. Returns a single
    effect (typically a ``MultiOutcomeEffect`` subclass); consume
    ``effect.candidates`` for the full outcome set.
    """
    from pyensembl import Transcript
    if not isinstance(transcript, Transcript):
        raise TypeError(
            "Expected pyensembl.Transcript, got %s" % type(transcript))

    sv_type = getattr(variant, "sv_type", None)
    if (not getattr(variant, "is_structural", False)
            or (sv_type != "BND" and sv_type not in _SPAN_EFFECTS)):
        return NotImplemented

    if not transcript.is_protein_coding:
        return NoncodingTranscript(variant, transcript)

    # A caller-resolved allele, when there is one, is preferred over
    # anything inferred from breakpoints; build it once per call.
    assembly = _build_alt_assembly_mutant_transcript(variant, transcript)

    if sv_type == "BND":
        effect = _annotate_breakend(variant, transcript, assembly)
    elif sv_type in _SPAN_EFFECTS:
        span_effect = _SPAN_EFFECTS[sv_type](
            variant, transcript, assembly)
        # A span with one end in this transcript and a sense-to-sense
        # partner at the other end is a fusion. What the span does to
        # this transcript's own exons stays on the fusion as a further
        # primary candidate. Insertions and CNVs have no junction, so
        # they never fuse.
        effect = _fusion_across_junction(
            variant, transcript, assembly)
        if effect is None:
            effect = span_effect
        elif isinstance(span_effect, StructuralVariantEffect):
            effect._attach_primary_effects((span_effect,))

    # Attach cryptic-exon candidates enumerated from flanking
    # sequence / long-read assembly (#337). They show up as
    # additional EffectCandidate entries with source="varcode_motif".
    _enumerate_and_attach_cryptics(variant, effect)
    # Attach splice-outcome candidates when any SV breakpoint
    # lands in a canonical splice window on the transcript
    # (#341). They show up as additional EffectCandidate entries with
    # source="varcode_splice".
    _enumerate_and_attach_splice_outcomes(variant, transcript, effect)
    return effect


def _enumerate_and_attach_splice_outcomes(variant, transcript, effect):
    """If any breakpoint of ``variant`` lands in a canonical
    splice window on ``transcript``, synthesize the matching
    splice-disrupting effect, feed it to
    :func:`~varcode.splice_outcomes.enumerate_splice_outcomes`,
    and attach the returned candidates (minus NormalSplicing,
    which the SV primary classification already covers) to
    ``effect`` (#341).

    No-op when ``effect`` isn't an
    :class:`StructuralVariantEffect` (Intronic/Intergenic
    fall-throughs) or when no breakpoint falls in a splice
    window.
    """
    from .effect_classes import (
        NormalSplicing, StructuralVariantEffect)
    from ..effect_candidates import EffectCandidate
    from ..splice_outcomes import enumerate_splice_outcomes
    if not isinstance(effect, StructuralVariantEffect):
        return
    affected_start, affected_end = _affected_span(variant)
    event_end = getattr(variant, "end", variant.start)
    if (affected_start, affected_end) != (variant.start, event_end):
        # Typed breakend pairs retain their original junction positions
        # separately from the bases affected between them.
        bp_positions = {
            position for _, position in variant.breakpoints}
    else:
        # Preserve the established symbolic-SV interpretation of POS
        # and END as the positions whose splice windows are inspected.
        bp_positions = {variant.start, event_end}
    sv_type = getattr(variant, "sv_type", None)
    sv_evidence = {"sv_type": sv_type} if sv_type is not None else {}
    attached = []
    seen_exons = set()
    for bp in bp_positions:
        window = _breakpoint_splice_window(transcript, bp)
        if window is None:
            continue
        splice_cls, exon, distance = window
        # If two endpoints land in the same exon's splice window,
        # only emit one set of splice candidates — the outcomes
        # are identical.
        if exon.exon_id in seen_exons:
            continue
        seen_exons.add(exon.exon_id)
        try:
            synthetic = splice_cls(
                variant=variant,
                transcript=transcript,
                nearest_exon=exon,
                distance_to_exon=distance)
            splice_set = enumerate_splice_outcomes(synthetic)
        except (AttributeError, KeyError, ValueError):
            continue
        for candidate in splice_set.candidates:
            if isinstance(candidate.effect, NormalSplicing):
                # Primary SV classification already covers the
                # "splicing proceeds normally" interpretation.
                continue
            attached.append(EffectCandidate(
                effect=candidate.effect,
                source="varcode_splice",
                evidence={**dict(candidate.evidence), **sv_evidence}))
    if attached:
        effect._attach_splice_outcomes(attached)


def _enumerate_and_attach_cryptics(variant, effect):
    """Enumerate cryptic-exon candidates around the SV and hand
    them to ``effect._attach_cryptic_candidates`` (#337). The
    annotator is the only legitimate caller of that method —
    separating enumeration (here) from storage (on the effect)
    keeps the motif-scoring dependency out of ``effect_classes``.

    No-op when ``effect`` isn't an
    :class:`StructuralVariantEffect` (e.g. intronic
    fall-through) or when the enumerator returns nothing.
    """
    from ..cryptic_exons import enumerate_from_structural_variant
    from .effect_classes import StructuralVariantEffect
    if not isinstance(effect, StructuralVariantEffect):
        return
    def enumerate_candidates():
        try:
            return enumerate_from_structural_variant(variant)
        except (AttributeError, KeyError, ValueError, OSError):
            # Genome sequence fetch failed (no FASTA cached, unknown
            # contig, etc.) — skip; consumers still get the primary
            # SV classification.
            return []

    # The candidates depend only on the variant, not the transcript.
    candidates = _cached(
        variant, "cryptic_candidates", enumerate_candidates)
    if candidates:
        effect._attach_cryptic_candidates(candidates)

# -- classification helpers -----------------------------------------


def _annotate_deletion(variant, transcript, assembly=None):
    """A deletion overlapping a transcript: :class:`LargeDeletion`
    when it covers one or more exons, :class:`Intronic` when it's
    purely intronic."""
    affected = _overlapping_exons(variant, transcript)
    if not affected:
        return _intronic_or_intergenic(variant, transcript)
    return LargeDeletion(
        variant=variant,
        transcript=transcript,
        affected_exons=affected,
        mutant_transcript=assembly or _build_deletion_mutant_transcript(
            variant, transcript))


def _annotate_duplication(variant, transcript, assembly=None):
    affected = _overlapping_exons(variant, transcript)
    if not affected:
        return _intronic_or_intergenic(variant, transcript)
    return LargeDuplication(
        variant=variant,
        transcript=transcript,
        affected_exons=affected,
        mutant_transcript=assembly or _build_duplication_mutant_transcript(
            variant, transcript))


def _annotate_inversion(variant, transcript, assembly=None):
    affected = _overlapping_exons(variant, transcript)
    if not affected:
        return _intronic_or_intergenic(variant, transcript)
    return Inversion(
        variant=variant,
        transcript=transcript,
        mutant_transcript=assembly or _build_inversion_mutant_transcript(
            variant, transcript))


def _annotate_insertion(variant, transcript, assembly=None):
    affected = _overlapping_exons(variant, transcript)
    if not affected:
        return _intronic_or_intergenic(variant, transcript)
    return LargeDuplication(
        variant=variant,
        transcript=transcript,
        affected_exons=affected,
        mutant_transcript=assembly or _build_duplication_mutant_transcript(
            variant, transcript))


def _annotate_breakend(variant, transcript, assembly=None):
    """A breakend whose breakpoint lies on ``transcript``. If the
    join links this transcript sense-to-sense with a protein-coding
    transcript at the mate, report :class:`GeneFusion`; otherwise
    :class:`TranslocationToIntergenic`.
    """
    if variant.junctions:
        # With a resolved alt_assembly the caller has already
        # settled orientation, so an unreadable ALT isn't worth a
        # warning.
        if (assembly is None
                and breakend_sides(variant.symbolic_alt) is None):
            _warn_on_unknown_breakend_orientation(variant)
        fusion = _fusion_across_junction(
            variant, transcript, assembly)
        if fusion is not None:
            return fusion
    return TranslocationToIntergenic(
        variant=variant,
        transcript=transcript,
        mutant_transcript=assembly or (
            _build_translocation_mutant_transcript(
                variant, transcript)))

# -- utilities ------------------------------------------------------


def _overlapping_exons(variant, transcript):
    """Return the exons (in transcript order) that overlap the SV
    span on this transcript's contig. Contig mismatch => empty."""
    if str(variant.contig) != str(transcript.contig):
        return []
    var_start, var_end = _affected_span(variant)
    overlapping = []
    for exon in transcript.exons:
        if exon.end < var_start or exon.start > var_end:
            continue
        overlapping.append(exon)
    return overlapping


def _intronic_or_intergenic(variant, transcript):
    """Classify an SV that doesn't overlap any exon — either
    intronic (breakpoints inside the transcript envelope) or
    intergenic."""
    var_start, var_end = _affected_span(variant)
    if (str(variant.contig) == str(transcript.contig)
            and var_start >= transcript.start
            and var_end <= transcript.end):
        # Inside the transcript envelope but outside every exon
        # => intronic.  Distance-to-exon detail is best-effort:
        # find the nearest exon start/end for diagnostic use.
        nearest_exon = min(
            transcript.exons,
            key=lambda e: min(
                abs(e.start - var_start),
                abs(e.end - var_start)))
        distance = min(
            abs(nearest_exon.start - var_start),
            abs(nearest_exon.end - var_start))
        return Intronic(
            variant=variant,
            transcript=transcript,
            nearest_exon=nearest_exon,
            distance_to_exon=distance)
    return Intergenic(variant=variant)


def _fusion_across_junction(variant, transcript, assembly):
    """The :class:`GeneFusion` ``transcript`` forms at one of
    ``variant``'s junctions, or ``None``.

    A junction end belongs to this transcript when the transcript
    spans it and the transcript's *gene* doesn't span the other
    end: an event with both ends in one gene is intragenic, however
    short the isoform being annotated. The kept side and the strand
    then give the transcript's role — keeping its 5' end makes it
    the 5' partner — and the partner at the other end has to keep
    the opposite end, which is what makes the join sense-to-sense.
    Junctions where this transcript is the 5' partner are tried
    first. A kept side of ``None`` (unreadable ALT) treats the
    transcript as the 5' partner and accepts any partner.
    """
    gene = transcript.gene

    def keeps_five_prime(breakend):
        return breakend.keeps is None or _retains_five_prime_end(
            transcript, breakend.keeps)

    ends = []
    for junction in getattr(variant, "junctions", ()):
        for near, far in (junction, junction[::-1]):
            if not _contains(transcript, near):
                continue
            if _contains(gene, far):
                continue
            ends.append((near, far))

    for near, far in sorted(
            ends, key=lambda pair: not keeps_five_prime(pair[0])):
        five_prime = keeps_five_prime(near)
        partner = _fusion_partner(
            variant, transcript, near, far, five_prime)
        if partner is None:
            continue
        if five_prime:
            five, five_position = transcript, near.position
            three, three_position = partner, far.position
        else:
            five, five_position = partner, far.position
            three, three_position = transcript, near.position
        return GeneFusion(
            variant=variant,
            transcript=transcript,
            partner_transcript=partner,
            five_prime_transcript=five,
            three_prime_transcript=three,
            mutant_transcript=assembly or _build_fusion_mutant_transcript(
                transcript, five, five_position, three, three_position))
    return None


def _fusion_partner(
        variant, transcript, near, far, transcript_is_five_prime):
    """The protein-coding transcript at ``far`` that fuses with
    ``transcript``, or ``None``.

    A partner is in another gene, keeps the end opposite to
    ``transcript``'s, and — like ``transcript`` — belongs to a gene
    that doesn't span both ends of the junction. An unknown kept
    side accepts any partner.
    """
    for candidate in _coding_transcripts_at(
            variant, far.contig, far.position):
        if candidate.gene_id == transcript.gene_id:
            continue
        if _contains(candidate.gene, near):
            continue
        if far.keeps is None or (
                _retains_five_prime_end(candidate, far.keeps)
                != transcript_is_five_prime):
            return candidate
    return None


def _coding_transcripts_at(variant, contig, position):
    """Protein-coding transcripts overlapping ``contig:position`` in
    the variant's genome, in pyensembl's order. Cached on the
    variant, since every transcript at one end of a junction asks
    the same question about the other end."""
    def lookup():
        try:
            transcripts = variant.genome.transcripts_at_locus(
                str(contig), int(position), int(position))
        except ValueError:
            # A contig the genome doesn't know, e.g. a mate on a
            # decoy or unplaced contig.
            return []
        return [t for t in transcripts if t.is_protein_coding]

    return _cached(
        variant, ("coding_transcripts", str(contig), int(position)),
        lookup)


# Reference-span events share the same implementation regardless of caller.
# Preserve the existing CNV-gain and symbolic-insertion interpretations.
_SPAN_EFFECTS = {
    "DEL": _annotate_deletion,
    "DUP": _annotate_duplication,
    "INV": _annotate_inversion,
    "CNV": _annotate_duplication,
    "INS": _annotate_insertion,
}
