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

"""``StructuralVariantAnnotator`` — minimal, pluggable SV effect
classifier (PR 10; #252 / #259 / #305).

Input: a :class:`~varcode.StructuralVariant` and a
``pyensembl.Transcript``.

Output: one of the SV effect classes from
:mod:`varcode.effects.effect_classes` (``LargeDeletion``,
``LargeDuplication``, ``Inversion``, ``GeneFusion``,
``TranslocationToIntergenic``), each of which is a
:class:`~varcode.effects.MultiOutcomeEffect` exposing
:attr:`outcomes` — a tuple of :class:`~varcode.outcomes.Outcome`
entries each carrying an effect + probability + source + evidence.

Scope
-----

This annotator is deliberately shallow:

* It identifies which transcript(s) an SV overlaps and classifies
  the top-level consequence (deletion of exons, fusion candidate,
  intergenic translocation, etc.).
* It does *not* compute the exact fused protein sequence, predict
  whether a cryptic exon will be retained, or rank outcomes by
  likelihood. Those require external tools (SpliceAI, RNA evidence
  via Isovar, or long-read assembly) that attach additional
  :class:`Outcome` entries with their own ``source`` and ``evidence``.

The point of shipping this annotator now, even with thin logic, is
to make the *shape* of SV output available to the rest of the
ecosystem: the ``MultiOutcomeEffect`` + ``Outcome`` contract is
what lets RNA and assembly tools integrate cleanly.

Integration hooks
-----------------

* **SV with ``alt_assembly`` populated** (long-read resolution):
  the annotator wraps the resolved allele as a single
  :class:`~varcode.ReferenceSegment` on the returned
  :class:`MutantTranscript`, so consumers that compute fused-protein
  sequences have a concrete assembled cDNA to read.
* **External splice predictor**: call the annotator, then wrap each
  returned effect in a new :class:`MultiOutcomeEffect` whose
  ``outcomes`` tuple includes a fresh ``Outcome(effect=cryptic,
  source="spliceai", probability=...)`` entry.
* **Short-read RNA evidence**: same pattern. Attach an
  ``Outcome(effect=existing, source="isovar", evidence={
  "junction_reads": N})`` entry alongside the varcode-nominated
  outcome.
"""

from ..effects.effect_classes import (
    GeneFusion,
    Intergenic,
    Intronic,
    Inversion,
    LargeDeletion,
    LargeDuplication,
    NoncodingTranscript,
    TranslocationToIntergenic,
)
from ..mutant_transcript import MutantTranscript, ReferenceSegment
from ..version import __version__ as _varcode_version


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
        length = exon.end - exon.start + 1
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


def _cdna_ranges_kept_after_deletion(variant, transcript):
    """Compute the cDNA ranges of ``transcript`` that survive after
    the genomic deletion described by ``variant``.

    Returns a merged list of ``(cdna_start, cdna_end)`` tuples in
    transcript order. Handles both forward and reverse strand
    transcripts and the case where the deletion cuts through the
    middle of an exon.
    """
    del_start = min(variant.start, variant.end)
    del_end = max(variant.start, variant.end)
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
    """cDNA ranges that fall INSIDE the SV span (used for DUP to
    derive the duplicated body and for INV to derive the flipped
    middle)."""
    sv_start = min(variant.start, variant.end)
    sv_end = max(variant.start, variant.end)
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


_REVERSE_COMPLEMENT_BND_ORIENTATIONS = frozenset({"[]", "]["})


def _warn_on_reverse_complement_orientation(variant):
    """Emit a single warning when the breakend's ``mate_orientation``
    indicates a reverse-complement join (``[]`` or ``][``) that
    ``_build_fusion_mutant_transcript`` doesn't yet handle correctly
    (#336). Canonical ``]]`` / ``[[`` pairings, or a missing
    orientation, pass silently.
    """
    import warnings
    orientation = getattr(variant, "mate_orientation", None)
    if orientation is None:
        return
    if orientation in _REVERSE_COMPLEMENT_BND_ORIENTATIONS:
        warnings.warn(
            "StructuralVariantAnnotator treats BND %r as a canonical "
            "5p→3p fusion; reverse-complement orientations (%s) are "
            "not yet interpreted. Pass a resolved sequence via "
            "StructuralVariant.alt_assembly to override." % (
                variant.short_description, orientation),
            stacklevel=3)


class _AssembledAllele:
    """Minimal :class:`ReferenceSegment.source` adapter for a
    caller-supplied assembled allele (long-read resolution, targeted
    local assembly, etc.). Exposes ``sequence`` so downstream segment
    consumers can slice into it the same way they slice a pyensembl
    ``Transcript``'s cDNA.
    """

    __slots__ = ("sequence",)

    def __init__(self, sequence):
        self.sequence = sequence

    def __repr__(self):
        length = len(self.sequence) if self.sequence else 0
        return "_AssembledAllele(len=%d)" % length


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
    source = _AssembledAllele(assembly)
    segments = (
        ReferenceSegment(
            source=source,
            start=0,
            end=len(assembly),
            strand="+",
            label="alt_assembly"),
    )
    return MutantTranscript(
        reference_transcript=transcript,
        reference_segments=segments,
        cdna_sequence=assembly,
        annotator_name="structural_variant",
        evidence={"source": "alt_assembly"})


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


def _cdna_offset_at_5p_breakpoint(transcript, breakpoint_pos):
    """Return the cDNA offset in the 5' partner at which fusion
    splits the transcript (#336).

    All cDNA bases at offsets ``[0, returned_offset)`` are retained
    in the fused product. Handles forward and reverse strand
    transcripts; for an intronic breakpoint the walker snaps to the
    end of the last transcript-order exon upstream of the
    breakpoint (the biologically expected behavior — the
    spliceosome completes the preceding exon before the junction).

    Boundary convention: a breakpoint exactly AT the 5' terminal
    base of an exon (``breakpoint_pos == exon.start`` on forward
    strand, ``== exon.end`` on reverse strand) excludes that exon
    from the retained 5' portion, so the 5p / 3p retained offsets
    partition the transcript cDNA without overlap.
    """
    offset = 0
    reverse = transcript.on_backward_strand
    for exon in transcript.exons:
        exon_len = exon.end - exon.start + 1
        if reverse:
            if exon.start > breakpoint_pos:
                offset += exon_len
            elif exon.end > breakpoint_pos:
                offset += exon.end - breakpoint_pos
                return offset
            else:
                return offset
        else:
            if exon.end < breakpoint_pos:
                offset += exon_len
            elif exon.start < breakpoint_pos:
                offset += breakpoint_pos - exon.start
                return offset
            else:
                return offset
    return offset


def _cdna_offset_at_3p_breakpoint(transcript, breakpoint_pos):
    """Return the cDNA offset in the 3' partner at which the fusion
    picks up (#336).

    All cDNA bases at offsets ``[returned_offset, end)`` are
    retained. For an intronic breakpoint the walker snaps to the
    start of the first transcript-order exon downstream of the
    breakpoint.

    Boundary convention mirrors :func:`_cdna_offset_at_5p_breakpoint`:
    a breakpoint AT the 5' terminal base of an exon keeps that exon
    in the 3' retained portion (its cDNA starts at the returned
    offset).
    """
    offset = 0
    reverse = transcript.on_backward_strand
    for exon in transcript.exons:
        exon_len = exon.end - exon.start + 1
        if reverse:
            if exon.start > breakpoint_pos:
                offset += exon_len
            elif exon.end > breakpoint_pos:
                return offset + (exon.end - breakpoint_pos)
            else:
                return offset
        else:
            if exon.end < breakpoint_pos:
                offset += exon_len
            elif exon.start < breakpoint_pos:
                return offset + (breakpoint_pos - exon.start)
            else:
                return offset
    return offset


def _translate_fused_cdna(fused_cdna, five_prime_transcript, five_prime_len):
    """Translate ``fused_cdna`` from the 5' partner's CDS start when
    that start codon lies in the retained 5' portion (#336).

    Returns the translated protein (stopping at the first stop
    codon) or ``None`` when the CDS start is past the breakpoint or
    translation otherwise fails.
    """
    from ..effects.codon_tables import (
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


def _build_fusion_mutant_transcript(variant, transcript, partner):
    """Build a :class:`MutantTranscript` for a BND classified as
    :class:`GeneFusion` (#336).

    The canonical case — intron:intron breakpoints joining a 5'
    partner's N-terminal exons to a 3' partner's C-terminal exons —
    yields a fused cDNA whose 5' portion is
    ``transcript.sequence[0:5p_cdna_offset]`` and whose 3' portion
    is ``partner.sequence[3p_cdna_offset:]``. When the 5' partner's
    start codon lies in the retained region, the fused protein is
    translated through the junction.

    BND ``mate_orientation`` semantics (which side of each
    breakpoint joins — ``]p]t`` vs ``t[p[`` etc.) are not yet
    consulted here; the builder assumes the canonical "5p
    N-terminus + 3p C-terminus" direction, which matches BRD4-NUTM1,
    BCR-ABL1, EWSR1-FLI1, and similar oncogenic fusions. Callers
    with reverse-complement breakend pairings can override by
    passing a pre-resolved ``alt_assembly`` on the variant.
    """
    five_prime_cdna = str(transcript.sequence)
    three_prime_cdna = str(partner.sequence)
    five_prime_offset = _cdna_offset_at_5p_breakpoint(
        transcript, variant.start)
    three_prime_offset = _cdna_offset_at_3p_breakpoint(
        partner, variant.mate_start)
    five_prime_retained = five_prime_cdna[:five_prime_offset]
    three_prime_retained = three_prime_cdna[three_prime_offset:]
    fused_cdna = five_prime_retained + three_prime_retained
    segments = (
        ReferenceSegment(
            source=transcript, start=0, end=five_prime_offset,
            strand="+", label="5p_partner"),
        ReferenceSegment(
            source=partner,
            start=three_prime_offset, end=len(three_prime_cdna),
            strand="+", label="3p_partner"),
    )
    mutant_protein = _translate_fused_cdna(
        fused_cdna, transcript, len(five_prime_retained))
    return MutantTranscript(
        reference_transcript=transcript,
        reference_segments=segments,
        cdna_sequence=fused_cdna,
        mutant_protein_sequence=mutant_protein,
        annotator_name="structural_variant")


def _build_translocation_mutant_transcript(variant, transcript):
    """Build a :class:`MutantTranscript` for a BND classified as
    :class:`TranslocationToIntergenic`.

    One segment covering the (this-side) transcript up to the
    breakpoint. Intergenic space beyond the breakpoint isn't
    representable as a transcript source; #338 and #341 extend this
    with caller-supplied genomic intervals / long-read assemblies.
    """
    full_len = len(str(transcript.sequence))
    segments = (
        ReferenceSegment(
            source=transcript, start=0, end=full_len,
            strand="+", label="translocation_5p"),)
    return MutantTranscript(
        reference_transcript=transcript,
        reference_segments=segments,
        cdna_sequence=None,
        annotator_name="structural_variant")


class StructuralVariantAnnotator:
    """Classify :class:`~varcode.StructuralVariant` consequences on
    a single transcript.

    ``supports`` advertises the SV kinds the annotator handles. The
    :class:`~varcode.FastEffectAnnotator` and
    :class:`~varcode.annotators.protein_diff.ProteinDiffEffectAnnotator`
    advertise ``{"snv", "indel", "mnv"}``; this one advertises the
    SV-type tokens used on :attr:`StructuralVariant.sv_type`.
    """

    name = "structural_variant"
    version = _varcode_version
    supports = frozenset({"DEL", "DUP", "INV", "INS", "CNV", "BND"})

    def __repr__(self):
        return "StructuralVariantAnnotator(name=%r, version=%r)" % (
            self.name, self.version)

    # -- public API -------------------------------------------------------

    def annotate_on_transcript(self, variant, transcript):
        """Classify ``variant`` on ``transcript``. Returns a single
        effect (typically a ``MultiOutcomeEffect`` subclass); consume
        ``effect.outcomes`` for the full outcome set.
        """
        from pyensembl import Transcript
        if not isinstance(transcript, Transcript):
            raise TypeError(
                "Expected pyensembl.Transcript, got %s" % type(transcript))

        if not transcript.is_protein_coding:
            return NoncodingTranscript(variant, transcript)

        sv_type = getattr(variant, "sv_type", None)

        if sv_type == "DEL":
            effect = self._annotate_deletion(variant, transcript)
        elif sv_type == "DUP":
            effect = self._annotate_duplication(variant, transcript)
        elif sv_type == "INV":
            effect = self._annotate_inversion(variant, transcript)
        elif sv_type in ("CNV",):
            # CNVs are ambiguous at the protein level (gain vs loss
            # depends on direction). Treat as duplication in the CNV-
            # gain case; the CN0 subtype is routed to deletion logic
            # at parse time (see sv_allele_parser).
            effect = self._annotate_duplication(variant, transcript)
        elif sv_type == "BND":
            effect = self._annotate_breakend(variant, transcript)
        elif sv_type == "INS":
            # A large symbolic insertion is *structurally* an SV but
            # functionally resembles an in-frame-or-frameshift
            # insertion if the sequence is known. For now we defer
            # to the deletion shape (reporting the anchor codon) —
            # callers with ``alt_assembly`` populated get a richer
            # result once the SV annotator materializes segments.
            effect = self._annotate_insertion(variant, transcript)
        else:
            # Unknown SV type — fall back to intergenic so the call
            # graph doesn't crash.
            return Intergenic(variant)

        # Attach cryptic-exon candidates enumerated from flanking
        # sequence / long-read assembly (#337). They show up as
        # additional Outcomes with source="varcode_motif".
        self._enumerate_and_attach_cryptics(variant, effect)
        return effect

    def _enumerate_and_attach_cryptics(self, variant, effect):
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
        from ..effects.effect_classes import StructuralVariantEffect
        if not isinstance(effect, StructuralVariantEffect):
            return
        try:
            candidates = enumerate_from_structural_variant(variant)
        except (AttributeError, KeyError, ValueError, OSError):
            # Genome sequence fetch failed (no FASTA cached, unknown
            # contig, etc.) — silently skip; consumers still get the
            # primary SV classification.
            return
        if candidates:
            effect._attach_cryptic_candidates(candidates)

    # -- classification helpers -----------------------------------------

    def _annotate_deletion(self, variant, transcript):
        """A deletion overlapping a transcript. If it covers one or
        more exons, report :class:`LargeDeletion`; if purely
        intronic, report :class:`Intronic`."""
        affected = self._overlapping_exons(variant, transcript)
        if not affected:
            return self._intronic_or_intergenic(variant, transcript)
        return LargeDeletion(
            variant=variant,
            transcript=transcript,
            affected_exons=affected,
            mutant_transcript=_build_alt_assembly_mutant_transcript(
                variant, transcript) or _build_deletion_mutant_transcript(
                variant, transcript))

    def _annotate_duplication(self, variant, transcript):
        affected = self._overlapping_exons(variant, transcript)
        if not affected:
            return self._intronic_or_intergenic(variant, transcript)
        return LargeDuplication(
            variant=variant,
            transcript=transcript,
            affected_exons=affected,
            mutant_transcript=_build_alt_assembly_mutant_transcript(
                variant, transcript) or _build_duplication_mutant_transcript(
                variant, transcript))

    def _annotate_inversion(self, variant, transcript):
        affected = self._overlapping_exons(variant, transcript)
        if not affected:
            return self._intronic_or_intergenic(variant, transcript)
        return Inversion(
            variant=variant,
            transcript=transcript,
            mutant_transcript=_build_alt_assembly_mutant_transcript(
                variant, transcript) or _build_inversion_mutant_transcript(
                variant, transcript))

    def _annotate_insertion(self, variant, transcript):
        # Large symbolic <INS>: treat like duplication-style disruption
        # for transcripts with overlapping exons.
        affected = self._overlapping_exons(variant, transcript)
        if not affected:
            return self._intronic_or_intergenic(variant, transcript)
        return LargeDuplication(
            variant=variant,
            transcript=transcript,
            affected_exons=affected,
            mutant_transcript=_build_alt_assembly_mutant_transcript(
                variant, transcript) or _build_duplication_mutant_transcript(
                variant, transcript))

    def _annotate_breakend(self, variant, transcript):
        """A breakend whose first breakpoint lies on ``transcript``.
        If the mate lands on another protein-coding transcript,
        report :class:`GeneFusion`; otherwise,
        :class:`TranslocationToIntergenic`.
        """
        mate_contig = getattr(variant, "mate_contig", None)
        mate_start = getattr(variant, "mate_start", None)
        assembly_mt = _build_alt_assembly_mutant_transcript(
            variant, transcript)
        # The reverse-complement warning only applies when we're
        # inferring the fused allele from reference + breakpoint; if
        # the caller supplied ``alt_assembly`` they already resolved
        # the orientation themselves.
        if assembly_mt is None:
            _warn_on_reverse_complement_orientation(variant)
        if mate_contig is None or mate_start is None:
            return TranslocationToIntergenic(
                variant=variant,
                transcript=transcript,
                mutant_transcript=assembly_mt or (
                    _build_translocation_mutant_transcript(
                        variant, transcript)))

        partner = self._coding_transcript_at(
            variant, mate_contig, mate_start)
        if partner is not None and partner.id != transcript.id:
            return GeneFusion(
                variant=variant,
                transcript=transcript,
                partner_transcript=partner,
                mutant_transcript=assembly_mt or (
                    _build_fusion_mutant_transcript(
                        variant, transcript, partner)))
        return TranslocationToIntergenic(
            variant=variant,
            transcript=transcript,
            mutant_transcript=assembly_mt or (
                _build_translocation_mutant_transcript(
                    variant, transcript)))

    # -- utilities ------------------------------------------------------

    def _overlapping_exons(self, variant, transcript):
        """Return the exons (in transcript order) that overlap the SV
        span on this transcript's contig. Contig mismatch => empty."""
        if str(variant.contig) != str(transcript.contig):
            return []
        var_start = variant.start
        var_end = getattr(variant, "end", None) or variant.start
        if var_end < var_start:
            var_start, var_end = var_end, var_start
        overlapping = []
        for exon in transcript.exons:
            if exon.end < var_start or exon.start > var_end:
                continue
            overlapping.append(exon)
        return overlapping

    def _intronic_or_intergenic(self, variant, transcript):
        """Classify an SV that doesn't overlap any exon — either
        intronic (breakpoints inside the transcript envelope) or
        intergenic."""
        if (str(variant.contig) == str(transcript.contig)
                and variant.start >= transcript.start
                and getattr(variant, "end", variant.start) <= transcript.end):
            # Inside the transcript envelope but outside every exon
            # => intronic.  Distance-to-exon detail is best-effort:
            # find the nearest exon start/end for diagnostic use.
            nearest_exon = min(
                transcript.exons,
                key=lambda e: min(
                    abs(e.start - variant.start),
                    abs(e.end - variant.start)))
            distance = min(
                abs(nearest_exon.start - variant.start),
                abs(nearest_exon.end - variant.start))
            return Intronic(
                variant=variant,
                transcript=transcript,
                nearest_exon=nearest_exon,
                distance_to_exon=distance)
        return Intergenic(variant=variant)

    def _coding_transcript_at(self, variant, contig, position):
        """Find the first protein-coding transcript overlapping
        ``contig:position`` in the variant's genome. Returns
        ``None`` if no coding transcript lives there.
        """
        genome = getattr(variant, "genome", None) or getattr(
            variant, "ensembl", None)
        if genome is None:
            return None
        try:
            contig_norm = str(contig)
            transcripts = genome.transcripts_at_locus(
                contig_norm, int(position), int(position))
        except Exception:
            return None
        for t in transcripts:
            try:
                if t.is_protein_coding:
                    return t
            except Exception:
                continue
        return None
