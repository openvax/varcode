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
    first_inside = inside[0][0]
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
    # Confirm the segments and the assembled cDNA agree — catches
    # any future drift between the two representations.
    assembled = "".join(full_cdna[s.start:s.end] for s in segments)
    assert assembled == joined, (
        "DUP segment layout inconsistent with assembled cDNA: "
        "segments=%r inside=%r" % (
            [(s.start, s.end) for s in segments], inside))
    # Reference `first_inside` in a way that documents the body
    # origin (used in `dup_body_copy` labels) without leaving a
    # dead variable behind.
    del first_inside
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


def _build_fusion_mutant_transcript(variant, transcript, partner):
    """Build a :class:`MutantTranscript` for a BND classified as
    :class:`GeneFusion`.

    Two segments: the 5' partner transcript (this transcript, up to
    its end) and the 3' partner transcript (from its start). The
    exact breakpoint-within-cDNA mapping and the fused-protein
    computation belong in #336 — here we fix the segment shape so
    consumers know a fusion has a 2-segment MutantTranscript.
    ``cdna_sequence`` is left None because the junction isn't
    resolved yet.
    """
    full_5p = len(str(transcript.sequence))
    full_3p = len(str(partner.sequence))
    segments = (
        ReferenceSegment(
            source=transcript, start=0, end=full_5p,
            strand="+", label="5p_partner"),
        ReferenceSegment(
            source=partner, start=0, end=full_3p,
            strand="+", label="3p_partner"),
    )
    return MutantTranscript(
        reference_transcript=transcript,
        reference_segments=segments,
        cdna_sequence=None,
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
            return self._annotate_deletion(variant, transcript)
        if sv_type == "DUP":
            return self._annotate_duplication(variant, transcript)
        if sv_type == "INV":
            return self._annotate_inversion(variant, transcript)
        if sv_type in ("CNV",):
            # CNVs are ambiguous at the protein level (gain vs loss
            # depends on direction). Treat as duplication in the CNV-
            # gain case; the CN0 subtype is routed to deletion logic
            # at parse time (see sv_allele_parser).
            return self._annotate_duplication(variant, transcript)
        if sv_type == "BND":
            return self._annotate_breakend(variant, transcript)
        if sv_type == "INS":
            # A large symbolic insertion is *structurally* an SV but
            # functionally resembles an in-frame-or-frameshift
            # insertion if the sequence is known. For now we defer
            # to the deletion shape (reporting the anchor codon) —
            # callers with ``alt_assembly`` populated get a richer
            # result once the SV annotator materializes segments.
            return self._annotate_insertion(variant, transcript)

        # Unknown SV type — fall back to intergenic so the call
        # graph doesn't crash.
        return Intergenic(variant)

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
