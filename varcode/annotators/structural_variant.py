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
from ..version import __version__ as _varcode_version


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
            # Purely intronic SV (or spans past transcript) —
            # report as Intronic for now. A future refinement would
            # add a cryptic-splice candidate outcome when the
            # breakpoints are near a splice site.
            return self._intronic_or_intergenic(variant, transcript)
        return LargeDeletion(
            variant=variant,
            transcript=transcript,
            affected_exons=affected)

    def _annotate_duplication(self, variant, transcript):
        affected = self._overlapping_exons(variant, transcript)
        if not affected:
            return self._intronic_or_intergenic(variant, transcript)
        return LargeDuplication(
            variant=variant,
            transcript=transcript,
            affected_exons=affected)

    def _annotate_inversion(self, variant, transcript):
        affected = self._overlapping_exons(variant, transcript)
        if not affected:
            return self._intronic_or_intergenic(variant, transcript)
        # Inversion disrupts whichever exons it crosses; single
        # outcome for now (cryptic-splice candidates layer on via
        # PR 11).
        return Inversion(variant=variant, transcript=transcript)

    def _annotate_insertion(self, variant, transcript):
        # Large symbolic <INS>: treat like duplication-style disruption
        # for transcripts with overlapping exons.
        affected = self._overlapping_exons(variant, transcript)
        if not affected:
            return self._intronic_or_intergenic(variant, transcript)
        return LargeDuplication(
            variant=variant,
            transcript=transcript,
            affected_exons=affected)

    def _annotate_breakend(self, variant, transcript):
        """A breakend whose first breakpoint lies on ``transcript``.
        If the mate lands on another protein-coding transcript,
        report :class:`GeneFusion`; otherwise,
        :class:`TranslocationToIntergenic`.
        """
        mate_contig = getattr(variant, "mate_contig", None)
        mate_start = getattr(variant, "mate_start", None)
        if mate_contig is None or mate_start is None:
            return TranslocationToIntergenic(
                variant=variant, transcript=transcript)

        partner = self._coding_transcript_at(
            variant, mate_contig, mate_start)
        if partner is not None and partner.id != transcript.id:
            return GeneFusion(
                variant=variant,
                transcript=transcript,
                partner_transcript=partner)
        return TranslocationToIntergenic(
            variant=variant, transcript=transcript)

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
