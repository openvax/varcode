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

"""Structural-variant, haplotype, and germline-phase consequence classes."""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

from .base import (
    MultiOutcomeEffect,
    MutationEffect,
    TranscriptMutationEffect,
    Unresolved,
)

if TYPE_CHECKING:
    from ..effect_candidates import EffectCandidate


# =====================================================================
# Structural-variant effects (#252 / PR 10)
#
# These describe consequences of SVs at the transcript level. A single
# SV commonly produces multiple plausible consequences (fusion vs exon
# loss vs readthrough), so each of these is a ``MultiOutcomeEffect``
# subclass exposing a :attr:`candidates` tuple per #299 / #382.
#
# The classes are deliberately thin wrappers — the important interface
# is :attr:`candidates`, which carries :class:`EffectCandidate` objects
# with source / evidence provenance so external tools (RNA evidence,
# SpliceAI, long-read assembly) can annotate them without subclassing.
# =====================================================================


class StructuralVariantEffect(TranscriptMutationEffect, MultiOutcomeEffect):
    """Base class for effects of a :class:`StructuralVariant` on a
    specific transcript. Local models contain consequence effects in
    :attr:`candidates`; their description follows the first prediction and
    their priority follows the highest-priority candidate. Gene fusions and
    translocations retain their specialized descriptions.

    :attr:`candidates` always returns a tuple of
    :class:`~varcode.effect_candidates.EffectCandidate` objects.
    Primary candidates carry ``source="varcode"``; cryptic-exon, splice,
    and RNA candidates are appended via :meth:`add_candidates` with their
    respective ``source`` tags (``"varcode_motif"``, ``"varcode_splice"``,
    ``"rna"``, etc.).
    """

    candidate_evidence = None

    @property
    def short_description(self):
        if self.primary_effects and self.primary_effects[0] is not self:
            return self.primary_effects[0].short_description
        return "unresolved-structural-transcript"

    def _format_short_description(self, deletion_notation):
        # Subclasses with dedicated labels (e.g. sv-deletion) retain them.
        if type(self).short_description is not StructuralVariantEffect.short_description:
            return self.short_description
        if self.primary_effects and self.primary_effects[0] is not self:
            return self.primary_effects[0].format_short_description(
                deletion_notation=deletion_notation)
        return self.short_description

    @property
    def priority_class(self):
        own = type(self) if type(self) is not StructuralVariantEffect else Unresolved
        if not self.primary_effects or self.primary_effects[0] is self:
            return own
        from .effect_ordering import effect_priority, transcript_effect_priority_dict

        def rank(candidate):
            # A candidate that is this set stands for its own class; ranking
            # it through this property would never terminate.
            if candidate.effect is self:
                return transcript_effect_priority_dict.get(own, -1)
            return effect_priority(candidate.effect)

        highest = max(self.candidates, key=rank).effect
        if highest is self:
            return own
        return getattr(highest, "priority_class", None) or type(highest)

    def __init__(
            self, variant, transcript,
            primary_effects=None, mutant_transcript=None):
        TranscriptMutationEffect.__init__(self, variant, transcript)
        self.primary_effects: tuple[MutationEffect, ...] = (
            tuple(primary_effects) if primary_effects is not None
            else (self,))
        self.mutant_transcript = mutant_transcript

    def add_primary_effects(self, effects: Iterable[MutationEffect]) -> None:
        """Append additional primary SV classifications (e.g. partner fusions or span effects)."""
        incoming = tuple(effects)
        if incoming:
            self.primary_effects = self.primary_effects + incoming

    def _get_mutant_protein_sequence(self) -> str | None:
        """Protein carried by the realized structural transcript, if any."""
        return getattr(
            self.mutant_transcript, "mutant_protein_sequence", None)

    def _get_modifies_coding_sequence(self) -> bool | None:
        """True if any candidate changes CDS, None if unresolved, else False.

        Protein equality alone does not establish coding-sequence equality.
        """
        from .sequence_change import modification_status
        return modification_status(self, "modifies_coding_sequence")

    def _get_modifies_protein_sequence(self) -> bool | None:
        """True if any candidate changes protein, None if unresolved, else False.

        A set with an unchanged primary and a changed alternative returns True.
        False requires every candidate to be known unchanged; None is not silent.
        """
        from .sequence_change import modification_status
        return modification_status(self, "modifies_protein_sequence")

    def _get_candidates(self) -> tuple[EffectCandidate, ...]:
        """Primary :class:`~varcode.EffectCandidate` entries for :attr:`primary_effects`.

        Additional motif, splice, and RNA candidates appended via
        :meth:`add_candidates` follow automatically.
        """
        from ..effect_candidates import EffectCandidate
        sv_type = getattr(self.variant, "sv_type", None)
        base_evidence = dict(self.candidate_evidence or {})
        if sv_type is not None:
            base_evidence["sv_type"] = sv_type
        return tuple(
            EffectCandidate(
                effect=effect,
                source="varcode",
                evidence=base_evidence)
            for effect in self.primary_effects)


class LargeDeletion(StructuralVariantEffect):
    """A deletion (``<DEL>`` / ``<CN0>``) that removes one or more
    exons — or an entire gene. Carries the list of affected exons
    for downstream analysis; the single ``outcome`` is this effect
    itself (callers that want to add RNA evidence or SpliceAI
    scoring construct additional outcomes and wrap the result)."""

    @property
    def short_description(self) -> str:
        return "sv-deletion"

    def __init__(
            self, variant, transcript, affected_exons,
            primary_effects=None, mutant_transcript=None):
        StructuralVariantEffect.__init__(
            self, variant, transcript,
            primary_effects=primary_effects,
            mutant_transcript=mutant_transcript)
        self.affected_exons = tuple(affected_exons)

    def __str__(self):
        return "LargeDeletion(%s, %s, exons=%d)" % (
            self.variant.short_description,
            self.transcript.id if self.transcript is not None else "?",
            len(self.affected_exons))


class LargeDuplication(StructuralVariantEffect):
    """A tandem duplication (``<DUP>``) overlapping exons.
    Biologically may produce a copy-number increase, a fused
    reading frame if junctions land in-frame, or a regulatory
    effect. Varcode reports the affected exons and leaves scoring
    to downstream tools."""

    @property
    def short_description(self) -> str:
        return "sv-duplication"

    def __init__(
            self, variant, transcript, affected_exons,
            primary_effects=None, mutant_transcript=None):
        StructuralVariantEffect.__init__(
            self, variant, transcript,
            primary_effects=primary_effects,
            mutant_transcript=mutant_transcript)
        self.affected_exons = tuple(affected_exons)


class Inversion(StructuralVariantEffect):
    """An inversion (``<INV>``) that flips a stretch of a transcript.
    Depending on whether breakpoints fall in exons or introns, the
    consequence is very different (exonic: likely disruptive; purely
    intronic: may or may not affect splicing). Reported as a single
    outcome here; subclasses / callers can enrich with cryptic-
    splice candidates per PR 11."""

    @property
    def short_description(self) -> str:
        return "sv-inversion"


class GeneFusion(StructuralVariantEffect):
    """A structural variant joining ``transcript`` to a sense-oriented
    transcript in another gene — the canonical fusion shape. Breakends
    produce it, and so do deletions, duplications and inversions with
    one end in each gene. A 3' partner that starts just past a junction
    end, rather than spanning it, is joined from its second exon and is
    only ever a further candidate behind the primary classification.

    The default annotator retains every compatible annotated partner isoform
    and junction in :attr:`candidates`, with no count cap. The first is the
    historical primary; producer order is not a probability ranking. Distinct
    transcript pairs remain separate even when their proteins match. Missing
    sequence leaves that candidate's transcript/protein unresolved.

    ``transcript`` is the transcript being annotated and
    ``partner_transcript`` the other one; ``five_prime_transcript``
    and ``three_prime_transcript`` say which is which (by default
    ``transcript`` is the 5' partner). When the annotator has enough
    context, :attr:`mutant_transcript` is built from
    :class:`ReferenceSegment` entries describing the fused
    allele. Predicting the exact fused-protein sequence requires
    knowing which exons are retained, which typically needs RNA
    evidence — outcomes beyond "this is a plausible fusion" are
    imported with :func:`varcode.make_fusion_outcome` or supplied by tools
    attaching :class:`EffectCandidate` objects with their own ``source`` tag.
    An RNA-observed gene fusion need not encode a complete protein.
    """

    @property
    def short_description(self) -> str:
        return "sv-gene-fusion"

    def __init__(
            self, variant, transcript, partner_transcript,
            mutant_transcript=None, primary_effects=None,
            five_prime_transcript=None, three_prime_transcript=None):
        StructuralVariantEffect.__init__(
            self, variant, transcript,
            primary_effects=primary_effects,
            mutant_transcript=mutant_transcript)
        self.partner_transcript = partner_transcript
        self.five_prime_transcript = (
            transcript if five_prime_transcript is None
            else five_prime_transcript)
        self.three_prime_transcript = (
            partner_transcript if three_prime_transcript is None
            else three_prime_transcript)


class TranslocationToIntergenic(StructuralVariantEffect):
    """A breakend that doesn't form a gene fusion: its mate lies in
    intergenic space, or the join meets a gene in an orientation that
    can't read sense-to-sense into it (e.g. two genes' 5' ends joined
    head to head). A junction with both ends in the transcript's gene is
    annotated as the local DEL, DUP or INV it describes instead. The downstream consequence depends on whether the
    sequence beyond the breakpoint contains cryptic splice / ORF
    signals. Fusions with a coding transcript that starts just past the
    mate follow as further candidates, and PR 11's cryptic-exon
    enumerator adds candidate outcomes when applicable."""

    @property
    def short_description(self) -> str:
        return "sv-translocation-intergenic"


class TranslocationToUnannotatedContig(TranslocationToIntergenic):
    """A breakend that forms no known gene fusion because its mate lies on
    a contig the annotation doesn't cover: a decoy, EBV, an unplaced or alt
    scaffold, or hg19's chrM under Ensembl's GRCh37 (whose MT is the rCRS,
    at other positions). Whether the mate is intergenic is unknown, so a
    fusion with it is not ruled out. ``partner_contig`` names the mate's
    contig. The local side is annotated as for
    :class:`TranslocationToIntergenic`, which this specializes."""

    @property
    def short_description(self) -> str:
        return "sv-translocation-unannotated-contig"

    def __init__(
            self, variant, transcript, partner_contig,
            mutant_transcript=None, primary_effects=None):
        TranslocationToIntergenic.__init__(
            self, variant, transcript,
            primary_effects=primary_effects,
            mutant_transcript=mutant_transcript)
        self.partner_contig = partner_contig


class CrypticExonCandidate(MutationEffect):
    """A region where an SV has brought novel sequence into range
    of the transcript, and motif scoring flags a plausible new
    splice acceptor / donor pair. Produced by PR 11's cryptic-exon
    enumerator; attached as additional :class:`EffectCandidate` entries on
    SV effects rather than standalone.

    Not a :class:`TranscriptMutationEffect` because the candidate
    region may not overlap any existing transcript — it's a *new*
    exon hypothesis. Carries the contig / interval and the motif
    scores as plain fields; external predictors (SpliceAI,
    Pangolin) attach their own scores via the enclosing
    :class:`EffectCandidate.evidence` dict.
    """

    @property
    def short_description(self) -> str:
        return "sv-cryptic-exon-candidate"

    _default_modifies_coding_sequence = None
    _default_modifies_protein_sequence = None

    def __init__(
            self, variant, contig, interval_start, interval_end,
            donor_score=None, acceptor_score=None):
        MutationEffect.__init__(self, variant)
        self.contig = contig
        self.interval_start = interval_start
        self.interval_end = interval_end
        self.donor_score = donor_score
        self.acceptor_score = acceptor_score


# =====================================================================
# Haplotype effects — joint consequence of multiple cis variants on
# the same transcript (#269).
#
# When a phase resolver (VCF PS, Isovar assembly, long-read haplotype
# caller) confirms two or more variants sit on the same physical copy
# of a transcript, their combined effect on the protein isn't just the
# sum of per-variant effects — they can share a codon, one can rescue
# another's frameshift, etc. HaplotypeEffect carries the joint
# MutantTranscript (via apply_variants_to_transcript) and coexists
# with the per-variant effects on the same EffectCollection so
# downstream consumers can pick the granularity they need.
# =====================================================================


class HaplotypeEffect(TranscriptMutationEffect, MultiOutcomeEffect):
    """Joint effect of two or more cis variants on the same
    transcript (#269).

    Emitted by :meth:`VariantCollection.effects` when a
    ``phase_resolver`` (VCF ``PS`` tag, Isovar assembly, or any
    other :class:`PhaseResolver`) groups cis variants together and
    :func:`varcode.mutant_transcript.apply_variants_to_transcript`
    can build a combined mutant cDNA. The per-variant effects stay
    on the collection — this is additive, not a replacement.

    Downstream consumers that want haplotype-level context (peptide
    prediction, ASE, compound het interpretation) read the
    :attr:`mutant_transcript` here. Consumers that want per-variant
    HGVS or annotation granularity iterate the per-variant effects
    as before.
    """

    def __init__(
            self, variants, transcript, mutant_transcript,
            phase_source=None):
        # The HaplotypeEffect is scoped to a set of variants — the
        # :class:`MutationEffect.variant` slot takes the first one
        # as a back-compat anchor so existing single-variant
        # consumers don't blow up, but real consumers read
        # :attr:`variants`.
        assert len(variants) >= 2, (
            "HaplotypeEffect requires ≥ 2 cis variants; got %d" % len(variants))
        TranscriptMutationEffect.__init__(self, variants[0], transcript)
        self.variants = tuple(variants)
        self.mutant_transcript = mutant_transcript
        # Which resolver produced the phase grouping (e.g. "vcf_ps"
        # for DNA-PS-tag phasing, "read_phasing" for an RNA / long-
        # read source). Useful for downstream filtering and for
        # auditing whether the cis call came from DNA-only phasing
        # or RNA-assembly evidence.
        self.phase_source = phase_source

    def _get_candidates(self) -> tuple[EffectCandidate, ...]:
        """Single-outcome wrapping: the haplotype IS the outcome,
        lifted into the uniform :class:`EffectCandidate` shape so
        consumers iterate :attr:`candidates` the same way across all
        :class:`MutationEffect` instances. Additional candidates appended
        via :meth:`add_candidates` (e.g. from RNA evidence) follow
        automatically.
        """
        from ..effect_candidates import EffectCandidate
        evidence = {}
        if self.phase_source is not None:
            evidence["phase_source"] = self.phase_source
        return (EffectCandidate(
            effect=self, source="varcode", evidence=evidence),)

    def _get_mutant_protein_sequence(self) -> str | None:
        """Joint protein translated from the haplotype's cDNA, or
        ``None`` if the edits don't land after the CDS start."""
        return getattr(
            self.mutant_transcript, "mutant_protein_sequence", None)

    def __str__(self):
        return "HaplotypeEffect(%d variants, transcript=%s, source=%s)" % (
            len(self.variants),
            self.transcript.name if self.transcript is not None else "?",
            self.phase_source or "unknown")

    @property
    def short_description(self):
        """HGVS haplotype notation: ``[c.1A>T;c.5G>C]`` for cis.

        Per-variant short descriptions are synthesized from each
        variant's genomic coords when the per-variant effect isn't
        available here — the transcript-level HGVS would require
        re-running the per-variant classifier, which
        ``HaplotypeEffect`` deliberately doesn't do (that work lives
        on the per-variant effects that coexist in the collection).
        """
        parts = [v.short_description for v in self.variants]
        return "[" + ";".join(parts) + "]"


# =====================================================================
# Phase candidate set — possibility set when somatic + germline
# share a window and phase between them is unknown (#268).
#
# Sibling of HaplotypeEffect: that one captures the *known-cis* case
# (multiple variants composed on the same haplotype), this one
# captures the *unknown-phase* case (multiple hypotheses, each with
# its own resulting effect). Both inherit the same MultiOutcomeEffect
# protocol so consumers iterate ``outcomes`` uniformly. This isn't a
# wrapper around a reference-relative effect — it's the primary
# effect class for the unknown-phase case, just like
# SpliceOutcomeSet is for splice ambiguity.
# =====================================================================


class PhaseCandidateSet(TranscriptMutationEffect, MultiOutcomeEffect):
    """Possibility set across phase hypotheses when a somatic variant
    and one or more germline variants share a window on a transcript
    and phase between them is unknown.

    The somatic effect at this locus depends on which haplotype the
    somatic landed on, and without phase data we can't say which.
    The honest output is the set of plausible effects, one per
    hypothesis. ``most_likely`` returns the highest-priority
    candidate; :attr:`candidates` exposes the full set with
    per-hypothesis evidence keys (``phase_state``, ``haplotype``,
    ``germline_variants``) so consumers — including RNA-evidence
    integrations from #259 — can align across axes.

    Emitted by :func:`varcode.germline.predict_germline_aware_effect`.
    Sibling of :class:`HaplotypeEffect`: that one captures the
    *known-cis* multi-variant case; this one captures the
    *unknown-phase* possibility set.
    """

    def __init__(
            self,
            variant,
            transcript,
            candidates,
            hypotheses,
            germline_variants):
        TranscriptMutationEffect.__init__(self, variant, transcript)
        if len(candidates) != len(hypotheses):
            raise ValueError(
                "PhaseCandidateSet needs one candidate per "
                "hypothesis; got %d candidates and %d hypotheses." % (
                    len(candidates), len(hypotheses)))
        if not candidates:
            raise ValueError(
                "PhaseCandidateSet requires at least one candidate.")
        self.candidate_effects = tuple(candidates)
        self.hypotheses = tuple(hypotheses)
        self.germline_variants = tuple(germline_variants)

    def _get_candidates(self) -> tuple[EffectCandidate, ...]:
        """One :class:`EffectCandidate` per hypothesis, sorted by
        underlying effect priority (most-severe first), carrying the
        phase metadata needed to align with RNA-evidence outcomes
        (#259, #382).

        ``evidence`` keys: ``phase_state`` (``"phased"`` /
        ``"implicit"`` / ``"unknown"``),
        ``haplotype`` (opaque tag), ``germline_variants`` (tuple of
        the cis germline variants on that hypothesis's haplotype).
        """
        from ..effect_candidates import EffectCandidate
        from .effect_ordering import effect_priority
        paired = sorted(
            zip(self.candidate_effects, self.hypotheses),
            key=lambda pair: -effect_priority(pair[0]))
        return tuple(
            EffectCandidate(
                effect=candidate,
                source="varcode_germline",
                evidence={
                    "phase_state": hypothesis.phase_state,
                    "haplotype": hypothesis.haplotype,
                    "germline_variants": tuple(hypothesis.cis),
                })
            for candidate, hypothesis in paired)

    # Note on accessor semantics for this class: the base-class
    # "producer order" contract is replaced here by
    # "sorted highest-impact-first" via
    # :func:`effect_priority`. Consequently
    # :attr:`most_likely_candidate` (== ``candidates[0]``) and
    # :attr:`highest_priority_candidate` coincide for this class —
    # both return the most protein-disruptive hypothesis, which is
    # also the conservative single-effect representation downstream
    # consumers see in ``short_description`` etc.

    @property
    def short_description(self):
        """``"?<most_likely_effect>"`` — the leading ``?`` flags the
        ambiguity. Consumers wanting the full possibility set read
        :attr:`candidates`."""
        top = self.most_likely_effect
        assert top is not None
        return "?" + top.short_description

    def _format_short_description(self, deletion_notation):
        top = self.most_likely_effect
        assert top is not None
        return "?" + top.format_short_description(
            deletion_notation=deletion_notation)

    @property
    def priority_class(self):
        """Rank as the most severe hypothesis, like other outcome sets."""
        highest = self.highest_priority_effect
        assert highest is not None
        return getattr(highest, "priority_class", None) or type(highest)

    def _get_mutant_protein_sequence(self) -> str | None:
        """Most-likely candidate's mutant protein. Consumers iterating
        over hypotheses pull per-candidate sequences from
        :attr:`candidates`."""
        return getattr(self.most_likely_effect, "mutant_protein_sequence", None)

    def _get_modifies_protein_sequence(self) -> bool | None:
        return getattr(
            self.most_likely_effect, "modifies_protein_sequence", False)

    def _get_modifies_coding_sequence(self) -> bool | None:
        return getattr(
            self.most_likely_effect, "modifies_coding_sequence", False)

    def __str__(self):
        return (
            "PhaseCandidateSet(variant=%s, transcript=%s, "
            "%d hypotheses, germline=%s)" % (
                self.variant,
                self.transcript.name if self.transcript is not None else "?",
                len(self.hypotheses),
                ",".join(g.short_description for g in self.germline_variants)))

__all__ = [
    "StructuralVariantEffect",
    "LargeDeletion",
    "LargeDuplication",
    "Inversion",
    "GeneFusion",
    "TranslocationToIntergenic",
    "TranslocationToUnannotatedContig",
    "CrypticExonCandidate",
    "HaplotypeEffect",
    "PhaseCandidateSet",
]
