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


"""Germline phase partitioning, hypothesis enumeration, and baseline classification."""

from __future__ import annotations

from dataclasses import dataclass
import logging
from typing import Optional, Tuple

from serializable import DataclassSerializable

from .effects.transcript_helpers import has_complete_coding_sequence
from .errors import HypothesisLimitError
from .germline import (
    Completeness,
    DEFAULT_MAX_PHASE_HYPOTHESES,
    GermlineContext,
    check_max_hypotheses,
    default_germline_window,
    detect_germline_overlap,
)
from .phasing import query_in_cis

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class PhaseHypothesis:
    """One way the somatic variant might be phased relative to the
    germline variants in its window.

    ``cis`` are germline variants on the same haplotype as the
    somatic; ``trans`` are on the other haplotype. The somatic's
    effect at this locus is determined entirely by the cis set —
    edits in trans live on the haplotype the somatic *didn't* hit and
    don't change the codon the somatic edits.

    ``haplotype`` is a stable opaque label (``"A"``, ``"B"``, etc.)
    that callers use to align outcomes across axes — e.g. an RNA
    resolver returning ``evidence={"haplotype": "B"}`` aligns with the
    hypothesis whose ``haplotype == "B"``.

    ``phase_state`` mirrors the table in #268: ``"phased"``,
    ``"implicit"`` (hemizygous), or ``"unknown"`` (enumerated).
    """
    cis: Tuple = ()
    trans: Tuple = ()
    haplotype: str = "unknown"
    phase_state: str = "unknown"




@dataclass(frozen=True)
class PhasePartition(DataclassSerializable):
    """Germline variants near a somatic variant, split by what is known
    about their phase relative to it.

    ``germline`` holds every nearby germline variant in input order;
    ``cis`` and ``trans`` are those whose phase relative to the somatic
    variant is known, and :attr:`unphased` is the rest. ``linked`` lists
    groups of unphased variants whose phase relative to each other is
    known, as ``(one_haplotype, other_haplotype)`` pairs with the group's
    first variant on ``one_haplotype``; each group changes sides as a unit. ``implicit`` is true on hemizygous contigs
    (chrM, chrY), where every germline variant is cis without evidence.
    Build one with :func:`partition_germline_by_phase`.
    """
    germline: Tuple = ()
    cis: Tuple = ()
    trans: Tuple = ()
    linked: Tuple = ()
    implicit: bool = False

    @property
    def unphased(self) -> Tuple:
        return tuple(
            v for v in self.germline if v not in self.cis and v not in self.trans)

    @property
    def blocks(self) -> Tuple:
        """Units of unknown phase as ``(one_haplotype, other_haplotype)``
        pairs: each ``linked`` group, then each other unphased variant."""
        grouped = tuple(v for block in self.linked for side in block for v in side)
        return self.linked + tuple(
            ((v,), ()) for v in self.unphased if v not in grouped)

    @property
    def hypothesis_count(self) -> int:
        """One hypothesis per assignment of each block to a haplotype."""
        return 2 ** len(self.blocks)

    def hypotheses(
            self, max_hypotheses: int = DEFAULT_MAX_PHASE_HYPOTHESES,
    ) -> Tuple[PhaseHypothesis, ...]:
        """Enumerate phase hypotheses, keeping known phase fixed in each.

        Known phase gives one hypothesis (``phase_state`` ``"phased"``,
        or ``"implicit"`` on hemizygous contigs). Otherwise there is one
        ``"unknown"`` hypothesis per assignment of the :attr:`blocks`.
        Haplotype labels depend only on the cis set: ``"A"`` when every
        germline variant is cis, ``"B"`` when none is, otherwise
        ``"A_mixed_<mask>"`` over :attr:`germline`.

        Raises
        ------
        HypothesisLimitError
            When :attr:`hypothesis_count` exceeds ``max_hypotheses``.
        """
        max_hypotheses = check_max_hypotheses(max_hypotheses)
        blocks = self.blocks
        if not blocks:
            return (self._hypothesis(
                self.cis, "implicit" if self.implicit else "phased"),)
        if self.hypothesis_count > max_hypotheses:
            raise HypothesisLimitError(max_hypotheses)
        return tuple(
            self._hypothesis(self.cis + tuple(
                v for i, (one, other) in enumerate(blocks)
                for v in (one if mask >> i & 1 else other)), "unknown")
            for mask in range(self.hypothesis_count))

    def _hypothesis(self, cis, phase_state):
        in_cis = [v in cis for v in self.germline]
        if all(in_cis):
            haplotype = "A"
        elif not any(in_cis):
            haplotype = "B"
        else:
            haplotype = "A_mixed_%d" % sum(
                1 << i for i, c in enumerate(in_cis) if c)
        return PhaseHypothesis(
            cis=tuple(v for v, c in zip(self.germline, in_cis) if c),
            trans=tuple(v for v, c in zip(self.germline, in_cis) if not c),
            haplotype=haplotype,
            phase_state=phase_state)


def partition_germline_by_phase(
        somatic_variant, germline_variants, phase_resolver=None,
        homozygous=(), cis_variants=()) -> PhasePartition:
    """Split ``germline_variants`` by their phase relative to ``somatic_variant``.

    Variants in ``homozygous`` sit on both haplotypes and are always cis,
    as is every germline variant on a hemizygous contig (chrM, chrY; chrX
    is treated as diploid, since Varcode has no sex information). The
    rest are placed by :func:`~varcode.phasing.query_in_cis` answers
    between every pair of variants, so phase can follow through another
    germline variant. Variants whose phase is known only relative to each
    other become ``linked`` groups. An answer that contradicts earlier
    ones is ignored with a warning.

    ``cis_variants`` supplies additional variants known to be in cis with
    the somatic variant. They anchor the phase graph before resolver answers,
    so germline phase can follow through any member of that group.
    Homozygous alleles do not link the two haplotypes through this graph.
    """
    germline = tuple(germline_variants)
    homozygous = tuple(homozygous)
    if str(somatic_variant.contig).lstrip("chr").upper() in ("M", "MT", "Y"):
        return PhasePartition(germline=germline, cis=germline, implicit=True)
    rest = [g for g in germline if g not in homozygous]
    # Union-find over the somatic variant (node 0), ``rest``, and any
    # additional cis-group members. Parity is 0 for the same haplotype
    # as the parent node and 1 for the other one.
    nodes = [somatic_variant] + rest
    cis_variants = tuple(cis_variants)
    for variant in cis_variants:
        if not detect_germline_overlap(variant, nodes + list(homozygous)):
            nodes.append(variant)
    parent = list(range(len(nodes)))
    parity = [0] * len(nodes)
    for i, variant in enumerate(nodes):
        if detect_germline_overlap(variant, cis_variants):
            parent[i] = 0

    def find(i):
        if parent[i] != i:
            root, root_parity = find(parent[i])
            parent[i], parity[i] = root, parity[i] ^ root_parity
        return parent[i], parity[i]

    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            answer = query_in_cis(phase_resolver, nodes[i], nodes[j])
            if answer is None:
                continue
            (root_i, parity_i), (root_j, parity_j) = find(i), find(j)
            relation = 0 if answer else 1
            if root_i != root_j:
                parent[root_i] = root_j
                parity[root_i] = parity_i ^ parity_j ^ relation
            elif parity_i ^ parity_j != relation:
                logger.warning(
                    "Ignoring phase answer for %s and %s that contradicts "
                    "earlier answers", nodes[i], nodes[j])
    somatic_root, somatic_parity = find(0)
    cis, trans, groups = [], [], {}
    for i, g in enumerate(rest, 1):
        root, node_parity = find(i)
        if root != somatic_root:
            groups.setdefault(root, []).append((g, node_parity))
        elif node_parity == somatic_parity:
            cis.append(g)
        else:
            trans.append(g)
    return PhasePartition(
        germline=germline,
        cis=tuple(g for g in germline if g in homozygous or g in cis),
        trans=tuple(trans),
        linked=tuple(
            (tuple(g for g, p in members if p == members[0][1]),
             tuple(g for g, p in members if p != members[0][1]))
            for members in groups.values() if len(members) > 1))


def enumerate_phase_hypotheses(
        somatic_variant,
        germline_in_window,
        phase_resolver=None,
        max_hypotheses: int = DEFAULT_MAX_PHASE_HYPOTHESES,
) -> Tuple[PhaseHypothesis, ...]:
    """Enumerate the phase configurations of ``somatic_variant``
    relative to ``germline_in_window``.

    Shorthand for :func:`partition_germline_by_phase` followed by
    :meth:`PhasePartition.hypotheses`. Raises
    :class:`~varcode.HypothesisLimitError` when the unphased variants
    need more than ``max_hypotheses`` hypotheses.
    """
    return partition_germline_by_phase(
        somatic_variant, germline_in_window, phase_resolver
    ).hypotheses(max_hypotheses)


def classify_against_patient_baseline(
        somatic_variant, transcript, hypothesis: PhaseHypothesis):
    """Build the patient-baseline and post-somatic haplotypes for a
    single phase hypothesis, then classify the somatic effect via
    the protein-diff classifier.

    Returns a :class:`MutationEffect` instance. The effect's class
    (Missense, FrameShift, etc.) reflects the diff against the
    patient's baseline — not the reference's — when ``hypothesis.cis``
    is non-empty. If those edits cannot be composed or the required proteins
    are unavailable, return ``Unresolved`` rather than discard the cis edits.
    Its ``reference_effect`` is diagnostic reference-only evidence, not a
    patient-relative prediction.

    For the all-trans case (``hypothesis.cis`` empty), the patient
    baseline is identical to the reference, so the diff degenerates
    to today's reference-relative classification. That's the
    structural back-compat: when no germline lands on the somatic's
    haplotype, this branch produces what the protein-diff annotator
    would produce on its own.
    """
    from .effects.classify import classify_from_protein_diff
    from .effects.codon_tables import codon_table_for_transcript, initiator_as_met
    from .mutant_transcript import apply_variants_to_transcript

    cis_germline = list(hypothesis.cis)
    baseline = None

    def unavailable(reason):
        from .effects import Unresolved, predict_variant_effect_on_transcript
        reference = predict_variant_effect_on_transcript(
            somatic_variant, transcript, annotator="fast")
        # With no cis edits there is no patient/reference distinction to lose.
        if not cis_germline:
            return reference
        return Unresolved(
            somatic_variant, transcript, mechanism="germline_composition",
            reason=reason, reference_effect=reference)

    # Patient baseline: just the cis germline applied. Empty cis →
    # baseline is the reference transcript (no edits).
    if cis_germline:
        baseline = apply_variants_to_transcript(cis_germline, transcript)
        if baseline is None:
            return unavailable(
                "Cis germline edits could not be applied to the reference transcript")
        baseline_protein = baseline.mutant_protein_sequence
        baseline_length_delta = baseline.total_length_delta
    else:
        # No cis germline → patient haplotype == reference at this
        # locus. Normalize annotated alternate initiators just as the
        # mutant builder does for germline and joint haplotypes (#634).
        baseline_protein = str(transcript.protein_sequence) if (
            getattr(transcript, "protein_sequence", None)) else None
        cds_start = min(transcript.start_codon_spliced_offsets)
        baseline_protein = initiator_as_met(
            baseline_protein, str(transcript.sequence)[cds_start:cds_start + 3],
            codon_table_for_transcript(transcript))
        baseline_length_delta = 0

    # Post-somatic: cis germline + somatic on the same haplotype.
    joint_variants = cis_germline + [somatic_variant]
    joint = apply_variants_to_transcript(joint_variants, transcript)
    if joint is None:
        return unavailable(
            "Cis germline and somatic edits could not be applied together to the reference transcript")
    if joint.mutant_protein_sequence is None:
        return unavailable(
            "Protein sequence is unavailable for the joint germline/somatic transcript")

    if baseline_protein is None:
        return unavailable(
            "Protein sequence is unavailable for the patient germline baseline")

    # Length delta of the somatic alone, after stripping the germline-
    # only baseline shift. classify_from_protein_diff expects the diff
    # delta, not the joint delta.
    somatic_length_delta = joint.total_length_delta - baseline_length_delta

    return classify_from_protein_diff(
        variant=somatic_variant,
        transcript=transcript,
        ref_protein=baseline_protein,
        mut_protein=joint.mutant_protein_sequence,
        length_delta=somatic_length_delta,
        mutant_transcript=joint,
        baseline_transcript=baseline)


def unresolved_hypothesis_limit(
    variant,
    transcript,
    max_hypotheses: int,
    phase: PhasePartition,
    reference_effect=None,
):
    """Build an ``Unresolved`` effect when outcomes exceed ``max_hypotheses``."""
    from .effects.effect_classes import Unresolved

    reason = "More than %d alternative outcomes" % max_hypotheses
    if phase.unphased:
        reason += " (%d unphased germline variants)" % len(phase.unphased)
    return Unresolved(
        variant,
        transcript,
        mechanism="hypothesis_limit",
        reason=reason,
        reference_effect=reference_effect,
        details={"max_hypotheses": max_hypotheses, "phase": phase},
    )


def predict_germline_aware_effect(
        somatic_variant,
        transcript,
        germline_ctx: GermlineContext,
        annotator,
        phase_resolver=None,
        window_fn=default_germline_window,
        max_hypotheses: Optional[int] = None):
    """Predict the effect of ``somatic_variant`` on ``transcript``
    against the patient's germline-applied transcript.

    Single entry point for germline-aware effect prediction.
    :func:`varcode.effects.predict_variant_effects` calls this whenever
    a non-empty :class:`GermlineContext` is supplied; otherwise it
    bypasses the germline path entirely and the existing annotator
    dispatch produces today's reference-relative output unchanged.

    Behaviour by case:

    * **Somatic allele matches a germline allele** — a
      ``GermlineAlleleOverlap``, without inferring zygosity, LOH, or a
      new allele sequence.
    * **Phase cannot matter** (no germline in the somatic's window, or a
      noncoding or incompletely annotated transcript) — delegate to
      ``annotator``. With no germline nearby, SPARSE / HOTSPOTS_ONLY
      contexts mark the result ``effect.germline_unknown = True``.
    * **Phase known** (resolver answers, homozygous germline, or a
      hemizygous contig) — classify
      against the single patient haplotype via
      :func:`classify_against_patient_baseline`.
    * **Phase unknown** — classify each hypothesis and wrap them in a
      :class:`~varcode.effects.effect_classes.PhaseCandidateSet`.
    * **Cis composition unavailable** — keep an ``Unresolved`` result for
      that hypothesis, with its reference-only annotation in ``reference_effect``.
    * **Too many hypotheses** — an :class:`~varcode.Unresolved` effect with
      ``mechanism="hypothesis_limit"`` recording the phase partition and cap in
      ``details`` plus the reference-relative effect in ``reference_effect``.

    ``max_hypotheses`` overrides ``germline_ctx.max_phase_hypotheses``.
    ``window_fn`` is the pluggable window selector — defaults to
    :func:`default_germline_window` (codon-level, with splice-signal
    expansion when the somatic is splice-adjacent). Callers that
    need different windows pass their own.
    """
    from pyensembl import Transcript
    from .effects.effect_classes import (
        GermlineAlleleOverlap,
        PhaseCandidateSet,
    )

    limit = check_max_hypotheses(
        germline_ctx.max_phase_hypotheses if max_hypotheses is None
        else max_hypotheses)
    if not isinstance(transcript, Transcript):
        # Mirrors annotator entry: SVs and other non-Transcript
        # consumers don't go through germline-aware prediction.
        return annotator.annotate_on_transcript(somatic_variant, transcript)

    contig, start, end = window_fn(somatic_variant, transcript)
    germline_in_window = germline_ctx.variants_in_window(contig, start, end)
    if detect_germline_overlap(somatic_variant, germline_in_window):
        return GermlineAlleleOverlap(somatic_variant, transcript)

    if not (germline_in_window and transcript.is_protein_coding
            and has_complete_coding_sequence(transcript)):
        effect = annotator.annotate_on_transcript(somatic_variant, transcript)
        if effect is NotImplemented:
            return effect
        if not germline_in_window and germline_ctx.completeness in (
                Completeness.SPARSE, Completeness.HOTSPOTS_ONLY):
            effect.germline_unknown = True
        return effect

    phase = partition_germline_by_phase(
        somatic_variant, germline_in_window, phase_resolver,
        homozygous=germline_ctx.homozygous(germline_in_window))
    try:
        hypotheses = phase.hypotheses(limit)
    except HypothesisLimitError:
        from .effects import predict_variant_effect_on_transcript
        return unresolved_hypothesis_limit(
            somatic_variant, transcript, limit, phase,
            predict_variant_effect_on_transcript(
                somatic_variant, transcript, annotator=annotator))

    if len(hypotheses) == 1:
        effect = classify_against_patient_baseline(
            somatic_variant, transcript, hypotheses[0])
        # Stash the phase metadata on the effect so consumers /
        # serializers can recover it. Single-hypothesis effects don't
        # need a possibility set, but the evidence is still useful.
        effect.germline_phase_state = hypotheses[0].phase_state
        effect.germline_variants_in_window = tuple(germline_in_window)
        return effect

    return PhaseCandidateSet(
        variant=somatic_variant,
        transcript=transcript,
        candidates=tuple(
            classify_against_patient_baseline(somatic_variant, transcript, h)
            for h in hypotheses),
        hypotheses=hypotheses,
        germline_variants=tuple(germline_in_window))


__all__ = [
    "PhaseHypothesis",
    "PhasePartition",
    "classify_against_patient_baseline",
    "enumerate_phase_hypotheses",
    "partition_germline_by_phase",
    "predict_germline_aware_effect",
]
