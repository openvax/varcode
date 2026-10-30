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

"""Shared transcript model: hypotheses → products → classification."""

from dataclasses import replace
from functools import lru_cache

from .effect_hypotheses import (
    ClassifiedOutcome,
    EffectHypothesis,
    RealizedTranscriptProduct,
    merge_classified_outcomes,
    order_realized_candidates,
)
from .effects.effect_classes import (
    GermlineAlleleOverlap,
    HypothesisLimit,
    IncompleteTranscript,
    NoncodingTranscript,
    Unresolved,
)
from .effects.effect_ordering import classification_priority
from .errors import HypothesisLimitError, ReferenceMismatchError
from .genome_sequence import _plus_strand_slice, reference_range
from .genomic_layout import GenomicLayout, NonlocalStructuralEdit
from .germline import (
    _check_max_hypotheses, detect_germline_overlap, partition_germline_by_phase)
from .nucleotides import reverse_complement
from .splice_graph import disrupted_splice_sites, enumerate_paired_splice_plans
from .transcript_layout import (
    build_exon_runs,
    classify_products,
    realize_splice_plan,
)
from .version import __version__ as _varcode_version


def _provider_for_genome(genome, transcript):
    """Reuse the selected annotation when DNA is absent; keep DNA priority."""
    has_dna = getattr(genome, "fasta", None) is not None

    @lru_cache(maxsize=1024)
    def provider(contig, start, end):
        sequence = ("" if has_dna or contig != transcript.contig else
                    _plus_strand_slice(transcript, start, end))
        if not sequence:
            sequence = reference_range(genome, contig, start, end)
        if not sequence:
            from .genomic_layout import SequenceUnavailable
            raise SequenceUnavailable(
                "No reference sequence for %s:%d-%d" % (
                    contig, start, end))
        return sequence
    return provider


def _phase_probability(phase_hypotheses, hypothesis):
    if hypothesis.phase_state == "unknown":
        return 1.0 / len(phase_hypotheses)
    return 1.0


def _phase_items(hypothesis):
    result = []
    result.extend((variant, "cis") for variant in hypothesis.cis)
    result.extend((variant, "trans") for variant in hypothesis.trans)
    return tuple(result)


def _status_map(transcript, layout):
    runs = build_exon_runs(transcript, layout)
    statuses = disrupted_splice_sites(transcript, layout, runs)
    return runs, {status.key: status for status in statuses}


def _change_cdna_offset(variant, transcript):
    start = (
        variant.affected_start
        if getattr(variant, "is_structural", False)
        else variant.trimmed_base1_start)
    try:
        return transcript.spliced_offset(start)
    except (KeyError, ValueError):
        return None


def _attach_candidate_set(ordered):
    top = ordered[0]
    top_effect = top.effect
    top_effect.candidates = ordered
    top_effect.realized_candidates = ordered
    top_effect.most_likely_candidate = top
    top_effect.highest_priority_candidate = max(
        ordered, key=lambda candidate: classification_priority(candidate.effect))
    top_effect.probability = top.probability
    top_effect.mechanism_rank = top.ordinal_key
    top_effect.hypotheses = top.hypotheses
    return top_effect


def _unique_alleles(variants):
    """Keep the first occurrence of each normalized allele, in input order."""
    unique = []
    for variant in variants:
        if not detect_germline_overlap(variant, unique):
            unique.append(variant)
    return tuple(unique)


def _refine_transcript_model_effect(effect, observed):
    """Reconcile RNA products into the existing phase/splice outcomes.

    Exact full-product matches support every compatible hypothesis, including
    distinct baselines that produce the same RNA. Partial or unmatched models
    remain explicit observed products with an unresolved comparison; their
    missing reference/phase mapping is never supplied by a DNA prediction.
    RNA support changes ordering, not probabilities or translation evidence.
    """
    from .mutant_transcript import MutantTranscript

    observed = tuple(observed)
    if not observed:
        return effect
    outcomes = [outcome for candidate in effect.realized_candidates
                for outcome in candidate.outcomes]
    changed = False
    for observation in observed:
        model = observation.effect.mutant_transcript
        model_evidence = model.evidence or {} if model is not None else {}
        record = {**model_evidence, **observation.evidence, "source": observation.source,
                  "variant": observation.effect.variant}
        partial = any(status in ("junction_fragment", "retained_reference_fragment")
                      or isinstance(status, str) and status.startswith("partial")
                      for status in (model_evidence.get("sequence_status"),
                                     observation.evidence.get("sequence_status")))
        completed = ("reference_completion" in record
                     or observation.source == "varcode_reference_completion")
        matches = [i for i, outcome in enumerate(outcomes)
                   if model is not None and model.cdna_sequence and not partial and not completed
                   and not outcome.hypothesis.evidence.get("rna_observation")
                   and model.cdna_sequence == outcome.mutant.cdna_sequence]
        if matches:
            for i in matches:
                outcome = outcomes[i]
                support = tuple(outcome.mutant.evidence.get("rna_support", ()))
                if record not in support:
                    mutant = replace(outcome.mutant, evidence={
                        **outcome.mutant.evidence, "rna_support": support + (record,)})
                    outcomes[i] = replace(outcome, mutant=mutant)
                    changed = True
            continue
        if any(outcome.hypothesis.evidence.get("rna_observation") == record
               and outcome.mutant.cdna_sequence == (model.cdna_sequence if model else None)
               for outcome in outcomes):
            continue
        if model is None:
            model = MutantTranscript(reference_transcript=effect.transcript,
                                     reference_segments=(), annotator_name=observation.source)
        mutant = RealizedTranscriptProduct(
            effect.transcript, model.cdna_sequence, model.mutant_protein_sequence,
            evidence={**(model.evidence or {}), "rna_observation": record,
                      "unresolved_mechanisms": ("rna_transcript_comparison",),
                      "unresolved_reason": "Supplied RNA model is not mapped to a supported phase/splice comparison"})
        outcomes.append(ClassifiedOutcome(
            EffectHypothesis(phase_probability=None, evidence={
                "phase_state": "unknown", "rna_observation": record},
                enumeration_index=len(outcomes)), observation.effect,
            RealizedTranscriptProduct(effect.transcript, None, None), mutant))
        changed = True

    if not changed:
        return effect

    # Reclassify through the same path as DNA products. Fresh effects also
    # keep refinement from mutating the original candidate graph.
    refined = []
    for outcome in outcomes:
        result = classify_products(
            effect.variant, effect.transcript, outcome.baseline, outcome.mutant)
        result.variants = getattr(outcome.effect, "variants", getattr(effect, "variants", (effect.variant,)))
        result.splice_signal = getattr(outcome.effect, "splice_signal", None)
        if outcome.hypothesis.evidence.get("rna_observation"):
            # Preserve the producer's partial structure and completeness as-is.
            if outcome.effect.mutant_transcript is not None:
                result.mutant_transcript = outcome.effect.mutant_transcript
        refined.append(replace(outcome, effect=result))
    result = _attach_candidate_set(order_realized_candidates(
        merge_classified_outcomes(refined), classification_priority))
    for name in ("germline_unknown", "phase_source", "annotator", "annotator_version"):
        if name in effect.__dict__:
            setattr(result, name, getattr(effect, name))
    return result


def predict_transcript_model_effect(
        variants, transcript, germline_variants=(), phase_resolver=None,
        sequence_provider=None, max_hypotheses=64, max_phase_hypotheses=None,
        homozygous_germline=()):
    """Predict one ordinary top effect with alternatives in ``.candidates``.

    Canonical and exon-skip paths resolve from transcript annotation alone.
    Intron retention and cryptic sites are also realized when genomic
    sequence is available, and remain explicit :class:`Unresolved`
    candidates at sequence-free tier 0. Rule order is never mislabeled as
    probability.

    More than ``max_phase_hypotheses`` phase hypotheses (default
    ``max_hypotheses``), or more than ``max_hypotheses`` combined
    phase/splice outcomes, gives a :class:`~varcode.HypothesisLimit`
    instead of a partial candidate set. ``homozygous_germline`` lists
    germline variants on both haplotypes, which are always cis.

    ``variants`` is a known-cis group. Members also present in the germline
    anchor the patient baseline; only novel alleles are applied to that
    baseline to obtain the mutant. Repeated normalized alleles are edited
    once. Other germline alleles retain their phase uncertainty.
    """
    max_hypotheses = _check_max_hypotheses(max_hypotheses)
    max_phase_hypotheses = max_hypotheses if max_phase_hypotheses is None else (
        _check_max_hypotheses(max_phase_hypotheses, "max_phase_hypotheses"))
    variants = tuple(variants)
    if not variants:
        raise ValueError("predict_transcript_model_effect requires a somatic variant")
    germline_variants = _unique_alleles(germline_variants)
    somatic_variants = tuple(
        v for v in _unique_alleles(variants)
        if not detect_germline_overlap(v, germline_variants))
    if not somatic_variants:
        result = GermlineAlleleOverlap(variants[0], transcript)
        result.variants = variants
        return result
    primary = somatic_variants[0]
    sv_type = getattr(primary, "sv_type", None)
    assembled = sv_type in ("DEL", "DUP", "INV") and primary.alt_assembly
    nonlocal_span = sv_type in ("DEL", "DUP", "INV") and (
        primary.affected_start < transcript.start
        or primary.affected_end > transcript.end)
    junction_record = sv_type is not None and any(
        bracket in getattr(primary, "symbolic_alt", "") for bracket in "[]")
    if sv_type in ("BND", "CNV", "INV") or assembled or nonlocal_span or junction_record:
        if len(variants) != 1 or germline_variants:
            return Unresolved(
                primary, transcript, mechanism="structural_haplotype",
                reason=("Composition across a supplied assembly or nonlocal junction "
                        "requires mapped patient-baseline coordinates"))
        # One established builder owns fusion partners, junction inserts and
        # supplied transcript assemblies. A finite single-transcript layout
        # must not replace those products with a clipped local span.
        from .effects.structural import predict_structural_variant_effect
        return predict_structural_variant_effect(primary, transcript)
    if not transcript.is_protein_coding:
        return NoncodingTranscript(primary, transcript)
    if not transcript.complete:
        return IncompleteTranscript(primary, transcript)
    if sequence_provider is None:
        sequence_provider = _provider_for_genome(primary.genome, transcript)

    phase = partition_germline_by_phase(
        primary, germline_variants, phase_resolver,
        homozygous=homozygous_germline, cis_variants=variants)
    try:
        outcomes = _realize_outcomes(
            variants, somatic_variants, transcript,
            phase.hypotheses(max_phase_hypotheses),
            sequence_provider, max_hypotheses,
            phase_source=getattr(phase_resolver, "phase_source", None))
    except ReferenceMismatchError as error:
        expected, observed = error.genomic_expected_ref, error.genomic_observed_ref
        if transcript.on_backward_strand:
            expected, observed = reverse_complement(expected), reverse_complement(observed)
        try:
            offset = transcript.spliced_offset(error.genome_start)
        except (KeyError, ValueError):
            offset = None
        raise ReferenceMismatchError(
            error.variant, transcript, expected, observed,
            transcript_offset=offset, genome_start=error.genome_start,
            genome_end=error.genome_end) from error
    except NonlocalStructuralEdit as error:
        return Unresolved(
            primary, transcript, mechanism="nonlocal_structural_variant",
            reason=str(error))
    except HypothesisLimitError as error:
        result = HypothesisLimit(
            primary, transcript, error.max_hypotheses, phase,
            _reference_effect(variants, transcript, germline_variants,
                              sequence_provider, max_hypotheses))
        result.variants = variants
        return result
    candidates = merge_classified_outcomes(outcomes)
    ordered = order_realized_candidates(candidates, classification_priority)
    return _attach_candidate_set(ordered)


def _realize_outcomes(
        variants, somatic_variants, transcript, phase_hypotheses,
        sequence_provider, max_hypotheses, phase_source=None):
    """Classify every phase hypothesis × splice plan.

    Raises :class:`~varcode.HypothesisLimitError` past ``max_hypotheses``
    outcomes and ``NonlocalStructuralEdit`` for edits the layout cannot hold.
    """
    primary = somatic_variants[0]
    outcomes = []
    enumeration_index = 0
    # Complete nucleotide alleles may cross a transcript boundary by more
    # than the splice scan flank. Keep their full genomic span before exon
    # projection; extending the layout does not extend the transcript.
    flank = 50
    for variant in (*somatic_variants, *(v for h in phase_hypotheses for v in h.cis)):
        if not getattr(variant, "is_structural", False):
            flank = max(flank, transcript.start - variant.trimmed_base1_start + 50,
                        variant.trimmed_base1_end - transcript.end + 50)
    reference = GenomicLayout.from_transcript(
        transcript, flank=flank, sequence_provider=sequence_provider)
    for phase_hypothesis in phase_hypotheses:
        baseline_layout = reference.apply_variants(phase_hypothesis.cis)
        mutant_layout = baseline_layout.apply_variants(somatic_variants)
        baseline_runs, baseline_statuses = _status_map(
            transcript, baseline_layout)
        mutant_runs, mutant_statuses = _status_map(transcript, mutant_layout)
        shared_keys = frozenset(
            key for key in baseline_statuses.keys() & mutant_statuses.keys()
            if baseline_statuses[key] == mutant_statuses[key])
        try:
            plans = enumerate_paired_splice_plans(
                transcript, baseline_statuses.values(), mutant_statuses.values(),
                tuple(run.key for run in baseline_runs),
                tuple(run.key for run in mutant_runs), shared_keys,
                max_plans=max_hypotheses - len(outcomes))
        except HypothesisLimitError as error:
            raise HypothesisLimitError(max_hypotheses) from error

        for baseline_plan, plan in plans:
            baseline_product, _ = realize_splice_plan(
                transcript, baseline_layout, baseline_plan,
                tuple(baseline_statuses.values()))
            mutant_product, _ = realize_splice_plan(
                transcript, mutant_layout, plan,
                tuple(mutant_statuses.values()))
            hypothesis = EffectHypothesis(
                phase=_phase_items(phase_hypothesis),
                splice_plan=plan,
                baseline_splice_plan=baseline_plan,
                phase_probability=_phase_probability(
                    phase_hypotheses, phase_hypothesis),
                evidence={
                    "phase_state": phase_hypothesis.phase_state,
                    "phase_source": phase_source,
                    "somatic_variants": somatic_variants,
                    "shared_splice_sites": tuple(sorted(shared_keys)),
                },
                enumeration_index=enumeration_index)
            result = classify_products(
                primary, transcript, baseline_product, mutant_product)
            result.variants = variants
            result.splice_signal = tuple(sorted(mutant_statuses)) or None
            outcome = ClassifiedOutcome(
                hypothesis=hypothesis,
                effect=result,
                baseline=baseline_product,
                mutant=mutant_product,
                change_aa_offset=getattr(
                    result, "aa_mutation_start_offset", None),
                change_cdna_offset=_change_cdna_offset(primary, transcript))
            outcomes.append(outcome)
            enumeration_index += 1
            if len(outcomes) > max_hypotheses:
                raise HypothesisLimitError(max_hypotheses)
    return outcomes


def _reference_effect(
        variants, transcript, germline_variants, sequence_provider,
        max_hypotheses):
    """The prediction without germline context, as ``effects()`` returns it,
    or ``None`` when that prediction is itself over the limit."""
    from .splice_outcomes import enumerate_splice_outcomes
    if not germline_variants:
        return None
    reference = predict_transcript_model_effect(
        variants, transcript, sequence_provider=sequence_provider,
        max_hypotheses=max_hypotheses)
    if isinstance(reference, HypothesisLimit):
        return None
    return enumerate_splice_outcomes(reference)


class TranscriptModelEffectAnnotator:
    """Default annotator backed by the shared transcript model."""

    name = "transcript_model"
    version = _varcode_version

    def annotate_on_transcript(self, variant, transcript):
        return self._predict(variant, transcript)

    @staticmethod
    def _predict(variant, transcript, **kwargs):
        return TranscriptModelEffectAnnotator._predict_variants(
            (variant,), transcript, **kwargs)

    @staticmethod
    def _predict_variants(variants, transcript, **kwargs):
        from .genomic_layout import SequenceUnavailable, UnsupportedLayoutEdit
        try:
            return predict_transcript_model_effect(variants, transcript, **kwargs)
        except UnsupportedLayoutEdit:
            return NotImplemented
        except SequenceUnavailable as error:
            return Unresolved(
                variants[0], transcript, mechanism="sequence_unavailable",
                reason=str(error))

    def annotate_with_context(
            self, variant, transcript, germline_ctx, phase_resolver=None):
        """Annotate through the same pipeline with patient germline edits."""
        return self._annotate_with_context(
            (variant,), transcript, germline_ctx, phase_resolver)

    def _annotate_with_context(
            self, variants, transcript, germline_ctx, phase_resolver):
        from .germline import Completeness

        # A full patient product must include supplied alleles across the
        # transcript, including sequence that a retention path may expose.
        germline = germline_ctx.variants_in_window(
            transcript.contig, transcript.start, transcript.end)
        result = self._predict_variants(
            variants, transcript,
            germline_variants=germline,
            phase_resolver=phase_resolver,
            max_phase_hypotheses=germline_ctx.max_phase_hypotheses,
            homozygous_germline=germline_ctx.homozygous(germline))
        if result is NotImplemented:
            return result
        if germline_ctx.completeness in (
                Completeness.SPARSE, Completeness.HOTSPOTS_ONLY):
            result.germline_unknown = True
        return result

    def annotate_haplotype(
            self, variants, transcript, germline_ctx=None, phase_resolver=None):
        """Realize the whole cis group with phase/splice alternatives.

        Reconcile supplied RNA models with the same phase/splice products;
        a per-variant observation need not cover every edit in the group.
        """
        if germline_ctx:
            result = self._annotate_with_context(
                variants, transcript, germline_ctx, phase_resolver)
        else:
            result = self._predict_variants(variants, transcript)
        if result is not NotImplemented:
            from .phasing import apply_phase_resolver_to_effects
            result = apply_phase_resolver_to_effects([result], phase_resolver)[0]
        return result
