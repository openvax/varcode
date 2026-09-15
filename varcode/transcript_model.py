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

"""Experimental transcript model: hypotheses → products → classification."""

from .effect_hypotheses import (
    ClassifiedOutcome,
    EffectHypothesis,
    SplicePlan,
    enumerate_splice_plans,
    merge_classified_outcomes,
    order_realized_candidates,
)
from .effects.effect_classes import (
    IncompleteTranscript,
    NoncodingTranscript,
    Unresolved,
)
from .effects.effect_ordering import effect_priority
from .genome_sequence import reference_range
from .genomic_layout import GenomicLayout
from .germline import enumerate_phase_hypotheses
from .splice_graph import disrupted_splice_sites, splice_axes
from .transcript_layout import (
    build_exon_runs,
    classify_products,
    realize_splice_plan,
)
from .version import __version__ as _varcode_version


def _provider_for_genome(genome):
    def provider(contig, start, end):
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
    if hypothesis.phase_state == "too_many_hypotheses":
        return None
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


def _combined_mutant_axes(
        baseline_statuses, mutant_statuses, transcript,
        baseline_run_keys, mutant_run_keys):
    """Return mutant axes, shared keys, and fixed baseline mechanisms."""
    axes = []
    shared = set()
    baseline_fixed = {}
    all_keys = sorted(set(baseline_statuses) | set(mutant_statuses))
    for key in all_keys:
        baseline = baseline_statuses.get(key)
        mutant = mutant_statuses.get(key)
        if (baseline is not None and mutant is not None
                and baseline.strength == mutant.strength
                and baseline.reason == mutant.reason):
            if mutant.strength == "lost":
                axes.extend(splice_axes(
                    transcript, (mutant,), mutant_run_keys))
                shared.add(key)
            # Weak germline-only changes remain canonical without a scorer.
            continue
        if mutant is not None:
            axes.extend(splice_axes(
                transcript, (mutant,), mutant_run_keys))
        if baseline is not None and baseline.strength == "lost":
            baseline_fixed[key] = splice_axes(
                transcript, (baseline,), baseline_run_keys)[0].options[0]
    return tuple(axes), frozenset(shared), baseline_fixed


def _selected_mechanisms(plan):
    return {key: option for key, option in plan.choices}


def _baseline_plan(
        transcript, baseline_runs, baseline_statuses, mutant_plan,
        shared_keys, baseline_fixed):
    """Map shared mutant choices and baseline-only disruptions to baseline."""
    kept = {run.key for run in baseline_runs}
    selected = _selected_mechanisms(mutant_plan)
    choices = []
    for key in shared_keys:
        option = selected.get(key)
        if option is None:
            # The corresponding mutant site disappeared behind another graph
            # rewrite. Baseline still has it, so use its declared first choice.
            status = baseline_statuses[key]
            option = _option_for_mechanism(
                status, selected=None, transcript=transcript,
                run_keys=tuple(run.key for run in baseline_runs))
        else:
            option = _option_for_mechanism(
                baseline_statuses[key], option.mechanism, transcript,
                tuple(run.key for run in baseline_runs))
        choices.append((key, option))
        if option.mechanism == "exon_skip":
            kept.discard(baseline_statuses[key].run_key)
    for key, option in baseline_fixed.items():
        choices.append((key, option))
        if option.mechanism == "exon_skip":
            kept.discard(baseline_statuses[key].run_key)
    return SplicePlan(
        choices=tuple(sorted(choices, key=lambda choice: choice[0])),
        kept_runs=tuple(run.key for run in baseline_runs if run.key in kept))


def _option_for_mechanism(status, selected, transcript, run_keys):
    """Return the first option, or the option matching ``selected``."""
    axis = splice_axes(transcript, (status,), run_keys)[0]
    if selected is None:
        return axis.options[0]
    return next(
        option for option in axis.options if option.mechanism == selected)


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
        ordered, key=lambda candidate: effect_priority(candidate.effect))
    top_effect.probability = top.probability
    top_effect.mechanism_rank = top.ordinal_key
    top_effect.hypotheses = top.hypotheses
    return top_effect


def predict_transcript_model_effect(
        variants, transcript, germline_variants=(), phase_resolver=None,
        sequence_provider=None, max_hypotheses=64):
    """Predict one ordinary top effect with alternatives in ``.candidates``.

    Canonical and exon-skip paths resolve from transcript annotation alone.
    Intron retention and cryptic sites are also realized when genomic
    sequence is available, and remain explicit :class:`Unresolved`
    candidates at sequence-free tier 0. Rule order is never mislabeled as
    probability.
    """
    variants = tuple(variants)
    if not variants:
        raise ValueError("predict_transcript_model_effect requires a somatic variant")
    germline_variants = tuple(germline_variants)
    primary = variants[0]
    if getattr(primary, "sv_type", None) == "BND":
        if len(variants) != 1 or germline_variants:
            return Unresolved(
                primary, transcript, mechanism="breakend_haplotype",
                reason=(
                    "BND composition with additional phased variants "
                    "requires an assembled allele"))
        # The established fusion builder already resolves partner transcript,
        # orientation, fused cDNA and translation from this same BND. Keep it
        # as the BND realization path while the layout engine handles local
        # span variants.
        from .annotators.structural_variant import StructuralVariantAnnotator
        return StructuralVariantAnnotator().annotate_on_transcript(
            primary, transcript)
    if not transcript.is_protein_coding:
        return NoncodingTranscript(primary, transcript)
    if not transcript.complete:
        return IncompleteTranscript(primary, transcript)
    if sequence_provider is None:
        sequence_provider = _provider_for_genome(primary.genome)

    phase_hypotheses = enumerate_phase_hypotheses(
        primary, germline_variants, phase_resolver=phase_resolver,
        max_hypotheses=max_hypotheses)
    outcomes = []
    enumeration_index = 0
    for phase_hypothesis in phase_hypotheses:
        reference = GenomicLayout.from_transcript(
            transcript, flank=50, sequence_provider=sequence_provider)
        baseline_layout = reference.apply_variants(phase_hypothesis.cis)
        mutant_layout = baseline_layout.apply_variants(variants)
        baseline_runs, baseline_statuses = _status_map(
            transcript, baseline_layout)
        mutant_runs, mutant_statuses = _status_map(transcript, mutant_layout)
        axes, shared_keys, baseline_fixed = _combined_mutant_axes(
            baseline_statuses, mutant_statuses, transcript,
            tuple(run.key for run in baseline_runs),
            tuple(run.key for run in mutant_runs))
        if axes:
            remaining = max_hypotheses - len(outcomes)
            if remaining < 1:
                raise ValueError(
                    "Combined phase/splice hypothesis count exceeds "
                    "max_hypotheses=%d" % max_hypotheses)
            try:
                plans = enumerate_splice_plans(
                    axes, tuple(run.key for run in mutant_runs),
                    max_plans=remaining)
            except ValueError as error:
                if "exceeds max_plans" not in str(error):
                    raise
                raise ValueError(
                    "Combined phase/splice hypothesis count exceeds "
                    "max_hypotheses=%d" % max_hypotheses) from error
        else:
            plans = (SplicePlan(
                choices=(),
                kept_runs=tuple(run.key for run in mutant_runs)),)

        for plan in plans:
            baseline_plan = _baseline_plan(
                transcript, baseline_runs, baseline_statuses, plan,
                shared_keys, baseline_fixed)
            baseline_product, baseline_unresolved = realize_splice_plan(
                transcript, baseline_layout, baseline_plan,
                tuple(baseline_statuses.values()))
            mutant_product, mutant_unresolved = realize_splice_plan(
                transcript, mutant_layout, plan,
                tuple(mutant_statuses.values()))
            hypothesis = EffectHypothesis(
                phase=_phase_items(phase_hypothesis),
                splice_plan=plan,
                phase_probability=_phase_probability(
                    phase_hypotheses, phase_hypothesis),
                evidence={
                    "phase_state": phase_hypothesis.phase_state,
                    "shared_splice_sites": tuple(sorted(shared_keys)),
                },
                enumeration_index=enumeration_index)
            unresolved = mutant_unresolved + baseline_unresolved
            if unresolved:
                mechanism = unresolved[0]
                result = Unresolved(
                    primary, transcript, mechanism=mechanism,
                    reason="genomic splice realization not available")
            else:
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
                raise ValueError(
                    "Combined phase/splice hypothesis count exceeds "
                    "max_hypotheses=%d" % max_hypotheses)

    candidates = merge_classified_outcomes(outcomes)
    ordered = order_realized_candidates(candidates, effect_priority)
    return _attach_candidate_set(ordered)


class TranscriptModelEffectAnnotator:
    """Experimental, opt-in annotator backed by the transcript model."""

    name = "transcript_model"
    version = _varcode_version

    def annotate_on_transcript(self, variant, transcript):
        return self._predict(variant, transcript)

    @staticmethod
    def _predict(variant, transcript, **kwargs):
        from .genomic_layout import SequenceUnavailable, UnsupportedLayoutEdit
        try:
            return predict_transcript_model_effect((variant,), transcript, **kwargs)
        except UnsupportedLayoutEdit:
            return NotImplemented
        except SequenceUnavailable as error:
            return Unresolved(
                variant, transcript, mechanism="sequence_unavailable",
                reason=str(error))

    def annotate_with_context(
            self, variant, transcript, germline_ctx, phase_resolver=None):
        """Annotate through the same pipeline with patient germline edits."""
        from .germline import Completeness, detect_loh

        if getattr(variant, "is_structural", False):
            start = variant.affected_start
            end = variant.affected_end
        else:
            start = variant.trimmed_base1_start
            end = variant.trimmed_base1_end
        start = max(transcript.start, start - 90)
        end = min(transcript.end, end + 90)
        germline = tuple(germline_ctx.variants_in_window(
            variant.contig, start, end))
        result = self._predict(
            variant, transcript,
            germline_variants=germline,
            phase_resolver=phase_resolver)
        if result is NotImplemented:
            return result
        if (not germline and germline_ctx.completeness in (
                Completeness.SPARSE, Completeness.HOTSPOTS_ONLY)):
            result.germline_unknown = True
        if detect_loh(variant, germline):
            result.is_loh = True
        return result
