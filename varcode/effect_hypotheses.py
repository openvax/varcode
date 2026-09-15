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

"""Composable hypothesis and realized-product data types.

These types deliberately keep three concepts separate:

* a hypothesis records uncertain biological choices;
* a transcript product records what sequence those choices produce;
* a classified outcome compares a mutant product with the patient's own
  baseline product.

Keeping the baseline in the merge identity is load-bearing.  Two phase
hypotheses can produce the same mutant protein while deleting different
patient-specific residues; those are different effects even though the
mutant molecule alone is identical.
"""

import itertools
from dataclasses import dataclass, field
from typing import Any, Mapping, Optional, Tuple


@dataclass(frozen=True, order=True)
class SpliceSiteKey:
    """Stable identity of a splice site in a rearranged layout."""

    side: str
    exon_number: int
    occurrence: int = 1

    def __post_init__(self):
        if self.side not in ("donor", "acceptor"):
            raise ValueError("side must be 'donor' or 'acceptor'")
        if self.exon_number < 1 or self.occurrence < 1:
            raise ValueError("exon_number and occurrence must be positive")


@dataclass(frozen=True)
class SpliceOption:
    """One graph rewrite offered at a disrupted splice site.

    ``ordinal_rank`` expresses preference only.  It is not a probability.
    A calibrated scorer can independently populate ``probability``.
    """

    mechanism: str
    ordinal_rank: int
    skipped_runs: Tuple[Any, ...] = field(default_factory=tuple)
    required_runs: Tuple[Any, ...] = field(default_factory=tuple)
    probability: Optional[float] = None
    evidence: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        if self.ordinal_rank < 0:
            raise ValueError("ordinal_rank must be non-negative")
        if self.probability is not None and not 0.0 <= self.probability <= 1.0:
            raise ValueError("probability must be between zero and one")


@dataclass(frozen=True)
class SpliceAxis:
    """Choices attached to one splice-graph site."""

    key: SpliceSiteKey
    anchor_runs: Tuple[Any, ...]
    options: Tuple[SpliceOption, ...]

    def __post_init__(self):
        if not self.anchor_runs:
            raise ValueError("A splice axis requires at least one anchor run")
        if not self.options:
            raise ValueError("A splice axis requires at least one option")


@dataclass(frozen=True)
class SplicePlan:
    """A graph-valid set of splice choices."""

    choices: Tuple[Tuple[SpliceSiteKey, SpliceOption], ...]
    kept_runs: Tuple[Any, ...]

    @property
    def ordinal_key(self):
        """Lexicographic mechanism preference; lower is preferred."""
        return tuple(option.ordinal_rank for _, option in self.choices)

    @property
    def probability(self):
        """Joint probability only when every selected option is calibrated."""
        probabilities = [option.probability for _, option in self.choices]
        if any(value is None for value in probabilities):
            return None
        result = 1.0
        for value in probabilities:
            result *= value
        return result


def enumerate_splice_plans(axes, run_keys, max_plans=64):
    """Enumerate graph-valid plans without multiplying dead splice sites.

    A Cartesian product treats every initially disrupted site as independent.
    That produces impossible combinations when one choice skips the exon that
    another choice needs.  Here every axis also has an implicit inactive
    choice.  It may be inactive only when another selected rewrite removes
    all of its anchor runs.  Conversely, a selected option is rejected when
    another rewrite removes one of its required runs.
    """
    axes = tuple(sorted(axes, key=lambda axis: axis.key))
    run_keys = tuple(run_keys)
    plans = []
    # Put real choices before the inactive sentinel so enumeration order
    # follows each axis's declared biological preference.
    selections = [(axis.options + (None,)) for axis in axes]
    for combination in itertools.product(*selections):
        skipped_by = {}
        for axis, option in zip(axes, combination):
            if option is None:
                continue
            for run in option.skipped_runs:
                skipped_by.setdefault(run, set()).add(axis.key)
        kept = tuple(run for run in run_keys if run not in skipped_by)
        kept_set = set(kept)
        valid = True
        choices = []
        for axis, option in zip(axes, combination):
            removed_by_other = any(
                any(owner != axis.key for owner in skipped_by.get(run, ()))
                for run in axis.anchor_runs)
            if option is None:
                if not all(run not in kept_set for run in axis.anchor_runs):
                    valid = False
                    break
                continue
            if removed_by_other:
                valid = False
                break
            if any(run not in kept_set for run in option.required_runs):
                valid = False
                break
            choices.append((axis.key, option))
        if not valid:
            continue
        plans.append(SplicePlan(tuple(choices), kept))
        if len(plans) > max_plans:
            raise ValueError(
                "Splice hypothesis count exceeds max_plans=%d" % max_plans)
    return tuple(plans)


@dataclass(frozen=True)
class EffectHypothesis:
    """One assignment of phase, splicing and structural resolution."""

    phase: Tuple[Tuple[Any, str], ...] = field(default_factory=tuple)
    splice_plan: Optional[SplicePlan] = None
    phase_probability: Optional[float] = 1.0
    evidence: Mapping[str, Any] = field(default_factory=dict)
    enumeration_index: int = 0


@dataclass(frozen=True)
class RealizedTranscriptProduct:
    """A materialized mRNA/protein plus its non-reference junctions."""

    transcript: Optional[Any]
    cdna_sequence: Optional[str]
    protein_sequence: Optional[str]
    junction_signature: Tuple[Tuple[Any, ...], ...] = field(
        default_factory=tuple)
    start_codon_present: Optional[bool] = None
    evidence: Mapping[str, Any] = field(default_factory=dict)

    def local_protein(self, offset, flank=30):
        """Protein window around a change, or the full product if unknown."""
        if self.protein_sequence is None:
            return None
        if offset is None:
            return self.protein_sequence
        start = max(0, offset - flank)
        end = min(len(self.protein_sequence), offset + flank + 1)
        return self.protein_sequence[start:end]

    def local_cdna(self, offset, flank=90):
        """cDNA window for products with no translated difference."""
        if self.cdna_sequence is None:
            return None
        if offset is None:
            return self.cdna_sequence
        start = max(0, offset - flank)
        end = min(len(self.cdna_sequence), offset + flank + 1)
        return self.cdna_sequence[start:end]


def _effect_identity(effect):
    """Stable semantic fields used in candidate merge identity."""
    return (
        type(effect).__name__,
        getattr(effect, "aa_mutation_start_offset", None),
        getattr(effect, "aa_ref", None),
        getattr(effect, "aa_alt", None),
        getattr(effect, "splice_signal", None),
    )


@dataclass(frozen=True)
class ClassifiedOutcome:
    """An effect classified between a baseline and mutant product."""

    hypothesis: EffectHypothesis
    effect: Any
    baseline: RealizedTranscriptProduct
    mutant: RealizedTranscriptProduct
    change_aa_offset: Optional[int] = None
    change_cdna_offset: Optional[int] = None

    def merge_key(self, protein_flank=30, cdna_flank=90):
        """Identity of the comparison, not merely the mutant product."""
        return (
            _effect_identity(self.effect),
            self.baseline.local_protein(self.change_aa_offset, protein_flank),
            self.mutant.local_protein(self.change_aa_offset, protein_flank),
            self.baseline.local_cdna(self.change_cdna_offset, cdna_flank),
            self.mutant.local_cdna(self.change_cdna_offset, cdna_flank),
            self.baseline.junction_signature,
            self.mutant.junction_signature,
        )


@dataclass(frozen=True)
class RealizedEffectCandidate:
    """Merged equivalent outcomes with all hypothesis provenance retained."""

    effect: Any
    outcomes: Tuple[ClassifiedOutcome, ...]
    ordinal_key: Tuple[int, ...]
    probability: Optional[float]
    phase_mass: Optional[float]
    first_enumeration_index: int

    @property
    def hypotheses(self):
        return tuple(outcome.hypothesis for outcome in self.outcomes)


def merge_classified_outcomes(outcomes, protein_flank=30, cdna_flank=90):
    """Merge only outcomes with identical baseline→mutant comparisons."""
    groups = {}
    order = []
    for outcome in outcomes:
        key = outcome.merge_key(protein_flank, cdna_flank)
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(outcome)

    candidates = []
    for key in order:
        group = tuple(groups[key])
        plans = [outcome.hypothesis.splice_plan for outcome in group]
        ordinal_keys = [
            plan.ordinal_key if plan is not None else () for plan in plans]
        probabilities = []
        for outcome, plan in zip(group, plans):
            splice_probability = (
                plan.probability if plan is not None else 1.0)
            phase_probability = outcome.hypothesis.phase_probability
            probabilities.append(
                splice_probability * phase_probability
                if splice_probability is not None
                and phase_probability is not None
                else None)
        probability = (
            sum(probabilities)
            if all(value is not None for value in probabilities) else None)
        phase_probabilities = [
            outcome.hypothesis.phase_probability for outcome in group]
        phase_mass = (
            sum(phase_probabilities)
            if all(value is not None for value in phase_probabilities)
            else None)
        candidates.append(RealizedEffectCandidate(
            effect=group[0].effect,
            outcomes=group,
            ordinal_key=min(ordinal_keys),
            probability=probability,
            phase_mass=phase_mass,
            first_enumeration_index=min(
                outcome.hypothesis.enumeration_index for outcome in group)))
    return tuple(candidates)


def order_realized_candidates(candidates, effect_priority_fn):
    """Order candidates without pretending ordinal rules are probabilities."""
    candidates = tuple(candidates)
    all_calibrated = bool(candidates) and all(
        candidate.probability is not None for candidate in candidates)

    def key(candidate):
        probability_key = (
            -candidate.probability if all_calibrated else 0.0)
        return (
            probability_key,
            candidate.ordinal_key,
            -effect_priority_fn(candidate.effect),
            candidate.first_enumeration_index,
        )

    return tuple(sorted(candidates, key=key))
