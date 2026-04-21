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

"""Unified ``Outcome`` type for multi-outcome effects (#299).

An ``Outcome`` is "one possible consequence of a variant" with
provenance attached. A single variant can have multiple outcomes —
splice disruption vs normal splicing (SNV at exon boundary), fusion
vs exon loss vs readthrough (translocation), missense vs UTR (variant
that spans a stop codon). Each outcome carries:

* the :class:`~varcode.effects.MutationEffect` it represents
* a probability / confidence score (or ``None`` if unscored)
* a string naming the source that produced the outcome
* an open-ended ``evidence`` dict for provenance

This shape is deliberately minimal. It's the interchange format
between:

* varcode's built-in classifiers (nominate candidates)
* external splice predictors (SpliceAI, Pangolin — score candidates)
* RNA-seq evidence callers (Isovar-style — upgrade / downgrade
  candidates with read support)
* long-read / assembly-based tools (resolve breakpoint ambiguity)

None of those integrations are wired up yet; the type exists so they
can plug in without churning the core classes.

Harmonization across effect kinds
---------------------------------

Today varcode has three overlapping multi-outcome shapes:

* :class:`~varcode.effects.effect_classes.ExonicSpliceSite` —
  ``alternate_effect`` is a single fallback.
* :mod:`varcode.splice_outcomes` — ``SpliceOutcomeSet`` is a richer
  list of splice-aware candidates.
* :class:`~varcode.effects.effect_classes.MultiOutcomeEffect` —
  marker base class with a ``candidates`` tuple.

All three converge on a single surface via the
:meth:`MultiOutcomeEffect.outcomes` property: ``tuple[Outcome, ...]``.
Existing ``candidates`` / ``alternate_effect`` accessors stay for
back-compat. New code should read ``outcomes``.
"""

from dataclasses import dataclass, field
from typing import Any, Mapping, Optional, Tuple

from serializable import DataclassSerializable


@dataclass(frozen=True)
class Outcome(DataclassSerializable):
    """One plausible consequence of a variant.

    Parameters
    ----------
    effect : MutationEffect
        The effect this outcome represents. Guaranteed to be a
        :class:`~varcode.effects.MutationEffect` instance — producers
        that can't compute a full coding effect use placeholder
        subclasses (e.g.
        :class:`~varcode.effects.effect_classes.PredictedIntronRetention`)
        so consumers can iterate ``outcome.effect.short_description``
        and ``outcome.effect.mutant_protein_sequence`` uniformly
        across SV, splice, and point-variant outcomes (#339).
    probability : float or None
        Estimated likelihood this outcome actually happens, in
        ``[0, 1]``. ``None`` means "not scored" — the outcome is in
        the candidate set but no tool has assigned a probability.
        Callers that treat ``None`` as "unknown" rather than "zero"
        will be robust to new outcomes added over time.
    source : str
        Name of the tool or annotator that produced this outcome.
        Defaults to ``"varcode"`` for built-in classifications.
        External integrations set their own (``"spliceai"``,
        ``"isovar"``, ``"longread_assembly"``, etc.) so downstream
        callers can filter by source.
    evidence : Mapping[str, Any]
        Open-ended provenance dict. Shape is source-specific; the
        convention is that keys match the source's native field names
        (e.g. SpliceAI scores under ``{"ds_ag": 0.12, "ds_al": ...}``,
        RNA read counts under ``{"junction_reads": 42}``). Consumers
        that need a particular shape should type-check at the call
        site rather than rely on a rigid schema here.
    description : str or None
        Optional human-readable sentence describing this specific
        outcome ("Exon 7 is skipped (in-frame, 15 aa removed)").
        Distinct from ``effect.short_description``, which is the
        effect's HGVS-style label. Producers that want a richer
        narrative attach it here rather than nesting it in
        ``evidence``. Declared last so positional calls of the form
        ``Outcome(effect, probability, source, evidence)`` continue
        to route the dict to ``evidence``.
    """

    effect: Any  # MutationEffect — typed loosely to avoid import cycle
    probability: Optional[float] = None
    source: str = "varcode"
    evidence: Mapping[str, Any] = field(default_factory=dict)
    description: Optional[str] = None

    def __post_init__(self):
        if self.probability is not None and not (
                0.0 <= self.probability <= 1.0):
            raise ValueError(
                "probability must be in [0, 1] or None, got %r"
                % (self.probability,))

    @property
    def short_description(self) -> str:
        """Convenience passthrough to the wrapped effect's
        ``short_description``. Lets callers build tables of outcomes
        without unpacking ``outcome.effect.short_description``
        everywhere."""
        return self.effect.short_description


def outcomes_from_candidates(
        candidates: Tuple[Any, ...],
        source: str = "varcode") -> Tuple[Outcome, ...]:
    """Wrap a tuple of :class:`MutationEffect` instances as
    :class:`Outcome` objects with a shared ``source`` string.

    Used by the default ``MultiOutcomeEffect.outcomes`` implementation
    to lift the existing ``candidates`` tuple into the new type without
    churning every subclass. No probabilities are assigned — callers
    that want scored outcomes construct :class:`Outcome` directly.
    """
    return tuple(Outcome(effect=c, source=source) for c in candidates)
