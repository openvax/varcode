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

"""``EffectOutcome``: one plausible effect with per-context provenance.

An :class:`EffectOutcome` pairs a :class:`MutationEffect` with a
``source`` / ``probability`` / ``evidence`` triple describing **how
this particular candidate came to be**. The point of the wrapper is
that the same :class:`MutationEffect` instance can appear in multiple
candidate sets with different provenance per context — the wrapper
carries the context-specific labels without forcing a copy of the
Effect itself.

Concrete example: a splice candidate produced by varcode's splice
classifier gets surfaced once inside its parent
:class:`SpliceOutcomeSet` (tagged ``source="varcode"`` with a
``splice_outcome`` enum in evidence). When an SV breakpoint sits in
the same splice window, the SV annotator re-uses that same Effect
object as an attached candidate on the SV's own
:class:`StructuralVariantEffect` — but tagged ``source="varcode_splice"``
and with the SV's ``sv_type`` merged into evidence. The Effect is
shared; the metadata diverges. That's what the wrapper exists for.

Without ``EffectOutcome``, the alternatives are:

* Copy the Effect at every re-tag site (breaks identity, adds
  overhead, has to be deep enough to detach evidence dicts).
* Put metadata on the Effect itself (forces all contexts to share
  one labeling — wrong).
* Carry parallel sidecar dicts keyed by ``id(effect)`` (fragile,
  awkward to serialize).

The wrapper is the cheapest answer.

Where ``EffectOutcome`` appears
---------------------------------

The :attr:`MultiOutcomeEffect.outcomes` property returns
``tuple[EffectOutcome, ...]``. ``candidates`` returns the raw
``tuple[MutationEffect, ...]`` — both accessors are first-class:

* Use ``outcomes`` when you need per-candidate provenance (filter by
  source, read evidence, score by probability).
* Use ``candidates`` when you only need the underlying Effects (e.g.
  iterate to render short descriptions, dispatch on effect type).

External integrations (SpliceAI/Pangolin scorers, RNA-evidence
callers, long-read assembly tools) construct
:class:`EffectOutcome` instances when they want to surface scored
or annotated effects through the same surface as varcode's
built-ins. None of those external integrations are wired up yet;
the type ships ahead of them so they can plug in without churning
core classes.
"""

from dataclasses import dataclass, field
from typing import Any, Mapping, Optional, Tuple

from serializable import DataclassSerializable


@dataclass(frozen=True)
class EffectOutcome(DataclassSerializable):
    """One plausible effect for a variant, paired with provenance.

    Parameters
    ----------
    effect : MutationEffect
        The effect this candidate represents. Guaranteed to be a
        :class:`~varcode.effects.MutationEffect` instance — producers
        that can't compute a full coding effect use placeholder
        subclasses (e.g.
        :class:`~varcode.effects.effect_classes.PredictedIntronRetention`)
        so consumers can read
        ``candidate.effect.short_description`` and
        ``candidate.effect.mutant_protein_sequence`` uniformly across
        SV, splice, and point-variant candidates.
    probability : float or None
        Estimated likelihood this candidate actually happens, in
        ``[0, 1]``. ``None`` means "not scored" — the candidate is in
        the set but no tool has assigned a probability. Callers that
        treat ``None`` as "unknown" rather than "zero" will be robust
        to new candidates added over time.
    source : str
        Name of the tool or annotator that produced this candidate.
        Defaults to ``"varcode"`` for built-in classifications.
        External integrations set their own opaque string. Varcode
        does not interpret this field beyond exact-string equality.
    evidence : Mapping[str, Any]
        Open-ended provenance dict. Shape is source-specific; the
        convention is that keys match the source's native field names.
        Consumers that need a particular shape should type-check at
        the call site rather than rely on a rigid schema here.
    """

    effect: Any  # MutationEffect — typed loosely to avoid import cycle
    probability: Optional[float] = None
    source: str = "varcode"
    evidence: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        if self.probability is not None and not (
                0.0 <= self.probability <= 1.0):
            raise ValueError(
                "probability must be in [0, 1] or None, got %r"
                % (self.probability,))

    @property
    def short_description(self) -> str:
        """Convenience passthrough to ``self.effect.short_description``.
        Lets callers build tables without unpacking
        ``candidate.effect.short_description`` everywhere."""
        return self.effect.short_description


def outcomes_from_effects(
        effects: Tuple[Any, ...],
        source: str = "varcode") -> Tuple[EffectOutcome, ...]:
    """Wrap a tuple of :class:`MutationEffect` instances as
    :class:`EffectOutcome` objects with a shared ``source`` string.

    Used by the default :attr:`MultiOutcomeEffect.outcomes`
    implementation to lift a raw ``candidates`` tuple into the wrapped
    form without churning every subclass. No probabilities are
    assigned — callers that want scored candidates construct
    :class:`EffectOutcome` directly.
    """
    return tuple(EffectOutcome(effect=c, source=source) for c in effects)
