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

"""``EffectCandidate``: one plausible effect with per-context provenance.

An :class:`EffectCandidate` pairs a :class:`MutationEffect` with a
``source`` / ``probability`` / ``evidence`` triple describing **how
this particular candidate came to be**. The wrapper exists because
the same :class:`MutationEffect` instance can appear in multiple
candidate sets with different provenance per context — the wrapper
carries the context-specific labels without forcing a copy of the
Effect itself.

Concrete example: a splice candidate produced by varcode's splice
classifier surfaces inside its parent
:class:`SpliceOutcomeSet` tagged ``source="varcode"`` with a
mechanism Effect. When an SV breakpoint sits in the same splice
window, the SV annotator re-uses that same Effect on the SV's own
:class:`StructuralVariantEffect` — but tagged ``source="varcode_splice"``
with the SV's ``sv_type`` merged into evidence. The Effect is shared;
the metadata diverges. That's why the wrapper exists.

Without ``EffectCandidate``, the alternatives are:

* Copy the Effect at every re-tag site (breaks identity, adds
  overhead, has to be deep enough to detach evidence dicts).
* Put metadata on the Effect itself (forces all contexts to share
  one labeling — wrong).
* Carry parallel sidecar dicts keyed by ``id(effect)`` (fragile,
  awkward to serialize).

The wrapper is the cheapest answer.

Where ``EffectCandidate`` appears
---------------------------------

:attr:`MultiOutcomeEffect.candidates` returns
``tuple[EffectCandidate, ...]`` on every subclass (#382). When
callers only need the underlying Effects (no provenance), the
:attr:`MultiOutcomeEffect.effects` convenience property unwraps:
``tuple(c.effect for c in self.candidates)``.

External integrations (splice predictors, RNA-evidence callers,
long-read assembly tools) construct :class:`EffectCandidate`
instances to surface scored or annotated effects through the same
surface as varcode's built-ins.
"""

from dataclasses import dataclass, field
from typing import Any, Mapping, Optional, Tuple

from serializable import DataclassSerializable


@dataclass(frozen=True)
class EffectCandidate(DataclassSerializable):
    """One plausible effect for a variant, paired with provenance.

    Parameters
    ----------
    effect : MutationEffect
        The effect this candidate represents. Guaranteed to be a
        :class:`~varcode.effects.MutationEffect` instance. For
        outcomes whose protein math isn't yet resolved (e.g. a
        splice-mechanism Effect built without a ``genomic_sequence``
        provider), the inner Effect still satisfies the interface —
        ``aa_ref`` / ``aa_alt`` / ``mutant_protein_sequence`` are
        simply ``None``, and consumers can read
        ``candidate.effect.short_description`` uniformly across SV,
        splice, and point-variant candidates.
    probability : float or None
        Optional source-scoped score in ``[0, 1]``. It answers only
        "how likely does this producer think this candidate is within
        this candidate set?" Varcode-generated splice candidates use
        hand-tuned DNA-only priors here so the set has a stable order;
        RNA/model integrations may supply empirical or calibrated
        estimates. Varcode stores the value unchanged and does not
        normalize across sources. ``None`` means "not scored", not
        impossible.
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


def candidates_from_effects(
        effects: Tuple[Any, ...],
        source: str = "varcode") -> Tuple[EffectCandidate, ...]:
    """Wrap a tuple of :class:`MutationEffect` instances as
    :class:`EffectCandidate` objects with a shared ``source`` string.

    Convenience for the common case where a producer has a tuple of
    Effects and wants to lift it into the wrapped form without
    setting probabilities. Callers that want scored candidates
    construct :class:`EffectCandidate` directly.
    """
    return tuple(EffectCandidate(effect=c, source=source) for c in effects)
