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

An :class:`EffectCandidate` pairs a :class:`MutationEffect` with
``source`` / ``evidence`` provenance describing **how this particular
candidate came to be**. The wrapper exists because
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
instances to surface annotated effects through the same surface as
varcode's built-ins.
"""

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Mapping, Tuple

from serializable import DataclassSerializable

if TYPE_CHECKING:
    from .effect_hypotheses import ClassifiedOutcome
    from .effects.effect_classes import MutationEffect


@dataclass(frozen=True)
class EffectCandidate(DataclassSerializable):
    """One plausible effect for a variant, paired with provenance and uncertainty.

    Parameters
    ----------
    effect : MutationEffect
        The effect this candidate represents. Guaranteed to be a
        :class:`~varcode.effects.MutationEffect` instance. Sequences live on
        ``effect`` (``mutant_protein_sequence``, ``original_protein_sequence``)
        and its optional ``mutant_transcript`` / ``baseline_transcript``.
    source : str
        Name of the tool or annotator that produced this candidate.
        Defaults to ``"varcode"`` for built-in classifications (``"varcode.dna"``
        from the default DNA annotator). External integrations set their own
        opaque string.
    evidence : Mapping[str, Any]
        Open-ended provenance dict.
    mechanisms : tuple of str, optional
        Splice mechanisms active for this candidate (e.g. ``("normal",)``,
        ``("exon_skip",)``, ``("intron_retention",)``, ``("cryptic",)``).
    phase : tuple of (Variant, str), optional
        Inherited germline phase assignments for this candidate, e.g.
        ``((germline_variant, "cis"),)``.
    uncertainty : tuple of str, optional
        Sorted axes of unresolved DNA uncertainty across sibling candidates
        (e.g. ``()``, ``("splicing",)``, ``("phasing",)``,
        ``("phasing", "splicing")``, ``("structural",)``).
    outcomes : tuple of ClassifiedOutcome, optional
        Paired baseline-vs-mutant transcript-model comparisons that collapsed
        into this candidate, when produced by the transcript-model / DNA
        composition pipeline. Preserved through serialization and excluded
        from ``repr()`` for readability. Empty tuple when unavailable.
    """

    effect: "MutationEffect"
    source: str = "varcode"
    evidence: Mapping[str, Any] = field(default_factory=dict)
    mechanisms: Tuple[str, ...] = ()
    phase: Tuple[Tuple[Any, str], ...] = ()
    uncertainty: Tuple[str, ...] = ()
    outcomes: Tuple["ClassifiedOutcome", ...] = field(default=(), repr=False)

    @property
    def short_description(self) -> str:
        """Convenience passthrough to ``self.effect.short_description``.
        Lets callers build tables without unpacking
        ``candidate.effect.short_description`` everywhere."""
        return self.effect.short_description


def candidates_from_effects(
        effects: Tuple["MutationEffect", ...],
        source: str = "varcode") -> Tuple[EffectCandidate, ...]:
    """Wrap a tuple of :class:`MutationEffect` instances as
    :class:`EffectCandidate` objects with a shared ``source`` string.

    Convenience for the common case where a producer has a tuple of
    Effects and wants to lift it into the wrapped form without
    setting evidence fields. Callers with producer-specific metadata
    construct :class:`EffectCandidate` directly.
    """
    return tuple(EffectCandidate(effect=c, source=source) for c in effects)
