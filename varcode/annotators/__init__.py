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

"""One built-in default annotator and optional partial implementations.

``transcript_model`` is the shared default. ``fast`` retains the compatibility
implementation, and ``protein_diff`` is an alternative for point edits. Third
parties register any object satisfying :class:`EffectAnnotator`; they need no
capability list.
"""

from __future__ import annotations

from types import NotImplementedType
from typing import TYPE_CHECKING, Protocol, runtime_checkable

if TYPE_CHECKING:
    from pyensembl import Transcript

    from ..effects.effect_classes import MutationEffect
    from ..variant import Variant

from .fast import FastEffectAnnotator
from .protein_diff import ProteinDiffEffectAnnotator
from ..realized_effects import RealizedEffectAnnotator
from ..transcript_model import TranscriptModelEffectAnnotator
from .registry import (
    get_annotator,
    get_default_annotator,
    register_annotator,
    resolve_annotator,
    set_default_annotator,
    use_annotator,
)


@runtime_checkable
class EffectAnnotator(Protocol):
    """Protocol for an object that annotates variant effects on
    transcripts.

    Conforming objects expose:

    * ``name`` — short identifier (e.g. ``"fast"``) used in the
      registry and in serialized provenance.
    * :meth:`annotate_on_transcript` — the per-transcript entry
      point, returning a ``MutationEffect`` or ``NotImplemented``.

    Return Python's ``NotImplemented`` singleton when this particular input
    is unsupported. Public prediction APIs expose it as an ``Unresolved``
    effect with a reason, retaining this annotator's provenance. They never
    silently replace a selected annotator's result with another's prediction.
    ``None`` is invalid; exceptions retain normal error-handling semantics.

    Optionally exposes ``version`` (string) — used in CSV provenance
    headers so readers can detect when a serialized collection came
    from a different annotator version. Built-in annotators track
    varcode's version; third-party annotators expose their own.

    Optionally implement ``annotate_with_context(variant, transcript,
    germline_ctx, phase_resolver=None)`` with the same return contract.
    Without it, nonempty germline context is unsupported. Empty context
    calls ``annotate_on_transcript`` as usual.

    Optionally implement ``annotate_haplotype(variants, transcript,
    germline_ctx=None, phase_resolver=None)`` with the same return contract.
    Collections delegate each known-cis group to this method. A missing hook
    or refusal leaves an explicit unresolved group alongside individual effects.

    The contract is duck-typed (``@runtime_checkable``) so third-party annotators
    don't need to inherit from varcode just to register.
    """

    name: str

    def annotate_on_transcript(
            self, variant: Variant, transcript: Transcript,
    ) -> MutationEffect | NotImplementedType:
        """Predict a transcript consequence or explicitly decline the input.

        Parameters
        ----------
        variant : Variant
            Genomic variant, including StructuralVariant subclasses.
        transcript : pyensembl.Transcript
            Reference transcript selected for annotation.

        Returns
        -------
        MutationEffect or types.NotImplementedType
            A consequence (which may retain multiple outcomes), or
            ``NotImplemented`` for an unsupported input. Never return None.
        """
        ...


__all__ = [
    "EffectAnnotator",
    "FastEffectAnnotator",
    "ProteinDiffEffectAnnotator",
    "RealizedEffectAnnotator",
    "TranscriptModelEffectAnnotator",
    "get_annotator",
    "get_default_annotator",
    "register_annotator",
    "resolve_annotator",
    "set_default_annotator",
    "use_annotator",
]
