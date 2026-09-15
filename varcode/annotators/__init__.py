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

``fast`` owns point-variant and structural routing. ``protein_diff`` and
``transcript_model`` remain experimental alternatives. Third parties register any
object satisfying :class:`EffectAnnotator`; they need no capability list.
"""

from typing import Protocol, runtime_checkable

from .fast import FastEffectAnnotator
from .protein_diff import ProteinDiffEffectAnnotator
from .structural_variant import StructuralVariantAnnotator
from ..realized_effects import RealizedEffectAnnotator
from ..transcript_model import TranscriptModelEffectAnnotator
from .registry import (
    UnsupportedVariantError,
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
    silently replace an experimental result with the default's prediction.
    ``None`` is invalid; exceptions retain normal error-handling semantics.

    Optionally exposes ``version`` (string) — used in CSV provenance
    headers so readers can detect when a serialized collection came
    from a different annotator version. Built-in annotators track
    varcode's version; third-party annotators expose their own.

    Optionally implement ``annotate_with_context(variant, transcript,
    germline_ctx, phase_resolver=None)`` with the same return contract.
    Without it, nonempty germline context is unsupported. Empty context
    calls ``annotate_on_transcript`` as usual.

    The contract is duck-typed (``@runtime_checkable``) so third-party annotators
    don't need to inherit from varcode just to register.
    """

    name: str

    def annotate_on_transcript(self, variant, transcript):
        ...


__all__ = [
    "EffectAnnotator",
    "FastEffectAnnotator",
    "ProteinDiffEffectAnnotator",
    "StructuralVariantAnnotator",
    "RealizedEffectAnnotator",
    "TranscriptModelEffectAnnotator",
    "UnsupportedVariantError",
    "get_annotator",
    "get_default_annotator",
    "register_annotator",
    "resolve_annotator",
    "set_default_annotator",
    "use_annotator",
]
