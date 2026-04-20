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

"""Effect annotator interface (openvax/varcode#271, stage 1).

An :class:`EffectAnnotator` takes a :class:`Variant` and a
:class:`Transcript` and returns a :class:`MutationEffect` (or a
:class:`MutantTranscript` consumed by one, depending on the
implementation). Annotators coexist behind a shared Protocol so
users can choose between:

* ``fast`` — the offset-based annotator that has shipped since
  2.0.0. Wraps :func:`varcode.effects.predict_variant_effect_on_transcript`.
* ``protein_diff`` — the coming annotator that materializes a
  :class:`MutantTranscript` and diffs its translated protein
  against the reference. Not in this stage; see #271.

Third parties (Isovar, Exacto) can register their own annotators by
implementing the Protocol and calling :func:`register_annotator`.

This stage 1 PR ships only the Protocol + registry + fast wrapper;
the protein-diff annotator, fast-path routing, per-call selection on
``Variant.effects()``, and ``EffectCollection`` provenance fields
land in follow-up PRs as outlined in #271.
"""

from typing import Protocol, runtime_checkable

from .fast import FastEffectAnnotator
from .protein_diff import ProteinDiffEffectAnnotator
from .structural_variant import StructuralVariantAnnotator
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
    * ``supports`` — set of variant-kind tags the annotator can
      handle (e.g. ``{"snv", "indel"}``). Callers that hand the
      annotator a variant outside this set get a clear
      :class:`UnsupportedVariantError` rather than silently wrong
      output.
    * :meth:`annotate_on_transcript` — the per-transcript entry
      point.

    Optionally exposes ``version`` (string) — used in CSV provenance
    headers so readers can detect when a serialized collection came
    from a different annotator version. Built-in annotators track
    varcode's version; third-party annotators expose their own.

    The protocol is intentionally narrow at this stage — additional
    methods (``annotate_collection``, ``annotate_with_context``) will
    be added as downstream work needs them. The contract is
    duck-typed (``@runtime_checkable``) so third-party annotators
    don't need to inherit from varcode just to register.
    """

    name: str
    supports: frozenset

    def annotate_on_transcript(self, variant, transcript):
        ...


__all__ = [
    "EffectAnnotator",
    "FastEffectAnnotator",
    "ProteinDiffEffectAnnotator",
    "StructuralVariantAnnotator",
    "UnsupportedVariantError",
    "get_annotator",
    "get_default_annotator",
    "register_annotator",
    "resolve_annotator",
    "set_default_annotator",
    "use_annotator",
]
