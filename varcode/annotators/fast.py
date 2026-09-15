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

"""The built-in default annotator, registered as ``fast`` for compatibility.

Routes point edits to the established offset predictor and rearrangements
to the structural implementation. Callers do not need to select by kind.
"""

from ..version import __version__ as _varcode_version


class FastEffectAnnotator:
    """Annotate point edits and structural variants through one interface."""

    name = "fast"

    version = _varcode_version
    """Built-in annotators track varcode's own version. Third-party
    annotators (isovar's plugin, exacto's plugin) expose their own
    version string here; CSV provenance headers and round-trip
    warnings read from this field. See #271."""

    def annotate_on_transcript(self, variant, transcript):
        """Delegate to the existing per-transcript prediction.

        Returns the raw effect class (``ExonicSpliceSite`` /
        ``SpliceDonor`` / etc. for splice disruptions), **not** wrapped
        in ``SpliceOutcomeSet``. The wrap is applied at the collection
        boundary in :func:`predict_variant_effects` so internal
        consumers (notably the ``protein_diff`` annotator's dual
        dispatch) can still pattern-match on the raw class.
        """
        if getattr(variant, "is_structural", False):
            from .structural_variant import StructuralVariantAnnotator
            return StructuralVariantAnnotator().annotate_on_transcript(
                variant, transcript)
        # Lazy import avoids a circular dep at package import time.
        from ..effects.effect_prediction import (
            _predict_variant_effect_on_transcript_raw,
        )
        return _predict_variant_effect_on_transcript_raw(variant, transcript)

    def annotate_with_context(
            self, variant, transcript, germline_ctx, phase_resolver=None):
        """Use the established patient-baseline path for point edits.

        Structural haplotype composition remains experimental in ``transcript_model``.
        Do not send an SV's placeholder alleles to the point-edit builder.
        """
        if getattr(variant, "is_structural", False):
            return NotImplemented
        from ..germline import predict_germline_aware_effect
        return predict_germline_aware_effect(
            variant, transcript, germline_ctx, annotator=self,
            phase_resolver=phase_resolver)
