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

"""The fast :class:`EffectAnnotator` — a thin wrapper around the
offset-based effect prediction that varcode has shipped since 2.0.0.

The default annotator (restored as the default in 7.0.0; see #397).
Produces byte-for-byte identical output to
``Variant.effect_on_transcript(transcript)``. The
:class:`~varcode.annotators.protein_diff.ProteinDiffEffectAnnotator`
offers an alternative protein-diff classification behind the same
:class:`EffectAnnotator` protocol (#271); as of 7.0.0 the two are
reconciled on SNVs / indels / MNVs.
"""

from ..version import __version__ as _varcode_version


class FastEffectAnnotator:
    """Wraps :func:`varcode.effects.predict_variant_effect_on_transcript`."""

    name = "fast"

    version = _varcode_version
    """Built-in annotators track varcode's own version. Third-party
    annotators (isovar's plugin, exacto's plugin) expose their own
    version string here; CSV provenance headers and round-trip
    warnings read from this field. See #271."""

    supports = frozenset({"snv", "indel", "mnv"})
    """Variant kinds this annotator handles. Splice-possibility
    sets, structural variants, and phased haplotypes fall outside
    the fast offset-based path and will be handled by the
    protein-diff annotator."""

    def annotate_on_transcript(self, variant, transcript):
        """Delegate to the existing per-transcript prediction.

        Returns the raw effect class (``ExonicSpliceSite`` /
        ``SpliceDonor`` / etc. for splice disruptions), **not** wrapped
        in ``SpliceOutcomeSet``. The wrap is applied at the collection
        boundary in :func:`predict_variant_effects` so internal
        consumers (notably the ``protein_diff`` annotator's dual
        dispatch) can still pattern-match on the raw class.
        """
        # Lazy import avoids a circular dep at package import time.
        from ..effects.effect_prediction import (
            _predict_variant_effect_on_transcript_raw,
        )
        return _predict_variant_effect_on_transcript_raw(variant, transcript)
