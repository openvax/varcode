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

"""The legacy :class:`EffectAnnotator` — a thin wrapper around the
offset-based effect prediction that varcode has shipped since 2.0.0.

Exists as an :class:`EffectAnnotator` Protocol implementation so
that the coming protein-diff annotator can be introduced behind
the same interface without churning callers (#271, stage 2). Until
then, this annotator is the default and produces byte-for-byte
identical output to ``Variant.effect_on_transcript(transcript)``.
"""

from ..version import __version__ as _varcode_version


class LegacyEffectAnnotator:
    """Wraps :func:`varcode.effects.predict_variant_effect_on_transcript`."""

    name = "legacy"

    version = _varcode_version
    """Built-in annotators track varcode's own version. Third-party
    annotators (isovar's plugin, exacto's plugin) expose their own
    version string here; CSV provenance headers and round-trip
    warnings read from this field. See #271."""

    supports = frozenset({"snv", "indel", "mnv"})
    """Variant kinds this annotator handles. Splice-possibility
    sets, structural variants, and phased haplotypes fall outside
    the legacy offset-based path and will be handled by the
    protein-diff annotator."""

    def annotate_on_transcript(self, variant, transcript):
        """Delegate to the existing per-transcript prediction.

        No fast-path / slow-path dispatch at this stage; that lives
        on the protein-diff annotator once it exists.
        """
        # Lazy import avoids a circular dep at package import time.
        from ..effects import predict_variant_effect_on_transcript
        return predict_variant_effect_on_transcript(variant, transcript)
