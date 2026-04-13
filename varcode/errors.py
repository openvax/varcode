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

"""
Exception types raised by varcode. ``ReferenceMismatchError`` subclasses
``ValueError`` for backwards compatibility with callers that catch
``ValueError`` already (including the internal
``predict_variant_effect_on_transcript_or_failure`` fallback).
"""


class SampleNotFoundError(KeyError):
    """Raised when genotype info is requested for a sample that isn't
    present in the VariantCollection's source VCF(s)."""


class ReferenceMismatchError(ValueError):
    """Raised when a variant's reported ref allele does not match the
    reference genome at the variant's position.

    This most often means one of:

    * The variant was called against a different reference build than
      the one being used for annotation (e.g. GRCh37 vs GRCh38).
    * The variant's ref field was populated with the patient's germline
      allele rather than the canonical reference. VCF requires the
      ref field to match the reference genome; germline variants at
      the same position should be encoded as separate variants.
    * Strand confusion: the variant is specified on the negative
      strand but varcode expects positive-strand coordinates.

    Callers who would rather continue past this error can pass
    ``raise_on_error=False`` to :meth:`Variant.effects` to receive
    ``Failure`` effects instead.
    """

    def __init__(self, variant, transcript, expected_ref, observed_ref,
                 transcript_offset=None, genome_start=None, genome_end=None):
        self.variant = variant
        self.transcript = transcript
        self.expected_ref = expected_ref
        self.observed_ref = observed_ref
        self.transcript_offset = transcript_offset
        self.genome_start = genome_start
        self.genome_end = genome_end

        location = ""
        if transcript_offset is not None:
            location = " at transcript offset %d" % transcript_offset
        if genome_start is not None and genome_end is not None:
            location += " (chromosome positions %d:%d)" % (
                genome_start, genome_end)

        message = (
            "Reference allele mismatch for %s on %s%s: variant reports "
            "ref=%r but the reference genome has %r at this position.\n"
            "This usually means the variant was called against a "
            "different genome build, the ref field was filled in with "
            "the patient's germline allele rather than the reference, "
            "or the variant is on the wrong strand. Pass "
            "raise_on_error=False to .effects() to receive a Failure "
            "effect instead of raising." % (
                variant, transcript, location, observed_ref, expected_ref)
        )
        super().__init__(message)
