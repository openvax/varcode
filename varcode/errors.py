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

from .nucleotides import reverse_complement


class SampleNotFoundError(KeyError):
    """Raised when genotype info is requested for a sample that isn't
    present in the VariantCollection's source VCF(s)."""


class InvalidContigError(ValueError):
    """The selected annotation dataset does not contain a variant's contig."""


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

    ``expected_ref`` and ``observed_ref`` retain the transcript-strand
    alleles supplied by the annotation engine. ``genomic_expected_ref``
    and ``genomic_observed_ref`` expose their forward-strand equivalents,
    which are used in the message alongside genomic coordinates.
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

        if getattr(transcript, "on_backward_strand", False):
            self.genomic_expected_ref = reverse_complement(expected_ref)
            self.genomic_observed_ref = reverse_complement(observed_ref)
        else:
            self.genomic_expected_ref = expected_ref
            self.genomic_observed_ref = observed_ref

        start = variant.trimmed_base1_start if genome_start is None else genome_start
        end = variant.trimmed_base1_end if genome_end is None else genome_end
        location = "%s:%d" % (variant.contig, start)
        if end != start:
            location += "-%d" % end

        context = []
        if transcript is not None:
            context.append(transcript.id)
            if transcript.name and transcript.name != transcript.id:
                context.append(transcript.name)
            context.append("strand %s" % transcript.strand)
        if transcript_offset is not None:
            context.append("transcript offset %d" % transcript_offset)

        message = (
            "Reference allele mismatch at %s: variant ref=%r, "
            "reference genome has %r%s. Check the genome build and allele strand." % (
                location, self.genomic_observed_ref, self.genomic_expected_ref,
                " (%s)" % ", ".join(context) if context else "")
        )
        super().__init__(message)


class HypothesisLimitError(ValueError):
    """Raised when a variant has more alternative outcomes than
    ``max_hypotheses`` allows, so they were not enumerated.

    Annotators catch this and return a :class:`~varcode.HypothesisLimit`
    effect instead of choosing among the outcomes.
    """

    def __init__(self, max_hypotheses):
        super().__init__(max_hypotheses)
        self.max_hypotheses = max_hypotheses

    def __str__(self):
        return (
            "More than max_hypotheses=%d alternative outcomes"
            % self.max_hypotheses)
