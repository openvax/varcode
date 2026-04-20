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

"""Tests for the :class:`ReferenceSegment` extension of
:class:`MutantTranscript` (PR 9; #252).

Coverage is intentionally narrow — this PR only adds the data shape.
Actual SV annotation that consumes ``reference_segments`` lands in
PR 10; cryptic-exon enumeration in PR 11.
"""

import pytest
from pyensembl import cached_release

from varcode import MutantTranscript, ReferenceSegment

ensembl_grch38 = cached_release(81)


def test_reference_segment_basic_fields():
    class _FakeSource:
        sequence = "ACGTACGTACGT"
    seg = ReferenceSegment(
        source=_FakeSource(), start=2, end=8, strand="+", label="test")
    assert seg.start == 2
    assert seg.end == 8
    assert seg.length == 6
    assert seg.strand == "+"
    assert seg.label == "test"


def test_reference_segment_bad_strand_rejected():
    with pytest.raises(ValueError):
        ReferenceSegment(source=object(), start=0, end=10, strand="?")


def test_reference_segment_end_before_start_rejected():
    with pytest.raises(ValueError):
        ReferenceSegment(source=object(), start=10, end=5)


def test_mutant_transcript_point_variant_shape_unchanged():
    """Existing callers that pass ``reference_transcript=`` keep
    working exactly as before."""
    transcript = ensembl_grch38.transcript_by_id("ENST00000003084")
    mt = MutantTranscript(
        reference_transcript=transcript,
        edits=(),
        cdna_sequence=None,
        mutant_protein_sequence=None,
    )
    assert mt.reference_transcript is transcript
    assert mt.reference_segments is None
    assert mt.is_structural is False
    assert mt.is_identical_to_reference is True


def test_mutant_transcript_sv_shape_from_segments():
    """A fusion-style MutantTranscript: two segments, no single
    reference transcript (though the 5' partner is typically kept
    in :attr:`reference_transcript` for ordering)."""
    class _Source:
        def __init__(self, seq): self.sequence = seq
    five_prime = _Source("AAAACCCCC")
    three_prime = _Source("GGGGGTTTTT")

    segments = (
        ReferenceSegment(
            source=five_prime, start=0, end=9,
            strand="+", label="5p_partner"),
        ReferenceSegment(
            source=three_prime, start=5, end=10,
            strand="+", label="3p_partner"),
    )
    mt = MutantTranscript(
        reference_transcript=None,
        reference_segments=segments,
        annotator_name="structural_variant",
    )
    assert mt.is_structural is True
    assert mt.is_identical_to_reference is False
    assert len(mt.reference_segments) == 2
    assert mt.reference_segments[0].label == "5p_partner"


def test_mutant_transcript_requires_ref_or_segments():
    """At least one of the two shapes must be populated."""
    with pytest.raises(ValueError):
        MutantTranscript(reference_transcript=None, reference_segments=None)


def test_mutant_transcript_sv_shape_length_delta_zero():
    """Structural-variant shape doesn't have a meaningful
    ``total_length_delta`` (the mutant isn't a delta on a single
    reference); it returns 0 so callers that sum it across variants
    don't crash."""
    class _Source:
        sequence = "AAAA"
    mt = MutantTranscript(
        reference_transcript=None,
        reference_segments=(
            ReferenceSegment(source=_Source(), start=0, end=4),
        ),
    )
    assert mt.total_length_delta == 0


def test_reference_segment_accepts_custom_source():
    """The ``source`` field is typed loosely so callers can plug in
    patient-specific contigs, long-read assemblies, or
    ``GenomicInterval`` stand-ins without varcode caring about the
    concrete type."""
    class _PatientContig:
        """Minimal custom source. Real pipelines would back this
        with a FASTA fragment from a personalized assembly."""
        sequence = "NNNNNNNNNNNN"

    seg = ReferenceSegment(
        source=_PatientContig(), start=0, end=12,
        label="patient_specific_contig_1")
    assert seg.length == 12
    assert seg.source.sequence == "NNNNNNNNNNNN"
