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

"""Tests for the MutantTranscript data model (openvax/varcode#271)."""

import pytest

import varcode
from varcode import MutantTranscript, TranscriptEdit


def test_transcript_edit_substitution_shape():
    edit = TranscriptEdit(cdna_start=10, cdna_end=13, alt_bases="TTT")
    assert not edit.is_insertion
    assert not edit.is_deletion
    assert edit.length_delta == 0


def test_transcript_edit_insertion_shape():
    edit = TranscriptEdit(cdna_start=10, cdna_end=10, alt_bases="ACG")
    assert edit.is_insertion
    assert not edit.is_deletion
    assert edit.length_delta == 3


def test_transcript_edit_deletion_shape():
    edit = TranscriptEdit(cdna_start=10, cdna_end=16, alt_bases="")
    assert not edit.is_insertion
    assert edit.is_deletion
    assert edit.length_delta == -6


def test_transcript_edit_rejects_negative_start():
    with pytest.raises(ValueError):
        TranscriptEdit(cdna_start=-1, cdna_end=0, alt_bases="A")


def test_transcript_edit_rejects_end_before_start():
    with pytest.raises(ValueError):
        TranscriptEdit(cdna_start=10, cdna_end=5, alt_bases="")


def test_transcript_edit_is_frozen():
    edit = TranscriptEdit(cdna_start=10, cdna_end=13, alt_bases="TTT")
    with pytest.raises(Exception):
        edit.cdna_start = 0


def test_mutant_transcript_no_edits_is_identity():
    mt = MutantTranscript(reference_transcript=object())
    assert mt.is_identical_to_reference
    assert mt.total_length_delta == 0
    assert mt.cdna_sequence is None
    assert mt.mutant_protein_sequence is None
    assert mt.annotator_name == "unknown"


def test_mutant_transcript_sums_length_deltas():
    mt = MutantTranscript(
        reference_transcript=object(),
        edits=(
            TranscriptEdit(10, 13, "T"),         # -2
            TranscriptEdit(100, 100, "AAAA"),    # +4
            TranscriptEdit(200, 210, ""),        # -10
        ),
    )
    assert mt.total_length_delta == -8
    assert not mt.is_identical_to_reference


def test_mutant_transcript_rejects_out_of_order_edits():
    with pytest.raises(ValueError):
        MutantTranscript(
            reference_transcript=object(),
            edits=(
                TranscriptEdit(100, 100, "A"),
                TranscriptEdit(10, 13, "T"),  # out of order
            ),
        )


def test_mutant_transcript_is_frozen():
    mt = MutantTranscript(reference_transcript=object())
    with pytest.raises(Exception):
        mt.annotator_name = "different"


def test_mutant_transcript_carries_sequences_when_producer_supplies_them():
    mt = MutantTranscript(
        reference_transcript=object(),
        edits=(TranscriptEdit(10, 13, "TTT"),),
        cdna_sequence="AAA...TTT...GGG",
        mutant_protein_sequence="MKL*",
        annotator_name="sequence_diff",
    )
    assert mt.cdna_sequence == "AAA...TTT...GGG"
    assert mt.mutant_protein_sequence == "MKL*"
    assert mt.annotator_name == "sequence_diff"


def test_mutant_transcript_is_exported_at_package_root():
    assert varcode.MutantTranscript is MutantTranscript
    assert varcode.TranscriptEdit is TranscriptEdit
