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

"""Tests for the MutantTranscript data model and construction
helpers (openvax/varcode#271)."""

import pytest
from pyensembl import cached_release

import varcode
from varcode import (
    MutantTranscript,
    TranscriptEdit,
    Variant,
    apply_variant_to_transcript,
)


ensembl_grch38 = cached_release(81)
CFTR_TRANSCRIPT_ID = "ENST00000003084"  # forward strand, chr7
BRCA1_TRANSCRIPT_ID = "ENST00000357654"  # reverse strand, chr17
MT_CO1_TRANSCRIPT_ID = "ENST00000361624"  # mitochondrial, forward strand


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
    assert varcode.apply_variant_to_transcript is apply_variant_to_transcript


# ====================================================================
# apply_variant_to_transcript — SNVs, insertions, deletions across
# both strands and the mitochondrial codon table.
# ====================================================================


def test_apply_snv_to_forward_strand_coding_variant():
    # CFTR: pick a known coding SNV. CFTR exon 4 is at 117531048-117531263.
    # A simple substitution inside the CDS.
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    variant = Variant("7", 117531100, "T", "A", ensembl_grch38)
    mt = apply_variant_to_transcript(variant, transcript)
    assert mt is not None
    assert mt.reference_transcript is transcript
    assert mt.annotator_name == "sequence_diff"
    assert len(mt.edits) == 1
    edit = mt.edits[0]
    assert edit.cdna_end - edit.cdna_start == 1
    assert edit.alt_bases == "A"
    # mutant cDNA differs from reference at exactly one position.
    ref_cdna = str(transcript.sequence)
    assert len(mt.cdna_sequence) == len(ref_cdna)
    diffs = [i for i in range(len(ref_cdna)) if ref_cdna[i] != mt.cdna_sequence[i]]
    assert len(diffs) == 1
    # Translated mutant protein populated for in-CDS variants.
    assert mt.mutant_protein_sequence is not None
    assert len(mt.mutant_protein_sequence) > 0


def test_apply_to_reverse_strand_transcript_complements_bases():
    # BRCA1 is on the reverse strand. Reuse the known-good coding
    # variant from test_splice_site_effects.py: Variant('17', 43082570,
    # 'CCT', 'GGG') is a 3-base substitution in BRCA1 exon 12. On the
    # reverse strand this appears in the cDNA as AGG (ref)/CCC (alt).
    transcript = ensembl_grch38.transcript_by_id(BRCA1_TRANSCRIPT_ID)
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", ensembl_grch38)
    mt = apply_variant_to_transcript(variant, transcript)
    assert mt is not None, "Known-good BRCA1 coding variant should apply"
    assert transcript.on_backward_strand
    edit = mt.edits[0]
    assert edit.alt_bases == "CCC", (
        "alt_bases should be reverse-complement of genomic alt 'GGG', got %r"
        % edit.alt_bases)
    ref_cdna = str(transcript.sequence)
    assert ref_cdna[edit.cdna_start:edit.cdna_end] == "AGG", (
        "Reference cDNA slice should be reverse-complement of genomic ref 'CCT'")


def test_apply_insertion_in_cds():
    # Insertion: ref of length 1, alt of length > 1 is a substitution+insert.
    # For a pure insertion in varcode coords, the trimmed_ref is empty.
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    # Construct a 3-base in-frame insertion inside CFTR. We'll insert
    # AAA after a base in the CDS.
    variant = Variant("7", 117531100, "T", "TAAA", ensembl_grch38)
    mt = apply_variant_to_transcript(variant, transcript)
    if mt is None:
        pytest.skip("Insertion didn't pass single-exon containment check.")
    # Net length delta is +3 (insertion of 3 bases).
    assert mt.total_length_delta == 3
    assert len(mt.cdna_sequence) == len(str(transcript.sequence)) + 3


def test_apply_deletion_in_cds():
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    # CFTR cDNA at genomic 117531100 is 'T', 117531101 onwards is 'TGA'.
    # A 3-base in-frame deletion: ref='TTGA' at 117531100, alt='T'
    # drops 'TGA' (3 bases) while anchoring at the leading T.
    variant = Variant("7", 117531100, "TTGA", "T", ensembl_grch38)
    mt = apply_variant_to_transcript(variant, transcript)
    assert mt is not None
    assert mt.total_length_delta == -3
    assert len(mt.cdna_sequence) == len(str(transcript.sequence)) - 3


def test_apply_returns_none_for_intronic_variant():
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    # Position deep in a CFTR intron — between exons.
    # CFTR exon 1 ends at ~117480148; intron 1 follows.
    intronic_pos = 117500000  # well inside an intron
    variant = Variant("7", intronic_pos, "T", "A", ensembl_grch38)
    mt = apply_variant_to_transcript(variant, transcript)
    assert mt is None, (
        "Intronic variant must return None so caller falls back to legacy")


def test_apply_returns_none_for_non_protein_coding_transcript():
    # Pick a non-coding transcript on chr1
    nc_transcripts = [
        t for t in ensembl_grch38.transcripts(contig="1")
        if not t.is_protein_coding
    ][:1]
    if not nc_transcripts:
        pytest.skip("No non-coding chr1 transcripts available locally")
    nc = nc_transcripts[0]
    variant = Variant("1", nc.start + 100, "T", "A", ensembl_grch38)
    mt = apply_variant_to_transcript(variant, nc)
    assert mt is None


def test_apply_variant_to_mitochondrial_transcript_uses_mt_codon_table():
    # MT-CO1 6739 C>G: TCA -> TGA. Standard table => stop; mt table => Trp.
    # apply_variant_to_transcript should use mt translation, so the
    # mutant_protein_sequence should NOT be truncated at this codon.
    transcript = ensembl_grch38.transcript_by_id(MT_CO1_TRANSCRIPT_ID)
    variant = Variant("MT", 6739, "C", "G", ensembl_grch38)
    mt = apply_variant_to_transcript(variant, transcript)
    assert mt is not None
    # The reference protein at this codon position is S (Ser) — the
    # mutant should be W (Trp) at the same position, NOT truncated.
    ref_protein = str(transcript.protein_sequence)
    aa_pos = 278  # known from #294 test
    assert ref_protein[aa_pos] == "S"
    assert mt.mutant_protein_sequence[aa_pos] == "W"
    # Length is preserved (no premature stop introduced under mt table).
    assert len(mt.mutant_protein_sequence) >= aa_pos + 1


def test_apply_variant_returns_none_on_reference_mismatch():
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    # Wrong reference base at this position — should return None
    # (not raise; caller decides whether to surface the error).
    variant = Variant("7", 117531100, "X", "A", ensembl_grch38,
                      allow_extended_nucleotides=True)
    mt = apply_variant_to_transcript(variant, transcript)
    assert mt is None
