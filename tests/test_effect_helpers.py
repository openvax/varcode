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

from varcode.effects.effect_helpers import (
    changes_exonic_splice_site,
    variant_overlaps_interval,
)


class DummyTranscript(object):
    def __init__(self, sequence, n_exons=2):
        self.sequence = sequence
        self.exons = [object()] * n_exons


def test_variant_overlaps_interval_single_base_at_interval_boundaries():
    assert variant_overlaps_interval(
        variant_start=3,
        n_ref_bases=1,
        interval_start=3,
        interval_end=5)
    assert variant_overlaps_interval(
        variant_start=5,
        n_ref_bases=1,
        interval_start=3,
        interval_end=5)


def test_variant_overlaps_interval_excludes_adjacent_bases():
    assert not variant_overlaps_interval(
        variant_start=2,
        n_ref_bases=1,
        interval_start=3,
        interval_end=5)
    assert not variant_overlaps_interval(
        variant_start=6,
        n_ref_bases=1,
        interval_start=3,
        interval_end=5)


def test_variant_overlaps_interval_spanning_variants():
    assert variant_overlaps_interval(
        variant_start=2,
        n_ref_bases=2,
        interval_start=3,
        interval_end=5)
    assert not variant_overlaps_interval(
        variant_start=1,
        n_ref_bases=2,
        interval_start=3,
        interval_end=5)


def test_variant_overlaps_interval_insertions():
    assert variant_overlaps_interval(
        variant_start=3,
        n_ref_bases=0,
        interval_start=3,
        interval_end=5)
    assert variant_overlaps_interval(
        variant_start=5,
        n_ref_bases=0,
        interval_start=3,
        interval_end=5)
    assert not variant_overlaps_interval(
        variant_start=2,
        n_ref_bases=0,
        interval_start=3,
        interval_end=5)


def test_changes_exonic_splice_site_breaking_substitution_in_motif():
    transcript = DummyTranscript("AATCAGAA")
    result = changes_exonic_splice_site(
        transcript_offset=4,
        transcript=transcript,
        transcript_ref="A",
        transcript_alt="C",
        exon_start_offset=0,
        exon_end_offset=5,
        exon_number=1)
    assert result is True


def test_changes_exonic_splice_site_preserving_substitution_in_motif():
    transcript = DummyTranscript("AATCAGAA")
    result = changes_exonic_splice_site(
        transcript_offset=3,
        transcript=transcript,
        transcript_ref="C",
        transcript_alt="A",
        exon_start_offset=0,
        exon_end_offset=5,
        exon_number=1)
    assert result is None


def test_changes_exonic_splice_site_ignores_noncanonical_reference_motif():
    # Exon-end motif is CTG (not MAG), so this path should not trigger.
    transcript = DummyTranscript("AATCTGAA")
    result = changes_exonic_splice_site(
        transcript_offset=4,
        transcript=transcript,
        transcript_ref="T",
        transcript_alt="A",
        exon_start_offset=0,
        exon_end_offset=5,
        exon_number=1)
    assert result is None


def test_changes_exonic_splice_site_skips_end_motif_check_for_last_exon():
    # Last exon has no downstream splice donor boundary in this model.
    transcript = DummyTranscript("AATCAGAA", n_exons=1)
    result = changes_exonic_splice_site(
        transcript_offset=4,
        transcript=transcript,
        transcript_ref="A",
        transcript_alt="C",
        exon_start_offset=0,
        exon_end_offset=5,
        exon_number=1)
    assert result is None


def test_changes_exonic_splice_site_breaking_deletion_in_motif():
    transcript = DummyTranscript("AATCAGAA")
    result = changes_exonic_splice_site(
        transcript_offset=4,
        transcript=transcript,
        transcript_ref="A",
        transcript_alt="",
        exon_start_offset=0,
        exon_end_offset=5,
        exon_number=1)
    assert result is True


def test_changes_exonic_splice_site_preserving_insertion_in_motif():
    transcript = DummyTranscript("AATCAGAA")
    result = changes_exonic_splice_site(
        transcript_offset=3,
        transcript=transcript,
        transcript_ref="",
        transcript_alt="A",
        exon_start_offset=0,
        exon_end_offset=5,
        exon_number=1)
    assert result is None


def test_changes_exonic_splice_site_breaking_insertion_in_motif():
    transcript = DummyTranscript("AATCAGAA")
    result = changes_exonic_splice_site(
        transcript_offset=3,
        transcript=transcript,
        transcript_ref="",
        transcript_alt="T",
        exon_start_offset=0,
        exon_end_offset=5,
        exon_number=1)
    assert result is True


def test_changes_exonic_splice_site_breaking_acceptor_purine_substitution():
    # Exon start base A->T for exon_number > 1 should break YAG|R.
    transcript = DummyTranscript("CCATGGTT", n_exons=2)
    result = changes_exonic_splice_site(
        transcript_offset=2,
        transcript=transcript,
        transcript_ref="A",
        transcript_alt="T",
        exon_start_offset=2,
        exon_end_offset=5,
        exon_number=2)
    assert result is True


def test_changes_exonic_splice_site_preserving_acceptor_purine_substitution():
    transcript = DummyTranscript("CCATGGTT", n_exons=2)
    result = changes_exonic_splice_site(
        transcript_offset=2,
        transcript=transcript,
        transcript_ref="A",
        transcript_alt="G",
        exon_start_offset=2,
        exon_end_offset=5,
        exon_number=2)
    assert result is None


def test_changes_exonic_splice_site_breaking_acceptor_deletion():
    # Deleting the first exon base (A) reveals C (non-purine) at exon start.
    transcript = DummyTranscript("CCACGGTT", n_exons=2)
    result = changes_exonic_splice_site(
        transcript_offset=2,
        transcript=transcript,
        transcript_ref="A",
        transcript_alt="",
        exon_start_offset=2,
        exon_end_offset=5,
        exon_number=2)
    assert result is True


def test_changes_exonic_splice_site_preserving_acceptor_deletion():
    # Deleting first exon base (A) reveals G (purine), preserving YAG|R.
    transcript = DummyTranscript("CCAGGGTT", n_exons=2)
    result = changes_exonic_splice_site(
        transcript_offset=2,
        transcript=transcript,
        transcript_ref="A",
        transcript_alt="",
        exon_start_offset=2,
        exon_end_offset=5,
        exon_number=2)
    assert result is None


def test_changes_exonic_splice_site_ignores_variant_before_splice_window():
    # The canonical exon-end splice motif is at offsets 3..5: CAG.
    # Variant is at offset 2, immediately before that window.
    transcript = DummyTranscript("AATCAGAA")
    result = changes_exonic_splice_site(
        transcript_offset=2,
        transcript=transcript,
        transcript_ref="T",
        transcript_alt="C",
        exon_start_offset=0,
        exon_end_offset=5,
        exon_number=1)
    assert result is None
