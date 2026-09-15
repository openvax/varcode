"""Tests for the coordinate-aware realized-genome layout."""

from types import SimpleNamespace

import pytest

from varcode.genomic_layout import (
    GenomicLayout,
    LayoutSegment,
    SequenceUnavailable,
    join_breakends,
    project_intervals,
)
from varcode.structural_variant import Breakend


SEQUENCES = {
    "a": "AACCGGTTAACC",
    "b": "TTGGCCAATTGG",
}


def provider(contig, start, end):
    return SEQUENCES[contig][start - 1:end]


def point(contig, start, ref, alt):
    return SimpleNamespace(
        contig=contig,
        trimmed_base1_start=start,
        trimmed_base1_end=start + len(ref) - 1 if ref else start,
        trimmed_ref=ref,
        trimmed_alt=alt,
        is_structural=False)


def structural(kind, start, end, affected_start=None, affected_end=None,
               alt_assembly=None):
    return SimpleNamespace(
        contig="a",
        start=start,
        end=end,
        affected_start=(start if affected_start is None else affected_start),
        affected_end=(end if affected_end is None else affected_end),
        sv_type=kind,
        alt_assembly=alt_assembly,
        is_structural=True)


def test_layout_is_lazy_until_materialized():
    layout = GenomicLayout.from_interval("a", 1, 12)
    assert layout.length == 12
    with pytest.raises(SequenceUnavailable):
        layout.materialize()


def test_point_variants_preserve_origins_on_both_strands():
    forward = GenomicLayout.from_interval(
        "a", 1, 12, sequence_provider=provider)
    reverse = GenomicLayout.from_interval(
        "a", 1, 12, strand="-", sequence_provider=provider)
    variant = point("a", 3, "C", "T")

    forward_mutant = forward.apply_point_variant(variant)
    reverse_mutant = reverse.apply_point_variant(variant)

    assert forward_mutant.materialize() == "AATCGGTTAACC"
    assert reverse_mutant.materialize() == "GGTTAACCGATT"
    assert forward_mutant.origins()[2] == ("a", 3, "alternate")
    assert reverse_mutant.origins()[-3] == ("a", 3, "alternate")


def test_indels_have_explicit_inserted_and_deleted_origins():
    layout = GenomicLayout.from_interval(
        "a", 1, 12, sequence_provider=provider)
    inserted = layout.apply_point_variant(point("a", 4, "", "TA"))
    deleted = layout.apply_point_variant(point("a", 3, "CC", ""))

    assert inserted.materialize() == "AACCTA GGTTAACC".replace(" ", "")
    assert [origin[1] for origin in inserted.origins()][4:6] == [None, None]
    assert deleted.materialize() == "AAGGTTAACC"
    assert 3 not in [origin[1] for origin in deleted.origins()]
    assert 4 not in [origin[1] for origin in deleted.origins()]


def test_typed_deletion_does_not_remove_retained_padding_anchor():
    layout = GenomicLayout.from_interval(
        "a", 1, 12, sequence_provider=provider)
    deletion = structural(
        "DEL", start=2, end=6, affected_start=3, affected_end=6)

    mutant = layout.apply_structural_variant(deletion)

    assert mutant.materialize() == "AATTAACC"
    assert [origin[1] for origin in mutant.origins()][:3] == [1, 2, 7]


def test_duplication_and_inversion_work_in_transcript_order():
    forward = GenomicLayout.from_interval(
        "a", 1, 12, sequence_provider=provider)
    reverse = GenomicLayout.from_interval(
        "a", 1, 12, strand="-", sequence_provider=provider)

    forward_dup = forward.apply_structural_variant(structural("DUP", 3, 6))
    reverse_dup = reverse.apply_structural_variant(structural("DUP", 3, 6))
    inverted = forward.apply_structural_variant(structural("INV", 2, 5))

    assert forward_dup.materialize() == "AACCGGCCGGTTAACC"
    assert reverse_dup.materialize() == "GGTTAACCGGCCGGTT"
    assert inverted.materialize() == "ACGGTGTTAACC"
    assert {kind for _, _, kind in inverted.origins()[1:5]} == {"inverted"}


def test_apply_variants_deletes_an_earlier_point_edit_on_same_haplotype():
    layout = GenomicLayout.from_interval(
        "a", 1, 12, sequence_provider=provider)
    variants = [
        structural("DEL", 2, 6, affected_start=3, affected_end=6),
        point("a", 4, "C", "T"),
    ]

    mutant = layout.apply_variants(variants)

    assert mutant.materialize() == "AATTAACC"
    assert all(kind != "alternate" for _, _, kind in mutant.origins())


def test_cross_contig_breakend_join_orients_retained_sides():
    left = GenomicLayout.from_interval(
        "a", 1, 12, sequence_provider=provider)
    right = GenomicLayout.from_interval(
        "b", 1, 12, strand="-", sequence_provider=provider)

    fusion = join_breakends(
        left,
        Breakend("a", 3, "left"),
        right,
        Breakend("b", 8, "right"))

    assert fusion.materialize() == "AACATTGG"
    assert fusion.endpoint_origin(first=False) == ("b", 12)
    assert fusion.origins()[2][:2] == ("a", 3)
    assert fusion.origins()[3][:2] == ("b", 8)


def test_project_intervals_keeps_duplicate_exon_occurrences():
    layout = GenomicLayout.from_interval(
        "a", 1, 12, sequence_provider=provider)
    duplicated = layout.apply_structural_variant(structural("DUP", 3, 6))

    runs = project_intervals(duplicated, (("a", 3, 6),))

    assert [run.materialize() for run in runs] == ["CCGG", "CCGG"]
    assert runs[0].segments[0].origin_kind == "reference"
    assert runs[1].segments[0].origin_kind == "duplicated"


def test_inserted_segment_requires_sequence():
    with pytest.raises(ValueError):
        LayoutSegment(None, None, None)
