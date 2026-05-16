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

import pytest

from varcode import (
    MolecularPhaseResolver,
    ReadPhaseResolver,
    ReadPhasingSource,
    RNAReadPhasingSource,
    Variant,
)

pysam = pytest.importorskip("pysam")


def _quality_array(length, low_at=None):
    qualities = [40] * length
    if low_at is not None:
        qualities[low_at] = 5
    return qualities


def _aligned_segment(
        name,
        start0,
        sequence,
        *,
        cigar=None,
        qualities=None,
        flag=0,
        mate_start0=0,
        template_length=0,
        mapping_quality=60):
    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = sequence
    read.flag = flag
    read.reference_id = 0
    read.reference_start = start0
    read.mapping_quality = mapping_quality
    read.cigar = cigar or [(0, len(sequence))]
    read.next_reference_id = 0
    read.next_reference_start = mate_start0
    read.template_length = template_length
    read.query_qualities = qualities or _quality_array(len(sequence))
    return read


def _write_bam(tmp_path, reads):
    bam_path = tmp_path / "reads.bam"
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": "chr1", "LN": 1000}],
    }
    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        for read in sorted(reads, key=lambda r: r.reference_start):
            bam.write(read)
    pysam.index(str(bam_path))
    return str(bam_path)


def _snv_read(name, start1, calls, *, low_quality_pos=None, flag=0):
    sequence = ["A"] * 80
    for pos1, base in calls.items():
        sequence[pos1 - start1] = base
    return _aligned_segment(
        name=name,
        start0=start1 - 1,
        sequence="".join(sequence),
        qualities=_quality_array(
            80,
            None if low_quality_pos is None else low_quality_pos - start1),
        flag=flag,
    )


def test_rna_read_phasing_source_satisfies_protocol(tmp_path):
    bam_path = _write_bam(tmp_path, [])
    source = RNAReadPhasingSource(bam_path)
    try:
        assert isinstance(source, ReadPhasingSource)
        assert MolecularPhaseResolver(source).phase_source == "rna_reads"
        assert ReadPhaseResolver(source).phase_source == "rna_reads"
    finally:
        source.close()


def test_supports_variant_counts_alt_fragments_with_contig_normalization(tmp_path):
    variant = Variant("1", 120, "A", "T")
    bam_path = _write_bam(tmp_path, [
        _snv_read("r1", 100, {120: "T"}),
        _snv_read("r2", 100, {120: "T"}),
        _snv_read("r3", 100, {120: "A"}),
    ])
    source = RNAReadPhasingSource(
        bam_path,
        max_distance_from_read_edge=0,
    )
    try:
        assert source.supports_variant(variant) == 2
        assert source.has_evidence(variant) is True
    finally:
        source.close()


def test_read_phase_resolver_reports_cis_from_same_read_alt_cooccurrence(tmp_path):
    v1 = Variant("1", 120, "A", "T")
    v2 = Variant("1", 140, "A", "G")
    bam_path = _write_bam(tmp_path, [
        _snv_read("r1", 100, {120: "T", 140: "G"}),
        _snv_read("r2", 100, {120: "T", 140: "G"}),
    ])
    source = RNAReadPhasingSource(
        bam_path,
        max_distance_from_read_edge=0,
    )
    try:
        resolver = MolecularPhaseResolver(source)
        assert resolver.in_cis(v1, v2) is True
        assert v2 in resolver.phased_partners(v1)
    finally:
        source.close()


def test_read_phase_resolver_reports_trans_from_mixed_reads(tmp_path):
    v1 = Variant("1", 120, "A", "T")
    v2 = Variant("1", 140, "A", "G")
    bam_path = _write_bam(tmp_path, [
        _snv_read("r1", 100, {120: "T", 140: "A"}),
        _snv_read("r2", 100, {120: "T", 140: "A"}),
    ])
    source = RNAReadPhasingSource(
        bam_path,
        max_distance_from_read_edge=0,
    )
    try:
        resolver = MolecularPhaseResolver(source)
        assert resolver.in_cis(v1, v2) is False
    finally:
        source.close()


def test_read_phase_resolver_returns_none_below_threshold(tmp_path):
    v1 = Variant("1", 120, "A", "T")
    v2 = Variant("1", 140, "A", "G")
    bam_path = _write_bam(tmp_path, [
        _snv_read("r1", 100, {120: "T", 140: "G"}),
    ])
    source = RNAReadPhasingSource(
        bam_path,
        max_distance_from_read_edge=0,
        min_alt_reads=2,
    )
    try:
        resolver = MolecularPhaseResolver(source)
        assert resolver.in_cis(v1, v2) is None
    finally:
        source.close()


def test_paired_end_fragment_can_support_cis(tmp_path):
    v1 = Variant("1", 120, "A", "T")
    v2 = Variant("1", 320, "A", "G")
    mate1 = _snv_read("fragment1", 100, {120: "T"}, flag=99)
    mate1.next_reference_start = 299
    mate1.template_length = 280
    mate2 = _snv_read("fragment1", 300, {320: "G"}, flag=147)
    mate2.next_reference_start = 99
    mate2.template_length = -280
    bam_path = _write_bam(tmp_path, [mate1, mate2])
    source = RNAReadPhasingSource(
        bam_path,
        max_distance_from_read_edge=0,
        min_alt_reads=1,
    )
    try:
        resolver = MolecularPhaseResolver(source)
        assert resolver.in_cis(v1, v2) is True
    finally:
        source.close()


def test_low_quality_alt_base_is_discounted(tmp_path):
    variant = Variant("1", 120, "A", "T")
    bam_path = _write_bam(tmp_path, [
        _snv_read("r1", 100, {120: "T"}, low_quality_pos=120),
        _snv_read("r2", 100, {120: "T"}),
    ])
    source = RNAReadPhasingSource(
        bam_path,
        max_distance_from_read_edge=0,
        min_alt_reads=2,
    )
    try:
        assert source.supports_variant(variant) == 1
        assert source.has_evidence(variant) is False
    finally:
        source.close()


def test_insertion_support_uses_inserted_cigar_bases(tmp_path):
    variant = Variant("1", 119, "", "T")
    sequence = "A" * 20 + "T" + "A" * 40
    bam_path = _write_bam(tmp_path, [
        _aligned_segment("r1", 99, sequence, cigar=[(0, 20), (1, 1), (0, 40)]),
        _aligned_segment("r2", 99, sequence, cigar=[(0, 20), (1, 1), (0, 40)]),
    ])
    source = RNAReadPhasingSource(
        bam_path,
        max_distance_from_read_edge=0,
    )
    try:
        assert source.supports_variant(variant) == 2
    finally:
        source.close()


def test_larger_insertion_does_not_support_smaller_insertion(tmp_path):
    variant = Variant("1", 119, "", "T")
    sequence = "A" * 20 + "TT" + "A" * 40
    bam_path = _write_bam(tmp_path, [
        _aligned_segment("r1", 99, sequence, cigar=[(0, 20), (1, 2), (0, 40)]),
        _aligned_segment("r2", 99, sequence, cigar=[(0, 20), (1, 2), (0, 40)]),
    ])
    source = RNAReadPhasingSource(
        bam_path,
        max_distance_from_read_edge=0,
        min_alt_reads=1,
    )
    try:
        assert source.supports_variant(variant) == 0
        assert source.has_evidence(variant) is False
    finally:
        source.close()


def test_shifted_insertion_does_not_support_requested_insertion(tmp_path):
    variant = Variant("1", 119, "", "T")
    sequence = "A" * 21 + "T" + "A" * 39
    bam_path = _write_bam(tmp_path, [
        _aligned_segment("r1", 99, sequence, cigar=[(0, 21), (1, 1), (0, 39)]),
        _aligned_segment("r2", 99, sequence, cigar=[(0, 21), (1, 1), (0, 39)]),
    ])
    source = RNAReadPhasingSource(
        bam_path,
        max_distance_from_read_edge=0,
        min_alt_reads=1,
    )
    try:
        assert source.supports_variant(variant) == 0
        assert source.has_evidence(variant) is False
    finally:
        source.close()


def test_deletion_support_uses_deleted_cigar_span(tmp_path):
    variant = Variant("1", 120, "A", "")
    sequence = "A" * 59
    bam_path = _write_bam(tmp_path, [
        _aligned_segment("r1", 99, sequence, cigar=[(0, 20), (2, 1), (0, 39)]),
        _aligned_segment("r2", 99, sequence, cigar=[(0, 20), (2, 1), (0, 39)]),
    ])
    source = RNAReadPhasingSource(
        bam_path,
        max_distance_from_read_edge=0,
    )
    try:
        assert source.supports_variant(variant) == 2
    finally:
        source.close()


def test_multibase_deletion_support_requires_exact_cigar_span(tmp_path):
    variant = Variant("1", 120, "AAA", "")
    sequence = "A" * 57
    bam_path = _write_bam(tmp_path, [
        _aligned_segment("r1", 99, sequence, cigar=[(0, 20), (2, 3), (0, 37)]),
        _aligned_segment("r2", 99, sequence, cigar=[(0, 20), (2, 3), (0, 37)]),
    ])
    source = RNAReadPhasingSource(
        bam_path,
        max_distance_from_read_edge=0,
    )
    try:
        assert source.supports_variant(variant) == 2
        assert source.has_evidence(variant) is True
    finally:
        source.close()


def test_larger_deletion_does_not_support_smaller_deletion(tmp_path):
    variant = Variant("1", 120, "A", "")
    sequence = "A" * 57
    bam_path = _write_bam(tmp_path, [
        _aligned_segment("r1", 99, sequence, cigar=[(0, 20), (2, 3), (0, 37)]),
        _aligned_segment("r2", 99, sequence, cigar=[(0, 20), (2, 3), (0, 37)]),
    ])
    source = RNAReadPhasingSource(
        bam_path,
        max_distance_from_read_edge=0,
        min_alt_reads=1,
    )
    try:
        assert source.supports_variant(variant) == 0
        assert source.has_evidence(variant) is False
    finally:
        source.close()


def test_smaller_deletion_does_not_support_larger_deletion(tmp_path):
    variant = Variant("1", 120, "AAA", "")
    sequence = "A" * 59
    bam_path = _write_bam(tmp_path, [
        _aligned_segment("r1", 99, sequence, cigar=[(0, 20), (2, 1), (0, 39)]),
        _aligned_segment("r2", 99, sequence, cigar=[(0, 20), (2, 1), (0, 39)]),
    ])
    source = RNAReadPhasingSource(
        bam_path,
        max_distance_from_read_edge=0,
        min_alt_reads=1,
    )
    try:
        assert source.supports_variant(variant) == 0
        assert source.has_evidence(variant) is False
    finally:
        source.close()


def test_shifted_deletion_does_not_support_requested_deletion(tmp_path):
    variant = Variant("1", 120, "A", "")
    sequence = "A" * 59
    bam_path = _write_bam(tmp_path, [
        _aligned_segment("r1", 99, sequence, cigar=[(0, 21), (2, 1), (0, 38)]),
        _aligned_segment("r2", 99, sequence, cigar=[(0, 21), (2, 1), (0, 38)]),
    ])
    source = RNAReadPhasingSource(
        bam_path,
        max_distance_from_read_edge=0,
        min_alt_reads=1,
    )
    try:
        assert source.supports_variant(variant) == 0
        assert source.has_evidence(variant) is False
    finally:
        source.close()


def test_deletion_phasing_ignores_larger_overlapping_deletion(tmp_path):
    deletion = Variant("1", 120, "A", "")
    snv = Variant("1", 150, "A", "T")
    sequence = list("A" * 87)
    sequence[47] = "T"
    bam_path = _write_bam(tmp_path, [
        _aligned_segment(
            "r1",
            99,
            "".join(sequence),
            cigar=[(0, 20), (2, 3), (0, 67)]),
        _aligned_segment(
            "r2",
            99,
            "".join(sequence),
            cigar=[(0, 20), (2, 3), (0, 67)]),
    ])
    source = RNAReadPhasingSource(
        bam_path,
        max_distance_from_read_edge=0,
        min_alt_reads=1,
    )
    try:
        resolver = MolecularPhaseResolver(source)
        assert source.supports_variant(deletion) == 0
        assert resolver.in_cis(deletion, snv) is None
    finally:
        source.close()


def test_reference_skip_does_not_support_deletion(tmp_path):
    variant = Variant("1", 120, "A", "")
    sequence = "A" * 40
    bam_path = _write_bam(tmp_path, [
        _aligned_segment("r1", 99, sequence, cigar=[(0, 20), (3, 50), (0, 20)]),
        _aligned_segment("r2", 99, sequence, cigar=[(0, 20), (3, 50), (0, 20)]),
    ])
    source = RNAReadPhasingSource(
        bam_path,
        max_distance_from_read_edge=0,
        min_alt_reads=1,
    )
    try:
        assert source.supports_variant(variant) == 0
        assert source.has_evidence(variant) is False
    finally:
        source.close()


def test_soft_clip_after_anchor_does_not_support_insertion(tmp_path):
    variant = Variant("1", 119, "", "TTTT")
    sequence = "A" * 20 + "TTTT"
    bam_path = _write_bam(tmp_path, [
        _aligned_segment("r1", 99, sequence, cigar=[(0, 20), (4, 4)]),
        _aligned_segment("r2", 99, sequence, cigar=[(0, 20), (4, 4)]),
    ])
    source = RNAReadPhasingSource(
        bam_path,
        max_distance_from_read_edge=0,
        min_alt_reads=1,
    )
    try:
        assert source.supports_variant(variant) == 0
        assert source.has_evidence(variant) is False
    finally:
        source.close()


def test_splice_junction_skip_does_not_support_intronic_variant(tmp_path):
    variant = Variant("1", 130, "A", "T")
    sequence = "A" * 40
    bam_path = _write_bam(tmp_path, [
        _aligned_segment("r1", 99, sequence, cigar=[(0, 20), (3, 50), (0, 20)]),
    ])
    source = RNAReadPhasingSource(
        bam_path,
        max_distance_from_read_edge=0,
        min_alt_reads=1,
    )
    try:
        assert source.supports_variant(variant) == 0
        assert source.has_evidence(variant) is False
    finally:
        source.close()
