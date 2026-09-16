"""Known MAP2 RNA haplotypes are sequence observations, not gap mechanisms."""

from types import SimpleNamespace

import pytest

from varcode import MolecularPhaseResolver, RNAReadPhasingSource, Variant

pysam = pytest.importorskip("pysam")

# GRCh38 MAP2, Ensembl ENST00000360351.8; public osteosarc audit (#441).
LEFT = 209694764
REFERENCE = "AAGACCTGGGCTACTGTGTGTTCAATAAGTACACAGTCCCATT"
COMBINED = "AAGACAGGGCCCATT"


@pytest.fixture
def variants(monkeypatch):
    genome = Variant("2", 209694769, "C", "A").genome
    values = [Variant("2", p, r, a, genome=genome) for p, r, a in [
        (209694769, "C", "A"), (209694770, "T", "G"),
        (209694773, "GCTACTGTGTGTTCAATAAGTACACAGT", "")]]
    monkeypatch.setattr("varcode.genome_sequence.reference_range",
                        lambda genome, contig, start, end: REFERENCE[start - LEFT:end - LEFT + 1])
    return values


def read(cigar, sequence=COMBINED, quality=40, name="molecule", start=LEFT):
    r = pysam.AlignedSegment()
    r.query_name, r.query_sequence = name, sequence
    r.reference_id, r.reference_start, r.mapping_quality = 0, start - 1, 60
    r.cigarstring = cigar
    if quality is not None:
        r.query_qualities = [quality] * len(sequence)
    return r


def source(tmp_path, reads):
    path = tmp_path / "map2.bam"
    with pysam.AlignmentFile(path, "wb", header={"SQ": [{"SN": "chr2", "LN": 243199373}]}) as bam:
        for r in sorted(reads, key=lambda r: r.reference_start):
            bam.write(r)
    pysam.index(str(path))
    return RNAReadPhasingSource(str(path), min_alt_reads=1, max_distance_from_read_edge=0)


@pytest.mark.parametrize("cigar", ["9M28D6M", "9M28N6M", "5M22D2M6D8M"])
@pytest.mark.parametrize("reverse", [False, True])
def test_map2_identical_sequence_across_gap_decompositions(tmp_path, variants, cigar, reverse):
    r = read(cigar)
    r.is_reverse = reverse  # BAM sequence remains reference-forward.
    s = source(tmp_path, [r])
    try:
        if "N" in cigar:
            assert s.supports_variant(variants[-1]) == 0
        s.register_haplotype(variants)
        assert all(s.supports_variant(v) == 1 for v in variants)
        resolver = MolecularPhaseResolver(s)
        assert all(resolver.in_cis(variants[-1], v) is True for v in variants[:2])
        assert set(resolver.phased_partners(variants[-1])) == set(variants[:2])
    finally:
        s.close()


@pytest.mark.parametrize("quality", [None, 0, 19])
@pytest.mark.parametrize("cigar", ["9M28D6M", "9M28N6M", "5M22D2M6D8M"])
def test_registered_haplotype_does_not_bypass_quality_filters(tmp_path, variants, quality, cigar):
    s = source(tmp_path, [read(cigar, quality=quality)])
    try:
        s.register_haplotype(variants)
        assert s.supports_variant(variants[-1]) == 0
        assert s.in_cis(variants[0], variants[-1]) is None
    finally:
        s.close()


@pytest.mark.parametrize("r", [
    read("5M19449N10M"),  # genuine large splice skip, not the known local allele
    read("9M28N6M", sequence="AAGACCTGGCCCATT"),  # deletion-only differs
    read("8M28N6M", sequence=COMBINED[1:], start=LEFT + 1),  # left anchor absent
    read("9M28N6M", sequence="TAGACAGGGCCCATT"),  # anchor mismatch
])
def test_no_anchor_or_different_sequence_is_unknown(tmp_path, variants, r):
    s = source(tmp_path, [r])
    try:
        s.register_haplotype(variants)
        assert s.supports_variant(variants[-1]) == 0
        assert s.in_cis(variants[0], variants[-1]) is None
    finally:
        s.close()


def test_no_assembly_of_partial_mates(tmp_path, variants):
    left = read("9M", COMBINED[:9])
    right = read("6M", COMBINED[9:], start=209694801)
    left.flag, right.flag = 99, 147
    left.next_reference_id = right.next_reference_id = 0
    left.next_reference_start, right.next_reference_start = right.reference_start, left.reference_start
    s = source(tmp_path, [left, right])
    try:
        s.register_haplotype(variants)
        assert s.supports_variant(variants[-1]) == 0
        assert s.in_cis(variants[0], variants[-1]) is None
    finally:
        s.close()


def test_read_edge_filter_applies_to_anchors(tmp_path, variants):
    s = source(tmp_path, [read("9M28N6M")])
    try:
        s.max_distance_from_read_edge = 1
        s.register_haplotype(variants)
        assert s.supports_variant(variants[-1]) == 0
    finally:
        s.close()


def test_registration_rejects_invalid_reference_and_edits(tmp_path, variants, monkeypatch):
    s = source(tmp_path, [])
    try:
        for group, flank in [([], 5), (variants, 0), (variants, True), ([variants[0]] * 2, 5)]:
            with pytest.raises(ValueError):
                s.register_haplotype(group, flanking_bases=flank)
        monkeypatch.setattr("varcode.genome_sequence.reference_range", lambda *args: "")
        with pytest.raises(ValueError, match="reference sequence"):
            s.register_haplotype(variants)
    finally:
        s.close()


def test_no_transitive_join_of_registered_haplotypes(tmp_path):
    # Same molecule name with complementary partial observations is insufficient
    # in anchored mode, even when individual matching hypotheses share an allele.
    s = source(tmp_path, [])
    a, b = Variant("2", 100, "A", "T"), Variant("2", 200, "A", "C")
    s._haplotypes = [(frozenset([s._variant_key(a)]), 90, 110, "A")]
    reads = [SimpleNamespace(query_name="same", has_tag=lambda tag: False, which=0),
             SimpleNamespace(query_name="same", has_tag=lambda tag: False, which=1)]
    s._fetch_reads_for_pair = lambda *args: reads
    s._passes_read_filters = lambda read: True
    s._read_allele = lambda read, v: "alt" if (v == a) == (read.which == 0) else None
    try:
        assert s.in_cis(a, b) is None
    finally:
        s.close()
