"""Public contig lookup against a caller's actual reference vocabulary (#645)."""

from dataclasses import FrozenInstanceError

import pytest

from varcode import (
    ContigResolver, StructuralVariant, Variant,
    normalize_contig_name, resolve_contig_name,
)


# UCSC hg19/hg38.chromAlias.txt (primary-source rows checked 2026-10-27).
@pytest.mark.parametrize("reference,ucsc,native", [
    ("hg19", "chrUn_gl000216", "GL000216.1"),
    ("GRCh37", "chr1_gl000191_random", "GL000191.1"),
    ("hg38", "chrUn_KI270438v1", "KI270438.1"),
    ("GRCh38", "chrUn_GL000216v2", "GL000216.2"),
    ("GRCh38", "chr12_GL383550v2_alt", "GL383550.2"),
    ("GRCh38", "chr2_KN538362v1_fix", "KN538362.1"),
])
def test_assembly_specific_aliases_work_in_both_directions(reference, ucsc, native):
    assert normalize_contig_name(ucsc, reference_name=reference) == native
    for supplied, target in [(ucsc, native), (native, ucsc)]:
        result = resolve_contig_name(supplied, [target], reference_name=reference)
        assert result.original_name == supplied
        assert result.normalized_name == native
        assert result.resolved_name == target
        assert result.reference_name == reference
        assert result.status == "alias"
        assert result.candidates == (target,)
        assert result.resolved
        assert result.diagnostic is None


@pytest.mark.parametrize("supplied", ["M", "MT", "chrM", "chrMT", "m", "mt", "chrm", "chrmt", "CHRM", "CHRMT"])
@pytest.mark.parametrize("target", ["M", "MT", "chrM", "chrMT", "m", "mt", "chrm", "chrmt", "CHRM", "CHRMT"])
def test_mitochondrial_spellings_preserve_the_actual_target(supplied, target):
    result = resolve_contig_name(supplied, [target], reference_name="GRCh38")
    assert result.normalized_name == "MT"
    assert result.resolved_name == target
    assert result.status == ("exact" if supplied == target else "alias")
    assert result.original_name == supplied


@pytest.mark.parametrize("name,target,key", [
    (1, "chr1", "1"), (1, "1", "1"),
    ("chrx", "X", "X"), ("X", "chrX", "X"),
    ("Mito", "MITO", "MITO"), ("MTDNA", "mtdna", "MTDNA"),
    ("contig_alpha", "contig_alpha", "contig_alpha"),
])
def test_case_integer_and_custom_names(name, target, key):
    result = resolve_contig_name(name, [target])
    assert result.original_name == name
    assert type(result.original_name) is type(name)
    assert result.normalized_name == key
    assert result.resolved_name == target


@pytest.mark.parametrize("name", ["M", "MT", "chrM"])
def test_exact_match_wins_over_multiple_aliases(name):
    result = resolve_contig_name(name, ["M", "MT", "chrM"])
    assert result.status == "exact"
    assert result.resolved_name == name
    assert result.candidates == (name,)


@pytest.mark.parametrize("names", [
    ["chrM", "MT", "M"], ["M", "MT", "chrM"], {"MT", "chrM", "M"},
])
def test_ambiguity_is_explicit_and_independent_of_iteration_order(names):
    result = resolve_contig_name("chrMT", names, reference_name="GRCh38")
    assert result.status == "ambiguous"
    assert not result.resolved
    assert result.resolved_name is None
    assert result.candidates == ("M", "MT", "chrM")
    assert "Ambiguous" in result.diagnostic
    assert "GRCh38" in result.diagnostic
    assert "chrMT" in result.diagnostic


def test_ambiguous_scaffold_and_case_aliases():
    result = resolve_contig_name("GL000216.1", ["chrUn_gl000216", "chrUn_GL000216v1"],
                                 reference_name="hg19")
    assert result.status == "ambiguous"
    assert result.resolved_name is None
    assert resolve_contig_name("mito", ["Mito", "MITO"]).status == "ambiguous"


@pytest.mark.parametrize("name,available,reference", [
    ("chrUn_gl000216", ["GL000216.2"], "hg19"),
    ("chrUn_gl000216", ["GL000216.1"], "GRCh38"),
    ("chrUn_gl000216", ["GL000216.1"], None),
    ("chrUn_KI270438v1", ["KI270438.2"], "GRCh38"),
    ("GL000216.1", ["GL000216.2"], "GRCh38"),
    ("GL000216", ["GL000216.1"], "hg19"),
    ("NC_012920.1", ["MT"], "GRCh38"),
    ("chr6_apd_hap1", ["6"], "hg19"),
    ("chrUn_unknown", ["unknown"], "GRCh38"),
    ("decoy_alpha", ["decoy_beta"], "custom"),
    ("2", ["1"], "GRCh38"),
])
def test_missing_names_do_not_infer_versions_assemblies_or_sequence_identity(name, available, reference):
    result = resolve_contig_name(name, available, reference_name=reference)
    assert result.status == "missing"
    assert result.resolved_name is None
    assert not result.resolved
    assert result.candidates == ()
    assert "absent" in result.diagnostic
    assert repr(name) in result.diagnostic


@pytest.mark.parametrize("reference", ["GRCh37", "hg19", "b37", "NCBI37", " HG19 "])
def test_reference_aliases_only_select_the_unversioned_spelling_rule(reference):
    assert normalize_contig_name("chrUn_gl000216", reference_name=reference) == "GL000216.1"


def test_unavailable_and_empty_reference_vocabularies_are_distinct():
    unavailable = resolve_contig_name("chr1", None, reference_name="GRCh38")
    missing = resolve_contig_name("chr1", [], reference_name="GRCh38")
    assert unavailable.status == "reference_unavailable"
    assert "unavailable" in unavailable.diagnostic
    assert missing.status == "missing"
    assert "absent" in missing.diagnostic
    for result in (unavailable, missing):
        assert result.normalized_name == "1"
        assert result.resolved_name is None
        assert not result.resolved


def test_resolver_consumes_a_vocabulary_once_and_does_not_cache_by_assembly():
    resolver = ContigResolver(iter(["chr1", "chrM", "chrM"]), reference_name="GRCh38")
    assert resolver.resolve("1").resolved_name == "chr1"
    assert resolver.resolve("MT").resolved_name == "chrM"
    assert resolver.resolve("MT").status == "alias"
    assert resolver.resolve("2").status == "missing"
    subset = ContigResolver(["2"], reference_name="GRCh38")
    assert subset.resolve("2").status == "exact"
    assert subset.resolve("1").status == "missing"
    assert resolver.available_contigs == frozenset(["chr1", "chrM"])


def test_resolver_snapshots_names_and_results_are_immutable():
    names = ["chrM"]
    resolver = ContigResolver(names)
    names.append("M")
    result = resolver.resolve("MT")
    assert result.status == "alias"
    assert result.resolved_name == "chrM"
    with pytest.raises(FrozenInstanceError):
        result.resolved_name = "M"


@pytest.mark.parametrize("name", [None, 1.2, b"1", [], ""])
def test_invalid_query_names_fail_explicitly(name):
    with pytest.raises((TypeError, ValueError)):
        resolve_contig_name(name, ["1"])


@pytest.mark.parametrize("names", ["chr1", b"chr1", [1], [None], [""]])
def test_invalid_vocabularies_fail_explicitly(names):
    with pytest.raises((TypeError, ValueError)):
        ContigResolver(names)


@pytest.mark.parametrize("reference", [1, [], "", " "])
def test_invalid_reference_context_fails_explicitly(reference):
    with pytest.raises((TypeError, ValueError)):
        ContigResolver(["1"], reference_name=reference)
    with pytest.raises((TypeError, ValueError)):
        normalize_contig_name("chr1", reference_name=reference)


def test_normalization_and_variant_flags_remain_compatible():
    assert normalize_contig_name("chrm", convert_ucsc=False) == "chrM"
    assert normalize_contig_name("chrm") == "MT"
    for normalize, convert, expected in [
        (True, True, "MT"), (True, False, "chrM"),
        (False, True, "m"), (False, False, "chrm"),
    ]:
        variant = Variant("chrm", 1, "A", "G", "GRCh38",
                          normalize_contig_names=normalize, convert_ucsc_contig_names=convert)
        assert variant.contig == expected
        assert variant.original_contig == "chrm"
        assert Variant.from_json(variant.to_json()).contig == expected
    assert Variant(1, 1, "A", "G", normalize_contig_names=False).contig == 1


def test_structural_breakend_mates_share_public_normalization():
    variant = StructuralVariant(
        "chr1", 10, "BND", ref="A", alt="A[chrUn_KI270438v1:20[",
        mate_contig="KI270438.1", mate_start=20,
        genome="GRCh38", convert_ucsc_contig_names=True)
    assert variant.contig == "1"
    assert variant.mate_contig == "KI270438.1"
    assert variant.original_contig == "chr1"
    assert "chrUn_KI270438v1" in variant.alt
    assert resolve_contig_name(variant.mate_contig, ["chrUn_KI270438v1"]).resolved


def test_resolved_names_can_fetch_from_an_actual_bam(tmp_path):
    pysam = pytest.importorskip("pysam")
    path = tmp_path / "aliases.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [
        {"SN": "chrM", "LN": 100}, {"SN": "chrUn_KI270438v1", "LN": 100},
    ]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as output:
        for reference_id, read_name in enumerate(["mitochondrial_read", "scaffold_read"]):
            read = pysam.AlignedSegment(output.header)
            read.query_name = read_name
            read.query_sequence = "ACGTACGTAC"
            read.reference_id = reference_id
            read.reference_start = 5
            read.mapping_quality = 60
            read.cigarstring = "10M"
            output.write(read)
    pysam.index(str(path))
    with pysam.AlignmentFile(str(path), "rb") as bam:
        resolver = ContigResolver(bam.references, reference_name="GRCh38")
        for query, expected in [("MT", "mitochondrial_read"), ("KI270438.1", "scaffold_read")]:
            result = resolver.resolve(query)
            assert result.resolved_name is not None
            assert [read.query_name for read in bam.fetch(result.resolved_name, 5, 15)] == [expected]
        assert resolver.resolve("KI270438.2").status == "missing"


@pytest.mark.parametrize("query,names,expected", [
    ("chr1", ["1", "CHR1"], ("1", "CHR1")),
    ("MT", ["CHRM", "chrMT"], ("CHRM", "chrMT")),
    ("x", ["X", "CHRX"], ("CHRX", "X")),
])
def test_canonical_case_aliases_preserve_ambiguity_and_exact_names(query, names, expected):
    result = ContigResolver(names).resolve(query)
    assert result.status == "ambiguous"
    assert result.candidates == expected
    for name in names:
        result = ContigResolver(names).resolve(name)
        assert result.status == "exact"
        assert result.resolved_name == name


def test_custom_reference_vocabulary_does_not_poison_unrelated_lookups(tmp_path):
    pysam = pytest.importorskip("pysam")
    names = ["0", "HPV16", "NC_007605.1", "hs37d5", "HLA-A*01:01:01:01", "chr6_apd_hap1"]
    path = tmp_path / "custom.bam"
    with pysam.AlignmentFile(path, "wb", header={"SQ": [{"SN": name, "LN": 100} for name in names]}) as bam:
        for index, name in enumerate(names):
            read = pysam.AlignedSegment(bam.header)
            read.query_name = "witness_" + str(index)
            read.query_sequence = "ACGT"
            read.reference_id = index
            read.reference_start = 10
            read.cigarstring = "4M"
            bam.write(read)
    pysam.index(str(path))
    with pysam.AlignmentFile(path) as bam:
        resolver = ContigResolver(bam.references, reference_name="custom")
        for index, name in enumerate(names):
            result = resolver.resolve(name)
            assert result.status == "exact"
            assert [read.query_name for read in bam.fetch(result.resolved_name, 10, 14)] == ["witness_" + str(index)]
        assert resolver.resolve(0).resolved_name == "0"
        assert resolver.resolve("EBV").status == "missing"
        assert resolver.resolve("NC_007605.2").status == "missing"
        assert resolver.resolve("6").status == "missing"
    assert ContigResolver(["1"]).resolve("0").status == "missing"
    with pytest.raises(ValueError, match="cannot be 0"):
        normalize_contig_name("0")
    with pytest.raises(ValueError, match="cannot be 0"):
        Variant("0", 1, "A", "G", "GRCh38")
