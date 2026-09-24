"""Synthetic VCF construction and independently specified sample identities.

No downloaded fixtures: the genotype patterns, LOH edits and tumor branches
used as ground truth are visible here alongside their expected comparisons.
"""

from dataclasses import replace
import gzip
import json
from pathlib import Path
import subprocess
import sys

import pytest

from varcode import SampleCheckConfig, SampleSpec, check_sample_identity, compare_samples, load_vcf_samples
from varcode.cli.main import main


def _vcf(tmp_path, name, sample_names, rows, *, reference="GRCh37", metadata=""):
    path = tmp_path / name
    text = (
        '##fileformat=VCFv4.2\n'
        + ('##reference=%s\n' % reference if reference else '')
        + metadata
        + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n'
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction">\n'
        '##FORMAT=<ID=FT,Number=1,Type=String,Description="Sample filter">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
        + '\t'.join(sample_names) + '\n'
        + ''.join('\t'.join(map(str, row)) + '\n' for row in rows))
    if name.endswith('.gz'):
        with gzip.open(path, 'wt') as out:
            out.write(text)
    else:
        path.write_text(text)
    return path


def _row(i, *calls, alt="C", fmt="GT:DP:GQ", info=".", chrom="1"):
    return [chrom, 1 + i * 200_000, '.', 'A', alt, '.', 'PASS', info, fmt, *calls]


def _donor(tmp_path, name="donor.vcf", *, swapped=False, loh=False, role="germline", **kwargs):
    genotypes = ["0/0", "0/1", "1/1"]
    if swapped:
        genotypes = ["1/1", "0/1", "0/0"]
    if loh:
        genotypes = ["0/0", "1/1", "1/1"]
    path = _vcf(tmp_path, name, ["S"], [_row(i, genotypes[i % 3] + ":30:60") for i in range(210)])
    return load_vcf_samples(path, samples=[SampleSpec("S", role=role, **kwargs)])[0]


def test_same_donor_and_swapped_labels(tmp_path):
    a = _donor(tmp_path, "normal-a.vcf", donor_id="A")
    b = _donor(tmp_path, "normal-b.vcf", swapped=True, donor_id="B")
    mislabeled = _donor(tmp_path, "tumor-b.vcf", role="tumor", donor_id="B")
    report = check_sample_identity([a, b, mislabeled])
    ab, at, bt = report["pairs"]
    assert ab["germline"]["status"] == "discordant"
    assert ab["germline"]["ibs0"] == 140
    assert at["germline"]["status"] == "compatible"
    assert at["flags"] == ["unexpected_donor_compatibility"]
    assert bt["flags"] == ["expected_donor_mismatch"]
    assert report["identity_conflicts"] == 1
    assert report["config"]["min_shared_snps"] == 100
    assert len(report["samples"][0]["sha256"]) == 64
    assert "snps" not in report["samples"][0]
    json.dumps(report, allow_nan=False)


def test_tumor_loh_is_tolerated_but_exact_normal_comparison_is_stricter(tmp_path):
    a = _donor(tmp_path)
    tumor = _donor(tmp_path, "loh.vcf", loh=True, role="tumor")
    result = compare_samples(a, tumor)["germline"]
    assert result["status"] == "compatible"
    assert result["genotype_concordance"] == pytest.approx(2 / 3)
    assert result["ibs0"] == 0
    tumor.role = "normal"
    assert compare_samples(a, tumor)["germline"]["status"] == "discordant"


def test_variant_only_tumors_cannot_match_donors_by_automatic_alt_sharing(tmp_path):
    path = _vcf(tmp_path, "ascertained.vcf", ["a", "b"], [
        _row(i, ["0/1:30:60", "1/1:30:60"][i % 2],
             ["1/1:30:60", "0/1:30:60"][i % 2]) for i in range(210)])
    a, b = load_vcf_samples(path, kind="somatic", samples=[SampleSpec("a", role="tumor"), SampleSpec("b", role="tumor")])
    result = compare_samples(a, b)["germline"]
    assert result["ibs0"] == 0
    assert result["status"] == "inconclusive"
    assert "too_few_explicit_tumor_comparison_reference_calls" in result["reasons"]


@pytest.mark.parametrize("pattern", ["0/0", "0/1", "1/1", "./.", "0/."])
def test_uninformative_or_partial_calls_are_inconclusive(tmp_path, pattern):
    path = _vcf(tmp_path, "sparse.vcf", ["a", "b"], [_row(i, pattern + ":30:60", pattern + ":30:60") for i in range(210)])
    a, b = load_vcf_samples(path)
    result = compare_samples(a, b)["germline"]
    assert result["status"] == "inconclusive"
    assert result["reasons"]


def test_two_alt_heterozygotes_do_not_fake_genotype_diversity(tmp_path):
    rows = [_row(i, ["0/1:30:60", "1/2:30:60"][i % 2], alt="C,G") for i in range(210)]
    sample = load_vcf_samples(_vcf(tmp_path, "heterozygotes.vcf", ["S"], rows))[0]
    result = compare_samples(sample, sample)["germline"]
    assert result["status"] == "inconclusive"
    assert result["reasons"] == ["insufficient_genotype_diversity"]


def test_absent_records_not_reference_and_linked_sites_not_independent(tmp_path):
    a = _donor(tmp_path)
    b = load_vcf_samples(_vcf(tmp_path, "few.vcf", ["S"], [_row(1, "0/1:30:60")]))[0]
    assert compare_samples(a, b)["germline"]["shared_snps"] == 1
    rows = [_row(i, ["0/0:30:60", "0/1:30:60", "1/1:30:60"][i % 3]) for i in range(210)]
    for i, row in enumerate(rows):
        row[1] = 100 + i
    b = load_vcf_samples(_vcf(tmp_path, "linked.vcf", ["S"], rows))[0]
    result = compare_samples(b, b)["germline"]
    assert result["status"] == "inconclusive"
    assert result["genomic_bins"] == 1


def test_alt_indexes_phase_and_chr_aliases_do_not_define_identity(tmp_path):
    one = _vcf(tmp_path, "one.vcf", ["S"], [_row(1, "0|2:30:60", alt="C,G")])
    two = _vcf(tmp_path, "two.vcf.gz", ["S"], [_row(1, "1/0:30:60", alt="G,C", chrom="chr1")], reference="hg19")
    a, b = load_vcf_samples(one)[0], load_vcf_samples(two)[0]
    assert a.snps == b.snps == {("1", 200001): ("A", ("A", "G"))}
    assert compare_samples(a, b)["germline"]["genotype_concordance"] == 1


def test_gvcf_blocks_and_unknown_alleles_are_not_imputed(tmp_path):
    rows = [_row(1, "0/1:30:60", alt="C,<NON_REF>"),
            _row(2, "0/2:30:60", alt="C,<NON_REF>"),
            _row(3, "0/0:30:60", alt="<NON_REF>", info="END=800000"),
            _row(4, "0/0:30:60", alt=".")]
    s = load_vcf_samples(_vcf(tmp_path, "gvcf.vcf", ["S"], rows))[0]
    assert s.snps == {("1", 200001): ("A", ("A", "C")), ("1", 800001): ("A", ("A", "A"))}


def test_filters_missing_quality_and_conflicting_duplicates(tmp_path):
    rows = [_row(1, "0/1:30:60"), _row(1, "1/1:30:60"), _row(1, "0/1:30:60"),
            _row(2, "0/1:2:60"), _row(3, "0/1:30:2"), _row(4, "0/1", fmt="GT"),
            _row(5, "0/1:LowQual", fmt="GT:FT")]
    rows.append(_row(6, "0/1:30:60"))
    rows[-1][6] = "LowQual"
    path = _vcf(tmp_path, "quality.vcf", ["S"], rows)
    a = load_vcf_samples(path)[0]
    assert a.snps == {("1", 200001): None, ("1", 800001): ("A", ("A", "C"))}
    assert a.qc["conflicting_snp_loci"] == 1
    assert a.qc["low_depth"] == a.qc["low_gq"] == a.qc["record_filter"] == a.qc["sample_filter"] == 1
    strict = load_vcf_samples(path, config=SampleCheckConfig(require_quality=True))[0]
    assert strict.qc["missing_quality"] == 1
    assert strict.snps == {("1", 200001): None}


def test_genotype_independent_bounded_sketch_is_order_independent(tmp_path):
    rows = [_row(i, ["0/0:30:60", "0/1:30:60", "1/1:30:60"][i % 3]) for i in range(300)]
    config = SampleCheckConfig(max_snps=100)
    a = load_vcf_samples(_vcf(tmp_path, "forward.vcf", ["S"], rows), config=config)[0]
    b = load_vcf_samples(_vcf(tmp_path, "reverse.vcf", ["S"], reversed(rows)), config=config)[0]
    assert a.snps == b.snps
    assert len(a.snps) == 100
    assert a.qc["snp_sketch_omissions"] == b.qc["snp_sketch_omissions"] == 200
    assert compare_samples(a, b)["germline"]["status"] == "compatible"


def test_reference_conflicts_and_assembly_mismatch(tmp_path):
    a = _donor(tmp_path)
    b = replace(a, assembly="GRCh38")
    pair = compare_samples(a, b)
    assert pair["germline"]["reasons"] == pair["somatic"]["reasons"] == ["assembly_mismatch"]
    b = replace(a, snps={k: ("G", v[1]) for k, v in a.snps.items()})
    result = compare_samples(a, b)["germline"]
    assert result["reference_conflicts"] == 210
    assert result["shared_snps"] == 0


def test_somatic_branch_overlap_is_separate_from_donor_identity(tmp_path):
    a = load_vcf_samples(_vcf(tmp_path, "a.vcf", ["T"], [_row(i, "0/1:30:60") for i in range(20)]), kind="somatic")[0]
    b = load_vcf_samples(_vcf(tmp_path, "b.vcf", ["T"], [_row(i, "0/1:30:60") for i in range(10, 80)]), kind="somatic")[0]
    result = compare_samples(a, b)
    assert result["germline"]["status"] == "inconclusive"
    somatic = result["somatic"]
    assert somatic["status"] == "shared_somatic_support"
    assert somatic["shared_variants"] == 10
    assert somatic["fraction_a"] == .5
    assert somatic["fraction_b"] == pytest.approx(1 / 7)
    assert somatic["jaccard"] == .125
    a.tumor_id = b.tumor_id = "same tumor"
    b.somatic.clear()
    b.somatic.update({("small", "2", i, "A", "C"): None for i in range(20)})
    result = compare_samples(a, b)
    assert result["somatic"]["status"] == "low_overlap"
    assert result["flags"] == ["expected_tumor_low_overlap"]
    assert not result["identity_conflict"]


def test_mixed_tumor_normal_uses_germline_and_somatic_evidence(tmp_path):
    rows = [_row(i, ["0/0:30:60", "0/1:30:60", "1/1:30:60"][i % 3],
                 ["0/0:30:60", "0/1:30:60", "1/1:30:60"][i % 3]) for i in range(210)]
    rows += [_row(i, "0/1:30:60", "0/0:30:60") for i in range(210, 230)]
    path = _vcf(tmp_path, "mixed.vcf", ["T", "N"], rows,
                metadata="##tumor_sample=T\n##normal_sample=N\n")
    tumor, normal = load_vcf_samples(path, kind="mixed")
    assert tumor.normal_sample == "N"
    assert len(tumor.somatic) == 20
    assert not normal.somatic
    assert tumor.qc["paired_normal_contrast"] == 20
    pair = compare_samples(tumor, normal)
    assert pair["paired_tumor_normal"]
    assert pair["germline"]["status"] == "compatible"
    assert pair["somatic"]["status"] == "inconclusive"
    # Selecting just the tumor still reads its matched normal column.
    selected = load_vcf_samples(path, kind="mixed", samples=[SampleSpec("T")])[0]
    assert selected.somatic == tumor.somatic


def test_no_gt_somatic_depths_missing_normal_and_vaf_correlation(tmp_path):
    rows = [_row(i, "%d,%d" % (30 - i, i), "30,0", fmt="AD") for i in range(3, 20)]
    rows += [_row(30, "20,10", ".", fmt="AD"), _row(31, "20,10", "20,10", fmt="AD")]
    path = _vcf(tmp_path, "depth.vcf", ["T", "N"], rows)
    specs = [SampleSpec("T", role="tumor"), SampleSpec("N", role="normal")]
    tumor, normal = load_vcf_samples(path, kind="mixed", samples=specs)
    assert len(tumor.somatic) == 17
    assert not tumor.snps and not normal.snps
    assert compare_samples(tumor, normal)["germline"]["status"] == "inconclusive"
    result = compare_samples(tumor, tumor)["somatic"]
    assert result["vaf_correlation"] == pytest.approx(1)
    assert result["vaf_pairs"] == 17


def test_mixed_unpaired_needs_somatic_flag_and_gt_takes_precedence(tmp_path):
    rows = [_row(1, "0/1:20,10", fmt="GT:AD"),
            _row(2, "0/1:20,10", fmt="GT:AD", info="SOMATIC"),
            _row(3, "0/0:20,10", fmt="GT:AD", info="SOMATIC")]
    sample = load_vcf_samples(_vcf(tmp_path, "mixed.vcf", ["T"], rows), kind="mixed",
                              samples=[SampleSpec("T", role="tumor")])[0]
    assert list(sample.somatic) == [("small", "1", 400001, "A", "C")]


def test_sv_normalization_deduplicates_reciprocal_breakends(tmp_path):
    rows = [["1", 101, "b1", "A", "A[2:201[", ".", "PASS", ".", "GT:DP:GQ", "0/1:30:60"],
            ["2", 201, "b2", "C", "]1:101]C", ".", "PASS", ".", "GT:DP:GQ", "0/1:30:60"],
            _row(2, "0/1:30:60", alt="<DEL>", info="END=400100;IMPRECISE"),
            _row(3, "0/1:30:60", alt="<INS>", info="SVLEN=200")]
    sample = load_vcf_samples(_vcf(tmp_path, "sv.vcf", ["T"], rows), kind="somatic")[0]
    assert len(sample.somatic) == 1
    assert sample.qc["duplicate_somatic_alleles"] == 1
    assert sample.qc["unsupported_or_imprecise_somatic_alleles"] == 2


@pytest.mark.parametrize("call,fmt", [("0/2", "GT"), ("9/.", "GT"), ("0/1:NaN", "GT:DP"),
                                      ("0/1:-1", "GT:DP"), ("10,2,3", "AD"),
                                      ("1.2", "AF"), ("0/1:30:extra", "GT:DP")])
def test_malformed_evidence_is_an_error_with_line_number(tmp_path, call, fmt):
    path = _vcf(tmp_path, "invalid.vcf", ["S"], [_row(1, call, fmt=fmt)])
    with pytest.raises(ValueError, match=r"invalid.vcf:\d+:"):
        load_vcf_samples(path)


def test_roles_builds_and_sample_selection_are_explicit(tmp_path):
    path = _vcf(tmp_path, "multi.vcf", ["tumor", "normal"], [])
    with pytest.raises(ValueError, match="Declare tumor/normal"):
        load_vcf_samples(path, kind="somatic")
    with pytest.raises(ValueError, match="contradicts"):
        load_vcf_samples(path, assembly="GRCh38")
    with pytest.raises(ValueError, match="Unknown VCF sample"):
        load_vcf_samples(path, samples=[SampleSpec("missing")])
    path = _vcf(tmp_path, "no-build.vcf", ["S"], [], reference=None)
    with pytest.raises(ValueError, match="Declare assembly"):
        load_vcf_samples(path)
    assert load_vcf_samples(path, assembly="hg38")[0].assembly == "GRCh38"
    path = _vcf(tmp_path, "duplicate.vcf", ["S", "S"], [])
    with pytest.raises(ValueError, match="unique sample"):
        load_vcf_samples(path)


def test_multiallelic_somatic_indexes_and_minimal_indel_representation(tmp_path):
    one = [_row(1, "0/2:10,0,5", alt="C,G", fmt="GT:AD"),
           ["1", 99, ".", "GAT", "GA", ".", "PASS", ".", "GT:DP:GQ", "0/1:30:60"]]
    two = [_row(1, "1/0:10,5,0", alt="G,C", fmt="GT:AD"),
           ["1", 100, ".", "AT", "A", ".", "PASS", ".", "GT:DP:GQ", "0/1:30:60"]]
    a = load_vcf_samples(_vcf(tmp_path, "one.vcf", ["T"], one), kind="somatic")[0]
    b = load_vcf_samples(_vcf(tmp_path, "two.vcf", ["T"], two), kind="somatic")[0]
    assert a.somatic == b.somatic == {("small", "1", 200001, "A", "G"): 1 / 3,
                                      ("small", "1", 101, "T", ""): None}
    assert compare_samples(a, b)["somatic"]["shared_variants"] == 2


def test_conflicting_vafs_and_missing_allele_depths_do_not_gain_evidence(tmp_path):
    rows = [_row(1, "0/1:20,10", fmt="GT:AD"), _row(1, "0/1:10,20", fmt="GT:AD"),
            _row(1, "0/1:20,10", fmt="GT:AD"), _row(2, ".:.,5", fmt="GT:AD")]
    a = load_vcf_samples(_vcf(tmp_path, "repeated.vcf", ["T"], rows), kind="somatic")[0]
    assert a.somatic == {("small", "1", 200001, "A", "C"): None}
    assert a.qc["duplicate_somatic_alleles"] == 2


def test_constant_vafs_have_no_correlation_even_with_floating_point_roundoff(tmp_path):
    rows = [_row(i, "0/1:20,10", fmt="GT:AD") for i in range(20)]
    sample = load_vcf_samples(_vcf(tmp_path, "constant.vcf", ["T"], rows), kind="somatic")[0]
    result = compare_samples(sample, sample)["somatic"]
    assert result["vaf_pairs"] == 20
    assert result["vaf_correlation"] is None


def test_multiple_normals_require_pairing_and_limits_are_explicit(tmp_path):
    path = _vcf(tmp_path, "multi.vcf", ["T", "N1", "N2"], [])
    specs = [SampleSpec("T", role="tumor"), SampleSpec("N1", role="normal"), SampleSpec("N2", role="normal")]
    with pytest.raises(ValueError, match="Multiple normals"):
        load_vcf_samples(path, kind="mixed", samples=specs)
    specs[0] = replace(specs[0], normal_sample="N2")
    assert load_vcf_samples(path, kind="mixed", samples=specs)[0].normal_sample == "N2"
    path = _vcf(tmp_path, "limit.vcf", ["T"], [_row(i, "0/1:30:60") for i in range(2)])
    with pytest.raises(ValueError, match="max_somatic_variants exceeded"):
        load_vcf_samples(path, kind="somatic", config=SampleCheckConfig(max_somatic_variants=1))


@pytest.mark.parametrize("options", [{"min_depth": -1}, {"max_snps": 10}, {"min_shared_snps": True},
                                     {"min_somatic_overlap": float("nan")}, {"require_quality": "yes"},
                                     {"compatible_concordance": .5}, {"max_compatible_ibs0": .5}])
def test_invalid_config(options):
    with pytest.raises(ValueError):
        SampleCheckConfig(**options)


def test_manifest_cli_json_tsv_and_mismatch_exit(tmp_path):
    _donor(tmp_path, "a.vcf")
    _donor(tmp_path, "b.vcf", swapped=True)
    manifest = tmp_path / "samples.tsv"
    manifest.write_text("path\tsample\trole\tdonor_id\tlabel\n"
                        "a.vcf\tS\tnormal\tpatient\ta\n"
                        "b.vcf\tS\tnormal\tpatient\tb\n")
    out, table = tmp_path / "report.json", tmp_path / "pairs.tsv"
    completed = subprocess.run([sys.executable, "-m", "varcode.cli.main", "check-samples", "--manifest", str(manifest),
                                "--json", str(out), "--tsv", str(table), "--fail-on-mismatch"],
                               text=True, capture_output=True)
    assert completed.returncode == 1, completed.stderr
    assert not completed.stdout
    report = json.loads(out.read_text())
    assert report["identity_conflicts"] == 1
    assert report["samples"][0]["path"] == str(tmp_path / "a.vcf")
    assert "expected_donor_mismatch" in table.read_text()


def test_cli_stdout_inconclusive_and_input_errors(tmp_path, capsys):
    path = _vcf(tmp_path, "empty.vcf", ["a", "b"], [])
    assert main(["check-samples", "--germline", str(path), "--fail-on-mismatch"]) == 0
    captured = capsys.readouterr()
    assert json.loads(captured.out)["pairs"][0]["germline"]["status"] == "inconclusive"
    with pytest.raises(SystemExit) as error:
        main(["check-samples", "--germline", str(path), "--json", str(path)])
    assert error.value.code == 2
    assert path.read_text().startswith("##fileformat")
    with pytest.raises(SystemExit) as error:
        main(["check-samples", "--germline", str(path), "--germline", str(path)])
    assert error.value.code == 2


def test_annotation_cli_dispatch_preserves_existing_arguments(monkeypatch):
    from varcode.cli import effects_script
    seen = []
    monkeypatch.setattr(effects_script, "main", lambda args: seen.extend(args))
    main(["--vcf", "original.vcf", "--genome", "GRCh38"])
    assert seen == ["--vcf", "original.vcf", "--genome", "GRCh38"]


def test_no_reference_data_access_for_small_variant_checks(tmp_path, monkeypatch):
    import varcode.reference
    def forbidden(*args, **kwargs):
        raise AssertionError("Sample checks must not resolve annotation releases")
    monkeypatch.setattr(varcode.reference, "infer_genome", forbidden)
    assert _donor(tmp_path).snps
