"""Generated regression data for #502; no source files or reference downloads."""

from io import StringIO
import os
from pathlib import Path
import subprocess
import sys

import pytest

from varcode import Variant, load_vcf
from varcode.vcf_output import variants_to_vcf
from varcode.vcf_parsing import VCFHeader


def _variant(pos=100, alt="C"):
    return Variant("1", pos, "A", alt, genome=81)


def _export(variants, metadata, **kwargs):
    out = StringIO()
    variants_to_vcf(variants, metadata, out=out, **kwargs)
    return out.getvalue()


def _reload(tmp_path, text):
    path = tmp_path / "export.vcf"
    path.write_text(text)
    return load_vcf(str(path), genome=81)


def test_samples_and_fields_are_selected_by_name(tmp_path):
    a, b = _variant(2), _variant(10)
    metadata = {
        a: {"sample_info": {"tumor": {"GT": "0/1", "DP": 11, "AD": [6, 5]},
                            "normal": {"AD": [22, 0], "DP": 22, "GT": "0/0"}}},
        b: {"sample_info": {"normal": {"GT": "0/0", "DP": 33, "AD": [33, 0]},
                            "tumor": {"DP": 44, "AD": [20, 24], "GT": "0/1"}}},
    }
    text = _export(iter([b, a]), metadata)
    rows = [line.split("\t") for line in text.splitlines() if not line.startswith("##")]
    assert rows[0][9:] == ["normal", "tumor"]
    assert [r[1] for r in rows[1:]] == ["2", "10"]
    assert all(r[8].split(":")[0] == "GT" for r in rows[1:])
    loaded = _reload(tmp_path, text)
    for variant in loaded:
        assert loaded.metadata[variant]["sample_info"] == metadata[variant]["sample_info"]


def test_union_of_samples_and_format_fields_preserves_missing_values(tmp_path):
    a, b, c = [_variant(p) for p in (2, 10, 20)]
    metadata = {
        a: {"sample_info": {"tumor": {"GT": "0|1", "AD": [None, 8]},
                            "normal": {"GT": "0/0", "DP": 30}}},
        b: {"sample_info": {"normal": {"GT": "0/1", "DP": 20}}},
        c: {},
    }
    loaded = _reload(tmp_path, _export([a, b, c], metadata))
    assert loaded.samples == ["normal", "tumor"]
    assert loaded.metadata[a]["sample_info"]["tumor"] == {"GT": "0|1", "AD": [None, 8], "DP": None}
    assert loaded.metadata[a]["sample_info"]["normal"] == {"GT": "0/0", "AD": None, "DP": 30}
    assert loaded.metadata[b]["sample_info"]["tumor"] == {"GT": ".", "DP": None}
    row = [r for r in _export([a, b, c], metadata).splitlines() if r.startswith("1\t20\t")][0]
    assert row.split("\t")[8:] == [".", ".", "."]


def _multiallelic(tmp_path, identifier="."):
    path = tmp_path / "source.vcf"
    path.write_text(
        '##fileformat=VCFv4.2\n'
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="Counts">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Depths">\n'
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor\tnormal\n'
        '1\t100\t' + identifier + '\tA\tG,C\t.\tPASS\tAC=1,1\tGT:AD:PL\t0|2:20,0,10:0,1,2,3,4,5\t0/1:15,8,0:5,4,3,2,1,0\n'
    )
    return load_vcf(str(path), genome=81), VCFHeader.from_path(str(path))


@pytest.mark.parametrize("identifier", [".", "rs123"])
def test_multiallelic_order_and_cardinality_survive_reordered_input(tmp_path, identifier):
    variants, header = _multiallelic(tmp_path, identifier)
    text = _export(sorted(variants, key=lambda v: v.alt), variants.metadata, header=header)
    row, = [r for r in text.splitlines() if not r.startswith("#")]
    assert row.split("\t")[4] == "G,C"
    assert 'ID=AC,Number=A,Type=Integer' in text
    loaded = _reload(tmp_path, text)
    for v in loaded:
        for sample in ("tumor", "normal"):
            assert dict(loaded.metadata[v]["sample_info"][sample]) == dict(
                variants.metadata[v]["sample_info"][sample])
        assert loaded.metadata[v]["info"] == variants.metadata[v]["info"]
    assert [v.alt for v in loaded.for_sample("tumor")] == ["C"]
    assert [v.alt for v in loaded.for_sample("normal")] == ["G"]


def test_filtered_multiallelic_group_is_rejected_before_writing(tmp_path):
    variants, _ = _multiallelic(tmp_path)
    out = StringIO()
    with pytest.raises(ValueError, match="incomplete source ALT list"):
        variants_to_vcf([variants[0]], variants.metadata, out=out)
    assert out.getvalue() == ""


def test_old_multiallelic_metadata_is_rejected():
    v = _variant()
    with pytest.raises(ValueError, match="source ALT list"):
        _export([v], {v: {"alt_allele_index": 1, "sample_info": {"tumor": {"GT": "0/2"}}}})
    with pytest.raises(ValueError, match="unavailable ALT"):
        _export([v], {v: {"sample_info": {"tumor": {"GT": "0/2"}}}})


def test_conflicting_multiallelic_metadata_is_rejected(tmp_path):
    variants, _ = _multiallelic(tmp_path)
    metadata = {v: dict(variants.metadata[v]) for v in variants}
    metadata[variants[0]]["qual"] = 42
    with pytest.raises(ValueError, match="Conflicting source metadata"):
        _export(variants, metadata)


def test_missing_info_and_sample_filters_are_valid_vcf(tmp_path):
    v = _variant()
    metadata = {v: {"info": {"X": [1, None], "SOMATIC": False}, "sample_info": {
        "tumor": {"GT": "0/1", "FT": ["LowDepth", "LowGQ"]},
        "normal": {"GT": "0/0", "FT": []}}}}
    text = _export([v], metadata)
    assert "X=1,." in text
    assert "LowDepth;LowGQ" in text
    loaded = _reload(tmp_path, text)
    assert loaded.metadata[v]["info"] == {"X": [1, None]}
    assert loaded.metadata[v]["sample_info"] == metadata[v]["sample_info"]


def test_empty_input_and_runtime_stdout(capsys):
    variants_to_vcf(iter([]), {})
    assert capsys.readouterr().out == (
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


def test_hash_seed_does_not_change_sample_columns():
    script = '''
from varcode import Variant
from varcode.vcf_output import variants_to_vcf
v = Variant("1", 100, "A", "C", genome=81)
values = {"tumor": "0/1", "normal": "0/0", "replicate": "1/1"}
metadata = {v: {"sample_info": {s: {"GT": values[s]} for s in set(values)}}}
variants_to_vcf([v], metadata)
'''
    outputs = [subprocess.check_output(
        [sys.executable, "-c", script], text=True,
        cwd=Path(__file__).resolve().parents[1],
        env=dict(os.environ, PYTHONHASHSEED=str(seed))) for seed in (1, 2, 42)]
    assert len(set(outputs)) == 1
    assert outputs[0].splitlines()[-1].split("\t")[9:] == ["0/0", "1/1", "0/1"]


@pytest.mark.parametrize("index", [0, 1])
def test_indexed_legacy_metadata_requires_source_alleles(index):
    v = _variant()
    with pytest.raises(ValueError, match="reload the VCF"):
        _export([v], {v: {"alt_allele_index": index, "sample_info": {"s": {"GT": "0/1"}}}})


def test_likelihoods_must_agree_with_allele_count_and_ploidy():
    v = _variant()
    with pytest.raises(ValueError, match="ALT/ploidy"):
        _export([v], {v: {"sample_info": {"s": {"GT": "0/1", "PL": [0, 1, 2, 3, 4, 5]}}}})
