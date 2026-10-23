"""Public path, comparison and validation contracts independent of annotation."""

import gzip
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest
from pyensembl import Genome

from varcode import Variant, load_maf, load_maf_dataframe, load_vcf
from varcode.nucleotides import normalize_nucleotide_string
from varcode.reference import infer_genome
from varcode.vcf import read_vcf_into_dataframe


DATA = Path(__file__).parent / "data"


class TextPath:
    def __init__(self, path):
        self.path = path

    def __fspath__(self):
        return str(self.path)


class BytesPath:
    def __fspath__(self):
        return b"example.vcf"


class InvalidPath:
    def __fspath__(self):
        return 123


@pytest.mark.parametrize("path_type", [Path, TextPath])
@pytest.mark.parametrize("compressed", [False, True])
def test_vcf_pathlike_preserves_alleles_and_source(tmp_path, path_type, compressed):
    path = tmp_path / ("variants.vcf.gz" if compressed else "variants.vcf")
    contents = (DATA / "simple.1.vcf").read_bytes()
    path.write_bytes(gzip.compress(contents) if compressed else contents)
    variants = load_vcf(path_type(path))
    variant, = variants
    assert (variant.contig, variant.start, variant.ref, variant.alt) == ("1", 13281, "C", "G")
    assert variants.path == str(path)
    frame = read_vcf_into_dataframe(path_type(path))
    assert list(frame.POS) == [13281]


@pytest.mark.parametrize("path_type", [Path, TextPath])
def test_maf_pathlike_preserves_alleles_and_source(path_type):
    path = DATA / "tcga_ov.head.maf"
    variants = load_maf(path_type(path))
    assert [v.start for v in variants] == [1650797, 23836447, 231401797, 124617502]
    assert variants.path == str(path)
    assert len(load_maf_dataframe(path_type(path))) == 4


@pytest.mark.parametrize("loader", [load_vcf, read_vcf_into_dataframe, load_maf, load_maf_dataframe])
@pytest.mark.parametrize("path", [None, 123, b"example.vcf", BytesPath(), InvalidPath()])
def test_loaders_reject_non_text_paths(loader, path):
    with pytest.raises(TypeError):
        loader(path)


def test_vcf_http_url_remains_a_string_and_retains_gzip_handling(tmp_path, monkeypatch):
    import urllib.request
    downloaded = tmp_path / "download"
    url = "https://example.test/variants.vcf.gz"

    def retrieve(actual_url):
        assert actual_url == url
        downloaded.write_bytes(gzip.compress((DATA / "simple.1.vcf").read_bytes()))
        return str(downloaded), {}

    monkeypatch.setattr(urllib.request, "urlretrieve", retrieve)
    variant, = load_vcf(url)
    assert (variant.start, variant.ref, variant.alt) == (13281, "C", "G")
    assert not downloaded.exists()
    assert not downloaded.with_suffix(".gz").exists()


@pytest.fixture
def variant():
    return Variant("1", 10, "A", "T", genome=Genome("test", "synthetic"))


@pytest.mark.parametrize("other", [None, 1, "variant", object(), SimpleNamespace(contig="1")])
def test_variant_comparison_with_other_types(variant, other):
    assert (variant == other) is False
    assert (other == variant) is False
    assert (variant != other) is True
    assert variant in [other, variant]
    with pytest.raises(TypeError):
        variant < other
    with pytest.raises(TypeError):
        other < variant


def test_variant_comparison_uses_reflected_dispatch(variant):
    class Reflected:
        def __eq__(self, other):
            return other is variant

        def __gt__(self, other):
            return other is variant

    assert variant == Reflected()
    assert variant < Reflected()


def test_variant_subclass_equality_and_locus_order_are_preserved(variant):
    class Subclass(Variant):
        pass

    same = Subclass("1", 10, "A", "T", genome=variant.genome)
    later = Variant("1", 11, "A", "T", genome=variant.genome)
    assert variant == same and same == variant
    assert variant != later
    assert sorted([later, variant]) == [variant, later]
    assert not variant < same


@pytest.mark.parametrize("value", [float("nan"), np.float64("nan")])
def test_nucleotide_nan_handling_precedes_string_checks(value):
    assert normalize_nucleotide_string(value) == ""
    with pytest.raises(TypeError, match="nucleotide string"):
        normalize_nucleotide_string(value, treat_nan_as_empty=False)


@pytest.mark.parametrize("value", [None, 5, 2.5, b"ACGT", ["A"]])
def test_invalid_nucleotide_types_raise_the_public_validation_error(value):
    with pytest.raises(TypeError, match="nucleotide string"):
        normalize_nucleotide_string(value)


def test_nucleotide_strings_keep_existing_normalization():
    assert normalize_nucleotide_string("acgt") == "ACGT"
    for value in ("", ".", "-"):
        assert normalize_nucleotide_string(value) == ""
    assert normalize_nucleotide_string("n", allow_extended_nucleotides=True) == "N"
    with pytest.raises(ValueError, match="Invalid character"):
        normalize_nucleotide_string("n")


@pytest.mark.parametrize("value", [75, True])
def test_reference_integer_dispatch_preserves_existing_types(monkeypatch, value):
    import varcode.reference as reference
    genome = Genome("test", "synthetic")

    def cached(release):
        assert release is value
        return genome

    monkeypatch.setattr(reference, "cached_ensembl_release", cached)
    assert infer_genome(value) == (genome, False)
    with pytest.raises(TypeError, match="Expected genome"):
        infer_genome(np.int64(75))


def test_varcode_import_does_not_require_typechecks_directly():
    # Dependencies may still use typechecks themselves. Block imports by
    # Varcode, which never declared it, in a fresh interpreter.
    code = '''
import builtins
original = builtins.__import__
def guarded(name, globals=None, *args, **kwargs):
    if name == "typechecks" and (globals or {}).get("__name__", "").startswith("varcode"):
        raise ModuleNotFoundError("Varcode must not import undeclared typechecks")
    return original(name, globals, *args, **kwargs)
builtins.__import__ = guarded
import varcode
'''
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
