"""Static consumer checks; checked by mypy, never executed as tests."""

from pathlib import Path
from typing import assert_type

from varcode import (
    ContigResolution, ContigResolver, Variant, VariantCollection,
    load_maf, load_maf_dataframe, load_vcf, normalize_contig_name, resolve_contig_name,
)
from varcode.nucleotides import normalize_nucleotide_string
from varcode.vcf import read_vcf_into_dataframe


def check_inputs(variant: Variant, other: object, path: Path) -> None:
    assert_type(load_vcf(path), VariantCollection)
    assert_type(load_maf(path), VariantCollection)
    load_maf_dataframe(path)
    read_vcf_into_dataframe(path)
    assert_type(normalize_nucleotide_string("ACGT"), str)
    assert_type(variant == other, bool)
    assert_type(variant < other, bool)
    assert_type(normalize_contig_name("chr1"), str)
    assert_type(resolve_contig_name("chr1", ["1"]), ContigResolution)
    assert_type(ContigResolver(["chrM"]).resolve("MT").resolved_name, str | None)

    # These ignores must remain necessary (warn_unused_ignores is enabled),
    # proving the public exports reject incompatible arguments statically.
    load_vcf(b"variants.vcf")  # type: ignore[arg-type]
    load_maf(42)  # type: ignore[arg-type]
    load_maf_dataframe(b"variants.maf")  # type: ignore[arg-type]
    read_vcf_into_dataframe(42)  # type: ignore[arg-type]
    normalize_nucleotide_string(None)  # type: ignore[arg-type]
    normalize_contig_name(None)  # type: ignore[arg-type]
    resolve_contig_name("chr1", [1])  # type: ignore[list-item]
