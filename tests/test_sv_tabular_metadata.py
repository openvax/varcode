"""Structural records must never masquerade as placeholder SNVs (#417)."""

import pickle

import pytest

from varcode import StructuralVariant, Variant, VariantCollection
from varcode.csv_helpers import STRUCTURAL_VARIANT_COLUMNS
from varcode.effects import EffectCollection, Intergenic


@pytest.mark.parametrize("ref", ["A", "C", "G", "T", "N"])
@pytest.mark.parametrize("kind,alt", [
    ("DEL", "<DEL>"), ("DUP", "<DUP>"), ("INV", "<INV>"),
    ("CNV", "<CN3>"), ("INS", "<INS:ME:ALU>"),
    ("BND", "G[2:200["), ("BND", ".ACGT"), ("INS", "ACGT"),
])
def test_actual_alleles_and_no_small_edit_flags(ref, kind, alt):
    sv = StructuralVariant("1", 100, kind, end=200, ref=ref, alt=alt)
    assert sv.start == sv.original_start == 100
    assert sv.ref == sv.original_ref == ref
    assert sv.alt == sv.original_alt == sv.symbolic_alt == alt
    for name in ["is_snv", "is_indel", "is_insertion", "is_deletion",
                 "is_transition", "is_transversion"]:
        assert getattr(sv, name) is False
    assert sv.is_structural
    for restored in [StructuralVariant.from_json(sv.to_json()), pickle.loads(pickle.dumps(sv))]:
        assert restored == sv
        assert (restored.ref, restored.alt) == (ref, alt)


def examples():
    snv = Variant("1", 50, "C", "T")
    sv = StructuralVariant("1", 100, "DEL", end=200, ref="T", affected_start=101,
                           mate_contig="2", mate_start=300)
    for variant in [snv, sv]:
        variant._gene_names = []
        variant._gene_ids = []
    return snv, sv


@pytest.mark.parametrize("effects", [False, True])
def test_mixed_tables_retain_structural_coordinates(effects, tmp_path):
    snv, sv = examples()
    collection = (EffectCollection([Intergenic(snv), Intergenic(sv)], sort_key=False)
                  if effects else VariantCollection([snv, sv]))
    df = collection.to_dataframe()
    assert list(df.columns) == list(collection._DATAFRAME_COLUMNS + STRUCTURAL_VARIANT_COLUMNS)
    row = df.iloc[1]
    assert (row["ref"], row["alt"], row["sv_type"]) == ("T", "<DEL>", "DEL")
    assert (row["start"], row["end"], row["affected_start"], row["affected_end"]) == (100, 200, 101, 200)
    assert (row["mate_contig"], row["mate_start"]) == ("2", 300)
    if effects:
        assert not row["is_snv"] and not row["is_transition"] and not row["is_transversion"]
        assert df.iloc[0]["is_snv"] and df.iloc[0]["is_transition"]
    path = tmp_path / "mixed.csv"
    collection.to_csv(str(path))
    with pytest.raises(ValueError, match="Structural rows cannot be reconstructed"):
        type(collection).from_csv(str(path))


@pytest.mark.parametrize("effects", [False, True])
@pytest.mark.parametrize("empty", [False, True])
def test_point_only_and_empty_schemas_unchanged(effects, empty):
    snv, _ = examples()
    variants = [] if empty else [snv]
    collection = (EffectCollection([Intergenic(v) for v in variants]) if effects
                  else VariantCollection(variants))
    assert tuple(collection.to_dataframe().columns) == collection._DATAFRAME_COLUMNS
