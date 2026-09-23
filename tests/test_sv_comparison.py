"""Small explicit VCF records are generated here; no external fixture cache."""
import json

import pytest

from varcode.sv_comparison import compare_sv_calls, main


def row(id, chrom="1", pos=100, ref="C", alt="C[2:200[", **kwargs):
    return dict(call_id=id, caller="caller-" + id, sample="sample-" + id,
                build="GRCh38", chrom=chrom, pos=str(pos), ref=ref, alt=alt, **kwargs)


def test_reciprocal_and_repeated_exports_preserve_all_evidence():
    a = row("a", record_id="vcf-a", source_url="https://example.org/a.vcf")
    b = row("b", chrom="2", pos=200, ref="G", alt="]1:100]G", record_id="vcf-b")
    c = dict(a, call_id="c", duplicate_export_of="a")
    result = compare_sv_calls([a, b, c])
    group, = result["groups"]
    assert group["status"] == "exact_reported_allele"
    assert group["call_count"] == 3
    assert len(group["exact_ids"]) == 1
    assert [m["raw"] for m in result["members"]] == [a, b, c]
    assert len(group["callers"]) == 2  # An exported copy is not another caller.


def test_inserted_sequence_disagreement_is_never_exact_equivalence():
    result = compare_sv_calls([row("a", alt="CA[2:200["), row("b", alt="CT[2:200[")])
    group, = result["groups"]
    assert group["status"] == "nearby_candidate"
    assert len(group["exact_ids"]) == 2
    assert group["insertion_disagreement"] is True


def test_breakpoint_disagreement_and_no_transitive_chain():
    rows = [row("a", pos=100), row("b", pos=180), row("c", pos=260)]
    result = compare_sv_calls(rows, max_distance=100)
    assert sorted(g["call_count"] for g in result["groups"]) == [1, 2]
    assert all(max(g["breakpoint_spread"]) <= 100 for g in result["groups"])
    assert result == compare_sv_calls(reversed(rows), max_distance=100)
    assert len(compare_sv_calls(rows, max_distance=0)["groups"]) == 3


@pytest.mark.parametrize("change", [dict(build="GRCh37"), dict(alt="C]2:200]"), dict(alt="C[3:200["), dict(alt="C[2:301[")])
def test_assembly_orientation_contig_and_second_breakpoint_are_required(change):
    assert len(compare_sv_calls([row("a"), dict(row("b"), **change)])["groups"]) == 2


@pytest.mark.parametrize("ref,alt,pos,bnd_alt,mate", [
    ("C", "CA", 100, "CA[1:101[", 101),
    ("CAAA", "C", 100, "C[1:104[", 104),
    ("CAAA", "CTT", 100, "CTT[1:104[", 104),
])
def test_explicit_indel_and_breakend_use_same_retained_base_coordinates(ref, alt, pos, bnd_alt, mate):
    result = compare_sv_calls([row("a", pos=pos, ref=ref, alt=alt), row("b", pos=pos, alt=bnd_alt)])
    group, = result["groups"]
    assert group["status"] == "exact_reported_allele"
    assert len(group["exact_ids"]) == 1
    assert result["members"][0]["positions"] == (100, mate)


def test_symbolic_and_explicit_deletion_agree_but_inversion_has_two_junctions():
    result = compare_sv_calls([row("a", alt="<DEL>", info="END=103"), row("b", ref="CAAA", alt="C")])
    assert len(result["groups"][0]["exact_ids"]) == 1
    inv = compare_sv_calls([row("a", alt="<INV>", info="END=103"), row("b", alt="C]1:103]")])
    assert len(inv["groups"]) == 2  # One adjacency doesn't prove the whole inversion.


def test_unknown_sequence_imprecision_and_unplaced_calls_stay_explicit():
    result = compare_sv_calls([row("a", alt="<INS>"), row("b", alt="CA[1:101["),
                              row("c", alt="C."), row("d", info="CIPOS=-5,5"),
                              dict(row("e"), build="unknown"), row("f", alt="bogus")])
    members = {m["raw"]["call_id"]: m for m in result["members"]}
    assert members["a"]["inserted_sequences"] == (None,)
    assert members["d"]["imprecise"] is True
    assert all(members[k]["status"] == "unresolved" for k in "cef")
    assert len(result["members"]) == 6
    assert all("reason" in members[k] for k in "cef")


def test_input_errors_and_empty_input():
    assert compare_sv_calls([])["groups"] == []
    with pytest.raises(ValueError, match="Duplicate"):
        compare_sv_calls([row("a"), row("a")])
    with pytest.raises(ValueError, match="required"):
        compare_sv_calls([{}])
    with pytest.raises(ValueError, match="non-negative"):
        compare_sv_calls([], -1)


def test_portable_cli_exports_every_input_and_hash(tmp_path):
    import csv
    source = tmp_path / "calls.csv"
    with source.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=row("a").keys())
        writer.writeheader()
        writer.writerows([row("a"), row("b", alt="C.")])
    destination = tmp_path / "result"
    main(["--calls", str(source), "--output", str(destination), "--max-distance", "10"])
    result = json.loads((destination / "comparison.json").read_text())
    assert result["provenance"]["call_count"] == 2
    assert len(result["provenance"]["input_sha256"]) == 64
    assert len(list(csv.DictReader((destination / "members.csv").open()))) == 2
