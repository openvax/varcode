"""Compare reported SV alleles without erasing source or breakpoint disagreements.

Run ``python -m varcode.sv_comparison --help`` for the portable CSV interface.
Coordinates follow VCF (one-based); no reference liftover or repeat left-alignment
is performed. Nearby groups nominate comparisons, not proven identical alleles.
"""

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from .nucleotides import reverse_complement
from .structural_variant import Breakend
from .sv_allele_parser import parse_symbolic_alt
from .variant import Variant

_REQUIRED = {"call_id", "caller", "sample", "build", "chrom", "pos", "ref", "alt"}


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()[:20]


def _info(text):
    result = {}
    for token in (text or "").split(";"):
        key, _, value = token.partition("=")
        if key in {"END", "SVLEN"} and value:
            result[key] = int(value.split(",")[0])
        elif key in {"CIPOS", "CIEND"} and value:
            result[key] = tuple(int(v) for v in value.split(","))
        elif key:
            result[key] = value or True
    return result


def _canonical(junction, inserted):
    first, second = junction
    if second < first:
        first, second = second, first
        if inserted is not None:
            inserted = reverse_complement(inserted)
    return (tuple(first), tuple(second)), inserted


def _normalize(row):
    """Normalize representation only; retain every raw field in the result."""
    assembly = row["build"].split(" (")[0]
    if assembly not in {"GRCh37", "GRCh38"}:
        raise ValueError("Unsupported/unknown assembly: %s" % row["build"])
    info = _info(row.get("info"))
    if row.get("end") and "END" not in info:
        info["END"] = int(row["end"])
    sv = parse_symbolic_alt(row["chrom"], int(row["pos"]), row["ref"], row["alt"],
                            info=info, genome=assembly, convert_ucsc_contig_names=True)
    adjacencies = []
    kind = "adjacency"
    if sv is not None:
        if sv.junctions:
            for junction in sv.junctions:
                if any(end.keeps is None for end in junction):
                    raise ValueError("Unknown breakend orientation")
                adjacencies.append(_canonical(junction, sv.junction_inserted_sequence(junction)))
        elif sv.sv_type == "INS":
            adjacencies.append(_canonical((Breakend(sv.contig, sv.start, "left"),
                                          Breakend(sv.contig, sv.start + 1, "right")),
                                         sv.alt_assembly))
        else:
            raise ValueError("No complete adjacency in %s allele" % sv.sv_type)
        imprecise = bool(info.get("IMPRECISE") or any(info.get(k) not in (None, (0, 0))
                                                    for k in ("CIPOS", "CIEND")))
    else:
        if not row["ref"] or not row["alt"] or set((row["ref"] + row["alt"]).upper()) - set("ACGTN"):
            raise ValueError("Unsupported explicit allele")
        v = Variant(row["chrom"], int(row["pos"]), row["ref"], row["alt"], genome=assembly,
                    convert_ucsc_contig_names=True)
        if not v.ref and not v.alt:
            raise ValueError("Reference-only allele")
        if len(v.ref) == len(v.alt):
            kind = "replacement"
        adjacencies.append(_canonical((Breakend(v.contig, v.start - 1 if v.ref else v.start, "left"),
                                      Breakend(v.contig, v.start + len(v.ref) if v.ref else v.start + 1, "right")), v.alt))
        imprecise = False
    adjacencies.sort(key=lambda item: item[0])
    loci = tuple(end for junction, _ in adjacencies for end in junction)
    inserts = tuple(insert for _, insert in adjacencies)
    topology = (assembly, kind, tuple((chrom, side) for chrom, _, side in loci))
    positions = tuple(pos for _, pos, _ in loci)
    # Explicit replacements need the reference allele too; their two boundaries
    # alone cannot establish the base identity on a nonstandard reference.
    key = (topology, positions, inserts, row["ref"].upper() if kind == "replacement" else "")
    return dict(exact_id=_digest(key), topology=topology, positions=positions,
                inserted_sequences=inserts, imprecise=imprecise, assembly=assembly,
                sequence_known=all(s is not None and "N" not in s.upper() for s in inserts),
                breakends=loci, raw=dict(row), status="normalized")


def compare_sv_calls(rows, max_distance=100):
    """Group complete reported adjacencies and retain all input records.

    Parameters
    ----------
    rows : iterable of mappings
        Fields ``call_id, caller, sample, build, chrom, pos, ref, alt`` are
        required. Optional VCF ``info`` and ``end`` aid symbolic parsing.
        Every input field is preserved. IDs must be unique within the input.
    max_distance : int
        Maximum separation at *every* corresponding breakpoint of *every*
        pair in a nearby group. Deterministic greedy complete-link grouping
        prevents chains from joining distant calls. Distinct inserted alleles
        can share a nearby group but always have distinct exact IDs.

    Returns
    -------
    dict
        Groups and exhaustive members. Malformed, unsupported and unplaced
        records remain singleton groups with an explicit reason. Counts are
        reported calls, not independent molecules or independent callers.
    """
    if not isinstance(max_distance, int) or isinstance(max_distance, bool) or max_distance < 0:
        raise ValueError("max_distance must be a non-negative integer")
    normalized, seen = [], set()
    for row in rows:
        missing = sorted(_REQUIRED - row.keys())
        if missing or not row.get("call_id"):
            raise ValueError("Missing required fields: %s" % (missing or ["call_id"]))
        if row["call_id"] in seen:
            raise ValueError("Duplicate call_id: %s" % row["call_id"])
        seen.add(row["call_id"])
        try:
            record = _normalize(row)
        except (ValueError, TypeError, KeyError) as error:
            record = dict(exact_id=_digest((row["call_id"], dict(row))), raw=dict(row),
                          status="unresolved", reason=str(error))
        normalized.append(record)
    buckets = defaultdict(list)
    for record in normalized:
        key = record.get("topology", ("unresolved", record["raw"]["call_id"]))
        buckets[key].append(record)
    groups = []
    for key in sorted(buckets, key=repr):
        clusters = []
        for record in sorted(buckets[key], key=lambda r: (r.get("positions", ()), r["exact_id"], r["raw"]["call_id"])):
            positions = record.get("positions", ())
            match = None
            # Earlier clusters cease to be eligible once the first breakpoint
            # is too far away; keeping only active clusters bounds dense scans.
            for cluster in reversed(clusters):
                if positions and positions[0] - cluster["minimum"][0] > max_distance:
                    break
                if positions and all(max(p, hi) - min(p, lo) <= max_distance
                                     for p, lo, hi in zip(positions, cluster["minimum"], cluster["maximum"])):
                    match = cluster
                    break
            if match is None:
                match = dict(members=[], minimum=positions, maximum=positions)
                clusters.append(match)
            match["members"].append(record)
            match["minimum"] = tuple(min(a, b) for a, b in zip(match["minimum"], positions))
            match["maximum"] = tuple(max(a, b) for a, b in zip(match["maximum"], positions))
        for cluster in clusters:
            members = cluster["members"]
            ids = sorted(r["raw"]["call_id"] for r in members)
            group_id = _digest(ids)
            exact = sorted({r["exact_id"] for r in members})
            unresolved = any(r["status"] == "unresolved" for r in members)
            sequence_known = all(r.get("sequence_known", False) for r in members)
            imprecise = any(r.get("imprecise", False) for r in members)
            status = ("unresolved" if unresolved else "nearby_candidate" if len(exact) > 1
                      else "same_reported_adjacency" if imprecise or not sequence_known
                      else "exact_reported_allele")
            for member in members:
                member["group_id"] = group_id
            groups.append(dict(group_id=group_id, status=status, call_count=len(members),
                               call_ids=ids, exact_ids=exact,
                               callers=sorted({r["raw"]["caller"] for r in members}),
                               samples=sorted({r["raw"]["sample"] for r in members}),
                               breakpoint_min=cluster["minimum"], breakpoint_max=cluster["maximum"],
                               breakpoint_spread=tuple(b - a for a, b in zip(cluster["minimum"], cluster["maximum"])),
                               insertion_disagreement=len({r.get("inserted_sequences") for r in members}) > 1,
                               imprecise=imprecise, sequence_known=sequence_known))
    return dict(schema_version=1, max_distance=max_distance,
                groups=sorted(groups, key=lambda g: g["group_id"]),
                members=sorted(normalized, key=lambda r: r["raw"]["call_id"]))


def _write_csv(path, rows):
    rows = list(rows)
    with path.open("w", newline="") as out:
        fields = sorted({key for row in rows for key in row})
        writer = csv.DictWriter(out, fieldnames=fields or ["group_id"])
        writer.writeheader()
        writer.writerows({key: json.dumps(value, sort_keys=True) if isinstance(value, (dict, tuple, list)) else value
                          for key, value in row.items()} for row in rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calls", required=True, type=Path, help="CSV with the documented VCF allele/provenance columns")
    parser.add_argument("--output", required=True, type=Path, help="New directory for exhaustive groups, members and manifest")
    parser.add_argument("--max-distance", default=100, type=int, help="Maximum pairwise breakpoint separation in bp (default: 100)")
    args = parser.parse_args(argv)
    try:
        with args.calls.open(newline="") as source:
            result = compare_sv_calls(csv.DictReader(source), args.max_distance)
        args.output.mkdir(parents=True, exist_ok=False)
        _write_csv(args.output / "groups.csv", result["groups"])
        _write_csv(args.output / "members.csv", result["members"])
        from .version import __version__
        result["provenance"] = dict(varcode_version=__version__, input_sha256=hashlib.sha256(args.calls.read_bytes()).hexdigest(),
                                    input=str(args.calls), call_count=len(result["members"]), group_count=len(result["groups"]))
        (args.output / "comparison.json").write_text(json.dumps(result, sort_keys=True) + "\n")
    except (OSError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
