"""Readers and independent expectations for the pinned public benchmarks."""

import csv
import gzip
import json
from pathlib import Path
import re

from varcode import StructuralVariant, Variant


DATA = Path(__file__).parent / "data/public_somatic"
MANIFEST = json.loads((DATA / "manifest.json").read_text())


def table_rows(name):
    with (DATA / name).open() as source:
        return list(csv.DictReader(source, delimiter="\t"))


def vcf_rows(name):
    path = DATA / name
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as source:
        return [line.rstrip("\n").split("\t")
                for line in source if not line.startswith("#")]


def hcc1395_structural_variant(row, genome):
    """Annotate the reported consensus interval, without inventing alleles.

    Table S4 contains clustered coordinates, sometimes in reverse order.
    It supplies neither insertion sequence nor a base-resolved adjacency.
    Preserve the entire original row; do not infer breakend orientation from
    the identifier's undocumented suffix. TRA and complex BND rows retain
    both reported endpoints with unknown orientation.
    """
    kind = row["SV_type"]
    start, end = int(row["Pos1"]), int(row["Pos2"])
    kwargs = {}
    if kind in ("BND", "TRA"):
        kind = "BND"
        kwargs.update(mate_contig=row["Chrom2"], mate_start=end)
        end = start
    else:
        assert row["Chrom1"] == row["Chrom2"]
        start, end = sorted((start, end))
        if kind == "INS":
            end = start
    return StructuralVariant(
        row["Chrom1"], start, kind, end=end, genome=genome,
        convert_ucsc_contig_names=True,
        info={"published_consensus": dict(row)}, **kwargs)


def canine_variant(row, genome):
    # The published table uses MAF coordinates: a "-" reference denotes an
    # insertion after Start_Position. Variant accepts that representation.
    return Variant(
        row["Chromosome"], int(row["Start_Position"]),
        row["Reference_Allele"], row["Tumor_Seq_Allele2"], genome=genome)


def assert_published_protein_change(hgvs, original, mutant):
    """Apply the paper's one-letter HGVS, independently of Varcode's formatter.

    Exact mutant sequences are checked for substitutions and in-frame edits.
    Frameshift HGVS specifies only the unchanged prefix, first new residue
    and (when known) termination distance; it cannot specify the whole tail.
    A terminal '*' belongs to HGVS but not to the stored peptide.
    """
    reference = original + "*"
    match = re.fullmatch(r"p\.([A-Z*])(\d+)([A-Z*])fs\*(\d+|\?)", hgvs)
    if match:
        old, position, new, stop = match.groups()
        position = int(position) - 1
        assert reference[position] == old
        assert mutant is not None
        assert mutant[:position] == original[:position]
        if new == "*":
            assert len(mutant) == position
        else:
            assert mutant[position:position + 1] == new
        if stop != "?":
            assert len(mutant) == position + int(stop) - 1
        return

    match = re.fullmatch(r"p\.([A-Z*])(\d+)([A-Z*=])", hgvs)
    if match:
        old, position, new = match.groups()
        position = int(position) - 1
        assert reference[position] == old
        expected = (original if new in (old, "=") else
                    original[:position] + new + original[position + 1:])
    elif (match := re.fullmatch(
            r"p\.([A-Z*])(\d+)(?:_([A-Z*])(\d+))?(delins|del|dup)([A-Z*]*)", hgvs)):
        first, start, last, end, operation, new = match.groups()
        start = int(start) - 1
        end = int(end) if end else start + 1
        assert reference[start] == first
        assert not last or reference[end - 1] == last
        expected = (original[:end] + original[start:end] + original[end:]
                    if operation == "dup" else original[:start] + new + original[end:])
    elif (match := re.fullmatch(
            r"p\.([A-Z*])(\d+)_([A-Z*])(\d+)ins([A-Z*]+)", hgvs)):
        first, start, last, end, new = match.groups()
        start, end = int(start), int(end)
        assert end == start + 1
        assert reference[start - 1] == first and reference[end - 1] == last
        expected = original[:start] + new + original[start:]
    else:
        raise AssertionError("Unrecognized published protein change: %s" % hgvs)
    assert mutant == expected.split("*", 1)[0], hgvs


def materialize_effects(effects):
    """Visit lazy alternatives too, so a successful wrapper is not enough."""
    pending, seen = list(effects), set()
    while pending:
        effect = pending.pop()
        if id(effect) in seen:
            continue
        seen.add(id(effect))
        assert effect.short_description
        getattr(effect, "mutant_protein_sequence", None)
        pending.extend(candidate.effect for candidate in getattr(effect, "candidates", ()))
