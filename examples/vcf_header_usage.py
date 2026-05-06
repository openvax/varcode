#!/usr/bin/env python3
# Licensed under the Apache License, Version 2.0 (the "License").
"""
Standalone usage example for ``varcode.vcf_parsing.VCFHeader``.

This demonstrates how to parse a VCF header and decode INFO / FORMAT cells
without going through ``load_vcf`` — useful when you want the typed dicts
but are streaming records yourself, building a custom data frame, etc.

Run from the repo root:

    python examples/vcf_header_usage.py

The fixture used here is the canonical VCF 4.2 spec example (samtools/hts-specs).
"""
import sys
from pathlib import Path

# Make the example runnable from anywhere — fall back to the local checkout
# even when a different varcode happens to be on sys.path first.
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from varcode.vcf_parsing import VCFHeader  # noqa: E402

TRIO_VCF = REPO_ROOT / "tests" / "data" / "spec_examples" / "vcf42_spec_trio.vcf"


def main() -> None:
    header = VCFHeader.from_path(str(TRIO_VCF))

    print("== Header ==")
    print(f"fileformat       : {header.metadata['fileformat']}")
    print(f"reference        : {header.metadata['reference']}")
    print(f"samples          : {header.samples}")
    print(f"INFO fields      : {list(header.info_fields)}")
    print(f"FORMAT fields    : {list(header.format_fields)}")

    # The typed declaration for an INFO field is a FieldDef dataclass.
    af = header.info_fields["AF"]
    print()
    print("== INFO[AF] declaration ==")
    print(f"id          : {af.id}")
    print(f"number      : {af.number}    # -1 means 'A' (one value per ALT)")
    print(f"type        : {af.type}")
    print(f"description : {af.description}")

    # Walk records and pretty-print decoded INFO + per-sample dicts. We
    # split rows ourselves to keep this example dependency-free.
    print()
    print("== Records ==")
    with open(TRIO_VCF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
            info = header.parse_info(cols[7])
            samples = header.parse_samples(cols[9:], cols[8])
            print(f"  {chrom}:{pos} {vid or '.'}  {ref} -> {alt}")
            print(f"    INFO    = {info}")
            for sample_name, fields in samples.items():
                print(f"    {sample_name:<8} = {dict(fields)}")
            print()


if __name__ == "__main__":
    main()
