"""Reproduce #407 with locally installed Ensembl 115 and public ESVEE VCFs.

Run this same script with PYTHONPATH pointing to each checkout being compared.
The output includes complete effect graphs for equality checks; only annotation
time and Varcode version are removed from those graphs.
"""

import argparse
import copy
import cProfile
import hashlib
import json
import pstats
import statistics
import time
from importlib.metadata import version
from pathlib import Path

from pyensembl import EnsemblRelease

import varcode
from varcode.transforms import pair_breakends


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("organoid_vcf", help="T1 organoid ESVEE somatic VCF")
    parser.add_argument("tumor_vcf", help="T2 tumor ESVEE somatic VCF")
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--profile", action="store_true")
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    genome = EnsemblRelease(115)
    if not genome.installed():
        parser.error("Install human Ensembl 115 before running this benchmark")
    organoid = varcode.load_vcf(
        args.organoid_vcf, genome=genome, parse_structural_variants=True)
    positions = {58909683, 54155308, 54153839, 36148058}
    cases = [v for v in organoid if v.start in positions]
    if len(cases) != 4:
        parser.error("Organoid VCF must contain the four reported breakends")
    tumor = pair_breakends(varcode.load_vcf(
        args.tumor_vcf, genome=genome, parse_structural_variants=True))
    cases.append(next(v for v in tumor if (
        v.contig == "2" and v.sv_type == "INV"
        and v.start == 25955847 and v.end == 63053516)))
    args.output_directory.mkdir(parents=True, exist_ok=True)
    report = {
        "versions": {"varcode": varcode.__version__, **{
            p: version(p) for p in ("pyensembl", "datacache", "numpy", "pandas")}},
        "source_sha256": {str(path): hashlib.sha256(Path(path).read_bytes()).hexdigest()
                          for path in (args.organoid_vcf, args.tumor_vcf)},
        "cases": [],
    }
    for variant in cases:
        times = []

        def annotate():
            # Retain the reference cache, but recompute variant-level candidates.
            fresh = copy.copy(variant)
            fresh._annotation_cache = {}
            return fresh.effects()

        for _ in range(args.repeats):
            started = time.perf_counter()
            effects = annotate()
            times.append(time.perf_counter() - started)
        data = json.loads(effects.to_json())
        data.pop("annotated_at", None)
        data.pop("annotator_version", None)
        encoded = json.dumps(data, sort_keys=True, separators=(",", ":")).encode()
        stem = args.output_directory / str(variant.start)
        stem.with_suffix(".json").write_bytes(encoded)
        row = {
            "case": variant.short_description, "seconds": times,
            "median_seconds": statistics.median(times),
            "effects": len(effects), "sha256": hashlib.sha256(encoded).hexdigest(),
        }
        report["cases"].append(row)
        print(json.dumps(row), flush=True)
        if args.profile:
            profiler = cProfile.Profile()
            profiler.runcall(annotate)
            profiler.dump_stats(str(stem.with_suffix(".prof")))
            with stem.with_suffix(".txt").open("w") as stream:
                pstats.Stats(profiler, stream=stream).sort_stats("cumulative").print_stats(50)
    (args.output_directory / "report.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
