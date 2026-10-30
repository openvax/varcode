"""Benchmark installed native annotations against frozen independent products.

Run from a checkout after installing the references described in
docs/native_references.md. Timings measure warm annotation of distinct variants,
including all overlapping transcripts; installation/index time is separate.
"""

import argparse
import json
import platform
from pathlib import Path
from time import perf_counter

from pyensembl import cached_release

from varcode import Variant, iter_candidate_effects
from varcode.reference_metadata import reference_metadata
from varcode.version import __version__


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", choices=("mouse", "dog", "chm13"), default="chm13")
    parser.add_argument("--chm13-cache")
    parser.add_argument("--genome-fasta")
    parser.add_argument("--repeat", type=int, default=3)
    args = parser.parse_args()
    if args.repeat < 1:
        parser.error("--repeat must be positive")
    if args.genome == "chm13":
        from orfdata import load_genome
        genome = load_genome("t2t-chm13", cache_root=args.chm13_cache,
                             genome_fasta_path_or_url=args.genome_fasta)
    else:
        genome = cached_release(115 if args.genome == "mouse" else "2020_10", species=args.genome)
        if args.genome_fasta:
            from varcode import Genome
            genome = Genome(genome, fasta=args.genome_fasta)
    path = Path(__file__).resolve().parents[1] / "tests/data/native_references/cases.json"
    cases = [c for c in json.loads(path.read_text()) if c["genome"] == args.genome]
    groups = {}
    for case in cases:
        key = (case["contig"], case["position"], case["ref"], case["alt"])
        groups.setdefault(key, []).append(case)
    def annotate(key, mode):
        # Use fresh variant objects and force lazy candidate protein sequences
        # for both implementations; neither gets a cached effect-only timing.
        effects = Variant(*key, genome=genome).effects(annotator=mode)
        for root in effects:
            for effect in iter_candidate_effects(root):
                _ = effect.mutant_protein_sequence
        return effects

    report = dict(varcode=__version__, python=platform.python_version(), platform=platform.platform(),
                  reference=reference_metadata(genome), chromosome_dna=bool(args.genome_fasta),
                  distinct_variants=len(groups), transcript_cases=len(cases),
                  timing_scope="fresh variants, all transcript candidates and protein sequences", runs=[])
    for mode in ("fast", "transcript_model"):
        start = perf_counter()
        products = {key: annotate(key, mode) for key in groups}
        first_pass = perf_counter() - start
        failures, checked = [], 0
        for key, rows in groups.items():
            effects = [e for root in products[key] for e in iter_candidate_effects(root)]
            for case in rows:
                if "expected_cdna" not in case or "protein_expectation" in case:
                    continue
                checked += 1
                if not any(e.transcript_id == case["transcript"]
                           and e.mutant_protein_sequence == case["expected_protein"] for e in effects):
                    failures.append(case["id"])
        elapsed = []
        for _ in range(args.repeat):
            start = perf_counter()
            for key in groups:
                annotate(key, mode)
            elapsed.append(perf_counter() - start)
        report["runs"].append(dict(annotator=mode, first_pass_seconds=first_pass,
                                   warm_batch_seconds=elapsed, checked_proteins=checked,
                                   failed_cases=failures))
    print(json.dumps(report, indent=2))
    if any(run["failed_cases"] for run in report["runs"]):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
