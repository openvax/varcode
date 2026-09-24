"""CLI for genotype concordance and somatic overlap screening."""

import argparse
import csv
import json
from pathlib import Path
import sys

from ..sample_identity import SampleCheckConfig, SampleSpec, check_sample_identity, load_vcf_samples


def _manifest(path, assembly, config):
    path = path.resolve()
    groups = {}
    allowed = {"path", "kind", "assembly", "sample", "role", "label", "donor_id", "tumor_id", "normal_sample"}
    with path.open(newline="") as source:
        reader = csv.DictReader(source, delimiter="\t" if path.suffix.lower() == ".tsv" else ",")
        if not reader.fieldnames or not {"path", "sample"} <= set(reader.fieldnames):
            raise ValueError("Manifest requires path and sample columns")
        if len(set(reader.fieldnames)) != len(reader.fieldnames) or set(reader.fieldnames) - allowed:
            raise ValueError("Duplicate or unknown manifest columns")
        for number, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()) or not row["path"] or not row["sample"]:
                raise ValueError("Malformed manifest row %d" % number)
            source_path = (path.parent / Path(row["path"]).expanduser()).resolve()
            key = (source_path, row.get("kind") or "germline", row.get("assembly") or assembly)
            spec = SampleSpec(**{k: v for k, v in row.items() if k not in {"path", "kind", "assembly"} and v})
            groups.setdefault(key, []).append(spec)
    result = []
    for (source_path, kind, build), specs in groups.items():
        result.extend(load_vcf_samples(source_path, kind=kind, assembly=build, samples=specs, config=config))
    return result


def _write_tsv(path, report):
    fields = ["sample_a", "sample_b", "germline_status", "shared_snps", "genotype_concordance",
              "ibs0_rate", "somatic_status", "shared_somatic", "somatic_overlap", "flags"]
    with path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for pair in report["pairs"]:
            g, s = pair["germline"], pair["somatic"]
            writer.writerow(dict(sample_a=pair["sample_a"], sample_b=pair["sample_b"],
                                 germline_status=g["status"], shared_snps=g.get("shared_snps"),
                                 genotype_concordance=g.get("genotype_concordance"), ibs0_rate=g.get("ibs0_rate"),
                                 somatic_status=s["status"], shared_somatic=s.get("shared_variants"),
                                 somatic_overlap=s.get("overlap_coefficient"), flags=";".join(pair["flags"])))


def main(argv=None):
    parser = argparse.ArgumentParser(prog="varcode check-samples", description=__doc__)
    for kind in ("germline", "somatic", "mixed"):
        parser.add_argument("--" + kind, action="append", default=[], type=Path, metavar="VCF",
                            help="%s VCF/VCF.gz; repeat for multiple files" % kind.capitalize())
    parser.add_argument("--manifest", type=Path, help="CSV/TSV sample selection, roles and expected identities")
    parser.add_argument("--assembly", help="GRCh37/GRCh38 when input reference headers are unrecognized")
    parser.add_argument("--config", type=Path, help="JSON object overriding SampleCheckConfig defaults")
    parser.add_argument("--json", type=Path, help="Full report path (default: stdout)")
    parser.add_argument("--tsv", type=Path, help="Optional compact table of every sample pair")
    parser.add_argument("--fail-on-mismatch", action="store_true",
                        help="Exit 1 for discordant expected donor pairs; inconclusive remains exit 0")
    args = parser.parse_args(argv)
    try:
        options = json.loads(args.config.read_text()) if args.config else {}
        if not isinstance(options, dict):
            raise ValueError("Config must be a JSON object")
        config = SampleCheckConfig(**options)
        samples = []
        for kind in ("germline", "somatic", "mixed"):
            for path in getattr(args, kind):
                samples.extend(load_vcf_samples(path, kind=kind, assembly=args.assembly, config=config))
        if args.manifest:
            samples.extend(_manifest(args.manifest, args.assembly, config))
        report = check_sample_identity(samples, config=config)
        # Never overwrite an input VCF or the manifest/config with a report.
        inputs = {Path(s.path) for s in samples} | {p.resolve() for p in (args.manifest, args.config) if p}
        outputs = [p.resolve() for p in (args.json, args.tsv) if p]
        if inputs & set(outputs) or len(set(outputs)) != len(outputs):
            raise ValueError("Output paths must differ from inputs and from one another")
        text = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
        if args.json:
            args.json.write_text(text)
        else:
            sys.stdout.write(text)
        if args.tsv:
            _write_tsv(args.tsv, report)
        print("Compared %d samples / %d pairs; %d expected donor conflicts." % (
            len(samples), len(report["pairs"]), report["identity_conflicts"]), file=sys.stderr)
        return int(args.fail_on_mismatch and report["identity_conflicts"] > 0)
    except (OSError, ValueError, TypeError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    sys.exit(main())
