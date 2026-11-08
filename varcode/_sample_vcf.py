"""Streaming VCF evidence extraction for :mod:`varcode.sample_identity`."""

from dataclasses import replace
from functools import partial
import gzip
import hashlib
import heapq
import math
from pathlib import Path

from .contigs import _annotation_contig_name
from .genotype import parse_gt_string
from .reference import infer_reference_name, ucsc_to_ensembl_reference_names
from .sample_identity import SampleFingerprint, SampleSpec
from .string_helpers import trim_shared_flanking_strings
from .vcf import _declares_rcrs_chrm
from .vcf_parsing import VCFHeader

_AUTOSOMES = frozenset(str(i) for i in range(1, 23))
_BUILDS = ("GRCh37", "GRCh38")


def _assembly(value):
    """The assembly ``value`` names, as :func:`load_vcf` reads it, or None."""
    if not value:
        return None
    try:
        name = infer_reference_name(value)
    except ValueError:
        return None
    return ucsc_to_ensembl_reference_names.get(name, name)


def _number(value):
    if isinstance(value, (list, tuple)):
        if len(value) != 1:
            raise ValueError("Expected a scalar quality/depth field")
        value = value[0]
    if value in (None, ".", ""):
        return None
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise ValueError("Quality, depth and allele fractions must be finite and nonnegative")
    return number


def _gt(value, nalleles):
    """Allele indexes of a fully called GT, or None for a whole or partial no-call."""
    alleles, _ = parse_gt_string(value)
    if any(a is not None and not 0 <= a < nalleles for a in alleles):
        raise ValueError("Invalid GT allele index: %s" % value)
    return None if None in alleles else alleles


def _evidence(row, nalleles, config):
    gt = _gt(row.get("GT"), nalleles)
    dp, gq = _number(row.get("DP")), _number(row.get("GQ"))
    ad = row.get("AD")
    if ad is not None:
        if not isinstance(ad, (list, tuple)) or len(ad) != nalleles:
            raise ValueError("AD cardinality does not match REF/ALT")
        ad = tuple(_number(n) for n in ad)
        if dp is None and all(n is not None for n in ad):
            dp = sum(ad)
    af = row.get("AF")
    if af is not None:
        af = af if isinstance(af, (list, tuple)) else [af]
        if len(af) != nalleles - 1:
            raise ValueError("AF cardinality does not match ALT")
        af = tuple(_number(n) for n in af)
        if any(n is not None and n > 1 for n in af):
            raise ValueError("AF must be in [0, 1]")
    reason = None
    if row.get("FT"):
        reason = "sample_filter"
    elif dp is not None and dp < config.min_depth:
        reason = "low_depth"
    elif gq is not None and gq < config.min_gq:
        reason = "low_gq"
    elif config.require_quality and (dp is None or gq is None):
        reason = "missing_quality"
    return dict(gt=gt, dp=dp, gq=gq, ad=ad, af=af, reason=reason)


def _vaf(evidence, index):
    ad = evidence["ad"]
    if ad is not None and all(n is not None for n in ad) and sum(ad):
        return ad[index] / sum(ad)
    af = evidence["af"]
    return af[index - 1] if af is not None else None


def _present(evidence, index, config):
    if evidence["reason"]:
        return False
    # A called genotype takes precedence over contradictory allele fractions.
    if evidence["gt"] is not None:
        return index in evidence["gt"]
    ad, vaf = evidence["ad"], _vaf(evidence, index)
    return (ad is not None and ad[index] is not None
            and ad[index] >= config.min_tumor_alt_reads
            and vaf is not None and vaf >= config.min_tumor_vaf
            and evidence["dp"] is not None and evidence["dp"] >= config.min_depth)


def _absent(evidence, index, config):
    if evidence is None or evidence["reason"]:
        return False
    gt, ad = evidence["gt"], evidence["ad"]
    if gt is not None and index in gt:
        return False
    if evidence["dp"] is None or evidence["dp"] < config.min_depth:
        return False
    if ad is not None:
        vaf = _vaf(evidence, index)
        return (ad[index] is not None and ad[index] <= config.max_normal_alt_reads
                and vaf is not None and vaf <= config.max_normal_vaf)
    return gt is not None


def _add_snp(sample, heap, locus, call):
    """Keep bottom-k hashed loci, including conflict tombstones.

    The cutoff only decreases, so an evicted locus can never reenter. A retained
    conflicting locus remains excluded regardless of subsequent record order.
    Hashing depends on assembly/locus, never genotype or input row order.
    """
    if locus in sample.snps:
        if sample.snps[locus] is not None and sample.snps[locus] != call:
            sample.snps[locus] = None
            sample.qc["conflicting_snp_loci"] += 1
        return
    priority = int.from_bytes(hashlib.sha256(
        ("%s:%s:%d" % (sample.assembly, *locus)).encode()).digest(), "big")
    if len(heap) == sample.config.max_snps:
        if priority >= -heap[0][0]:
            sample.qc["snp_sketch_omissions"] += 1
            return
        _, removed = heapq.heappop(heap)
        del sample.snps[removed]
        sample.qc["snp_sketch_omissions"] += 1
    heapq.heappush(heap, (-priority, locus))
    sample.snps[locus] = call


def _allele_key(contig, chrom, pos, ref, alt, info, assembly):
    if set(ref + alt) <= set("ACGT"):
        # Minimal representation, trimming shared bases as Variant does; no
        # repeat left-alignment without FASTA.
        ref, alt, prefix, _ = trim_shared_flanking_strings(ref, alt)
        return ("small", contig, pos + len(prefix), ref, alt) if ref or alt else None
    # Reuse the public adjacency normalizer, which names the contigs of both
    # sides itself; do not count reciprocal BND rows twice or conflate
    # inserted sequences. Imprecise/incomplete SVs are excluded.
    from .sv_comparison import compare_sv_calls
    result = compare_sv_calls([dict(call_id="record", caller="VCF", sample="sample",
                                   build=assembly, chrom=chrom, pos=pos, ref=ref, alt=alt,
                                   info=info)], max_distance=0)
    if result["groups"][0]["status"] == "exact_reported_allele":
        return ("sv", result["members"][0]["exact_id"])
    return None


def _select(header, specs, kind):
    if kind not in {"germline", "somatic", "mixed"}:
        raise ValueError("kind must be germline, somatic or mixed")
    names = header.samples
    if not names or len(set(names)) != len(names) or any(not n.strip() for n in names):
        raise ValueError("VCF must have nonempty, unique sample columns")
    declared = {}
    for role in ("tumor", "normal"):
        name = header.get_metadata(role + "_sample")
        if name:
            if name not in names or name in declared:
                raise ValueError("Invalid or conflicting declared tumor/normal sample: %s" % name)
            declared[name] = role
    specs = list(specs) if specs is not None else [SampleSpec(n) for n in names]
    if not specs or len({s.sample for s in specs}) != len(specs):
        raise ValueError("Select each sample at most once")
    selected = []
    for spec in specs:
        if spec.sample not in names:
            raise ValueError("Unknown VCF sample: %s" % spec.sample)
        role = spec.role or declared.get(spec.sample)
        if role is None:
            role = "germline" if kind == "germline" else "tumor" if len(names) == 1 and kind == "somatic" else None
        if role not in {"germline", "normal", "tumor"}:
            raise ValueError("Declare tumor/normal roles for sample %s" % spec.sample)
        if spec.sample in declared and role != declared[spec.sample]:
            raise ValueError("Role contradicts VCF header for sample %s" % spec.sample)
        if role != "tumor" and (spec.normal_sample or spec.tumor_id):
            raise ValueError("normal_sample and tumor_id apply only to tumors")
        selected.append(replace(spec, role=role))
    normals = {s.sample for s in selected if s.role == "normal"} | {n for n, r in declared.items() if r == "normal"}
    for i, spec in enumerate(selected):
        if spec.role != "tumor":
            continue
        normal = spec.normal_sample
        if normal is None:
            if len(normals) > 1:
                raise ValueError("Multiple normals: set normal_sample for %s" % spec.sample)
            normal = next(iter(normals), None)
        if normal is not None:
            if normal not in names or normal == spec.sample:
                raise ValueError("Invalid matched normal: %s" % normal)
            other_role = next((s.role for s in selected if s.sample == normal), declared.get(normal))
            if other_role not in (None, "normal", "germline"):
                raise ValueError("Matched normal is declared as a tumor: %s" % normal)
        selected[i] = replace(spec, normal_sample=normal)
    return selected


def _consume(fields, header, contig_name, samples, heaps, config):
    chrom, pos, _, ref, alt, _, filt, info_text, fmt = fields[:9]
    if int(pos) < 1:
        raise ValueError("VCF POS must be positive")
    pos = int(pos)
    contig = contig_name(chrom)
    ref = ref.upper()
    alts = [] if alt == "." else alt.split(",")
    alts = [a.upper() if not any(c in a for c in "[]<>") else a for a in alts]
    alleles = [ref] + alts
    if not ref or any(not a for a in alts) or len(set(alleles)) != len(alleles):
        raise ValueError("Empty, repeated or reference-equal ALT allele")
    for sample in samples:
        sample.qc["records"] += 1
    if filt not in {"PASS", "."}:
        for sample in samples:
            sample.qc["record_filter"] += 1
        return
    format_keys = fmt.split(":")
    if len(set(format_keys)) != len(format_keys) or any(len(s.split(":")) > len(format_keys) for s in fields[9:]):
        raise ValueError("Duplicate FORMAT keys or too many FORMAT values")
    rows = header.parse_samples(fields[9:], fmt)
    needed = {s.sample for s in samples} | {s.normal_sample for s in samples if s.normal_sample}
    evidence = {name: _evidence(rows[name], len(alleles), config) for name in needed}
    info = header.parse_info(info_text)
    # gVCF reference blocks are not hundreds of explicit SNP observations.
    end = _number(info.get("END"))
    snp_locus = (contig in _AUTOSOMES and len(ref) == 1 and ref in "ACGT"
           and (end is None or end == pos))
    keys = {}
    for sample, heap in zip(samples, heaps):
        ev = evidence[sample.sample]
        for key in ("dp", "gq"):
            if ev[key] is None:
                sample.qc["missing_" + key] += 1
        if ev["reason"]:
            sample.qc[ev["reason"]] += 1
            continue
        gt = ev["gt"]
        if (snp_locus and gt is not None and len(gt) == 2
                and all(len(alleles[i]) == 1 and alleles[i] in "ACGT" for i in gt)):
            sample.qc["eligible_snp_records"] += 1
            _add_snp(sample, heap, (contig, pos), (ref, tuple(sorted(alleles[i] for i in gt))))
        else:
            sample.qc["not_callable_diploid_autosomal_snp"] += 1
        if sample.role != "tumor" or sample.kind == "germline":
            continue
        normal = evidence.get(sample.normal_sample)
        for index, allele in enumerate(alts, 1):
            if not _present(ev, index, config):
                continue
            if sample.normal_sample:
                if not _absent(normal, index, config):
                    sample.qc["somatic_without_normal_absence"] += 1
                    continue
                basis = "paired_normal_contrast"
            elif sample.kind == "somatic" or info.get("SOMATIC") is True:
                basis = "caller_reported_somatic"
            else:
                sample.qc["nonreference_without_somatic_evidence"] += 1
                continue
            if index not in keys:
                keys[index] = _allele_key(contig, chrom, pos, ref, allele, info_text, sample.assembly)
            key = keys[index]
            if key is None:
                sample.qc["unsupported_or_imprecise_somatic_alleles"] += 1
                continue
            value = _vaf(ev, index)
            if key in sample.somatic:
                sample.qc["duplicate_somatic_alleles"] += 1
                # Repeated records do not provide a uniquely defined VAF.
                if sample.somatic[key] != value:
                    sample.somatic[key] = None
            else:
                if len(sample.somatic) >= config.max_somatic_variants:
                    raise ValueError("max_somatic_variants exceeded; increase the explicit config limit")
                sample.somatic[key] = value
            sample.qc[basis] += 1


def read_samples(path, kind, assembly, specs, config):
    path = Path(path).expanduser().resolve()
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as source:
        lines = []
        for line in source:
            lines.append(line)
            if line.startswith("#CHROM\t"):
                break
            if not line.startswith("##"):
                raise ValueError("Expected a VCF header before records: %s" % path)
        if not lines or not lines[-1].startswith("#CHROM\t"):
            raise ValueError("Missing #CHROM header: %s" % path)
        columns = lines[-1].rstrip("\r\n").split("\t")
        if columns[:9] != ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]:
            raise ValueError("Expected standard VCF columns including FORMAT")
        header = VCFHeader.from_lines(lines)
        declared = _assembly(header.get_metadata("reference"))
        explicit = _assembly(assembly)
        if assembly and explicit not in _BUILDS:
            raise ValueError("Unsupported assembly: %s (use GRCh37 or GRCh38)" % assembly)
        if explicit and declared and explicit != declared:
            raise ValueError("Explicit assembly contradicts VCF reference")
        build = explicit or declared
        if build not in _BUILDS:
            raise ValueError("Declare assembly=GRCh37 or GRCh38; VCF reference is %s"
                             % ("unrecognized" if build is None else "unsupported: " + build))
        # Name contigs as load_vcf does: chr1 and 1 are one contig, as are
        # hg38's chrM and MT, but hg19's chrM is not the rCRS that GRCh37
        # names MT unless the header gives chrM the rCRS length.
        contig_name = partial(_annotation_contig_name, reference_name=build,
                              chrm_is_rcrs=build == "GRCh37" and _declares_rcrs_chrm(header))
        selected = _select(header, specs, kind)
        samples = [SampleFingerprint(label=s.label or "%s::%s" % (path, s.sample),
                                     sample=s.sample, path=str(path), sha256=digest.hexdigest(),
                                     assembly=build, kind=kind, role=s.role, config=config,
                                     donor_id=s.donor_id, tumor_id=s.tumor_id,
                                     normal_sample=s.normal_sample) for s in selected]
        heaps = [[] for _ in samples]
        for number, line in enumerate(source, len(lines) + 1):
            if not line.strip():
                continue
            try:
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != len(columns):
                    raise ValueError("Record column count differs from header")
                _consume(fields, header, contig_name, samples, heaps, config)
            except (ValueError, TypeError, IndexError) as error:
                raise ValueError("%s:%d: %s" % (path, number, error)) from error
    return samples
