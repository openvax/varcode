"""Heuristic donor concordance and somatic overlap, without absent-call imputation.

See ``docs/sample_identity.md`` for assumptions, thresholds and limitations.
These checks nominate mix-ups for review; they do not certify identity.
"""

from collections import Counter
from dataclasses import asdict, dataclass, field
from itertools import combinations
import math
from typing import Optional


@dataclass(frozen=True)
class SampleCheckConfig:
    """Configurable, uncalibrated screening thresholds (all reported in output).

    Missing DP/GQ are counted but accepted unless ``require_quality`` is true.
    SNPs use a deterministic, genotype-independent bottom-k locus sketch.
    The somatic limit raises an error instead of silently truncating evidence.
    """

    min_depth: int = 10
    min_gq: int = 20
    require_quality: bool = False
    min_shared_snps: int = 100
    min_nonreference_snps: int = 20
    min_homozygous_snps: int = 20
    min_tumor_hom_reference_snps: int = 20
    min_genomic_bins: int = 20
    genomic_bin_size: int = 100_000
    max_snps: int = 50_000
    max_somatic_variants: int = 100_000
    compatible_concordance: float = 0.98
    discordant_concordance: float = 0.90
    max_compatible_ibs0: float = 0.01
    min_discordant_ibs0: float = 0.05
    min_somatic_variants: int = 10
    min_shared_somatic: int = 5
    min_somatic_overlap: float = 0.20
    min_tumor_alt_reads: int = 3
    min_tumor_vaf: float = 0.03
    max_normal_alt_reads: int = 1
    max_normal_vaf: float = 0.02

    def __post_init__(self):
        for name, value in asdict(self).items():
            if name == "require_quality":
                if not isinstance(value, bool):
                    raise ValueError("require_quality must be boolean")
            elif isinstance(getattr(type(self), name), int):
                minimum = 0 if name in {"min_depth", "min_gq", "max_normal_alt_reads"} else 1
                if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
                    raise ValueError("%s must be an integer >= %d" % (name, minimum))
            elif isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 <= value <= 1:
                raise ValueError("%s must be a finite fraction in [0, 1]" % name)
        if self.discordant_concordance >= self.compatible_concordance:
            raise ValueError("discordant_concordance must be below compatible_concordance")
        if self.max_compatible_ibs0 >= self.min_discordant_ibs0:
            raise ValueError("max_compatible_ibs0 must be below min_discordant_ibs0")
        if self.max_snps < self.min_shared_snps:
            raise ValueError("max_snps must be >= min_shared_snps")


@dataclass(frozen=True)
class SampleSpec:
    """Select a VCF column and optionally declare its role and expected identity.

    Parameters
    ----------
    sample : str
        Exact VCF column name. Names are never interpreted as biological roles.
    role : str, optional
        ``germline``, ``normal`` or ``tumor``. Ambiguous somatic/mixed inputs
        require roles here or in tumor_sample/normal_sample header metadata.
    label : str, optional
        Unique report identifier; defaults to absolute path plus sample name.
    donor_id, tumor_id : str, optional
        Expected identities used to flag contradictions, not scoring inputs.
    normal_sample : str, optional
        Matched normal column in the same VCF (may be unselected for output).
    """

    sample: str
    role: Optional[str] = None
    label: Optional[str] = None
    donor_id: Optional[str] = None
    tumor_id: Optional[str] = None
    normal_sample: Optional[str] = None


@dataclass
class SampleFingerprint:
    """Loaded evidence with input provenance and exclusion counts.

    Construct with :func:`load_vcf_samples`. ``snps`` maps (chromosome,
    position) to (reference base, sorted diploid base genotype), or None for a
    conflicting locus. ``somatic`` maps normalized allele keys to VAF or None.
    Report serialization omits these potentially large, identifying maps.
    """

    label: str
    sample: str
    path: str
    sha256: str
    assembly: str
    kind: str
    role: str
    config: SampleCheckConfig
    donor_id: Optional[str] = None
    tumor_id: Optional[str] = None
    normal_sample: Optional[str] = None
    snps: dict = field(default_factory=dict, repr=False)
    somatic: dict = field(default_factory=dict, repr=False)
    qc: Counter = field(default_factory=Counter)

    def summary(self):
        """Return JSON-compatible metadata and QC, excluding genotypes."""
        result = {key: getattr(self, key) for key in (
            "label", "sample", "path", "sha256", "assembly", "kind", "role",
            "donor_id", "tumor_id", "normal_sample")}
        result.update(qc=dict(self.qc), retained_snps=sum(v is not None for v in self.snps.values()),
                      somatic_variants=len(self.somatic))
        return result


def load_vcf_samples(path, *, kind="germline", assembly=None, samples=None, config=None):
    """Stream a local VCF/VCF.gz into selected sample fingerprints.

    Parameters
    ----------
    path : str or pathlib.Path
        Local VCF, optionally gzip compressed. No annotation data is needed.
    kind : str
        ``germline``, ``somatic`` (caller-reported somatic records), or ``mixed``.
        Mixed tumor calls require SOMATIC or an explicit matched-normal contrast.
    assembly : str, optional
        GRCh37/hg19/b37 or GRCh38/hg38. Required if the header is unrecognized.
        An explicit value may not contradict a recognized header reference.
    samples : iterable of SampleSpec, optional
        Selection, roles and expectations. Default: all VCF sample columns.
    config : SampleCheckConfig, optional
        Shared loading and scoring thresholds for the cohort.

    Returns
    -------
    list of SampleFingerprint
        One entry per selected column, including samples with no usable calls.
    """
    from ._sample_vcf import read_samples
    return read_samples(path, kind, assembly, samples, config or SampleCheckConfig())


def _fraction(numerator, denominator):
    return numerator / denominator if denominator else None


def _germline(a, b, config):
    shared = a.snps.keys() & b.snps.keys()
    pairs = [(locus, a.snps[locus], b.snps[locus]) for locus in shared
             if a.snps[locus] is not None and b.snps[locus] is not None]
    ref_conflicts = sum(x[0] != y[0] for _, x, y in pairs)
    pairs = [(locus, x, y) for locus, x, y in pairs if x[0] == y[0]]
    total = len(pairs)
    exact = ibs0 = homo = homo_exact = nonref_a = nonref_b = homref_a = homref_b = 0
    patterns_a, patterns_b, bins = set(), set(), set()
    for (chrom, pos), (ref, x), (_, y) in pairs:
        exact += x == y
        ibs0 += not (set(x) & set(y))
        both_homo = x[0] == x[1] and y[0] == y[1]
        homo += both_homo
        homo_exact += both_homo and x == y
        nonref_a += x != (ref, ref)
        nonref_b += y != (ref, ref)
        homref_a += x == (ref, ref)
        homref_b += y == (ref, ref)
        patterns_a.add(0 if x == (ref, ref) else 2 if x[0] == x[1] else 1)
        patterns_b.add(0 if y == (ref, ref) else 2 if y[0] == y[1] else 1)
        bins.add((chrom, (pos - 1) // config.genomic_bin_size))
    concordance, opposite = _fraction(exact, total), _fraction(ibs0, total)
    tumor = "tumor" in {a.role, b.role}
    result = dict(status="inconclusive", shared_snps=total, exact_matches=exact,
                  genotype_concordance=concordance, ibs0=ibs0, ibs0_rate=opposite,
                  both_homozygous=homo, homozygous_concordance=_fraction(homo_exact, homo),
                  nonreference_a=nonref_a, nonreference_b=nonref_b, genomic_bins=len(bins),
                  homozygous_reference_a=homref_a, homozygous_reference_b=homref_b,
                  reference_conflicts=ref_conflicts, tumor_tolerant=tumor)
    reasons = []
    if total < config.min_shared_snps:
        reasons.append("too_few_shared_snps")
    if min(nonref_a, nonref_b) < config.min_nonreference_snps:
        reasons.append("too_few_nonreference_snps")
    if min(len(patterns_a), len(patterns_b)) < 2:
        reasons.append("insufficient_genotype_diversity")
    if len(bins) < config.min_genomic_bins:
        reasons.append("insufficient_genomic_spread")
    if tumor and homo < config.min_homozygous_snps:
        reasons.append("too_few_shared_homozygous_snps")
    if tumor and min(homref_a, homref_b) < config.min_tumor_hom_reference_snps:
        # Two nonreference-only biallelic callsets necessarily share ALT at
        # every overlapping site, so allele-sharing cannot establish identity.
        reasons.append("too_few_explicit_tumor_comparison_reference_calls")
    if not reasons:
        # Allele sharing tolerates a heterozygote becoming homozygous in a tumor.
        # It cannot distinguish some relatives; 'compatible' is not identity.
        score = 1 - opposite if tumor else concordance
        if opposite >= config.min_discordant_ibs0 or score < config.discordant_concordance:
            result["status"] = "discordant"
        elif opposite <= config.max_compatible_ibs0 and score >= config.compatible_concordance:
            result["status"] = "compatible"
        else:
            reasons.append("intermediate_concordance")
    result["reasons"] = reasons
    return result


def _somatic(a, b, config):
    shared = a.somatic.keys() & b.somatic.keys()
    na, nb, ns = len(a.somatic), len(b.somatic), len(shared)
    overlap = _fraction(ns, min(na, nb))
    result = dict(status="inconclusive", variants_a=na, variants_b=nb, shared_variants=ns,
                  fraction_a=_fraction(ns, na), fraction_b=_fraction(ns, nb),
                  jaccard=_fraction(ns, na + nb - ns), overlap_coefficient=overlap)
    values = [(a.somatic[key], b.somatic[key]) for key in sorted(shared)
              if a.somatic[key] is not None and b.somatic[key] is not None]
    correlation = None
    if (len(values) >= config.min_shared_somatic
            and len({x for x, _ in values}) > 1 and len({y for _, y in values}) > 1):
        mx, my = (sum(v[i] for v in values) / len(values) for i in (0, 1))
        dx, dy = [x - mx for x, _ in values], [y - my for _, y in values]
        denominator = math.sqrt(sum(x * x for x in dx) * sum(y * y for y in dy))
        if denominator:
            correlation = max(-1.0, min(1.0, sum(x * y for x, y in zip(dx, dy)) / denominator))
    result.update(vaf_pairs=len(values), vaf_correlation=correlation)
    if a.role != "tumor" or b.role != "tumor":
        result["reasons"] = ["requires_two_tumors"]
    elif min(na, nb) < config.min_somatic_variants:
        result["reasons"] = ["too_few_somatic_variants"]
    else:
        result["status"] = ("shared_somatic_support" if ns >= config.min_shared_somatic
                            and overlap >= config.min_somatic_overlap else "low_overlap")
        result["reasons"] = []
    return result


def compare_samples(a, b, *, config=None):
    """Return pairwise donor/tumor evidence and expectation flags as a dict.

    Both fingerprints must have been loaded with the same config. Different
    assemblies yield inconclusive results. ``flags`` are review prompts;
    only ``expected_donor_mismatch`` sets ``identity_conflict`` to true.
    """
    config = config or a.config
    if a.config != config or b.config != config:
        raise ValueError("Load all samples with the comparison config")
    if a.assembly != b.assembly:
        germline = dict(status="inconclusive", reasons=["assembly_mismatch"])
        somatic = dict(status="inconclusive", reasons=["assembly_mismatch"])
    else:
        germline, somatic = _germline(a, b, config), _somatic(a, b, config)
    paired = a.path == b.path and (a.normal_sample == b.sample or b.normal_sample == a.sample)
    same_donor = a.donor_id == b.donor_id if a.donor_id and b.donor_id else None
    flags = []
    if paired and same_donor is False:
        flags.append("paired_normal_has_different_donor_label")
    if (same_donor is True or paired) and germline["status"] == "discordant":
        flags.append("expected_donor_mismatch")
    if same_donor is False and germline["status"] == "compatible":
        flags.append("unexpected_donor_compatibility")
    same_tumor = a.tumor_id == b.tumor_id if a.tumor_id and b.tumor_id else None
    if same_tumor is True and somatic["status"] == "low_overlap":
        flags.append("expected_tumor_low_overlap")
    return dict(sample_a=a.label, sample_b=b.label, expected_same_donor=same_donor,
                paired_tumor_normal=paired, expected_same_tumor=same_tumor,
                germline=germline, somatic=somatic, flags=flags,
                identity_conflict="expected_donor_mismatch" in flags)


def check_sample_identity(samples, *, config=None):
    """Check every unordered pair and return a JSON-compatible cohort report.

    At least two uniquely labelled samples are required. No automatic renaming
    or definitive swap assignment is performed. Config, sample QC, provenance,
    all pairs and heuristic limitations are included for reproducibility.
    """
    from .version import __version__
    samples = list(samples)
    if len(samples) < 2:
        raise ValueError("At least two samples are required")
    if len({s.label for s in samples}) != len(samples):
        raise ValueError("Sample labels must be unique across the cohort")
    if len({(s.path, s.sample) for s in samples}) != len(samples):
        raise ValueError("The same VCF sample was supplied more than once")
    config = config or samples[0].config
    pairs = [compare_samples(a, b, config=config) for a, b in combinations(samples, 2)]
    return dict(schema_version=1, varcode_version=__version__, config=asdict(config),
                samples=[s.summary() for s in samples], pairs=pairs,
                identity_conflicts=sum(p["identity_conflict"] for p in pairs),
                limitations=["Uncalibrated screening heuristics, not identity probabilities or proof.",
                             "Absent variants are unknown, not reference calls.",
                             "Relatives, twins, LOH, contamination and assay differences can confound results.",
                             "Low somatic overlap does not establish a different donor or tumor."])
