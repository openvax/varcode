# Check sample identities

`varcode check-samples` screens a cohort for possible sample mix-ups. It reports
two separate kinds of evidence for every pair:

- **Donor:** compatibility of shared, explicitly called germline SNP genotypes.
- **Tumor:** overlap of reported somatic alleles, including precise SV adjacencies.

This is an uncalibrated heuristic screen, not a probability of identity, a kinship
test, a contamination estimate or a definitive swap assignment. “Compatible”
means the available evidence did not distinguish the samples. Relatives and
identical twins can remain compatible. Check suspected mix-ups with an informative
SNP panel and a validated identity workflow.

## CLI

```bash
varcode check-samples \
  --germline normal1.vcf.gz --germline normal2.vcf.gz \
  --somatic tumor1.vcf.gz --somatic tumor2.vcf.gz \
  --assembly GRCh38 --json checks.json --tsv pairs.tsv
```

Each file contributes all its sample columns. Repeat `--germline`, `--somatic`,
or `--mixed` as needed. `--somatic` declares that records are caller-reported
somatic calls. Use `--mixed` for a combined germline/somatic callset. These flags
do not imply a tumor/normal role from a column's name. Multi-sample somatic or
mixed VCFs need `##tumor_sample` / `##normal_sample` header declarations or roles
in a manifest. A single-sample somatic VCF defaults to a tumor.

JSON goes to stdout unless `--json` is supplied; the short summary goes to stderr.
The optional TSV contains every pair, statuses, counts, concordance, somatic
overlap and flags. JSON also contains reasons for inconclusive results, QC counts,
absolute source paths, input SHA-256 hashes, sample roles, Varcode version and all
thresholds. It does not export individual genotypes or mutation lists.

### Declare expected identities

A CSV or `.tsv` manifest selects exact sample columns and records expectations.
Paths are relative to the manifest, not the shell's working directory:

```csv
path,sample,kind,role,label,donor_id,tumor_id,normal_sample,assembly
normal.vcf.gz,N,germline,normal,normal-A,patient-A,,,GRCh38
paired.vcf.gz,N,mixed,normal,paired-normal-A,patient-A,,,GRCh38
paired.vcf.gz,T,mixed,tumor,tumor-A,patient-A,tumor-A,N,GRCh38
recurrence.vcf.gz,T,somatic,tumor,recurrence-A,patient-A,tumor-A,,GRCh38
```

```bash
varcode check-samples --manifest samples.csv \
  --json checks.json --tsv pairs.tsv --fail-on-mismatch
```

Only `path` and `sample` are required columns. `kind` defaults to `germline`.
Labels must be unique; otherwise they default to `absolute-path::sample`.
The same file/sample cannot be supplied twice. A tumor's `normal_sample` names a
column in the same VCF, even if that normal is not selected for output. A single
declared normal is paired automatically; multiple normals require an explicit
selection. Names are never guessed from words such as `TUMOR` or `NORMAL`.

Expected identities only affect flags, never the scores:

| Flag | Meaning |
| --- | --- |
| `expected_donor_mismatch` | Discordant donor evidence despite equal donor IDs or a declared tumor/normal pairing. |
| `unexpected_donor_compatibility` | Compatible evidence despite different donor IDs; review potential swaps, relatedness and evidence limitations. |
| `expected_tumor_low_overlap` | Equal tumor IDs but low somatic overlap; review coverage, purity, calling and tumor evolution. |
| `paired_normal_has_different_donor_label` | The declared pairing and donor labels disagree. |

`--fail-on-mismatch` returns exit 1 only for `expected_donor_mismatch`. Successful
inconclusive analyses return 0; invalid inputs return 2. A low somatic overlap is
not sufficient to call a different donor or tumor. Without expectations, the
report provides pairwise evidence rather than declaring an automatic swap.

## Python API

```python
from varcode import (
    SampleCheckConfig, SampleSpec, load_vcf_samples,
    compare_samples, check_sample_identity,
)

config = SampleCheckConfig()
normals = load_vcf_samples("normal.vcf.gz", assembly="GRCh38", config=config)
paired = load_vcf_samples(
    "paired.vcf.gz", kind="mixed", assembly="GRCh38", config=config,
    samples=[
        SampleSpec("N", role="normal", donor_id="patient-A"),
        SampleSpec("T", role="tumor", donor_id="patient-A", normal_sample="N"),
    ],
)
pair = compare_samples(paired[0], paired[1])
report = check_sample_identity(normals + paired)
```

`load_vcf_samples` returns `SampleFingerprint` objects. `compare_samples` and
`check_sample_identity` return JSON-compatible dictionaries. Use the same config
when loading every sample. No transcript annotation or reference sequence download
is needed. Local VCF and gzip-compressed VCF are supported, not BCF or remote URLs.
Assemblies must be GRCh37/hg19/b37 or GRCh38/hg38/b38. Recognized `##reference`
values are used unless explicitly supplied; contradictions are errors. Different
builds produce inconclusive comparisons. There is no liftover.

## Evidence and thresholds

### Donor checks

Only complete diploid GTs at individually represented autosomal SNP positions
(chromosomes 1–22) are used. Comparison uses actual base genotypes, so phase and
ALT order do not change identity. `chr1` and `1` are equivalent. Missing records,
partial/haploid GTs, gVCF blocks, indels and SVs do not supply donor genotypes.
Individually called SNPs remain usable when an uncalled `<NON_REF>` or indel ALT
is also listed; a GT that actually calls such an allele is excluded.
An absent call is **never** imputed as homozygous reference.

Record FILTER and sample FT must be PASS or missing. Available depth must be at
least 10 and available GQ at least 20. Complete AD can supply depth when DP is
absent. Missing quality is counted and accepted by default; set
`require_quality=True` to require both depth and GQ. Malformed GT indexes,
nonfinite/negative quality and invalid AD/AF cardinality stop analysis with a file
and line number. Conflicting passing SNP records at a locus exclude that locus;
reference-base disagreements between samples are separately counted and excluded.

Default evidence requirements are 100 shared callable SNPs, 20 nonreference SNPs
in each sample, at least two reference/heterozygous/homozygous-ALT genotype categories
in each, and 20 distinct 100 kb genomic bins. These are screening safeguards,
not a claim that the sites are independent or population-informative.

For two non-tumor samples, exact genotype concordance of at least 0.98 and IBS0
(no shared allele) at most 0.01 is compatible. Concordance below 0.90 or IBS0 at
least 0.05 is discordant; intermediate results are inconclusive.

For a comparison involving a tumor, heterozygous-to-homozygous changes can reflect
loss of heterozygosity. The score uses allele sharing (`1 - IBS0 rate`) instead
of exact genotype concordance, while reporting both. It additionally requires
20 jointly homozygous sites and 20 explicit homozygous-reference calls in **each**
sample at shared sites. Without reference calls, biallelic variant-only VCFs
automatically share ALT at every overlapping site, making this tolerant test
uninformative. Such data can still contribute somatic evidence. Normal columns
containing only `0/0` calls at somatic sites also remain inconclusive for identity.

Up to 50,000 SNP loci per sample are retained by a deterministic SHA-256 bottom-k
sketch keyed only by assembly and position. This bounds memory and is independent
of genotype and record order. Comparisons report the actual retained intersection,
not estimated whole-genome counts; omissions and conflicting loci are in QC.

### Somatic checks

Only tumor samples contribute somatic alleles. The allele must be present in a
complete GT; without a called GT it needs AD support of at least 3 ALT reads,
VAF at least 0.03 and depth at least 10. GT takes precedence over contradictory AD.
AF is retained for VAF comparison but AF alone does not establish a call.

When a matched normal is declared, that allele must be absent from its GT, with
adequate normal depth. If AD is available, at most 1 ALT read and VAF at most 0.02
are required. Missing or filtered normal evidence cannot establish absence.
Otherwise, a `somatic` input declaration or a VCF `SOMATIC` flag is required.
Unlabelled nonreference alleles in a mixed tumor VCF are not assumed somatic.
QC distinguishes paired-normal contrasts from caller-reported evidence.

Small alleles are minimally trimmed, without repeat left-alignment. Precise,
sequence-resolved SVs use the existing [SV comparison normalizer](sv_comparison.md),
including reciprocal-BND deduplication. Imprecise, unknown-insertion and unsupported
SVs are excluded with counts; nearby breakpoints are not treated as identical.
Equivalent repeats or caller representations may consequently underlap.
Explicit sequence alleles and symbolic SVs use separate keys; comparisons across
those representations require the dedicated SV comparison workflow.

The report includes shared/total counts, both directional shared fractions,
Jaccard similarity and the overlap coefficient (`shared / min(count_a, count_b)`).
At least 10 alleles in each tumor, 5 shared alleles and an overlap coefficient of
0.20 are required for `shared_somatic_support`. Adequately sized sets below those
overlap thresholds are `low_overlap`; smaller sets are inconclusive. VAF Pearson
correlation is descriptive only, requires at least 5 shared VAFs and nonzero
variance, and does not determine status. Repeated alleles count once; conflicting
duplicate VAFs are omitted from correlation. A 100,000-allele limit raises an
explicit error rather than silently truncating sets.

Somatic overlap tolerates branching evolution and differing callset sizes; it
cannot prove that two samples came from the same tumor. Recurrent drivers,
germline leakage, assay overlap, purity, copy number and caller sensitivity can
confound both high and low overlap. Coverage/callable-region modeling and
population allele frequencies are not included.

### Configure and reproduce

Every threshold is a `SampleCheckConfig` field and appears in report `config`.
The CLI accepts the same names in a JSON object, for example:

```json
{"require_quality": true, "min_shared_snps": 200, "max_snps": 100000}
```

Pass it with `--config thresholds.json`. Changing thresholds changes the evidence
requirements; the defaults have **not** been calibrated on an independent clinical
cohort. Synthetic ground-truth VCF generators and regression cases are checked in
as `tests/test_sample_identity.py` alongside the tests.

The distinction between genotype concordance and tumor-specific allele-fraction
changes is informed by the primary [Somalier paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7362544/)
and its [cancer concordance guidance](https://github.com/brentp/somalier/blob/master/cancer-concordance-contamination.md).
[NGSCheckMate](https://pmc.ncbi.nlm.nih.gov/articles/PMC5499645/) illustrates a
validated SNP-panel allele-fraction approach. This implementation does not reproduce
their calibrated models, panels or accuracy claims. Input semantics follow the
[VCF specification](https://samtools.github.io/hts-specs/VCFv4.5.pdf).

## API reference

::: varcode.sample_identity
    options:
      members:
        - SampleCheckConfig
        - SampleSpec
        - SampleFingerprint
        - load_vcf_samples
        - compare_samples
        - check_sample_identity
