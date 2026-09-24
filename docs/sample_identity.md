# Check sample identities

`varcode check-samples` screens a cohort of VCFs for possible sample mix-ups,
such as a tumor labelled with the wrong patient or two swapped samples. For
every pair of samples it reports two separate kinds of evidence:

- **Donor:** do the samples' germline SNP genotypes look like they come from
  the same person?
- **Tumor:** do two tumor samples share somatic mutations, including precise
  SV junctions?

This is a screening heuristic with uncalibrated thresholds. It does not give a
probability of identity, test kinship, estimate contamination, or decide which
samples were swapped. "Compatible" means only that the available evidence did
not distinguish the samples; relatives and identical twins can be compatible.
Confirm a suspected mix-up with an informative SNP panel and a validated
identity workflow.

<a id="cli"></a>

## Run a check

```bash
varcode check-samples \
  --germline normal1.vcf.gz --germline normal2.vcf.gz \
  --somatic tumor1.vcf.gz --somatic tumor2.vcf.gz \
  --assembly GRCh38 --json checks.json --tsv pairs.tsv
```

Each file contributes all of its sample columns. Declare what each file
contains, repeating a flag as often as needed:

| Flag | Use for |
|---|---|
| `--germline` | Germline calls, such as a normal sample's VCF. |
| `--somatic` | Caller-reported somatic calls. A single-sample somatic VCF is treated as a tumor. |
| `--mixed` | A combined germline/somatic callset. |

A sample's role is never guessed from its column name. Multi-sample somatic or
mixed VCFs need `##tumor_sample` / `##normal_sample` header lines, or roles
in a [manifest](#declare-expected-identities).

Outputs:

- A short summary on stderr.
- A JSON report, written to stdout unless you pass `--json`. It includes every
  pair's statuses and counts, reasons for inconclusive results, QC counts,
  absolute input paths and SHA-256 hashes, sample roles, the Varcode version,
  and all thresholds. It does not export individual genotypes or mutation lists.
- With `--tsv`, a table with one row per pair: statuses, counts, concordance,
  somatic overlap, and flags.

## Declare expected identities

To check samples against what you expect, list them in a CSV or `.tsv`
manifest. It selects exact sample columns and records which donor and tumor
each should belong to. Paths are relative to the manifest, not the shell's
working directory:

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

- Only `path` and `sample` are required. `kind` defaults to `germline`.
- `label` names the sample in reports and must be unique. It defaults to
  `absolute-path::sample`.
- The same file and sample can't be listed twice.
- A tumor's `normal_sample` names a column in the same VCF, even if that
  normal isn't selected for output. A single declared normal is paired
  automatically; with several normals, name one explicitly.
- Roles are never guessed from words such as `TUMOR` or `NORMAL`.

Expected identities never change the scores. They only add
[flags](#flags) where the results contradict your expectations.

## Reading the results

Each pair gets a donor status and a tumor status. An inconclusive status comes
with `reasons` in the JSON report.

| Donor status | Meaning |
|---|---|
| `compatible` | Genotypes agree closely enough that the samples may come from the same person. Relatives and identical twins can also be compatible. |
| `discordant` | Genotypes disagree enough to suggest different people. |
| `inconclusive` | Too little usable shared evidence, a result between the compatible and discordant thresholds, or different genome builds. |

| Tumor status | Meaning |
|---|---|
| `shared_somatic_support` | Both samples are tumors, and they share enough somatic alleles to suggest a common tumor. |
| `low_overlap` | Both tumors have enough somatic alleles, but few are shared. On its own, this does not show a different tumor or donor. |
| `inconclusive` | One sample isn't a tumor, a tumor has too few somatic alleles, or the genome builds differ. |

Somatic overlap cannot prove that two samples came from the same tumor.
Recurrent drivers, germline leakage, assay overlap, purity, copy number, and
caller sensitivity can all raise or lower it.

<a id="flags"></a>

When you've declared expected identities, these flags mark results to review:

| Flag | Meaning |
| --- | --- |
| `expected_donor_mismatch` | Discordant donor evidence despite equal donor or tumor IDs, or a declared tumor/normal pairing. |
| `unexpected_donor_compatibility` | Compatible evidence despite different donor IDs; review potential swaps, relatedness and evidence limitations. |
| `expected_tumor_low_overlap` | Equal tumor IDs but low somatic overlap; review coverage, purity, calling and tumor evolution. |
| `paired_normal_has_different_donor_label` | The declared pairing and donor labels disagree. |
| `unexpected_tumor_overlap` | Shared somatic support despite different tumor IDs; review shared origins, recurrent mutations or possible duplicate/mislabelled samples. |
| `same_tumor_has_different_donor_labels` | Equal tumor IDs conflict with different donor IDs. |

Exit codes: with `--fail-on-mismatch`, the command returns 1 only for
`expected_donor_mismatch`. A completed analysis otherwise returns 0, even when
results are inconclusive; invalid inputs return 2. Without declared
expectations, the report gives pairwise evidence and never declares a swap.

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
`check_sample_identity` return JSON-compatible dictionaries. Use the same
config when loading every sample.

Requirements and limits:

- No transcript annotation or reference sequence download is needed.
- Local VCF and gzip-compressed VCF files are supported, not BCF or remote URLs.
- Assemblies must be GRCh37/hg19/b37 or GRCh38/hg38/b38. A recognized
  `##reference` header is used unless you supply the assembly explicitly;
  contradictions are errors.
- Samples on different builds give inconclusive comparisons. There is no liftover.

<a id="evidence-and-thresholds"></a>

## How the evidence is judged

The numbers below are defaults; each can be [changed](#configure-and-reproduce).

### Donor checks

**Which genotypes count.** Only complete diploid genotypes at individually
represented SNPs on autosomes (chromosomes 1–22; `chr1` and `1` are
equivalent). Genotypes are compared as actual bases, so phase and ALT order
don't matter.

- The record FILTER and sample FT must be `PASS` or missing.
- Depth must be at least 10 and GQ at least 20 where available. Complete AD
  can supply depth when DP is absent. Missing values are counted and accepted
  by default; set `require_quality=True` to require both.
- Missing records, partial or haploid GTs, gVCF blocks, indels, and SVs don't
  count. A missing call is **never** assumed to be homozygous reference.
- An individually called SNP stays usable when an uncalled `<NON_REF>` or
  indel ALT is also listed; a GT that actually calls such an allele is excluded.
- Loci with conflicting passing SNP records are excluded. Reference-base
  disagreements between samples are counted separately and excluded.
- Malformed GT indexes, non-finite or negative quality values, and AD/AF with
  the wrong number of values stop the analysis with a file and line number.

**Minimum evidence.** Without all of these, the status is `inconclusive`:

- 100 shared callable SNPs;
- 20 non-reference SNPs in each sample;
- at least two of the three genotype categories (homozygous reference,
  heterozygous, homozygous ALT) in each sample;
- SNPs spread over at least 20 distinct 100 kb genomic bins.

These are screening safeguards, not a claim that the sites are independent or
population-informative.

**Decision.** IBS0 is the fraction of sites where the two genotypes share no
allele.

| Comparison | Score | Compatible | Discordant |
|---|---|---|---|
| Two non-tumor samples | Exact genotype concordance | score ≥ 0.98 and IBS0 ≤ 0.01 | score < 0.90 or IBS0 ≥ 0.05 |
| Either sample is a tumor | Allele sharing (`1 - IBS0`) | score ≥ 0.98 and IBS0 ≤ 0.01 | score < 0.90 or IBS0 ≥ 0.05 |

Anything in between is `inconclusive`.

**Comparisons involving a tumor.** Loss of heterozygosity can turn a
heterozygous genotype into a homozygous one, so exact concordance would
penalize a true match. Allele sharing tolerates that change; both measures are
reported. These comparisons also require 20 sites that are homozygous in both
samples, and 20 explicit homozygous-reference calls in **each** sample at
shared sites. Without reference calls, two variant-only VCFs automatically
share an ALT allele at every overlapping site, which would make the test
meaningless. Such data can still contribute somatic evidence. A normal column
containing only `0/0` calls at somatic sites also stays inconclusive for
identity.

**Sampling.** Up to 50,000 SNP loci per sample are kept, chosen by a
deterministic SHA-256 bottom-k sketch keyed only by assembly and position.
This bounds memory and doesn't depend on genotypes or record order.
Comparisons report counts for the loci actually retained, not whole-genome
estimates; omitted and conflicting loci are listed in QC.

### Somatic checks

**Which alleles count.** Only tumor samples contribute somatic alleles.

- The allele must be present in a complete GT. Without a called GT, it needs
  at least 3 ALT reads in AD, VAF at least 0.03, and depth at least 10. GT
  takes precedence over contradictory AD. AF is kept for VAF comparison, but
  AF alone does not establish a call.
- The allele must also be shown to be somatic:
    - With a declared matched normal, the allele must be absent from the
      normal's GT with adequate normal depth and, if AD is available, at most
      1 ALT read and VAF at most 0.02. Missing or filtered normal evidence
      cannot establish absence.
    - Otherwise, the input must be declared `somatic` or the record must carry
      a VCF `SOMATIC` flag. Unlabelled non-reference alleles in a mixed tumor
      VCF are not assumed to be somatic.
- QC distinguishes paired-normal contrasts from caller-reported evidence.

**Matching alleles.** Small alleles are minimally trimmed, without repeat
left-alignment. Precise, sequence-resolved SVs use the
[SV comparison normalizer](sv_comparison.md), including reciprocal-BND
deduplication. Imprecise, unknown-insertion, and unsupported SVs are excluded
with counts, and nearby breakpoints are not treated as identical. As a result,
equivalent repeats or caller representations may undercount overlap. Explicit
sequence alleles and symbolic SVs use separate keys; comparing across those
representations requires the dedicated SV comparison workflow.

**Decision.** The report includes shared and total counts, the shared fraction
in each direction, Jaccard similarity, and the overlap coefficient
(`shared / min(count_a, count_b)`). The overlap coefficient tolerates
branching evolution and callsets of different sizes.

| Tumor status | Requirement |
|---|---|
| `shared_somatic_support` | At least 10 alleles in each tumor, at least 5 shared, and an overlap coefficient of at least 0.20 |
| `low_overlap` | At least 10 alleles in each tumor, but fewer than 5 shared or an overlap coefficient below 0.20 |
| `inconclusive` | Fewer than 10 alleles in either tumor |

VAF Pearson correlation is descriptive only and never determines status. It
requires at least 5 shared VAFs with nonzero variance. Repeated alleles count
once, and conflicting duplicate VAFs are left out of the correlation. A sample
with more than 100,000 somatic alleles raises an explicit error rather than
being silently truncated. Coverage and callable regions are not modeled, and
population allele frequencies are not used.

### Configure and reproduce

Every threshold is a `SampleCheckConfig` field and appears in the report's
`config`. The CLI accepts the same names in a JSON object, for example:

```json
{"require_quality": true, "min_shared_snps": 200, "max_snps": 100000}
```

Pass it with `--config thresholds.json`. Changing thresholds changes the
evidence requirements; the defaults have **not** been calibrated on an
independent clinical cohort. Synthetic ground-truth VCF generators and
regression cases are in `tests/test_sample_identity.py`.

The distinction between genotype concordance and tumor-specific allele-fraction
changes is informed by the primary [Somalier paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7362544/)
and its [cancer concordance guidance](https://github.com/brentp/somalier/blob/master/cancer-concordance-contamination.md).
[NGSCheckMate](https://pmc.ncbi.nlm.nih.gov/articles/PMC5499645/) illustrates a
validated SNP-panel allele-fraction approach. This implementation does not
reproduce their calibrated models, panels, or accuracy claims. Input semantics
follow the [VCF specification](https://samtools.github.io/hts-specs/VCFv4.5.pdf).

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
