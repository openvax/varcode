# Variants and files API

For examples, see [file loading](getting_started.md#load-a-file),
[sample queries](genotype.md), and [transforms](transforms.md).

## Reference genomes

### Pin an annotation release

Use the same reference string in Python's `genome=` argument and the CLI's
`--genome` option:

| Reference string | Assembly and annotation |
|---|---|
| `GRCh38:93` | Human GRCh38, numbered Ensembl release 93 |
| `GRCh38:2026_04` | Human GRCh38, annotation dated April 2026 |
| `GRCh37:2013_09` | Human GRCh37, historical annotation dated September 2013 |
| `ROS_Cfam_1.0:2020_10` | Dog ROS_Cfam_1.0, annotation dated October 2020 |

Selection is offline. Install the matching annotation before using it, for example:

```bash
pyensembl install --release 2013_09 --assembly GRCh37
varcode --genome GRCh37:2013_09 --vcf variants.vcf --output-csv effects.csv
```

Dates use `YYYY_MM`. An incompatible assembly/date combination or malformed
explicit date raises an error. An assembly name without a release selects the
most recent installed **numbered** release, or the newest numbered release
when none is installed. Selecting a release does not convert coordinates.

Historical dated datasets keep PyEnsembl's accession-based reference names
and cache directories. For example, `GRCh37:2013_09` has `genome.assembly ==
"GRCh37"` and `genome.reference_name == "GCA_000001405.14"`.

### Resolve chromosome names across readers

Use the names actually present in an annotation, BAM or FASTA. An exact match
wins; otherwise the resolver requires one unambiguous alias:

```python
from varcode import ContigResolver, resolve_contig_name, normalize_contig_name

result = resolve_contig_name("MT", ["chr1", "chrM"], reference_name="GRCh38")
assert result.status == "alias"
assert result.original_name == "MT"
assert result.resolved_name == "chrM"  # spelling to pass to the target reader

assert resolve_contig_name("MT", ["M", "chrM"]).status == "ambiguous"
assert resolve_contig_name("MT", ["MT", "chrM"]).resolved_name == "MT"
assert normalize_contig_name("chrUn_gl000216", reference_name="hg19") == "GL000216.1"
```

For many queries, build the index once with
`resolver = ContigResolver(bam.references, reference_name="GRCh38")`, then call
`resolver.resolve(variant.contig)`. For annotation lookup, supply
`genome.contigs()` and prefer `genome.assembly` as the reference name when it is
available; historical annotations can use accession-based `genome.reference_name`
values. The resolver does not open or download reference data itself.
When fetching alignments for an input variant, prefer `variant.original_contig`
so an exact source spelling wins even if another alias is also in the header.
The reference context selects spelling rules; it does not validate that the BAM
and annotation use identical sequences or coordinates.

| Status | Meaning |
| --- | --- |
| `exact` | The supplied spelling occurs in the target vocabulary. |
| `alias` | Exactly one target name shares the normalized spelling. |
| `ambiguous` | Multiple aliases match; `resolved_name` is `None`. |
| `missing` | No target name matches, including an empty target vocabulary. |
| `reference_unavailable` | Target names were supplied as `None`. |

`candidates` retains the selected name or sorted ambiguous names.
`normalized_name` is a comparison key, not necessarily a name in the target file.
Unresolved results provide `diagnostic` text without emitting warnings; aggregate
these by input source/contig when processing a cohort. Missing reference names
must not be interpreted as intergenic annotation or absence of RNA evidence.

Canonical chromosome case variants (`CHR1`, `CHRX`, `CHRM`, `CHRMT`) participate
in alias resolution. If both `1` and `CHR1` occur, querying `chr1` is ambiguous;
querying either exact name retains it. Custom/viral/decoy names remain available
as supplied, including the valid SAM reference name `0`. Accepting `0` in an
external vocabulary does not relax `Variant`/annotation validation. No `EBV` ↔
GenBank-accession synonym is inferred, and alternative loci are not mapped onto
their primary chromosomes. Genuine accession versions remain distinct.

Supported spelling rules include ordinary `chr` prefixes, `M`/`MT`/`chrM`/`chrMT`,
PyEnsembl case normalization (including yeast `Mito`), and accession-bearing
UCSC scaffold names. Explicit versions remain distinct: `KI270438.1` does not
match `KI270438.2`. Unversioned UCSC GL aliases require explicit GRCh37/hg19
context (also `b37`/`NCBI37`); the resolver never infers that context from the
available names. Unknown names remain eligible for exact matches.

These are spelling rules, **not sequence or coordinate equivalence checks**.
Callers must establish compatible references before using matched names.
The assembly-specific [UCSC hg19 alias table](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chromAlias.txt)
and [hg38 alias table](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromAlias.txt)
even assign different sequence accessions to `chrM`. The
[NCBI sequence report](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/genome-sequence/)
records names and accessions within a specific assembly. This API performs no
liftover and does not infer RefSeq/GenBank accession equivalence.

::: varcode.normalize_contig_name

::: varcode.resolve_contig_name

::: varcode.ContigResolver

::: varcode.ContigResolution

### Chromosome DNA

**Provide chromosome DNA once on your reference genome.** Varcode then uses it
for SV breakpoint windows, cryptic-exon searches, intron-retention products,
ncORF consequences and indel left alignment. Without chromosome DNA, exonic
consequences still work; outcomes needing unavailable intronic/intergenic bases
remain unresolved or unavailable.

Use the same assembly for the annotation, DNA and variant coordinates:

```python
from pyensembl import EnsemblRelease
from varcode import load_vcf

genome = EnsemblRelease(81, genome_fasta="/data/GRCh38.fa")
variants = load_vcf("variants.vcf", genome=genome, parse_structural_variants=True)
effects = variants.effects(extra_orfs="all")  # optional: requires varcode[orf]
```

Construction and rewrapping do not open or download native DNA. The first
read can index a local FASTA or decompress it into PyEnsembl's cache. For a
remote source, install DNA explicitly before querying:

```python
genome = EnsemblRelease(81, genome_fasta=True)
genome.download_genome_fasta()  # explicit, potentially large download
genome.index_genome_fasta()
```

Varcode requires PyEnsembl 2.30 or later. For a pre-opened FASTA reader, use
`varcode.Genome(native, fasta="/data/GRCh38.fa")`.
An explicit `Genome(native, fasta="/data/other.fa")` overrides native DNA;
the supplied reference must match the annotation assembly.

Lookups behave as follows:

- Without DNA, `reference_base` and `reference_range` fall back to transcript
  cDNA and return `""` for uncovered positions.
- `sequence` reads only chromosome DNA. Configured native DNA preserves
  PyEnsembl's missing-DNA and invalid-interval errors, while unconfigured
  genomes return `""`.
- Explicit FASTA overrides retain their legacy permissive lookup behavior.
- Coordinates are 1-based inclusive, and returned bases are uppercase on the
  genomic plus strand.

The native PyEnsembl genome owns its reader. Calling `genome.close()` leaves
native DNA and caller-provided reader objects open; it closes only readers
that this wrapper opened from an explicit `fasta=` path. Rewrapping borrows an
explicit reader, so its owning wrapper must stay open. Close the native genome
or a caller-provided reader only after all its borrowers have finished.

The command-line equivalent is:

```bash
varcode --genome GRCh38:81 --genome-fasta /data/GRCh38.fa \
  --vcf variants.vcf --extra-orfs all --output-csv effects.csv
```

`--genome-fasta` works for VCF, MAF, individual variants and variant JSON.
All combined inputs must use one assembly. `--extra-orfs` is optional;
repeat an evidence selection (for example `--extra-orfs proteomics
--extra-orfs immunopeptidomics`) to include ORFs supported by either.

### Current dog reference

For variants called on **ROS_Cfam_1.0** (GCA_014441545.1), select the matching
dog annotation explicitly. Ensembl's current dog annotation is dated
**2020_10**, as verified on 2026-10-20; this is its build date, not the date
of the website release.

```bash
pyensembl install --release 2020_10 --species dog
varcode --genome ROS_Cfam_1.0:2020_10 --vcf ros_cfam_variants.vcf --output-csv effects.csv
```

```python
variants = load_vcf("ros_cfam_variants.vcf", genome="ROS_Cfam_1.0:2020_10")
```

For SV or intronic sequence analysis, install with `--with-genome-fasta`
and pass `EnsemblRelease("2020_10", species="dog", genome_fasta=True)` as
`genome=`, or supply a matching local FASTA. The complete genomic DNA download
requires additional disk space.

Keep CanFam3.1 calls on their original reference unless their coordinates and
alleles have been explicitly converted and validated. Selecting a new assembly
does not convert variants. Our [native dog benchmark](https://github.com/openvax/varcode/tree/main/tests/data/public_somatic)
uses published ROS_Cfam_1.0 calls; the historical CanFam3.1 comparison stays pinned.

::: varcode.Genome

### Reference sequence access

These helpers use one-based inclusive genomic coordinates and return
forward-strand sequence. An empty string means sequence is unavailable;
it does not establish that the queried region is absent.

::: varcode.reference_base

::: varcode.reference_range

## Variants

::: varcode.Variant

::: varcode.VariantCollection

::: varcode.StructuralVariant

::: varcode.parse_symbolic_alt

SV parsing raises `ValueError` for malformed consumed INFO fields, including
Boolean string fields, noninteger END values and incomplete or reversed
CIPOS/CIEND pairs. Valid signed confidence offsets are preserved on symbolic,
paired-breakend and single-breakend objects. Missing intervals remain `None`;
they are not converted to a claim of exact localization.

Direct parsing accepts the confidence pair for one ALT. The VCF loader selects
the corresponding pair from VCF 4.5 per-ALT vectors while retaining the complete
original INFO mapping and ALT indexes in metadata. Older shared two-value
intervals remain accepted. No source mapping is modified in place.

::: varcode.SV_TYPES

## Genotypes

::: varcode.Genotype

::: varcode.Zygosity

## VariantCollection transforms

::: varcode.transforms.pair_breakends

::: varcode.transforms.left_align_indels

## File loading

::: varcode.vcf.load_vcf

::: varcode.load_vcf_fast

::: varcode.load_maf

::: varcode.load_maf_dataframe

## Exceptions

::: varcode.InvalidContigError

::: varcode.ReferenceMismatchError

::: varcode.SampleNotFoundError

::: varcode.GenomeBuildMismatchError

::: varcode.HypothesisLimitError
