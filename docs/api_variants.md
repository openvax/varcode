# Variants and files API

For examples, see [file loading](getting_started.md#load-a-file),
[sample queries](genotype.md), and [transforms](transforms.md).

## Reference genomes

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

The `genome_fasta=` PyEnsembl spelling requires PyEnsembl 2.28 or later.
For older supported versions or a pre-opened FASTA reader, use
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
```

```python
genome = EnsemblRelease("2020_10", species="dog")
variants = load_vcf("ros_cfam_variants.vcf", genome=genome)
```

For SV or intronic sequence analysis, install with `--with-genome-fasta`
and construct the dataset with `genome_fasta=True`, or supply a matching
local FASTA. The complete genomic DNA download requires additional disk space.

Keep CanFam3.1 calls on their original reference unless their coordinates and
alleles have been explicitly converted and validated. Selecting a new assembly
does not convert variants. Our [native dog benchmark](https://github.com/openvax/varcode/tree/main/tests/data/public_somatic)
uses published ROS_Cfam_1.0 calls; the historical CanFam3.1 comparison stays pinned.

::: varcode.Genome

## Variants

::: varcode.Variant

::: varcode.VariantCollection

::: varcode.StructuralVariant

::: varcode.parse_symbolic_alt

::: varcode.SV_TYPES

## Genotypes

::: varcode.Genotype

::: varcode.Zygosity

## VariantCollection transforms

::: varcode.transforms.pair_breakends

::: varcode.transforms.left_align_indels

## File loading

::: varcode.vcf.load_vcf

::: varcode.load_maf

::: varcode.load_maf_dataframe

## Exceptions

::: varcode.ReferenceMismatchError

::: varcode.SampleNotFoundError

::: varcode.GenomeBuildMismatchError

::: varcode.HypothesisLimitError
