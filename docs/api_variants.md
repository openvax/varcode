# Variants and files API

For examples, see [file loading](getting_started.md#load-a-file),
[sample queries](genotype.md), and [transforms](transforms.md).

## Reference genomes

`Genome` inherits reference DNA configured on a PyEnsembl genome (PyEnsembl
2.11 or later). This makes intronic and intergenic bases available to Varcode's
reference lookups as well as features such as indel left alignment.

```python
from pyensembl import EnsemblRelease
from varcode import Genome

native = EnsemblRelease(81, genome_fasta_path="/data/GRCh38.fa")
genome = Genome(native)
bases = genome.reference_range("7", 117_480_000, 117_480_050)
```

Construction and rewrapping do not open or download native DNA. The first
read can index a local FASTA or decompress it into PyEnsembl's cache. For a
remote source, install DNA explicitly before querying:

```python
native = EnsemblRelease(81, download_genome_fasta=True)
native.download_genome_fasta()  # explicit, potentially large download
native.index_genome_fasta()
genome = Genome(native)
```

An explicit `Genome(native, fasta="/data/other.fa")` overrides native DNA;
the supplied reference must match the annotation assembly. Without DNA,
`reference_base` and `reference_range` fall back to transcript cDNA and return
`""` for uncovered positions. `sequence` reads only chromosome DNA: configured
native DNA preserves PyEnsembl's missing-DNA and invalid-interval errors,
while unconfigured genomes return `""`. Explicit FASTA overrides retain their
legacy permissive lookup behavior. Coordinates are 1-based inclusive, and
returned bases are uppercase on the genomic plus strand.

The native PyEnsembl genome owns its reader. Calling `genome.close()` leaves
native DNA and caller-provided reader objects open; it closes only readers
that this wrapper opened from an explicit `fasta=` path. Rewrapping borrows an
explicit reader, so its owning wrapper must stay open. Close the native genome
or a caller-provided reader only after all its borrowers have finished.

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
