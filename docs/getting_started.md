# Getting started

This guide takes you from installation to a predicted effect and protein
sequence, then shows how to annotate a file and read the results. The examples
use Python 3.9 or later.

## Install

```bash
pip install varcode
```

## Reference data

Varcode uses [PyEnsembl](https://github.com/openvax/pyensembl) for gene annotations
and transcript/protein sequences. Download the release used in these examples:

```bash
pyensembl install --release 81 --species human
```

This is a one-time download. The examples pin Ensembl 81 (GRCh38) so their
coordinates and transcript IDs are reproducible; you don't need to use that
release for your own data.

For your own data, use an annotation that matches the genome build your
variants were called against. GRCh37 input, for example, needs a GRCh37
annotation such as Ensembl 75, installed separately. Pass a release number
(`genome=81`), an assembly with a release (`genome="GRCh38:93"` or
`"GRCh38.93"`), or a
PyEnsembl genome object to pin the annotation exactly. An assembly name alone,
such as `genome="GRCh38"`, uses the most recent installed release of that
assembly. A release that does not provide the assembly, such as `"GRCh38:75"`,
is an error.

## Annotate one variant

```python
from varcode import Variant

variant = Variant("7", 117_531_100, "T", "A", genome=81)
transcript = variant.genome.transcript_by_id("ENST00000003084")
effect = variant.effect_on_transcript(transcript)

print(effect.gene.name)          # CFTR
print(effect.short_description)  # p.L159M
```

This predicts a leucine-to-methionine substitution at amino acid 159 of the
selected transcript. To annotate all overlapping transcripts instead, use
`variant.effects()`. A variant may have different consequences on different
transcripts, so keep the transcript identity with each result.

## Read the predicted protein

```python
protein = effect.mutant_protein_sequence
if protein is not None:
    print(protein[158])  # M: Python offsets start at zero
```

`None` means the sequence could not be determined for that effect, not that
the protein is unchanged. Some effects contain several possible outcomes, each
with its own protein; see [alternative outcomes](effect_annotation.md#alternative-outcomes).

## Load a file

Replace the path with your own VCF called against GRCh38:

```python
from varcode import load_vcf

variants = load_vcf("variants.vcf", genome=81)
effects = variants.effects()
```

`load_vcf` loads passing records by default. Use `only_passing=False` if you
also want records marked by the caller's FILTER field. If your input contains
structural variants, add `parse_structural_variants=True`; otherwise symbolic
and breakend records are skipped with a warning. See the
[structural variant guide](structural_variants.md#basic-usage).

For a MAF, use `variants = varcode.load_maf("variants.maf")` after importing
`varcode`; the MAF's build information supplies the reference. See
[file-loading parameters](api_variants.md#file-loading) for format-specific options.

## Summarize and save results

Choose one effect per variant when you need a compact report:

```python
for variant, effect in effects.top_priority_effect_per_variant().items():
    print(variant.short_description, effect.short_description)

effects.to_csv("effects.csv")
```

Keep the full `effects` collection when you need every transcript's prediction
or the alternative outcomes. Note that `effects.top_priority_effect()` (without
`_per_variant`) returns a single effect for the entire collection.

CSV is a report format, not a complete archive: it does not keep every
alternative outcome or the evidence behind it. Keep your input files too. See
[saving and reloading results](csv.md) for what survives a round trip.

## How to read results

Varcode's main objects:

| Object | What it is |
|---|---|
| `Variant` | One alternate allele at a genomic position. A VCF row with two ALT alleles becomes two variants. |
| `VariantCollection` | The variants loaded from your files, with their sample genotypes and source metadata. |
| Effect | The predicted consequence of one variant on one transcript, such as `Substitution` (`p.L159M`), `FrameShift`, or `Intronic`. All effects subclass `MutationEffect`. |
| `EffectCollection` | The effects for a collection of variants. Filter, summarize, or save it. |
| Multi-outcome effect | An effect with several possible `candidates`, used for splice variants, structural variants, and unknown phase. See [alternative outcomes](effect_annotation.md#alternative-outcomes). |
| `Unresolved` | Varcode could not determine the consequence from the available sequence and annotation. |

When interpreting them:

- **Predictions are per transcript.** Keep the transcript with each result;
  the same variant can be missense on one isoform and noncoding on another.
- **Priority means severity.** `top_priority_effect_per_variant()` uses
  Varcode's consequence ordering. It is not the likeliest outcome, a
  pathogenicity score, or a clinical classification.
- **`None` means unknown.** A missing protein sequence, or a `None` answer to
  "does this change the protein?", is not the same as "unchanged".
- **Predictions are not observations.** A predicted protein is not evidence
  that the transcript or protein is expressed.

## Command line

The `varcode` command annotates files without writing Python, and
`varcode-genes` lists the genes each variant overlaps:

```bash
varcode --genome GRCh38 --vcf variants.vcf --output-csv effects.csv
varcode-genes --genome GRCh38 --variant chr12 25245350 C T --output-csv genes.csv
```

Inputs can be combined and repeated: `--vcf`, `--maf`, `--variant`, and
`--json-variants`. `--genome GRCh38:93` pins Ensembl release 93, as in Python. The commands differ from `load_vcf` in two ways:

- They load structural variants automatically.
- They skip records whose FILTER is not `PASS` or `.`, and report how many were
  skipped. Use `--include-filtered` to keep them.

An annotation error stops the command. `--skip-errors` continues and records
each failure in the output; see [continuing past errors](errors.md#continuing-past-errors).
For sample mix-up screening, see [`varcode check-samples`](sample_identity.md).

## Coordinate conventions

Genomic positions are 1-based and ranges are inclusive. Specify variant alleles
on the reference genome's forward strand, even for genes on the reverse strand.
Protein/cDNA offsets exposed for Python slicing are 0-based. Structural records
also distinguish [junction coordinates from affected spans](sv_reference.md#alleles-coordinates-and-exports).

## Next steps

- [Effect annotation](effect_annotation.md): inspect results, sequences, and alternatives.
- [Sample-aware queries](genotype.md): filter a multi-sample VCF.
- [Sample identity checks](sample_identity.md): compare donor genotypes and somatic tumor calls.
- [Structural variants](structural_variants.md): SVs, fusions, and observed RNA.
- [Germline and phasing](germline.md): use a patient's baseline and linked variants.
- [Troubleshooting](errors.md): check genome build, alleles, and sample names.

Detailed signatures and options are in the [API reference](api.md).
