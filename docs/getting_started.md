# Getting started

This guide takes you from installation to a predicted effect and protein sequence,
then shows how to annotate a file. The examples use Python 3.9 or later.

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

This is a one-time download and indexing step. These examples pin Ensembl 81
(GRCh38) so the transcript and coordinates are reproducible; it is not a
recommendation to use that historical release for every analysis.

For your own data, match the annotation to the reference assembly used for
variant calling. Pass an explicit release number or PyEnsembl genome object
when reproducibility matters. An assembly name such as `genome="GRCh38"` is
also accepted, but does not itself pin an annotation release. GRCh37 input
needs a GRCh37 annotation, such as Ensembl 75, installed separately.

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

Protein sequence is available through the ordinary interface; no special
annotator selection is needed. `None` means no sequence is available for that
effect, not that the protein is unchanged. Some effects contain several possible
outcomes; see [alternative outcomes](effect_annotation.md#alternative-outcomes).
These are predictions, not evidence that a protein was expressed.

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

The summary uses Varcode's consequence priority, not the likelihood of an outcome
or a pathogenicity assessment. Keep the full `effects` collection when you need
all transcript predictions or alternatives. `effects.top_priority_effect()`
returns just one effect across the entire collection, not one per variant.

CSV is an inspection/report format, not a complete archive of every prediction.
Keep your input variants and evidence. See [saving and reloading results](csv.md)
for round-trip limits, especially for structural and multi-outcome effects.

## Coordinate conventions

Genomic positions are 1-based and ranges are inclusive. Specify variant alleles
on the reference genome's forward strand, even for genes on the reverse strand.
Protein/cDNA offsets exposed for Python slicing are 0-based. Structural records
also distinguish [junction coordinates from affected spans](sv_reference.md#alleles-coordinates-and-exports).

## Next steps

- [Effect annotation](effect_annotation.md): inspect results, sequences, and alternatives.
- [Sample-aware queries](genotype.md): filter a multi-sample VCF.
- [Structural variants](structural_variants.md): SVs, fusions, and observed RNA.
- [Germline and phasing](germline.md): use a patient's baseline and linked variants.
- [Troubleshooting](errors.md): check genome build, alleles, and sample names.

Detailed signatures and options are in the [API reference](api.md).
