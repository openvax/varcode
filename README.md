[![Tests](https://github.com/openvax/varcode/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/varcode/actions/workflows/tests.yml)
[![Coverage Status](https://coveralls.io/repos/openvax/varcode/badge.svg?branch=main&service=github)](https://coveralls.io/github/openvax/varcode?branch=main)
[![PyPI](https://img.shields.io/pypi/v/varcode.svg?maxAge=1000)](https://pypi.org/project/varcode/)
[![PyPI downloads](https://img.shields.io/pypi/dm/varcode.svg)](https://pypistats.org/packages/varcode)

# Varcode

Varcode is a Python library that predicts what genomic variants do to genes,
transcripts, and proteins. Load a VCF or MAF file, and for each variant on each
overlapping transcript, Varcode reports the predicted consequence (for example
a substitution such as `p.L159M`, a frameshift, or a gene fusion). Where the
sequence can be determined, it also gives the mutant protein.

## What it does

- **Small variants:** SNVs, indels, and multi-base substitutions are classified
  as coding changes (substitution, frameshift, premature stop, and others) or
  by where they fall (UTR, intron, splice site, intergenic).
- **Splice variants:** returns the possible outcomes, such as exon skipping,
  intron retention, or a cryptic splice site, each with its own predicted
  protein when it can be determined.
- **Structural variants:** deletions, duplications, inversions, and breakends
  get transcript consequences and gene-fusion candidates.
- **Patient context (optional):** germline variants, phase from a phased VCF
  or RNA reads, and observed RNA structures refine predictions.
- **Sample checks:** screens a cohort for possible sample mix-ups using donor
  genotypes and shared somatic mutations.

## Installation

Requires Python 3.9 or later:

```bash
pip install varcode
pyensembl install --release 81 --species human
```

The second command downloads the gene annotation and transcript sequences used
in the example below (Ensembl 81, GRCh38). For your own data, choose an
annotation that matches your input's genome build; see
[reference setup](https://openvax.github.io/varcode/getting_started/#reference-data).

<a id="example"></a>

## Quick start

Predict the effect of a single GRCh38 variant on a CFTR transcript:

```python
from varcode import Variant

variant = Variant("7", 117_531_100, "T", "A", genome=81)
transcript = variant.genome.transcript_by_id("ENST00000003084")
effect = variant.effect_on_transcript(transcript)

print(effect.short_description)  # p.L159M
protein = effect.mutant_protein_sequence
```

Annotate every variant in a VCF called against the same genome build, and
print the most severe effect for each:

```python
from varcode import load_vcf

variants = load_vcf("variants.vcf", genome=81)
effects = variants.effects()
for variant, effect in effects.top_priority_effect_per_variant().items():
    print(variant.short_description, effect.short_description)
```

Or from the command line:

```bash
varcode --genome GRCh38 --vcf variants.vcf --output-csv effects.csv
```

<a id="reading-results"></a>

## Reading the results

- **One variant, many transcripts.** Each prediction belongs to one
  transcript, and a variant's effect can differ between transcripts.
  `top_priority_effect_per_variant()` picks one per variant for a summary.
- **"Top priority" means most severe, not most likely.** The ranking is
  Varcode's ordering of consequence severity. It is not a probability, a
  pathogenicity score, or a clinical classification.
- **Some effects have several possible outcomes.** Splice variants,
  structural variants, and variants with unknown phase keep all their
  candidate outcomes rather than guessing one.
- **`None` means unknown, not unchanged.** A missing protein sequence means
  Varcode could not determine it.
- **Predictions are not observations.** A predicted protein doesn't show that
  the transcript or protein is expressed.

<a id="further-reading"></a>

## Learn more

The [documentation](https://openvax.github.io/varcode/) starts with a
[getting-started guide](https://openvax.github.io/varcode/getting_started/)
and then covers sample filtering, splice and structural variants, germline
and phasing, RNA evidence, and saving results.

<a id="effect-types"></a>
<a id="coordinate-system"></a>

Reference pages cover [effect types](https://openvax.github.io/varcode/effect_types/),
[coordinate conventions](https://openvax.github.io/varcode/getting_started/#coordinate-conventions),
and the [API](https://openvax.github.io/varcode/api/).

Varcode is part of the [OpenVax](https://github.com/openvax) tools.
[Isovar](https://github.com/openvax/isovar) reconstructs variant sequences from
RNA reads, and [Vaxrank](https://github.com/openvax/vaxrank) evaluates the resulting
protein/peptide candidates; see
[how the libraries fit together](https://openvax.github.io/varcode/library_roles/).

For bugs or questions, [open an issue](https://github.com/openvax/varcode/issues).
Contributions are welcome; see
[CONTRIBUTING.md](https://github.com/openvax/varcode/blob/main/CONTRIBUTING.md)
and the [changelog](https://github.com/openvax/varcode/blob/main/CHANGELOG.md).
