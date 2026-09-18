[![Tests](https://github.com/openvax/varcode/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/varcode/actions/workflows/tests.yml)
[![Coverage Status](https://coveralls.io/repos/openvax/varcode/badge.svg?branch=main&service=github)](https://coveralls.io/github/openvax/varcode?branch=main)
[![PyPI](https://img.shields.io/pypi/v/varcode.svg?maxAge=1000)](https://pypi.python.org/pypi/varcode/)
[![PyPI downloads](https://img.shields.io/pypi/dm/varcode.svg)](https://pypistats.org/packages/varcode)

# Varcode

Varcode helps you work with genomic variants in Python and predict their effects
on genes, transcripts, and protein sequences.

Load variants from VCF or MAF files, annotate them, and inspect the results.
The same interface handles small variants and structural variants; optional
germline, phasing, and RNA evidence can refine predictions.

## Installation

Requires Python 3.9 or later:

```bash
pip install varcode
pyensembl install --release 81 --species human
```

The second command downloads the annotation and transcript sequences for the
GRCh38 example below. Release 81 is pinned for reproducibility, not a requirement
to use that release for your own data. Choose an annotation matching your
input's genome build; see [reference setup](https://openvax.github.io/varcode/getting_started/#reference-data).

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

For a VCF called against the same genome build:

```python
from varcode import load_vcf

variants = load_vcf("variants.vcf", genome=81)
effects = variants.effects()
for variant, effect in effects.top_priority_effect_per_variant().items():
    print(variant.short_description, effect.short_description)
```

A variant can affect several transcripts. Priority is a summary of predicted
consequence, not proof of expression or pathogenicity. A protein sequence may be
unavailable for unresolved effects.

<a id="further-reading"></a>

## Learn more

Varcode predicts transcript/protein consequences;
[Isovar](https://github.com/openvax/isovar) handles RNA reconstruction and evidence;
[Vaxrank](https://github.com/openvax/vaxrank) evaluates protein/peptide candidates.
See [how the libraries fit together](https://openvax.github.io/varcode/library_roles/)
for their responsibilities and current integration limits.

Start with the [getting-started guide](https://openvax.github.io/varcode/getting_started/)
for file loading, result access, and saving a table. Then follow the
[task guides](https://openvax.github.io/varcode/#find-your-next-task) for sample
filtering, structural variants, phasing, and RNA evidence.

<a id="effect-types"></a>
<a id="coordinate-system"></a>

Detailed [effect types](https://openvax.github.io/varcode/effect_types/),
[coordinate conventions](https://openvax.github.io/varcode/getting_started/#coordinate-conventions),
and the [API reference](https://openvax.github.io/varcode/api/) live in the docs.

For bugs or questions, [open an issue](https://github.com/openvax/varcode/issues).
Contributions are welcome; see [CONTRIBUTING.md](CONTRIBUTING.md) and the
[changelog](CHANGELOG.md).
