[![Tests](https://github.com/openvax/varcode/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/varcode/actions/workflows/tests.yml)
[![Coverage Status](https://coveralls.io/repos/openvax/varcode/badge.svg?branch=main&service=github)](https://coveralls.io/github/openvax/varcode?branch=main)
[![PyPI](https://img.shields.io/pypi/v/varcode.svg?maxAge=1000)](https://pypi.org/project/varcode/)
[![PyPI downloads](https://img.shields.io/pypi/dm/varcode.svg)](https://pypistats.org/packages/varcode)
[![Docs](https://img.shields.io/badge/docs-openvax.github.io%2Fvarcode-blue)](https://openvax.github.io/varcode/)

# Varcode

Varcode is a Python library that predicts what genomic variants do to genes,
transcripts, and proteins. Load a VCF or MAF file, and for each variant on each
overlapping transcript, Varcode reports the predicted consequence (for example
a substitution such as `p.L159M`, a frameshift, or a gene fusion). Where the
sequence can be determined, it also gives the mutant protein.

In the [OpenVax](https://github.com/openvax) tools, Varcode reads gene
annotations through [PyEnsembl](https://github.com/openvax/pyensembl).
[Isovar](https://github.com/openvax/isovar) reconstructs variant sequences from
RNA reads, and [Topiary](https://github.com/openvax/topiary) and
[Vaxrank](https://github.com/openvax/vaxrank) use Varcode's predicted proteins
to choose and rank candidate peptides; see
[how the libraries fit together](https://openvax.github.io/varcode/library_roles/).

## Install

<a id="installation"></a>Requires Python 3.11 or later:

```bash
pip install varcode
pyensembl install --release 81 --species human
```

The second command downloads the gene annotation and transcript sequences used
in the example below (Ensembl 81, GRCh38). For your own data, choose an
annotation that matches your input's genome build; see
[reference setup](https://openvax.github.io/varcode/getting_started/#reference-data).

Varcode installs [ORFdata](https://github.com/openvax/orfdata), which supplies
the ORFs for [noncanonical ORF consequences](https://openvax.github.io/varcode/uorf_annotations/).
The optional `varcode[rna]` extra adds pysam for
[phasing variants from RNA-seq reads](https://openvax.github.io/varcode/phasing/#rna-evidence).

<a id="example"></a>

## Quick start

Predict the effect of a single GRCh38 variant on a CFTR transcript:

```python
from varcode import Variant

variant = Variant("7", 117_531_100, "T", "A", genome=81)
transcript = variant.genome.transcript_by_id("ENST00000003084")
effect = variant.effect_on_transcript(transcript)

print(effect.gene.name, effect.short_description)
protein = effect.mutant_protein_sequence
print(protein[158])  # Python offsets start at zero
```

```text
CFTR p.L159M
M
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

For a VCF containing only the variant above, this prints:

```text
chr7 g.117531100T>A p.L159M
```

<a id="reading-results"></a>
<a id="reading-the-results"></a>

Before relying on the output, read
[how to read results](https://openvax.github.io/varcode/getting_started/#how-to-read-results).
In short: each prediction belongs to one transcript; "top priority" means most
severe, not most likely; `None` means unknown, not unchanged; and a predicted
protein is not evidence that it is expressed.

## Command line

The `varcode` command annotates variant files without writing Python and saves
one row per variant and transcript:

```bash
varcode effects --genome GRCh38 --vcf variants.vcf --output-csv effects.csv
```

`varcode genes` lists the genes each variant overlaps, and
`varcode check-samples` screens a cohort of VCFs for possible sample mix-ups.
`varcode --help` lists the subcommands. See the [command-line guide](https://openvax.github.io/varcode/getting_started/#command-line)
and [sample identity checks](https://openvax.github.io/varcode/sample_identity/).

<a id="further-reading"></a>
<a id="learn-more"></a>
<a id="effect-types"></a>
<a id="coordinate-system"></a>

## Documentation

Full documentation: [openvax.github.io/varcode](https://openvax.github.io/varcode/).

- [Getting started](https://openvax.github.io/varcode/getting_started/): reference
  data, one variant, a VCF, saving results, and
  [coordinate conventions](https://openvax.github.io/varcode/getting_started/#coordinate-conventions).
- [Effect annotation](https://openvax.github.io/varcode/effect_annotation/) and
  [saving results](https://openvax.github.io/varcode/csv/) (CSV, JSON, VCF).
- [Sample queries](https://openvax.github.io/varcode/genotype/) for multi-sample VCFs.
- [Splice variants](https://openvax.github.io/varcode/splice_variants/) and
  [structural variants](https://openvax.github.io/varcode/structural_variants/),
  including gene-fusion candidates.
- [Germline](https://openvax.github.io/varcode/germline/),
  [phasing](https://openvax.github.io/varcode/phasing/), and
  [observed RNA](https://openvax.github.io/varcode/rna_structures/): optional
  patient context that refines predictions.
- [Effect types](https://openvax.github.io/varcode/effect_types/) and the
  [API reference](https://openvax.github.io/varcode/api/).

## Development

```bash
./develop.sh   # editable install into the active environment
./lint.sh      # Ruff and mypy
./test.sh      # pytest; needs the Ensembl releases listed in CONTRIBUTING.md
mkdocs build --strict
```

See [CONTRIBUTING.md](https://github.com/openvax/varcode/blob/main/CONTRIBUTING.md)
for the development environment and reference data, and the
[changelog](https://github.com/openvax/varcode/blob/main/CHANGELOG.md) for
release history. For bugs or questions,
[open an issue](https://github.com/openvax/varcode/issues).

## License

Varcode is released under the
[Apache 2.0 license](https://github.com/openvax/varcode/blob/main/LICENSE).
