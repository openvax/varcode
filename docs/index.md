# varcode

Varcode is a Python library for manipulating genomic variants and
predicting their effects on protein sequences.

## Install

```bash
pip install varcode
```

Reference genome data is managed through
[PyEnsembl](https://github.com/openvax/pyensembl):

```bash
pyensembl install --release 75 76
```

## Quick example

```python
import varcode

variants = varcode.load_maf("tcga-ovarian-cancer-variants.maf")
TP53_effects = variants.groupby_gene_name()["TP53"].effects()
print(TP53_effects.top_priority_effect())
```

See the [project README](https://github.com/openvax/varcode#readme)
for a longer walkthrough and the effect-class table.

## Feature guides

- [**Effect annotation**](effect_annotation.md) — how variants
  become effects, splice outcome representations, pluggable
  annotators, possibility sets for ambiguous outcomes, and the
  SV-extensibility roadmap. Start here.
- [**Genotypes and sample-aware queries**](genotype.md) — per-sample
  zygosity on multi-sample VCFs. New in 2.3.
- [**CSV round-trip and metadata headers**](csv.md) — `to_csv` /
  `from_csv`, genome recovered from the header. New in 2.1, refined in 2.2.
- [**Error handling**](errors.md) — `ReferenceMismatchError`,
  `SampleNotFoundError`, and `raise_on_error=False`.

## API reference

The [API reference](api.md) is auto-generated from docstrings in the
source.

## Change log

See the [changelog](changelog.md) for release history.
