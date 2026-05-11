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

- [**Effect annotation**](effect_annotation.md) — how variants become
  effects, splice outcome representations, pluggable annotators,
  possibility sets for ambiguous outcomes, and structural variants.
  Start here.
- [**Germline-aware annotation**](germline.md) — classify somatic
  variants against the patient's germline-applied transcript;
  possibility sets when phase is unknown; LOH detection. New in 4.19.
- [**Genotypes and sample-aware queries**](genotype.md) — per-sample
  zygosity on multi-sample VCFs.
- [**VariantCollection transforms**](transforms.md) — pure
  `VC -> VC` refinements; ships `pair_breakends` for collapsing
  MATEID-paired BND rows into one combined variant.
- [**CSV round-trip and metadata headers**](csv.md) — `to_csv` /
  `from_csv` with genome recovered from the header.
- [**Error handling**](errors.md) — `ReferenceMismatchError`,
  `GenomeBuildMismatchError`, `SampleNotFoundError`, and
  `raise_on_error=False`.

## API reference

The [API reference](api.md) is auto-generated from docstrings in the
source.

## Change log

See the [changelog](changelog.md) for release history.
