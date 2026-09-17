# Varcode

Varcode is a Python library for working with genomic variants and predicting
their effects on transcripts and proteins.

**New here?** Follow [Getting started](getting_started.md) to install reference
data, annotate a variant or a VCF, inspect the predicted protein, and save results.
You do not need to choose an annotation implementation.

## The usual workflow

```python
import varcode

# Use an annotation release matching your VCF's genome build.
variants = varcode.load_vcf("variants.vcf", genome=81)  # GRCh38
effects = variants.effects()
for variant, effect in effects.top_priority_effect_per_variant().items():
    print(variant.short_description, effect.short_description)
```

The [setup instructions](getting_started.md#reference-data) install the reference
data for this example. Predictions are transcript-specific; a top-priority
effect is a useful summary, not a measure of likelihood or clinical significance.

## Find your next task

- [Read effects and protein sequences](effect_annotation.md).
- [Select variants for a sample or compare tumor and normal](genotype.md).
- [Save and reload tables](csv.md).
- [Diagnose reference and sample errors](errors.md).
- [Load structural variants and inspect fusion predictions](structural_variants.md).
- [Interpret splice alternatives](effect_annotation.md#splice-disrupting-variants).
- [Include patient germline and phasing](germline.md).
- [Attach observed RNA structures or Exacto protein predictions](structural_variants.md#importing-observed-rna-structures).
- [Pair breakends or left-align indels](transforms.md).

## Go deeper

The guides begin with usage, then explain evidence, limitations, and detailed
behavior. [Experimental annotators](effect_annotation.md#annotator-selection)
are optional; most users can skip them.

Use the [effect-type reference](effect_types.md) for class definitions and the
[API reference](api.md) for signatures and parameters. The
[annotator contract](annotator_contract.md) is for integration authors, and the
[changelog](changelog.md) records release history.
