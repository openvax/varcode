# Varcode

Varcode is a Python library for working with genomic variants and predicting
their effects on transcripts and proteins.

**New here?** Follow [Getting started](getting_started.md) to install reference
data, annotate a variant or a VCF, inspect the predicted protein, and save results.
You do not need to choose an annotation implementation.

## Annotate a VCF

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

<a id="find-your-next-task"></a>

## Guides

- [Read effects and protein sequences](effect_annotation.md).
- [Select variants for a sample or compare tumor and normal](genotype.md).
- [Save and reload tables](csv.md).
- [Diagnose reference and sample errors](errors.md).
- [Load structural variants and inspect fusion predictions](structural_variants.md).
- [Interpret splice alternatives](splice_variants.md).
- [Include patient germline](germline.md) or [phase linked variants](phasing.md).
- [Attach observed RNA structures or Exacto protein predictions](rna_structures.md).
- [Pair breakends or left-align indels](transforms.md).

<a id="go-deeper"></a>

## Reference and extensions

- [Varcode, Isovar, and Vaxrank](library_roles.md): responsibilities and evidence handoffs.
- [Transcript models](transcript_models.md): structures and sequence completeness.
- [Effect types](effect_types.md) and [API reference](api.md): classes and parameters.
- [Experimental annotators](experimental_annotators.md): alternative implementations.
- [Writing an annotator](annotator_contract.md): custom integrations.
- [Changelog](changelog.md): release history.
