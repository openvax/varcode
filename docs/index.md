# Varcode

Varcode is a Python library that predicts what genomic variants do to genes,
transcripts, and proteins. For each variant on each overlapping transcript, it
reports the predicted consequence and, where it can be determined, the mutant
protein sequence.

**New here?** [Getting started](getting_started.md) walks through installing
reference data, annotating a variant and a VCF, reading the predicted protein,
and saving results.

## Annotate a VCF

```python
import varcode

# Use an annotation release matching your VCF's genome build.
variants = varcode.load_vcf("variants.vcf", genome=81)  # GRCh38
effects = variants.effects()
for variant, effect in effects.top_priority_effect_per_variant().items():
    print(variant.short_description, effect.short_description)
```

The [setup instructions](getting_started.md#reference-data) install the
reference data for this example. Before relying on the output, read
[how to read results](getting_started.md#how-to-read-results): predictions are
per transcript, and "top priority" means most severe, not most likely.

<a id="find-your-next-task"></a>

## Guides

**Everyday tasks**

- [Read effects and protein sequences](effect_annotation.md)
- [Select variants for a sample or compare tumor and normal](genotype.md)
- [Check a cohort for sample mix-ups](sample_identity.md)
- [Save and reload tables](csv.md)
- [Fix reference, allele, and sample-name errors](errors.md)

**Specific variant types and evidence**

- [Splice variants and their possible outcomes](splice_variants.md)
- [Structural variants and gene fusions](structural_variants.md)
- [Compare SV calls across callers and samples](sv_comparison.md)
- [Include patient germline variants](germline.md) or [phase linked variants](phasing.md)
- [Attach observed RNA structures or Exacto protein predictions](rna_structures.md)
- [Pair breakends or left-align indels](transforms.md)

<a id="go-deeper"></a>

## Reference and extensions

- [Varcode, Isovar, and Vaxrank](library_roles.md): which library does what.
- [Transcript models](transcript_models.md): cDNA, partial structures, and missing sequence.
- [Effect types](effect_types.md) and [API reference](api.md): classes and parameters.
- [Experimental annotators](experimental_annotators.md): optional alternative implementations.
- [Writing an annotator](annotator_contract.md): plugging in your own model.
- [Changelog](changelog.md): release history.
