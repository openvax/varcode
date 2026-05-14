# API reference

Auto-generated from in-source docstrings via
[mkdocstrings](https://mkdocstrings.github.io/).

## Variants

### `varcode.Variant`

::: varcode.Variant

### `varcode.VariantCollection`

::: varcode.VariantCollection

### `varcode.StructuralVariant`

::: varcode.StructuralVariant

### `varcode.parse_symbolic_alt`

::: varcode.parse_symbolic_alt

### `varcode.SV_TYPES`

::: varcode.SV_TYPES

## Genotypes

### `varcode.Genotype`

::: varcode.Genotype

### `varcode.Zygosity`

::: varcode.Zygosity

## Effects

### `varcode.MutationEffect`

::: varcode.MutationEffect

### `varcode.NonsilentCodingMutation`

::: varcode.NonsilentCodingMutation

### `varcode.MultiOutcomeEffect`

::: varcode.MultiOutcomeEffect

### `varcode.EffectCollection`

::: varcode.EffectCollection

### `varcode.EffectCandidate`

::: varcode.EffectCandidate

### Priority ordering

::: varcode.effect_priority

::: varcode.top_priority_effect

## Mutant transcripts

### `varcode.MutantTranscript`

::: varcode.MutantTranscript

### `varcode.apply_variant_to_transcript`

::: varcode.apply_variant_to_transcript

### `varcode.apply_variants_to_transcript`

::: varcode.apply_variants_to_transcript

## Annotators

### `varcode.EffectAnnotator`

::: varcode.EffectAnnotator

### `varcode.FastEffectAnnotator`

::: varcode.FastEffectAnnotator

### `varcode.ProteinDiffEffectAnnotator`

::: varcode.ProteinDiffEffectAnnotator

### `varcode.StructuralVariantAnnotator`

::: varcode.StructuralVariantAnnotator

### Registry

::: varcode.register_annotator

::: varcode.get_annotator

::: varcode.get_default_annotator

::: varcode.set_default_annotator

::: varcode.use_annotator

## Phasing

### `varcode.ReadPhasingSource`

::: varcode.ReadPhasingSource

### `varcode.MutantTranscriptSource`

::: varcode.MutantTranscriptSource

### `varcode.ReadPhaseResolver`

::: varcode.ReadPhaseResolver

### `varcode.VCFPhaseResolver`

::: varcode.VCFPhaseResolver

### `varcode.apply_phase_resolver_to_effects`

::: varcode.apply_phase_resolver_to_effects

## RNA evidence

### `varcode.RNAEvidenceResolver`

::: varcode.RNAEvidenceResolver

### `varcode.NullRNAEvidenceResolver`

::: varcode.NullRNAEvidenceResolver

### `varcode.apply_rna_evidence_to_effects`

::: varcode.apply_rna_evidence_to_effects

### `varcode.make_rna_candidate`

::: varcode.make_rna_candidate

## Germline-aware annotation

### `varcode.GermlineContext`

::: varcode.GermlineContext

### `varcode.Completeness`

::: varcode.Completeness

### `varcode.predict_germline_aware_effect`

::: varcode.predict_germline_aware_effect

### `varcode.apply_germline_to_transcript`

::: varcode.apply_germline_to_transcript

### `varcode.enumerate_phase_hypotheses`

::: varcode.enumerate_phase_hypotheses

### `varcode.detect_loh`

::: varcode.detect_loh

### `varcode.default_germline_window`

::: varcode.default_germline_window

## VariantCollection transforms

### `varcode.transforms.pair_breakends`

::: varcode.transforms.pair_breakends

### `varcode.transforms.left_align_indels`

::: varcode.transforms.left_align_indels

## File loading

### `varcode.load_vcf`

::: varcode.vcf.load_vcf

### `varcode.load_maf`

::: varcode.load_maf

::: varcode.load_maf_dataframe

## Exceptions

### `varcode.ReferenceMismatchError`

::: varcode.ReferenceMismatchError

### `varcode.SampleNotFoundError`

::: varcode.SampleNotFoundError

### `varcode.GenomeBuildMismatchError`

::: varcode.GenomeBuildMismatchError

### `varcode.UnsupportedVariantError`

::: varcode.UnsupportedVariantError
