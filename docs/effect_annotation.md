# Effect annotation

Call `effects()` to predict the consequences of your variants. If you have not
installed reference data yet, start with [Getting started](getting_started.md).

<a id="basic-usage"></a>

## Annotate variants

```python
import varcode

variants = varcode.load_vcf("variants.vcf", genome=81)  # GRCh38
effects = variants.effects()
for variant, effect in effects.top_priority_effect_per_variant().items():
    print(variant.short_description, effect.short_description)
```

The collection holds one prediction per variant per overlapping transcript,
and a variant's effect can differ between transcripts. The loop above prints
one summary per variant, picking the most severe effect by Varcode's
consequence ordering, not by likelihood or clinical significance.
`effects.top_priority_effect()` instead selects one effect from the entire
collection. See [how to read results](getting_started.md#how-to-read-results).

The same interface handles structural variants. To include them when loading
a VCF, pass `parse_structural_variants=True`; see [SV loading](structural_variants.md#basic-usage).

<a id="read-an-effect"></a>

## Protein sequences

Given an effect from the collection:

```python
protein = effect.mutant_protein_sequence  # may be None
```

`None` means the sequence could not be determined, not that the protein is
unchanged. The prediction belongs to `effect.transcript`; intergenic effects
have no transcript. For cDNA, partial structures, and sequence evidence, see
[Transcript models](transcript_models.md).

To keep only effects that may change the protein, use
`effects.drop_silent_and_noncoding()`. It removes known silent and noncoding
predictions but keeps unresolved ones by default; see
[filtering by protein change](structural_variants.md#filtering-by-protein-change).

<a id="reading-alternatives"></a>

## Alternative outcomes

When Varcode can't decide between several possible consequences, it returns a
`MultiOutcomeEffect` holding every candidate instead of guessing. This happens
for splice variants, structural variants, and variants whose phase relative to
a nearby germline variant is unknown:

```python
from varcode import MultiOutcomeEffect

if isinstance(effect, MultiOutcomeEffect):
    for candidate in effect.candidates:
        print(candidate.effect.short_description)
        print(candidate.effect.mutant_protein_sequence)
        print(candidate.source, candidate.evidence)
```

Two shortcuts pick a single candidate, and they answer different questions:

- `most_likely_effect` returns the first candidate in the order the predictor
  listed them. That order reflects the predictor's preference; it is not a
  calibrated probability.
- `highest_priority_effect` returns the most severe candidate.

When you report a candidate's sequence, keep its `candidate.source` (which
predictor or evidence produced it) and `candidate.evidence` alongside it.

<a id="choose-a-deeper-topic"></a>

## Related guides

- <a id="splice-disrupting-variants"></a><a id="when-splice-disruption-is-in-play"></a><a id="splice-and-coding-effects-can-co-occur"></a><a id="the-spliceoutcomeset-shape"></a><a id="common-questions"></a><a id="rna-evidence-reconciliation"></a><a id="candidate-provenance"></a><a id="picking-a-single-candidate"></a><a id="limitations"></a>[Splice variants](splice_variants.md): affected signals, candidate mechanisms, and RNA evidence.
- <a id="structural-variants"></a>[Structural variants](structural_variants.md): deletions, duplications, inversions, and fusions.
- [Germline](germline.md) and [phasing](phasing.md): patient-specific baselines and linked variants.
- <a id="how-it-composes"></a><a id="the-four-primitives"></a>[Transcript models](transcript_models.md): cDNA, partial structures, and missing sequences.
- <a id="provenance"></a>[Saving results](csv.md#annotation-provenance): tables, provenance, and round-trip limits.
- <a id="annotator-selection"></a><a id="advanced-annotators-and-implementation-limits"></a><a id="what-the-optional-implementations-change"></a><a id="reading-sequences-and-alternatives"></a><a id="current-boundaries-and-legacy-names"></a>[Experimental annotators](experimental_annotators.md): optional implementations and their limits.
- <a id="writing-an-annotator"></a><a id="downstream-consumers"></a>[Writing an annotator](annotator_contract.md): custom implementations and unsupported inputs.
- [Effect types](effect_types.md) and [API reference](api.md): class definitions and parameters.
