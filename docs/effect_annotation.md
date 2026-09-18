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

The collection contains predictions for the relevant transcripts of each variant.
A variant can have different effects on different transcripts; the loop above
prints one summary per variant. `effects.top_priority_effect()` instead selects
one effect from the entire collection.

Priority is Varcode's consequence ordering, not a probability or a clinical
classification. Keep the full collection when you need all transcript predictions.

The same interface handles structural variants. To include them when loading
a VCF, pass `parse_structural_variants=True`; see [SV loading](structural_variants.md#basic-usage).

<a id="read-an-effect"></a>

## Protein sequences

Given an effect from the collection:

```python
protein = effect.mutant_protein_sequence  # may be None
```

You do not need a different annotator to get a protein sequence. `None` means
the sequence is unavailable, not that the protein is unchanged. The prediction
belongs to `effect.transcript` when a transcript applies; intergenic effects
have none. Predicted sequence is not evidence of expression.

For cDNA, partial structures, and sequence evidence, see [Transcript models](transcript_models.md).

<a id="reading-alternatives"></a>

## Alternative outcomes

Splice, structural, and phase-dependent effects may contain several candidates:

```python
from varcode import MultiOutcomeEffect

if isinstance(effect, MultiOutcomeEffect):
    for candidate in effect.candidates:
        print(candidate.effect.short_description)
        print(candidate.effect.mutant_protein_sequence)
        print(candidate.source, candidate.evidence)
```

`most_likely_effect` returns the producer's first candidate;
`highest_priority_effect` selects by consequence severity. The first candidate
is not necessarily the most disruptive, and its position is not a calibrated
probability. Keep `candidate.source` and `candidate.evidence` with any sequence
you report.

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
