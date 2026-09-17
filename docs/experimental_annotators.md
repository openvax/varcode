# Experimental annotators

Varcode includes two opt-in implementations for comparing predictions:
`protein_diff` and `transcript_model`. They do not support every input handled
by the default. For ordinary annotation, use [effects()](effect_annotation.md)
without selecting an implementation.

## Supported inputs

| Selection | Method | Inputs |
|---|---|---|
| Default (omit `annotator=`) | Predicts small-edit and structural consequences; compares small edits against a patient baseline when supplied | Small variants and SVs. Germline/phase context for small edits, but not general SV-plus-germline combinations. |
| `annotator="protein_diff"` | Edits and translates transcripts, then compares proteins; shares the default's splice, location, and germline helpers | Small variants. SVs are unsupported. |
| `annotator="transcript_model"` | Builds transcripts for phase/splice hypotheses, compares each with the patient baseline, and merges equivalent results | Small variants and local DEL/DUP/INV. Insertions need `alt_assembly`; CNVs are unsupported. |

The transcript-model experiment delegates BNDs to the structural helper; it does
not resolve BND-plus-haplotype combinations. SVs extending beyond the selected
transcript also have a known classification problem; see [limitations](#known-limitations).

The default's registry name is `fast`; it appears in provenance but need not be
passed explicitly. It already predicts protein sequences and handles structural
variants. `protein_diff` is an alternative implementation, not a protein-output
option.

## Comparing annotators

Given a `variant` and one of its `transcript` objects, selection is explicit:

```python
effect = variant.effect_on_transcript(transcript)  # ordinary use
comparison = variant.effect_on_transcript(transcript, annotator="protein_diff")
experimental = variant.effect_on_transcript(transcript, annotator="transcript_model")
```

The transcript model can use `germline=` and `phase_resolver=` context. Canonical
and exon-skip paths need only transcript annotation; a genome with reference
FASTA additionally supplies sequence for intron retention and cryptic splice
sites. Without a calibrated scorer, mechanism preference is an ordering rule,
not a probability. Selecting this experiment does not guarantee that every
structural or combined input is supported.

## Transcript-model results

Use the [sequence accessors](transcript_models.md#sequence-access) for a single
effect, including the same checks for missing sequences and partial structures.

The ordinary and experimental candidate wrappers are not interchangeable:

| Candidate type | Shared access | Additional information |
|---|---|---|
| `EffectCandidate` (ordinary splice/SV/RNA outcome sets) | `candidate.effect` | `source`, `evidence`; sequences are on the effect or its optional `mutant_transcript` |
| `RealizedEffectCandidate` (transcript-model hypothesis pipeline) | `candidate.effect` | `outcomes`, `hypotheses`, `probability`, `ordinal_key`; each outcome has `baseline` and `mutant` products with `cdna_sequence`, `protein_sequence` and `evidence` |

The experiment returns an ordinary top effect with candidates attached, not
necessarily a `MultiOutcomeEffect`. Its BND delegation uses ordinary candidates;
early noncoding, incomplete or unresolved results may have no candidates at all.
Inspect the returned shape rather than assuming it from the annotator name:

```python
from varcode import EffectCandidate
from varcode.effect_hypotheses import RealizedEffectCandidate

print(experimental.short_description)
for candidate in getattr(experimental, "candidates", ()):
    print(candidate.effect.short_description)
    if isinstance(candidate, RealizedEffectCandidate):
        # A merged candidate can retain several hypotheses and their products.
        for outcome in candidate.outcomes:
            print(outcome.baseline.protein_sequence, outcome.mutant.protein_sequence)
            print(outcome.mutant.cdna_sequence, outcome.hypothesis.evidence)
    elif isinstance(candidate, EffectCandidate):
        print(candidate.source, candidate.evidence)
        print(candidate.effect.mutant_protein_sequence)
```

## Known limitations

- Combined variants: `VariantCollection.effects()` still builds joint haplotype
  effects outside the selected annotator when given a `phase_resolver`
  ([#437](https://github.com/openvax/varcode/issues/437)). The separate
  `predict_transcript_model_effect(variants, transcript, ...)` function accepts
  known-cis variants, but selecting the experiment does not replace that
  collection-level step.
- Cross-boundary SVs: the experiment can clip a DUP spanning CFTR's coding
  sequence to a local model and misclassify it as `FivePrimeUTR`
  ([#449](https://github.com/openvax/varcode/issues/449)). This is a bug, not
  evidence that the event is harmless.
- SV filtering: `modifies_protein_sequence=False` does not reliably mean an
  unchanged SV protein, and `drop_silent_and_noncoding()` can discard unresolved
  SVs ([#418](https://github.com/openvax/varcode/issues/418)). Retain them for review.

## Previous names

`realized`, `RealizedEffectAnnotator`, and `predict_realized_effect` are aliases
for the transcript-model implementation. Use the `transcript_model` names in
new code. The `RealizedEffectCandidate` wrapper still has a different interface
from ordinary `EffectCandidate`, as shown above.

To implement your own, see [Writing an annotator](annotator_contract.md).
