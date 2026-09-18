# Transcript models

An effect may carry a `MutantTranscript` describing edits, retained reference
segments, or an observed RNA structure. It can contain cDNA and protein sequence,
but it may also be partial or have no resolved sequence.

You usually want the effect's `mutant_protein_sequence`; a transcript model is
useful when you also need the RNA structure or its provenance. Not every effect
has a model, even when a predicted protein is available.

## Sequence access

Given an annotated `effect`:

```python
protein = effect.mutant_protein_sequence
model = effect.mutant_transcript
cdna = model.cdna_sequence if model is not None else None
evidence = model.evidence if model is not None else None
```

`None` means unavailable. It does not establish an unchanged protein, an absent
transcript, or a harmless variant.

## Selenocysteine

Ensembl marks selenocysteine (Sec) as `U` in the reference protein, encoded by
an in-frame UGA. UGA is decoded as Sec only with a SECIS element in the same
mRNA's 3′ UTR, whose position isn't annotated. Predicted proteins read an
annotated Sec codon as `U` when the model still maps it (no edit touches it)
and keeps some of a selenoprotein's 3′ UTR. Where no selenoprotein 3′ UTR
remains, for example a fusion downstream of Sec, UGA is read as a stop.
A sequence without reference coordinates, such as an imported RNA assembly,
is translated literally. A partly kept 3′ UTR still gets the Sec reading; the
[SV change flags](structural_variants.md#filtering-by-protein-change) report
such cases as unresolved.

## Partial structures

A model with reference segments is not necessarily a complete mutant RNA.
For example, a breakend can retain only a transcript's local 5′ prefix or 3′
suffix. That model has `evidence["sequence_status"] = "retained_reference_fragment"`,
while full cDNA and protein remain unknown. Concatenating its segments does
not reconstruct the missing partner or establish the full allele.

Likewise, a sequence supplied by an RNA assembler may be only a junction
fragment. Keep its completeness and source information with the sequence.
The [RNA import guide](rna_structures.md)
explains how observed structures and sequence-predicted proteins are represented;
neither alone proves translation.

## Related objects

| Object | Contents |
|---|---|
| `MutationEffect` | A predicted consequence or an unresolved result |
| `MutantTranscript` | Edits or transcript structure, provenance, and optional sequences |
| `MultiOutcomeEffect` | Several candidate effects, each with its own consequence and evidence |

The [API reference](api_rna.md#mutant-transcripts) lists model fields and construction
helpers. The experimental transcript-model annotator uses additional baseline
and mutant product wrappers; see [its result format](experimental_annotators.md#transcript-model-results).
