# Annotators API

Ordinary annotation uses the default without selecting an implementation.
See [annotation choices](experimental_annotators.md) for alternatives and
[Writing an annotator](annotator_contract.md) for extensions.

## Annotators

::: varcode.EffectAnnotator

::: varcode.FastEffectAnnotator

::: varcode.ProteinDiffEffectAnnotator

## Registry

::: varcode.register_annotator

::: varcode.registered_annotators

::: varcode.get_annotator

::: varcode.get_default_annotator

::: varcode.resolve_annotator

::: varcode.set_default_annotator

::: varcode.use_annotator

<a id="experimental-transcript-model"></a>

## Shared transcript model

Select `annotator="transcript_model"` for nested `RealizedEffectCandidate`
outcome graphs, or `annotator="fast"` for the Varcode 10 compatibility
implementation.
See [supported inputs and results](experimental_annotators.md) before using it.

::: varcode.TranscriptModelEffectAnnotator

::: varcode.predict_transcript_model_effect

## Default DNA annotator

`annotator="dna"` (the default) combines a guarded interior-exon coding
fast-path with the full transcript/splice/phase model, returning standard
`MutationEffect`, `EffectCandidate`, and `MutantTranscript` objects across
point, splice, germline/phase, and structural variants.

::: varcode.DNAEffectAnnotator
