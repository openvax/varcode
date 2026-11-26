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

`annotator="dna"` (the default) combines a safe interior-exon coding fast-path
with the full transcript/splice/phase model, flat `DNAHypothesis` candidates,
ascending genomic intron coordinates (`GenomicSpliceJunction`), and
deterministic `DNASearchAnchor` search geometry for large insertions and
structural rearrangements. `UnifiedDNAEffectAnnotator` and `annotator="unified"`
remain as deprecated aliases for 12.3.x callers.

::: varcode.DNAEffectAnnotator

::: varcode.UnifiedDNAEffectAnnotator

::: varcode.DNAHypothesis

::: varcode.GenomicSpliceJunction

::: varcode.DNASearchAnchor

::: varcode.build_dna_search_anchor
