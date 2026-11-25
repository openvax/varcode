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

This is the default annotator; `annotator="fast"` selects the compatibility
implementation.
See [supported inputs and results](experimental_annotators.md) before using it.

::: varcode.TranscriptModelEffectAnnotator

::: varcode.predict_transcript_model_effect

## Unified DNA annotator (opt-in)

`annotator="unified"` combines a safe interior-exon SNV fast-path with flat
`DNAHypothesis` candidates and deterministic `DNASearchAnchor` search geometry
for large insertions and structural rearrangements.

::: varcode.UnifiedDNAEffectAnnotator

::: varcode.DNAHypothesis

::: varcode.GenomicSpliceJunction

::: varcode.DNASearchAnchor

::: varcode.build_dna_search_anchor
