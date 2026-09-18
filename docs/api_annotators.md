# Annotators API

Ordinary annotation uses the default without selecting an implementation.
See [experimental annotators](experimental_annotators.md) for alternatives and
[Writing an annotator](annotator_contract.md) for extensions.

## Annotators

::: varcode.EffectAnnotator

::: varcode.FastEffectAnnotator

::: varcode.ProteinDiffEffectAnnotator

## Registry

::: varcode.register_annotator

::: varcode.get_annotator

::: varcode.get_default_annotator

::: varcode.set_default_annotator

::: varcode.use_annotator

## Experimental transcript model

This annotator is opt-in with `annotator="transcript_model"`; `fast` remains
the default. The former `realized` name and public imports remain aliases.
See [supported inputs and results](experimental_annotators.md) before using it.

::: varcode.TranscriptModelEffectAnnotator

::: varcode.predict_transcript_model_effect
