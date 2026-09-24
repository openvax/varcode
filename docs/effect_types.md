# Effect types

Use this page to look up a class after reading the
[effect annotation guide](effect_annotation.md). Definitions below are generated
from the source docstrings; this is reference material, not a required first read.

<a id="find-an-effect-family"></a>

## Effect families

| Question | Examples |
|---|---|
| [Coding changes](#coding-effects) | `Substitution`, `Insertion`, `Deletion`, `FrameShift`, `PrematureStop`, `Silent` |
| [Noncoding and untranslated regions](#noncoding-effects) | `Intronic`, `FivePrimeUTR`, `ThreePrimeUTR`, `NoncodingTranscript`, `Intergenic` |
| [Splice effects](#splice-effects) | `SpliceOutcomeSet`, with `NormalSplicing`, `ExonSkipping`, `IntronRetention`, or cryptic-site candidates |
| [Structural effects](#structural-effects) | `StructuralVariantEffect`, `GeneFusion`, `TranslocationToIntergenic` |
| [Linked variants and phase](#haplotype-effects) | `HaplotypeEffect`, `PhaseCandidateSet` |
| [Inherited allele overlap](germline.md#loss-of-heterozygosity-loh) | `GermlineAlleleOverlap`: no new allele sequence; LOH unassessed |
| [Unknown results and failures](#unknown-results) | `Unresolved`, `IncompleteTranscript`, `Failure` |

These labels describe predictions and their limits, not clinical significance.

Two distinctions are worth keeping in mind while reading:

- **Effects that carry a protein consequence vs. location-only
  effects.** Coding effects (`Substitution`, `FrameShift`,
  `PrematureStop`, …) and the splice *mechanism* effects
  (`NormalSplicing`, `ExonSkipping`, …) describe a change to the
  protein. The splice-signal *disruption* effects (`SpliceDonor`,
  `SpliceAcceptor`, `IntronicSpliceSite`, `ExonicSpliceSite`, all
  sharing the `SpliceSite` base) and the region effects (`Intronic`,
  `FivePrimeUTR`, …) describe *where* a variant landed and carry no
  protein consequence on their own.
- **Single effects vs. multi-outcome containers.** Most effects are a
  single answer. `MultiOutcomeEffect` subclasses (`SpliceOutcomeSet`,
  `ExonicSpliceSite`, the structural-variant effects, `HaplotypeEffect`,
  `PhaseCandidateSet`) bundle several candidate effects when the
  protein-level outcome isn't deterministic; each exposes
  `.candidates`, `.most_likely_effect`, and `.highest_priority_effect`.

Severity ordering across types is set by
[`effect_priority`](api_effects.md#varcode.effect_priority).

<a id="class-reference"></a>
<a id="varcode.effects.effect_classes"></a>

## Coding effects

::: varcode.effects.effect_classes.Silent

::: varcode.effects.effect_classes.AlternateStartCodon

::: varcode.effects.effect_classes.StartLoss

::: varcode.effects.effect_classes.Substitution

::: varcode.effects.effect_classes.ComplexSubstitution

::: varcode.effects.effect_classes.Insertion

::: varcode.effects.effect_classes.Deletion

::: varcode.effects.effect_classes.PrematureStop

::: varcode.effects.effect_classes.StopLoss

::: varcode.effects.effect_classes.FrameShift

::: varcode.effects.effect_classes.FrameShiftTruncation

::: varcode.effects.effect_classes.ExonLoss

## Noncoding effects

::: varcode.effects.effect_classes.Intergenic

::: varcode.effects.effect_classes.NoncodingTranscript

::: varcode.effects.effect_classes.FivePrimeUTR

::: varcode.effects.effect_classes.ThreePrimeUTR

::: varcode.effects.effect_classes.Intronic

## Splice effects

<a id="splice-outcome-container"></a>

See the [splice guide](splice_variants.md) for candidate access.

::: varcode.SpliceOutcomeSet

::: varcode.effects.effect_classes.SpliceSite

::: varcode.effects.effect_classes.IntronicSpliceSite

::: varcode.effects.effect_classes.SpliceDonor

::: varcode.effects.effect_classes.SpliceAcceptor

::: varcode.effects.effect_classes.ExonicSpliceSite

::: varcode.effects.effect_classes.SpliceMechanismEffect

::: varcode.effects.effect_classes.NormalSplicing

::: varcode.effects.effect_classes.ExonSkipping

::: varcode.effects.effect_classes.IntronRetention

::: varcode.effects.effect_classes.CrypticSpliceSiteEffect

::: varcode.effects.effect_classes.CrypticDonor

::: varcode.effects.effect_classes.CrypticAcceptor

::: varcode.effects.effect_classes.CrypticExonCandidate

## Structural effects

::: varcode.effects.effect_classes.StructuralVariantEffect

Legacy event-type classes remain importable; the default annotator emits
[transcript consequences](structural_variants.md#reading-sv-results).

::: varcode.effects.effect_classes.LargeDeletion

::: varcode.effects.effect_classes.LargeDuplication

::: varcode.effects.effect_classes.Inversion

::: varcode.effects.effect_classes.GeneFusion

::: varcode.effects.effect_classes.TranslocationToIntergenic

## Haplotype effects

::: varcode.effects.effect_classes.HaplotypeEffect

::: varcode.effects.effect_classes.PhaseCandidateSet

::: varcode.effects.effect_classes.GermlineAlleleOverlap

## Unknown results

::: varcode.effects.effect_classes.Unresolved

::: varcode.effects.effect_classes.IncompleteTranscript

::: varcode.effects.effect_classes.Failure

## Base classes

::: varcode.effects.effect_classes.MutationEffect

::: varcode.effects.effect_classes.MultiOutcomeEffect

::: varcode.effects.effect_classes.Intragenic

::: varcode.effects.effect_classes.TranscriptMutationEffect

::: varcode.effects.effect_classes.Exonic

::: varcode.effects.effect_classes.CodingMutation

::: varcode.effects.effect_classes.NonsilentCodingMutation

::: varcode.effects.effect_classes.KnownAminoAcidChange
