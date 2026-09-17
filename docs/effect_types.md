# Effect types

Use this page to look up a class after reading the
[effect annotation guide](effect_annotation.md). Definitions below are generated
from the source docstrings; this is reference material, not a required first read.

## Find an effect family

| Question | Examples |
|---|---|
| Did the coding sequence change? | `Substitution`, `Insertion`, `Deletion`, `FrameShift`, `PrematureStop`, `Silent` |
| Is the location noncoding or untranslated? | `Intronic`, `FivePrimeUTR`, `ThreePrimeUTR`, `NoncodingTranscript`, `Intergenic` |
| Could splicing change? | `SpliceOutcomeSet`, with `NormalSplicing`, `ExonSkipping`, `IntronRetention`, or cryptic-site candidates |
| Is this a structural event? | `LargeDeletion`, `LargeDuplication`, `Inversion`, `GeneFusion`, `TranslocationToIntergenic` |
| Do linked variants or unknown phase matter? | `HaplotypeEffect`, `PhaseCandidateSet` |
| Is the result unknown or failed? | `Unresolved`, `IncompleteTranscript`, `Failure` |

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
[`effect_priority`](api.md#varcode.effect_priority).

## Class reference

::: varcode.effects.effect_classes

## Splice outcome container

`SpliceOutcomeSet` lives in a separate module. See the
[splice guide](effect_annotation.md#splice-disrupting-variants) for usage.

::: varcode.SpliceOutcomeSet
