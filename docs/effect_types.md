# Effect types

Every effect class varcode can attach to a `(variant, transcript)`
pair, auto-generated from the source docstrings in
`varcode.effects.effect_classes` — so this list never drifts from the
code.

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

For a grouped quick-reference, see the
[Effect Types table in the README](https://github.com/openvax/varcode#effect-types).
Severity ordering across types is set by
[`effect_priority`](api.md#varcode.effect_priority).

::: varcode.effects.effect_classes
