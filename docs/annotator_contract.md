# Unified default annotator

## Scope

Keep `fast` as the compatible registry name for the built-in default. It owns
point-variant prediction, structural-variant routing, and the existing
germline-aware point-variant path. Both explicit and implicit selection must
use that same routing. Keep the structural implementation in its own module.
Do not promote the experimental transcript model or change its ranking.

Structural prediction lives in `varcode.effects.structural` as plain internal
helpers. The default and transcript model's fusion path call those helpers
directly, not another annotator. The registry selects an explicitly requested
implementation, never an implementation based on variant kind. There is no
separate structural annotator or implicit fallback between annotators.

## Partial annotators

The protocol requires only `name` and `annotate_on_transcript`. An annotator
returns a `MutationEffect`, or Python's `NotImplemented` singleton when it
cannot annotate this particular input. No capability list is required.
An optional `annotate_with_context` method accepts patient germline context;
without that method a nonempty context is unsupported, not silently ignored.

At public prediction boundaries, convert `NotImplemented` to
`Unresolved(mechanism="unsupported_annotation", reason=...)`. Preserve the
selected annotator's provenance. Never substitute another annotator after an
experimental annotator declines. Actual exceptions retain the existing
`raise_on_error` behavior; `None` is not a valid annotation result.

`protein_diff` and `transcript_model` remain optional experimental implementations.
Varcode 9 removes the old structural annotator class/module/registry entry and
`UnsupportedVariantError`. Structural callers use the default; partial plugins
use the return-value contract above. Structural mutant transcripts retain their
historical builder provenance (`structural_variant`); effect collections record
the selected annotator, including `fast` or `transcript_model`. A stored builder
label is not a selectable annotator name.

Joint haplotype construction still occurs in `VariantCollection.effects()`
outside the selected annotator. Moving that active behavior into the annotator
is tracked separately in [#437](https://github.com/openvax/varcode/issues/437),
not part of the structural-helper cleanup.

## Verification

Test default and explicit routing, direct transcript calls, partial plugins
without capability lists, context-dependent refusal, unresolved serialization
and provenance, and genuine errors. Keep protein-diff parity tests restricted
to the implementations' shared domain; unsupported SV inputs must be tested
as explicit unknowns. Run lint, the full suite, GitHub CI, then merge and deploy.
