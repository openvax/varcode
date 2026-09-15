# Unified default annotator

## Scope

Keep `fast` as the compatible registry name for the built-in default. It owns
point-variant prediction, structural-variant routing, and the existing
germline-aware point-variant path. Both explicit and implicit selection must
use that same routing. Keep the structural implementation in its own module.
Do not promote the experimental realized-effect model or change its ranking.

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

`protein_diff` and `realized` remain optional experimental implementations.
The structural-only entry point remains available for compatibility. Remove
the unused `supports` metadata; retain the exported `UnsupportedVariantError`
as a compatibility import, not the new contract. Use a major version bump
because capability metadata and scoped SV-selection behavior change.

## Verification

Test default and explicit routing, direct transcript calls, partial plugins
without capability lists, context-dependent refusal, unresolved serialization
and provenance, and genuine errors. Keep protein-diff parity tests restricted
to the implementations' shared domain; unsupported SV inputs must be tested
as explicit unknowns. Run lint, the full suite, GitHub CI, then merge and deploy.
