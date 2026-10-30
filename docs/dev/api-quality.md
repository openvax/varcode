# Maintaining the public API and documentation

The common traversal helpers (`candidate_children`, `iter_effects`, and
`iter_effect_paths`) cover both candidate shapes with typed effect access.
Collection iteration, indexing, filtering, cloning and generic grouping are
typed; slicing preserves Sercol's underlying-sequence return rather than
silently creating a new collection. Subclass-preserving clones/groups/filters
retain the existing metadata behavior. Broader model protocols, specialized
grouping/set operations, and subclass field annotations remain gradual (#653).

The public reference and type annotations should describe the same behavior
users get from a normal installation. Preserve established callers while
improving clarity; naming consistency alone is not a reason to change a return
shape, coordinate convention, or uncertainty policy.

## API quality audit: October 2026

The audit at 10.22.0 inspected package-root exports, core annotation and
collection methods, extension contracts, generated reference coverage, error
guidance, and installation instructions. It was a public-surface audit, not a
line-by-line review of every internal algorithm.

The baseline had 114 declared root exports: 111 classes/functions, including
51 functions. Only seven of those functions had complete parameter and return
annotations. Five additional intended imports were omitted from `__all__`.
Thirty-six declared classes/functions lacked a direct generated reference entry;
some appeared in task guides or were aliases, so this was a discoverability
measure rather than a count of wholly undocumented features.

The 26 public methods/functions in `variant.py`, 29 in
`variant_collection.py`, and 29 in `effects/effect_collection.py` had no return
annotations. Among the 78 public methods in `effect_classes.py`, only one had a
return annotation. These counts exclude underscore-prefixed methods and include
property getters; they are a snapshot, not a coverage target.

### Addressed in the first pass

- Explicit root exports for the chromosome resolver and InvalidContigError,
  enforced by consumer checks that disallow implicit Varcode re-exports.
- Return contracts for ordinary annotation, filtering, and ranking, plus typed
  annotator registration/selection and explicit unsupported return values.
- Nullable base-effect sequence/offset fields and three-state change flags.
  Previously, type checkers inferred None-only sequence and bool-only change
  fields from base-class defaults, overlooking valid subclass results.
- Reference entries for previously omitted root APIs, with separate ORF and
  sample-identity pages and explanations of compatibility return shapes.
- Docstrings for legacy uORF wrappers, corrected empty-collection ranking
  behavior, and an explicit process-global override contract.
- Current reference-mismatch diagnostics and ORF extra dependency requirements.

See [#651](https://github.com/openvax/varcode/issues/651) for implementation
and validation. Simultaneous chromosome-resolver fixes are tracked in
[#650](https://github.com/openvax/varcode/issues/650).

### Follow-up priorities

Collection traversal and remaining model contracts are tracked in
[#653](https://github.com/openvax/varcode/issues/653). Present transcript/gene
adapter objects remain explicitly gradual until their shared protocols are
defined; their base-class None defaults do not mean they are always absent.

| Priority | Work | Acceptance evidence |
|---|---|---|
| First | Type collection iteration, indexing, grouping, and clone operations without losing subclass/metadata behavior | A typed file-to-effects loop retains concrete element types; slicing and grouped results preserve existing behavior |
| First | Type effect sequence/coordinate properties and nullable/tri-state results | Consumer checks distinguish unavailable sequence, false change flags, and unresolved alternatives |
| Next | Type optional ORF, phase, RNA and germline inputs and complete their parameter/return docstrings | Supported duck-typed integrations remain valid; invalid shapes are rejected without requiring unnecessary optional imports |
| Next | Reconcile result-model availability and composition | Address [cDNA access #179](https://github.com/openvax/varcode/issues/179) and [uncertainty composition #421](https://github.com/openvax/varcode/issues/421) with behavior tests, not aliases that hide differences |
| Next | Complete documentation style and integration-guide upkeep | Follow [#570](https://github.com/openvax/varcode/issues/570), check examples/links and rendered pages, and date downstream compatibility claims |

Dependency work in [#640](https://github.com/openvax/varcode/issues/640) remains
separate. Correcting an installation paragraph does not consolidate optional
dependencies. This pass does not make every method statically checked.

## Adding or changing an API

1. Define input/return shapes, coordinate conventions, unavailable states,
   exceptions, and resource ownership before choosing names.
2. Use an existing entry point when it already expresses the operation. Keep
   compatibility wrappers where callers rely on their selection or export
   behavior, and document the canonical alternative.
3. Add intended root imports to `__all__` and the relevant reference page.
   Mark experimental and compatibility APIs beside their entries.
4. Add useful parameter and return annotations. Use protocols for actual
   extension contracts and type-only imports to avoid runtime dependency
   cycles. Do not replace unknown contracts with Any to claim full coverage.
5. Add static consumer checks that exercise a real workflow, assert inferred
   types, and reject incompatible inputs. `warn_unused_ignores` ensures
   negative checks fail if an API silently becomes untyped.

## Writing docstrings and guides

Use NumPy-style Parameters, Returns/Yields, Raises, and Notes sections where
they add information. Explain units, interval endpoints, None/empty behavior,
and material exceptions. A type name alone rarely explains the contract.

Start task guides with a small runnable example and interpret its output.
Keep complete signatures in generated reference pages. Put uncertainty beside
the result it qualifies. Use inline code for literal symbols, values, and
commands; use ordinary prose for concepts and library names.

Check version claims against package metadata and executable behavior. Prefer
links to maintained dependency documentation over copying changing compatibility
tables. Preserve existing page URLs and useful anchors when reorganizing.

## Validation

Install `.[lint]` for type/lint checks and MkDocs, Material for MkDocs, and
`mkdocstrings[python]` for documentation builds. Run:

```bash
./lint.sh
./test.sh
mkdocs build --strict
```

Keep user-facing example checks in `tests/test_documentation_examples.py`,
which reads the actual Markdown. Review generated signatures and navigation in
the built site as well as checking the build exit status. Follow AGENTS.md for
releases: version and release notes in the PR, green CI, then deployment from
clean synchronized main and verification of the published version.
