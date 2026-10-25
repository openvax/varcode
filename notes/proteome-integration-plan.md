# Varcode Proteome integration (#620)

## Scope and dependencies

Expose existing Varcode results through ORFdata's `Proteome`, without changing
annotation or translation. #195 is independent: partial and unresolved results
must be represented faithfully by this adapter today. The user is handling
ORFdata #124 (source identity independent of exporter) and #125 (unknown source
taxon) separately. Do not modify that checkout or implement substitutes here.
Release against the published shared contract once those fixes are available.
ORFdata #126 additionally tracks available producer metadata omitted from its
shared FASTA header. Keep explicit contract tests for #124/#125/#126; do not
merge this PR until all three pass against a published release. The fresh
environment also reproduced Varcode's existing #628 mypy failure; include the
small typing correction, with no runtime or dependency downgrade.

## Public workflow

`varcode.to_proteome(effects, sample_id=None, run_id=None)` accepts an
`EffectCollection` or an iterable of native effects/candidates. An explicit
`orf_annotations=` input retains standalone ORF annotations, including failed
bindings that cannot appear as ordinary effects. Callers then use the shared
`merge`, `by_sequence`, `write_fasta`, and `save`/`load` methods. ORFdata imports
stay at the optional adapter boundary; no new protein hierarchy is introduced.

## Mapping

- Retain each supplied event, candidate context, and every classified outcome.
  A representative effect is never a substitute for a realized candidate's
  full `outcomes` collection. Traverse nested multi-outcome effects without
  losing RNA evidence candidates or creating cycles.
- Emit supplied mutant amino-acid sequences as source-specific external
  entries. Keep baseline products and comparisons in versioned producer
  records, even when the baseline is absent from the FASTA collection.
- Retain focal/background variant identities, phase and splice assumptions,
  source models, ordered sequence segments, affected intervals with coordinate
  conventions, evidence, probabilities, phase mass, and ordinal rank. Native
  identifiers and a versioned Varcode record type link entries to their source.
- Distinguish unresolved, annotation failure, failed binding, hypothesis limit,
  known no-product, unchanged, and available/partial product states. Preserve
  unknown boundary states and omitted-sequence reasons. Do not fabricate an
  origin, taxon, gene, generating tool, or complete sequence.
- Source identity is independent of iteration order and sequence grouping.
  Use readable native/semantic labels with deterministic disambiguation when
  necessary. Shared entry/record identity validates conflicting content;
  exporter provenance must not split one source record (#124).
- Reuse the shared readable FASTA token/manifest contract. Default export
  retains distinct sources; explicit exact-sequence deduplication preserves
  every contributing association. Record the supplied sequence's terminal-stop
  convention and compute length consistently with that sequence.
- A namespaced, versioned JSON producer graph must survive offline load without
  fetching genomes/catalogues or invoking arbitrary Python constructors.
  Unknown or unsupported metadata must be retained or explicitly rejected.
- Preserve known reference/taxon/strain identities and validate composition
  through ORFdata; unknown source taxon stays unknown (#125).

## Validation and delivery

Use native synthetic results with independently specified expected sequences
and relationships. Cover coding and ncORF changes, alternate splice outcomes,
locally merged candidates with distinct full proteins, identical mutant
sequences with different baselines, empty deletion intervals, fusion/multigene
products, unknown labels/taxa, no-product/failed/capped states, and source-ID
escaping/conflicts. Round-trip FASTA plus the manifest offline, merge with a
reference and an Isovar-shaped producer graph, and verify optional imports.

Add API documentation and a composition example. Bump a minor version for the
new public adapter. Run Ruff/mypy, the full test suite and documentation build;
review the PR; require green CI with the released ORFdata contract before merge.
Deploy from clean main, verify PyPI artifacts/tag, and update the shared Varcode
wheel while preserving editable checkouts and unrelated packages.
