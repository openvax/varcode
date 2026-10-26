# Proteome identity and provenance corrections (#631)

Reproduce the three review cases on merged Varcode 10.19.0 with the published
ORFdata 0.36 API before editing implementation. Keep changes in Varcode.

1. Include non-sequence mechanism coordinates and model evidence in event
   identity, so standalone alternatives retain distinct source records. Keep
   predicted products out of event identity so conflicting content still raises.
2. Retain single-hypothesis germline phase state and in-window variants in the
   native archive and source identity. Preserve the same evidence through native
   effect serialization and proteome save/load without changing annotation math.
3. Deduplicate transcript aliases by stable transcript/reference identity, not
   Python object identity. Preserve distinct reference scopes and fusion roles;
   original and native-JSON-restored predictions must compose idempotently.

Test direct and parent exports, known-trans/cis versus absent phase context,
reference-scoped aliases, actual Ensembl fusion round trips, stable repeated
export, and real conflicting-output rejection. Run lint, full tests, docs and
CI, then merge and deploy a patch version from clean synchronized main.

Self-review refinement: source-model identity must also exclude sequence bytes
embedded in assembled ReferenceSegment providers. Project only declared model
evidence, reference identities and coordinates; keep complete sequence content
in the conflict-checked archive. Add a same-model/different-assembled-sequence
regression before the final full-suite run.
