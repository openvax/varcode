# Reference DNA, species and ORF integration

## Scope

One Varcode PR addresses #610, #611, CLI reference/ORF access, and #593.
Keep ordinary annotations and existing input formats compatible.

1. Use the variant genome's reference DNA as the default fallback for ORF
   calculations, preserving bound transcript sequence and explicit providers.
2. Select translation codes from explicit source metadata and authoritative
   species/compartment assignments. Support the needed NCBI codes; reject
   unknown assignments instead of silently assuming vertebrate mitochondria.
3. Expose chromosome FASTA and ORF evidence selection in the CLI with early
   validation. Explain one reference setup path and sequence limitations.
4. Add versioned, lossless ORF effect/collection serialization, including
   failed bindings, alternative products, reference identity and evidence.
   Deserialization must not require a catalogue download or reannotation.
   Preserve ordinary JSON/CSV compatibility.

## Verification

- Regression tests for native and overridden DNA, both strands, missing DNA,
  splice/SV products, and collection/haplotype calls.
- Species/code tests against NCBI definitions, including human/mouse,
  invertebrate and yeast mitochondria, explicit overrides and unknown codes.
- CLI tests exercise actual parsing, annotation and exported ORF identity.
- ORF round trips cover mixed collections, failures, source ID collisions,
  phased outcomes and portable reference reconstruction.
- Run lint and the full suite, obtain green CI, merge and deploy from clean main.

## Dependency scope

Inspect PyEnsembl/ORFdata APIs and report missing capabilities upstream.
Implementation and release work remain in Varcode.

## Implementation and checks

- The ORF adapter now defaults to the existing genome sequence provider.
  Tiny GTF/RNA/DNA fixtures exercise native and overridden DNA on both strands,
  intronic edits, donor disruption, local SVs and known-cis groups.
- NCBI gc.prt 4.6 plus Taxonomy EUtils supply the packaged, hashed offline
  snapshot. Explicit models override species/compartment defaults. Unknown
  assignments and context-dependent stop/sense tables remain explicit errors.
  PyEnsembl's missing metadata contract is filed as openvax/pyensembl#481.
  Broader ORF binding-code coverage remains upstream in openvax/orfdata#57.
- Reuse the existing effect graph for ORF dataclasses and validated public
  ORFdata record constructors. Snapshot full RNAs so reference reconstruction
  does not consult a catalogue or annotation DB. Incomplete fusion partners
  retain only defined start/stop features.
- ORF-containing CSVs carry each effect graph beside the summary. Import
  rejects missing payloads and edited identifying columns.
- First full run: 3,512 passed; one documentation test stopped before the
  export block because the new setup snippet changed the fixed block count.
  Corrected the count; the documentation regression passes independently.
- Focused reference/ORF/CLI/species tests: 73 passed before the final scoped
  evidence/multiple-placement addition; that addition also passes independently.
- Source and wheel build successfully; wheel import verifies that the NCBI
  data is packaged. Lint and shared-environment `pip check` pass.
- Final full suite: **3,520 passed**, 64 existing warnings, 93% coverage.
  This includes exact-assembly species resolution for mirrored human/yeast
  references and explicit rejection of an unregistered mitochondrial reference.
