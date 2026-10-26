# Transcript selection (#195)

## Current behavior

`effects()` already includes overlapping noncoding and incomplete transcripts.
The missing selection API forces callers to annotate transcripts they later
discard. Ordinary protein prediction still uses the `protein_coding` biotype;
selection and eligibility for protein prediction are distinct decisions.

## Proposed contract

- Add an optional `transcript_filter(transcript)` predicate to
  `Variant.effects()`, `VariantCollection.effects()` and
  `predict_variant_effects()`. Omission preserves current results.
- Evaluate selection before invoking the annotator. Keep genuine intergenic
  and intragenic results; rejecting an overlapping transcript must never turn
  it into an intergenic/intragenic result. An entirely rejected transcript
  set may produce an empty collection.
- Preserve selection through germline, phase and RNA processing. Additional
  ORF placements are selected using their original reference transcript;
  unresolved records without a transcript remain visible as failures.
- A non-callable filter is a configuration error, including for empty inputs.
  Predicate exceptions follow `raise_on_error`. Document the distinction between selecting a
  transcript and supplying an additional coding model through ORFdata.
- Do not silently change the coding biotype or fabricate complete coding
  boundaries. This PR implements the selection request from #195; broader
  coding-biotype eligibility remains separate follow-up work. The user was
  offered both scopes; absent a reply, selection first is the stated default.
- Select primary annotation targets; fusion partners discovered inside an
  annotator remain part of the selected target's prediction. Do not hide
  mechanism alternatives by post-filtering an effect graph.

## Validation

Use pinned references and small local GTF fixtures to cover coding, noncoding,
pseudogene and incomplete transcript visibility; selection by ID and biotype;
skipped annotator calls; false intergenic/intragenic prevention; collections,
structural inputs, additional ORFs and phased joint effects. Run focused tests,
lint, the full suite and docs, then open the versioned PR and verify CI.

## Primary references

- [Ensembl biotypes](https://www.ensembl.org/info/genome/genebuild/biotypes.html)
- [GENCODE biotypes](https://www.gencodegenes.org/pages/biotypes.html)

These distinguish transcript classification from an available coding model;
NMD, pseudogene and IG/TR labels alone do not establish a complete translated
product or its expression in a sample.

## Upstream finding

[PyEnsembl #499](https://github.com/openvax/pyensembl/issues/499) records an
absent-feature-table error found while constructing the incomplete-transcript
fixture. This PR tests an explicitly truncated start codon and makes no
PyEnsembl changes.
