# General ORF reference and evidence layer

Status: architectural correction to the unmerged oncoref/Varcode/Isovar PRs.
Working distribution name: orfdata. No release is authorized by this document.

## Ownership

PyEnsembl remains responsible for reference genomes, annotations, transcripts,
exons and canonical proteins. A small separate library adds source-defined ORFs
on those transcripts and their experimental evidence. oncoref remains a consumer
for cancer-specific analysis. Neither Varcode nor Isovar owns reference datasets.
PepData's immune-epitope/property scope is not a fit.

The new library has a dependency-light data layer and a PyEnsembl-compatible
binding API. It must not import Varcode, Isovar or oncoref. No new PyEnsembl core
API is needed to bind a Transcript; avoid reverse dependencies or monkey-patching.

## Public objects

- ORFModel: immutable source ORF identity, assembly, species (NCBI taxon ID and
  scientific name), translation table, parent transcript, spliced coding blocks,
  protein sequence, ORF class. ORF class is distinct from transcript biotype.
- ORFEvidence: ribosome observations and peptide-level MS observations with their
  source, accession, mapping ambiguity, assay and quality. Native evidence does
  not become evidence for a modified product; orthology is not direct evidence.
- ORFRecord: model plus evidence, retaining legacy flat accessors for Isovar.
- ReferenceORF: validated binding of a record to a PyEnsembl Transcript, retaining
  the original transcript (including noncoding biotype), cDNA coding interval,
  coding sequence and actual resolved annotation/version. Multiple ORFs per RNA
  are normal. Annotation lookup never relabels the transcript protein-coding.
- ORFCatalog: explicitly one species/assembly per catalogue. Query by transcript,
  gene, locus, class or evidence. for_transcript(transcript) validates and returns
  ReferenceORFs. Load custom versioned catalogues offline. The initial bundled
  human data is identified as human-only; no implied mouse evidence.
- Varcode ORFAnnotation: reference ORF plus effect and realized products. Exports
  separately identify reference model/evidence and modifications, preserving
  phase/splice alternatives. A private transcript-engine adapter permits reuse
  of existing variant logic without exposing a fake transcript as the base ORF.
- Isovar: attach those same annotations to RNA products with candidate-specific
  RNA support, maintaining conditional translation and unresolved RNA polarity.

## Generality and compatibility

Support uORF, overlapping uORF, internal/out-of-frame, downstream and noncoding-
transcript ORFs with the same model and evidence schema. Preserve source classes
verbatim rather than equating them with Ensembl transcript biotypes. Complete
stopped ORFs are the initial binding contract; partial models are rejected
explicitly rather than given invented boundaries.

A generic v1 schema requires species, assembly, translation table and source
provenance. Legacy Isovar human data may be read only through its explicit known
human/GRCh38 provenance. Keep legacy Isovar uorf_* imports/call shapes working.
New generic APIs use orf terminology; defaults remain explicit about dataset scope.

Multi-species path: one versioned catalogue per organism/assembly, importer input
using the same model/evidence schema, optional organism-specific source adapters
(e.g. OpenProt). Test a mouse noncoding transcript with real PyEnsembl objects,
plus human coding transcripts, different frames, minus strand, splice junctions,
assembly/species mismatches and unknown translation tables. Synthetic examples
are identified as test fixtures, never scientific observations. OpenProt already
publishes multiple species with Ensembl/RefSeq mappings and evidence filters.

## Review corrections

Remove the oncoref runtime dependency from both consumers and withdraw its ORF
PR when the replacement is reviewable. Keep unrelated lint fixes independently.
Review evidence retention for duplicate reference codons (retain all copies),
explicit translation tables and noncoding parents; preserve all Sid alternatives.
Document unsupported de novo initiation discovery and partial ORFs honestly.

## Verification

New library: existing evidence/context/builder tests, generic data round-trip,
PyEnsembl binding and multi-species isolation, independently installed wheel.
Varcode/Isovar: focused tests then required full lint/test scripts. Regenerate
Sid outputs from the original packaged reads. Existing human evidence counts and
native peptide attribution must remain unchanged. Do not claim release until CI,
merge and PyPI publication actually complete; publication remains unapproved.

Primary sources:
- https://www.nature.com/articles/s41587-022-01369-0
- https://www.nature.com/articles/s41586-026-10459-x
- https://openprot.org/p/help
- https://github.com/openvax/pyensembl
