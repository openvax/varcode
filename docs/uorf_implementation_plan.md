# Shared uORF evidence and variant annotation

The released Isovar catalogue contains 3,771 upstream ORF models with model-specific ribosome and MS observations. Varcode currently annotates canonical UTRs only. Implement the existing Varcode #583 pathway with additive annotations.

Ownership: oncoref owns the unchanged checksum-pinned empirical catalogue, records, source manifests/licenses and reproducible builder. Varcode owns sequence-context and variant consequence calculations using an explicit validated alternative transcript frame and existing transcript-model realization/classification. Isovar retains compatible catalogue/context imports and its RNA reconstruction duties. No Varcode -> Isovar dependency. Varcode can accept supplied model records without oncoref; an optional uorf extra loads the shared catalogue.

Deliverables:
- Move data/APIs to oncoref, keep legacy data schema and provenance byte-for-byte, test offline packaged access and source builder. Preserve Isovar compatibility and CLI.
- Expose per-model annotations for variants and known-cis variant groups: validate assembly, entire spliced model path and reference protein; preserve canonical effects; include contexts, exact ORF identity, coding consequences, reference vs mutant sequences and native peptide retention, missingness and phase/splice alternatives. Use full reference transcript sequence for stop loss and frameshifts.
- Annotate mapped realized RNA products through an explicit reference start/frame. Exercise TPST1-CRCP with the existing checksum-pinned Sid RNA fixture and retained native/altered peptide intervals. Unresolved starts/sequence stay explicit; never promote reference MS into mutant evidence.
- Document usage and bounded discovery limits. Regression tests: minus strand, split codon, REF mismatch, dual +4 roles, start/stop loss, indels, phase/splice alternatives, Sid altered uORF and lossless JSON.
- Version/release notes in every PR; lint/full suites/CI green, merge then deploy each dependency in order; verify live PyPI and relevant issue triage.

Primary sources read: Deutsch et al. 2026 (doi:10.1038/s41586-026-10459-x); catalogue-pinned Mudge et al. 2022 reference and author exports; Kozak context conventions documented in released Isovar and primary short-leader experiments (PMID 1820208). Evidence coverage remains limited to the catalogue sources; Isovar #466 tracks missing complete 2026 curated shotgun data.
