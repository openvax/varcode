# Public somatic variant validation

## Scope

Add complete, version-pinned public HCC1395 and COLO829 benchmark call sets
(small variants and published structural variants), and published canine
somatic variants with independent consequence annotations. Keep original
source identity, confidence, assembly, transcript release and licensing.

1. Identify primary releases and record checksums, row counts and variant types.
   Use complete truth call sets, including variants outside high-confidence
   intervals; preserve confidence categories without treating candidate calls
   as validated truth. Explain any unavailable source explicitly.
2. Add reproducible data provisioning and fixtures appropriate for routine CI.
   Exercise every loaded allele in the complete datasets, accounting for source
   rows, duplicates and paired breakends; no silent dropped records.
3. Annotate every benchmark variant against the correct reference. Report
   errors and unresolved annotations by reason; fix Varcode defects found.
4. Compare canine consequences to published transcript-specific annotations
   on matching reference versions. Retain and explain all disagreements;
   never generate the expected labels with Varcode itself.
5. Run focused and full tests, lint, CI, merge and release a patch after 10.18.0.

## Boundaries

Varcode implementation only. File upstream reference or dataset problems in
the appropriate repository. Public truth sets are a reproducible definition
of the requested benchmark coverage, not every historical caller output for
these cell lines. Primary sources and pinned versions define completeness.
