# Varcode, Isovar, and Vaxrank

The libraries answer different questions about the same variant. Use Varcode
for predicted consequences, Isovar for RNA reconstruction and evidence, and
Vaxrank for evaluating protein/peptide candidates.

| Library | Responsibility |
|---|---|
| **Varcode** | Generate structural/transcript hypotheses and predict their coding consequences. |
| **Isovar** | Reconstruct RNA-supported sequences, compare them with those hypotheses, and preserve unresolved alternatives. |
| **Vaxrank** | Evaluate the resulting protein/peptide candidates, retaining their evidence. |

This is the shared responsibility split, not a claim that every path is already
connected. The implementation limits below matter when building a pipeline.

## Varcode's part

Varcode models possible effects on transcripts and proteins. Its inputs can
include DNA variants, germline/phase context, and externally reconstructed RNA
structures. It should not independently collect BAM reads or decide which RNA
assembly the patient expresses.

Use the ordinary `effects()` interface; no special annotator selection is needed.
Keep candidate sets when a variant has several plausible consequences. Candidate
order and effect severity are not RNA support or calibrated probabilities.

Isovar can return a structure that changes the original hypothesis, including a
previously unmodeled junction. Varcode's job is to represent that structure and
predict its coding consequences, not force it into the first reference isoform.
[RNA imports](rna_structures.md) describe the existing entry points.

## What must survive the handoff

- Nucleotide sequence and structure, coordinates, reference/annotation identity,
  and compatible transcript IDs—not only a protein string.
- Partial versus complete sequence, reading-frame evidence, and unresolved regions.
- Sample/library and read/fragment provenance, support, conflicts, and assumptions.
- Separate alternatives when evidence cannot distinguish them. Grouping identical
  proteins must preserve the contributing structures and their evidence.

RNA support for a local junction does not establish a full-length transcript or
prove translation. Missing coverage is not evidence against a hypothesis. A
predicted protein change and evidence that the RNA exists are separate facts.

## Available today and remaining work

- **Varcode:** structural partner candidates and supplied RNA imports exist.
  Splice/phase/SV combinations are not exhaustively composed
  ([#423](https://github.com/openvax/varcode/issues/423)); combined-haplotype
  ownership still needs cleanup ([#437](https://github.com/openvax/varcode/issues/437)).
  Imported partial proteins also need completeness-aware change flags
  ([#462](https://github.com/openvax/varcode/issues/462)); do not interpret an
  unobserved suffix as a demonstrated protein deletion.
- **Isovar:** the small-variant path reconstructs RNA context; the separate
  supplied-fusion path validates RNA and retains alternative frame hypotheses.
  Automated collection, alternative-path assembly, and competitive reconciliation
  for nominated SVs are tracked in
  [Isovar #305](https://github.com/openvax/isovar/issues/305).
- **Vaxrank:** the ordinary RNA path selects Isovar's top protein; a separate
  supplied-fusion adapter retains coding hypotheses. The opt-in DNA fallback
  still collapses fusion candidates too early
  ([Vaxrank #482](https://github.com/openvax/vaxrank/issues/482)).

A DNA-only fallback remains a prediction without RNA confirmation. It must not
silently discard alternatives or relabel missing RNA as evidence of absence.

## Other library guides

- [Isovar: RNA reconstruction and reconciliation](https://github.com/openvax/isovar/blob/master/docs/library-responsibilities.md).
- [Vaxrank: candidate evaluation and evidence](https://openvax.github.io/vaxrank/library-responsibilities/).
