# Common effects, noncanonical ORFs, and reference validation

## Scope and order

1. Ship Varcode's shared candidate traversal first (#653). Expose typed traversal
   of ordinary, composed, legacy, structural, and ORF effects. Preserve selected
   effects, alternative paths, patient baselines, unknown sequence, and provenance.
   Keep downstream Isovar/Vaxrank changes with their owning agents.
2. Audit noncanonical ORF consequence fields and defaults. Annotate supplied ORFs
   alongside host transcripts; cover changes to their coding paths independently
   of host biotype. Implement bounded new-start discovery where appropriate (#583),
   with explicit sequence-only evidence and unavailable/discovery status.
3. Add publicly sourced native GRCm39 mouse and ROS_Cfam dog variants, including
   both strands and available indel/splice/mitochondrial/phase/rearrangement cases.
   Preserve source coordinates, annotations, versions, and acquisition checksums.
   Freeze independent expected consequences/sequences and state the limits of
   each source; never generate expectations using Varcode or relabel lifted calls
   as native calls.
4. Complete and validate CHM13 setup: pinned annotation and chromosome DNA,
   ordinary/intronic/structural annotation, reference identity in exports, and
   native benchmark variants including duplicated/newly resolved sequence.
   Make conversion explicit; reference aliases must not perform liftover.

## First implementation: candidate traversal

Use existing effect/candidate classes and a shared child-enumeration helper.
Provide an identity-deduplicated effect iterator for filtering and a path iterator
that retains distinct provenance. Terminate self-references/cycles per path without
discarding candidate records. Reuse this in existing sequence-status, ranking,
and proteome consumers only where their aggregation contracts agree. Document
candidate-specific baseline/mutant sequences and unknown flags. Add meaningful
static consumer checks and runtime regressions using the public CFTR example,
nested wrappers, shared effects under distinct evidence, unresolved products,
and supplied noncanonical ORFs.

## Validation and releases

Read primary scientific sources before changing consequence semantics. Run focused
regressions, ./lint.sh, ./test.sh, and dependency checks in an isolated release
environment. Each PR includes a version bump and release notes; wait for green CI,
merge, then run ./deploy.sh from clean synchronized main. Verify published wheels.
Preserve existing checkouts and the shared Python environment. File new problems
on the owning repository and link them from the relevant PR.
