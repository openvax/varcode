# Effects API

See [effect annotation](effect_annotation.md) for result access and the
[effect catalog](effect_types.md) for individual consequence classes.

## Effects

::: varcode.MutationEffect
    options:
      skip_local_inventory: true

::: varcode.NonsilentCodingMutation
    options:
      skip_local_inventory: true

::: varcode.MultiOutcomeEffect
    options:
      skip_local_inventory: true

::: varcode.EffectCollection

::: varcode.EffectCandidate

::: varcode.candidates_from_effects

## Candidate traversal

See [traversing candidates](experimental_annotators.md#traversing-candidates)
for terminal comparisons, provenance paths and unresolved states.

::: varcode.get_effect_candidates

::: varcode.iter_candidate_paths

::: varcode.iter_candidate_effects

## Splice candidates

See [splice variants](splice_variants.md) for the distinction between signal
disruption and candidate consequences. The scores below are sequence heuristics,
not probabilities that a candidate is expressed.

::: varcode.SpliceCandidateRNAEvidence

::: varcode.enumerate_cryptic_exon_candidates

::: varcode.enumerate_from_structural_variant

::: varcode.score_donor

::: varcode.score_acceptor

## Protein collections

See [protein collections and FASTA](proteome.md) for `varcode.to_proteome()`,
offline provenance, alternatives and reference composition.

## Priority ordering

::: varcode.effect_priority

::: varcode.classification_priority

::: varcode.top_priority_effect
