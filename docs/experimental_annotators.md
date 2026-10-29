# Annotation pipeline and compatibility

Since Varcode 11, the default `transcript_model` annotator uses one path for
reference edits, patient germline baselines, known-cis groups and splice
alternatives: realize the baseline and mutant transcripts, then classify their
sequence difference. Read collection, alignment and RNA assembly remain in
upstream tools such as Isovar. Varcode consumes their supplied models.

## Supported inputs

| Selection | Method | Inputs |
|---|---|---|
| Default, or `annotator="transcript_model"` | Composes phase/splice hypotheses and compares each product with its patient baseline | Small variants and local DEL/DUP; sequence insertions need `alt_assembly`. Reuses existing structural builders for junctions, inversions and supplied assemblies. |
| `annotator="fast"` | Compatibility implementation from Varcode 10 | Small variants and SVs; local germline lookup; limited joint composition. |
| `annotator="protein_diff"` | Alternative protein-difference implementation | Small variants; shares legacy splice, location and germline helpers; no SV support. |

Typed junctions, BNDs, inversions, nonlocal rearrangements and supplied
transcript assemblies reuse the existing structural builders. Combining these
with germline or other cis edits requires mapped patient-baseline coordinates
and remains unresolved. CNVs retain dosage uncertainty. See [limitations](#known-limitations).

## Comparing annotators

Given a `variant` and one of its `transcript` objects, selection is explicit:

```python
effect = variant.effect_on_transcript(transcript)  # ordinary use
comparison = variant.effect_on_transcript(transcript, annotator="protein_diff")
legacy = variant.effect_on_transcript(transcript, annotator="fast")
```

The transcript model can use `germline=` and `phase_resolver=` context. Canonical
and exon-skip paths need only transcript annotation; a genome with reference
FASTA additionally supplies sequence for intron retention and cryptic splice
sites. Without a calibrated scorer, mechanism preference is an ordering rule,
not a probability. No implementation resolves every structural or combined input.

Use the same entry point as more evidence becomes available:

```python
effects = variant.effects(
    germline=germline_ctx,       # or None when unavailable
    phase_resolver=phaser,       # or None when phase is unknown
    rna_resolver=rna_evidence,   # or None when no RNA models are available
)
```

The germline lookup covers the whole transcript span, including introns that
a splice hypothesis may retain. This can expose more phase uncertainty than
the legacy annotator's local lookup. Sparse/hotspot contexts retain
`germline_unknown=True` even when some inherited alleles are supplied.

Through `effects(germline=...)`, the transcript model uses the context phase cap (`GermlineContext.max_phase_hypotheses`, default 8). It
also caps the combined phase × baseline-splice × mutant-splice outcomes at `max_hypotheses=64`, which is
the phase cap too when `predict_transcript_model_effect` is called directly.
Exceeding either returns a [`HypothesisLimit`](germline.md#phase-enumeration-limit)
rather than a partial set of candidates.

### Mixed inherited and somatic haplotypes

`predict_transcript_model_effect(variants, transcript, germline_variants=...)`
and `annotate_haplotype` take a known-cis group. Members matching the supplied
germline are included once in the patient baseline. Only novel alleles are
applied to that baseline to obtain the mutant; `effect.variant` is the first
novel member and `effect.variants` retains the supplied group. Matching uses
reference assembly, contig, normalized position, REF and ALT. A group with no
novel allele returns `GermlineAlleleOverlap`, without inferring LOH.

The group itself supplies cis constraints. Phase can propagate from any member
to other germline alleles; unlinked alleles still produce separate hypotheses.
Homozygous alleles are on both haplotypes and do not bridge their phase.
As with other phase constraints, conflicting resolver answers are logged and
ignored. Each hypothesis retains its germline cis/trans assignments, phase
state, resolver source and novel `somatic_variants` in its evidence.

`VariantCollection.effects` constructs these groups only when the resolver
supports cis. Unknown or trans relationships retain individual predictions;
listing alleles in a collection alone does not assert that they are in cis.

## Transcript-model results

Use `effect.baseline_transcript` and `effect.mutant_transcript` for the products
used in classification, and the usual `original_protein_sequence` and
`mutant_protein_sequence` properties for proteins. Candidate effects use these
same accessors. See [sequence availability](transcript_models.md#sequence-access).

The shared pipeline and legacy candidate wrappers expose different details:

| Candidate type | Shared access | Additional information |
|---|---|---|
| `EffectCandidate` (legacy splice/SV sets and imported RNA) | `candidate.effect` | `source`, `evidence`; sequences are on the effect or its optional `mutant_transcript` |
| `RealizedEffectCandidate` (transcript-model hypothesis pipeline) | `candidate.effect` | `outcomes`, `hypotheses`, `probability`, `ordinal_key`; each outcome has `baseline` and `mutant` products with `cdna_sequence`, `protein_sequence` and `evidence` |

The default returns an ordinary top effect with candidates attached, not
necessarily a `MultiOutcomeEffect`. Its structural builders use ordinary candidates;
early noncoding, incomplete or unresolved results may have no candidates at all.
Inspect the returned shape rather than assuming it from the annotator name:

```python
from varcode import EffectCandidate
from varcode.effect_hypotheses import RealizedEffectCandidate

print(effect.short_description)
for candidate in getattr(effect, "candidates", ()):
    print(candidate.effect.short_description)
    if isinstance(candidate, RealizedEffectCandidate):
        # A merged candidate can retain several hypotheses and their products.
        for outcome in candidate.outcomes:
            print(outcome.baseline.protein_sequence, outcome.mutant.protein_sequence)
            print(outcome.mutant.cdna_sequence, outcome.hypothesis.evidence)
    elif isinstance(candidate, EffectCandidate):
        print(candidate.source, candidate.evidence)
        print(candidate.effect.mutant_protein_sequence)
```

### Baseline splice uncertainty

Each hypothesis records `baseline_splice_plan` and `splice_plan`. Unchanged
active sites share a mechanism between baseline and mutant. If the site changes
or disappears on one side, all compatible baseline choices remain; Varcode does
not select an arbitrary patient splice product. Both plans are graph-valid.
Unchanged annotated noncanonical sites are not called variant-induced loss.

### Migrating from Varcode 10

Ordinary calls now return a classified effect such as `Deletion`, `Substitution`
or `Unresolved`, with alternatives in `candidates`. Code requiring
`SpliceOutcomeSet`, `LargeDeletion`, `PhaseCandidateSet` or `HaplotypeEffect`
results can explicitly select `annotator="fast"`. Existing classes, imports and
serialized results remain supported. A candidate rank is not a probability.

The shared default builds full baseline/mutant products and retains uncertainty;
this costs more CPU and memory than the legacy local classifier, especially
across large call sets or many transcripts. For workflows requiring the older
throughput and result contract, select `annotator="fast"` explicitly.

Collection filtering and effect priority consider every candidate, so a
preferred silent comparison cannot hide a protein-changing alternative. The
selected effect's own sequence properties describe that comparison; use
`highest_priority_candidate.effect` when you need the most severe alternative.

## Known limitations

RNA evidence refines the existing hypothesis candidates. An exact match to a
full predicted cDNA adds `rna_support` to that product's evidence and ranks
supported comparisons first. All compatible phase baselines remain, and
unobserved alternatives are retained. Read counts do not become probabilities,
and RNA sequence agreement does not prove translation.

Partial or unmatched RNA models retain their actual supplied sequence as an
`Unresolved` candidate with an unavailable baseline. No reference sequence is
appended, and no phase or CDS mapping is inferred. Such observations follow the
DNA candidates in ordering. The same behavior applies to RNA models supplied
by a phase resolver for a known-cis group; no separate observation list is
needed on this path. Legacy delegated BND/assembly results keep their existing
contracts. Read-based transcript reconstruction belongs upstream; broader mapped
rearrangement consequence support is tracked in
[#466](https://github.com/openvax/varcode/issues/466).
Models containing assumed reference completion do not count as full observed
RNA support, even when their completed cDNA matches a prediction.

- Combined variants: the selected annotator makes known-cis joint predictions
  through `VariantCollection.effects(phase_resolver=...)`. The default composes
  germline, splice and local structural edits. `fast` and `protein_diff` retain
  their legacy point-edit haplotype behavior. Unsupported groups remain
  explicit unresolved results. See the [joint contract](annotator_contract.md#joint-haplotypes).
- Rearrangements beyond the layout: existing structural builders retain their
  fusion candidates and supplied transcript sequences. A clipped DUP/INV cannot
  establish an unchanged protein or a complete rearranged transcript
  ([#449](https://github.com/openvax/varcode/issues/449)); unavailable sequence
  and unmapped CDS/translation remain uncertain. Local symbolic DEL/DUP events
  use the shared genomic layout.

Unresolved results report `None`, not `False`, for sequence-change flags.
`drop_silent_and_noncoding()` retains them by default; see
[SV filtering](structural_variants.md#filtering-by-protein-change).

## Previous names

`realized`, `RealizedEffectAnnotator`, and `predict_realized_effect` are aliases
for the transcript-model implementation. Use the `transcript_model` names in
new code. The `RealizedEffectCandidate` wrapper still has a different interface
from ordinary `EffectCandidate`, as shown above.

To implement your own, see [Writing an annotator](annotator_contract.md).
