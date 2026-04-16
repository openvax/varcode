# Effect annotation

How varcode turns a variant into one or more `MutationEffect`
objects, and how the pieces fit together across the basic case,
splice-disrupting variants, and (eventually) structural variants.

## The pipeline

```
Variant ─(Annotator)─> MutantTranscript(s) ─(Classifier)─> Effect(s)
                              │
                              └── wrapped in MultiOutcomeEffect
                                   when >1 plausible outcome
                                       │
                                       └── narrowed by RNA evidence
                                           (isovar / Exacto plug-ins)
```

Single DNA event can produce one or more plausible mutant
proteins. varcode represents each concrete mutant as a
`MutantTranscript`. When the DNA alone is ambiguous, the
possibility set is wrapped in a `MultiOutcomeEffect`. RNA
evidence (when available) narrows the set.

## The four primitives

| Primitive | What it represents | Module |
|---|---|---|
| `MutationEffect` (and subclasses) | One deterministic consequence: `Substitution`, `Silent`, `FrameShift`, `PrematureStop`, ... | `varcode.effects.effect_classes` |
| `MutantTranscript` | One concrete mutant protein with edit provenance and the annotator that produced it | `varcode.mutant_transcript` |
| `MultiOutcomeEffect` | Possibility set: a sequence of candidate outcomes ordered by a prior (most likely first) | `varcode.effects.effect_classes` |
| `EffectAnnotator` | How a variant becomes effects or mutant transcripts | `varcode.annotators` |

Everything in effect annotation is an implementation or consumer
of one of those four.

## Basic usage

```python
import varcode

variants = varcode.load_maf("my_variants.maf")

# Simplest path: get an EffectCollection.
effects = variants.effects()
effects.top_priority_effect()

# Filter by transcript consequence.
nonsilent = effects.drop_silent_and_noncoding()
```

`variants.effects()` calls the current default annotator
(`"fast"`) on every `(variant, transcript)` pair and returns
an `EffectCollection`. Each element is a `MutationEffect`
subclass — `Substitution`, `Silent`, `PrematureStop`, and so
on.

## Splice-disrupting variants: two representations

When a variant sits in the canonical splice window (last 3
exonic bases, first 3–6 intronic, canonical donor/acceptor),
varcode recognizes it as splice-disrupting. Two ways the
effect is expressed, at different richness levels.

### Default: lightweight 2-outcome form

```python
variant = Variant("17", 43082575 - 5, "C", "T", "GRCh38")
effect = variant.effect_on_transcript(transcript)
# ExonicSpliceSite(...)
#   .alternate_effect -> Substitution(...)  # if splicing proceeds
```

`ExonicSpliceSite` carries `alternate_effect`: the coding
consequence that applies *if splicing still works*. Exactly
two outcomes, represented as a primary effect + one
alternate field. Cheap. Ships unconditionally.

`SpliceDonor`, `SpliceAcceptor`, and `IntronicSpliceSite`
don't expose `alternate_effect` today because the variant
is intronic — there's no coding consequence to attach.

### Opt-in: full possibility set

```python
effects = variant.effects(splice_outcomes=True)
# SpliceOutcomeSet(...) replaces the splice effect
#   .candidates ordered most-plausible-first:
#     SpliceCandidate(NORMAL_SPLICING, plausibility=0.1,
#                     coding_effect=Substitution(...))
#     SpliceCandidate(EXON_SKIPPING, plausibility=0.5,
#                     coding_effect=Deletion(...))
#     SpliceCandidate(INTRON_RETENTION, plausibility=0.3)
#     SpliceCandidate(CRYPTIC_DONOR, plausibility=0.1)
```

`SpliceOutcomeSet` replaces the splice effect with a set of
candidate outcomes, each carrying a plausibility score
(hand-tuned heuristic, not a probability) and — where
computable from cDNA — a concrete `coding_effect`. The
`NORMAL_SPLICING` candidate carries the same information as
`alternate_effect` in the default form.

When you opt in, `SpliceDonor` / `SpliceAcceptor` /
`IntronicSpliceSite` also get wrapped, so every splice-
disrupting variant produces a `SpliceOutcomeSet`.

### Relationship between the two

`SpliceOutcomeSet` is the N-outcome case of the same idea
`alternate_effect` expresses with 2 outcomes. The
`NORMAL_SPLICING` candidate is `alternate_effect` promoted to
first-class; the other candidates express alternatives you
don't see from the default form. Both are progressions of
the same pattern:

| # candidates | Class |
|---|---|
| 1 | plain `Substitution` / `Silent` / etc. — not wrapped |
| 2 | `ExonicSpliceSite` with `alternate_effect` |
| N (≥3) | `SpliceOutcomeSet` (opt-in via `splice_outcomes=True`) |

`ExonicSpliceSite` is a `MultiOutcomeEffect` subclass since
[#299][i299] Part 1. Downstream code can write
`isinstance(e, MultiOutcomeEffect)` and iterate `.candidates`
uniformly, regardless of richness level. `alternate_effect` is
available on both `ExonicSpliceSite` (the 2-outcome case, as a
first-class instance attribute) and `SpliceOutcomeSet` (where it
resolves to the `NORMAL_SPLICING` candidate's `coding_effect`).
Same field, same meaning, works on both shapes — no parallel
mechanism.

```python
# Both now work:
if isinstance(effect, MultiOutcomeEffect):
    for candidate in effect.candidates:
        ...
if effect.alternate_effect is not None:
    coding_consequence = effect.alternate_effect
```

The element types in `.candidates` differ by shape —
`ExonicSpliceSite.candidates` yields `MutationEffect` instances
directly (self + alternate_effect); `SpliceOutcomeSet.candidates`
yields `SpliceCandidate` objects with `.outcome` / `.plausibility`
/ `.coding_effect` fields. Downstream code that needs the richer
per-candidate metadata can branch on element type; code that just
wants "is there a coding consequence" reads `.alternate_effect`
and gets a `MutationEffect`-or-None regardless of which class
produced it.

### Variants in the CDS and in the splice window simultaneously

A missense at the last exonic base is biologically *both* a
splice effect (disrupts the splice signal) AND a protein-level
coding effect (changes a residue). This is exactly what
`ExonicSpliceSite` exists for:

```python
# Default: 2 outcomes
effect = ExonicSpliceSite
effect.alternate_effect      # the coding change, if splicing proceeds

# Opt-in: N outcomes including the coding change as NORMAL_SPLICING
for c in SpliceOutcomeSet.candidates:
    if c.outcome is SpliceOutcome.NORMAL_SPLICING:
        coding_change = c.coding_effect
```

Both representations encode the "splice + coding" case; the
richer `SpliceOutcomeSet` just adds the other candidates
(exon skipping, intron retention, cryptic splice).

### What varcode does NOT classify as splice-affecting today

The splice classifier is **position-based** — it fires on the
canonical window and nothing else. Variants that affect
splicing through *sequence-based* signals aren't flagged:

- **Exonic splicing enhancer / silencer (ESE / ESS) disruption**
  mid-exon: ~6–10 nt motifs that recruit or repel SR proteins.
  A missense that disrupts an ESE can cause exon skipping but
  gets classified as plain `Substitution` today.
- **Branch points** ~20–50 nt upstream of the acceptor.
- **Deep intronic cryptic splice sites** far from canonical
  boundaries.

Detecting these needs ML predictors (SpliceAI, Pangolin,
MMSplice, SpliceTransformer) trained on RNA-seq, or direct RNA
evidence. Tracked as [#297][i297].

## Annotator selection

Two annotators coexist behind the `EffectAnnotator` protocol:

| Annotator | Algorithm | Status |
|---|---|---|
| `FastEffectAnnotator` | Offset arithmetic against the reference CDS | Default |
| `ProteinDiffEffectAnnotator` | Materializes `MutantTranscript`, translates, diffs against reference protein | [#309][i309] — WIP |

Both emit the same `MutationEffect` classes, share a fast path
for trivial SNVs, and are interchangeable at the output level.
The difference is internal: protein-diff catches
boundary-codon cases and frameshift realignments that offset
arithmetic can miss.

```python
# Default (fast):
effects = variant.effects()

# Opt into protein-diff (once available):
effects = variant.effects(annotator="protein_diff")

# Scoped default swap:
with varcode.use_annotator("protein_diff"):
    effects = variant_collection.effects()
```

Third-party annotators (isovar, Exacto) register via the
registry:

```python
varcode.register_annotator(my_annotator)
variant.effects(annotator=my_annotator.name)
```

Any object exposing `name` / `supports` / `version` /
`annotate_on_transcript` satisfies the protocol.

## Provenance

Every `EffectCollection` produced by `predict_variant_effects`
records:

- `annotator` — name of the annotator that ran (`"fast"`,
  `"protein_diff"`, etc.)
- `annotator_version` — version string
- `annotated_at` — ISO-8601 UTC timestamp

Fields are preserved through `clone_with_new_elements`
(so `filter` / `groupby` keep them), written to CSV headers
(`# annotator=fast`, etc.), and recovered by `from_csv`
verbatim — restored collections remember *when* they were
originally produced.

A mismatch between the CSV's annotator and the current default
raises a warning on load; wrap `from_csv` in
`use_annotator(<csv's annotator>)` if you need the original
annotator's output specifically.

## Structural variants (roadmap)

The current `MutantTranscript` handles point variants and
indels. Fusions, translocations, inversions, and large
duplications need a multi-source generalization — a
`reference_segments` alternative to the single
`reference_transcript` field, where a fusion is two segments
from two transcripts and a translocation-to-intergenic is one
transcript segment plus a genomic-interval segment.

One DNA event can produce many candidate ORFs that only RNA
evidence can resolve (long-read transcript assembly from
Exacto, assembled contigs from isovar). SV annotators return
`List[MutantTranscript]` wrapped in a `MultiOutcomeEffect`;
each class gets a prior over outcomes so DNA-only callers can
still read `.most_likely`:

- Splice outcomes → existing plausibility tables.
- Fusions → canonical-transcript preferred, known recurrent
  junctions short-circuit to possibility-set size 1.
- Translocations → rank by `(in-frame, ORF length,
  canonical-transcript involvement)`.

`classify_from_protein_diff` (the 3d classifier) is unchanged
for SVs — once the mutant protein is materialized,
classification is the same. The difference is upstream: how
you build the `MutantTranscript`.

Tracked across [#252][i252] (fusions), [#257][i257] (SV
types), [#259][i259] (RNA evidence), [#299][i299]
(MultiOutcomeEffect formalization), [#305][i305]
(splice_outcomes rewrite on MutantTranscript).

## Downstream consumers

- **topiary** consumes mutant proteins as `ProteinFragment`.
  The `MutantTranscript` → `ProteinFragment` conversion is 1:1
  at the prediction boundary: `fragment.sequence =
  mt.mutant_protein_sequence`; provenance fields carry through
  as fragment metadata. Use `MutantTranscript` for
  transcript-level reasoning (splicing, UTR readthrough,
  codon tables); convert at the prediction boundary.
- **vaxrank** consumes the `EffectCollection` +
  `ProteinFragment` pair to score neoantigens. After
  [#305][i305], splice outcomes deliver `MutantTranscript`s
  directly; after fusion support, SV events deliver
  `MultiOutcomeEffect[List[MutantTranscript]]` consumable the
  same way.
- **isovar** / **Exacto** plug in as registered annotators,
  producing `MutantTranscript`s from assembled RNA evidence
  rather than inferred from DNA alone.

[i252]: https://github.com/openvax/varcode/issues/252
[i257]: https://github.com/openvax/varcode/issues/257
[i259]: https://github.com/openvax/varcode/issues/259
[i271]: https://github.com/openvax/varcode/issues/271
[i297]: https://github.com/openvax/varcode/issues/297
[i299]: https://github.com/openvax/varcode/issues/299
[i305]: https://github.com/openvax/varcode/issues/305
[i309]: https://github.com/openvax/varcode/issues/309
