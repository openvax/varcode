# Effect annotation

How varcode turns a variant into one or more `MutationEffect` objects.

## How it composes

A single DNA event can produce one or more plausible mutant proteins.
varcode represents each concrete mutant as a `MutantTranscript`. The
annotator turns a variant into one or more of these; the classifier
turns each `MutantTranscript` into a typed `MutationEffect`. When the
DNA alone admits multiple plausible outcomes — splice ambiguity, SV
breakpoint resolution, unphased germline-overlapping codons — the
results are packaged in a `MultiOutcomeEffect` whose `outcomes`
property exposes the set. An optional RNA-evidence resolver narrows
the set to observed isoforms or appends observed-only outcomes.

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

| # candidates | Class |
|---|---|
| 1 | plain `Substitution` / `Silent` / etc. — not wrapped |
| 2 | `ExonicSpliceSite` with `alternate_effect` |
| N | `SpliceOutcomeSet` (opt-in via `splice_outcomes=True`) |

Both `ExonicSpliceSite` and `SpliceOutcomeSet` are `MultiOutcomeEffect`
subclasses, so consumers iterate `.outcomes` uniformly without caring
about which form they're holding. `alternate_effect` works on both:
on `ExonicSpliceSite` it's the splicing-proceeds outcome directly; on
`SpliceOutcomeSet` it resolves to the `NORMAL_SPLICING` candidate's
`coding_effect`. The element types inside `.outcomes` differ
(`MutationEffect` vs `SpliceCandidate`), but `outcome.effect.short_description`
is uniform.

### Limitations

The splice classifier is **position-based** — it fires on the
canonical window (last 3 exonic, first 3-6 intronic, donor/acceptor)
and nothing else. Sequence-based signals are not flagged: exonic
splicing enhancer/silencer disruption mid-exon (~6-10nt SR-protein
motifs), branch points (~20-50nt upstream of the acceptor), deep
intronic cryptic sites. Detecting these needs ML predictors (SpliceAI,
Pangolin, MMSplice, SpliceTransformer) or direct RNA evidence;
tracked in [#297][i297].

## Annotator selection

Three annotators ship behind the `EffectAnnotator` protocol:

| Annotator | Algorithm | Used for |
|---|---|---|
| `ProteinDiffEffectAnnotator` | Builds a `MutantTranscript`, translates, diffs against the reference protein | Default for SNVs / indels / MNVs |
| `FastEffectAnnotator` | Offset arithmetic against the reference CDS | Opt-in for byte-for-byte 2.x parity or perf-sensitive paths |
| `StructuralVariantAnnotator` | Reassembles SV outcomes (deletions, duplications, inversions, fusions, translocations) | Routed automatically when the variant is a `StructuralVariant` |

All three emit the same `MutationEffect` hierarchy. `protein_diff`
catches boundary-codon and frameshift-realignment cases that
offset-arithmetic can miss; for trivial SNVs the two produce
identical output. The SV annotator dispatches on `variant.is_structural`
and isn't user-selectable for point variants.

```python
# Default (protein_diff for point variants, structural_variant for SVs):
effects = variant.effects()

# Opt into the legacy fast path:
effects = variant.effects(annotator="fast")

# Scoped swap:
with varcode.use_annotator("fast"):
    effects = variant_collection.effects()
```

Third-party annotators (isovar, Exacto) register via the registry:

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

## Structural variants

`StructuralVariant` (a `Variant` subclass) carries SV-specific fields:
`sv_type` (one of `DEL`, `DUP`, `INV`, `INS`, `CNV`, `BND`), `end`,
breakend mate fields, confidence intervals, and an open-ended `info`
dict. Pass `parse_structural_variants=True` to `load_vcf` to load
symbolic ALTs (`<DEL>`, `<INS:ME:ALU>`, `<CN0>`, breakends) as
`StructuralVariant` objects rather than dropping them.

```python
from varcode import load_vcf

vc = load_vcf("manta.vcf", parse_structural_variants=True)
sv_effects = [
    e for e in vc.effects()
    if e.variant.__class__.__name__ == "StructuralVariant"
]
```

SV effects (`LargeDeletion`, `LargeDuplication`, `Inversion`,
`GeneFusion`, `TranslocationToIntergenic`) are `MultiOutcomeEffect`
subclasses — `e.outcomes` exposes the candidate ORFs / cryptic-splice
outcomes the annotator generated, ordered by a per-class prior.
External scorers (RNA evidence, long-read assembly) plug in via
`apply_rna_evidence_to_effects` to narrow the set or append observed
outcomes; see [Germline-aware annotation](germline.md) for the same
composition pattern applied to germline.

Limitations:

- Mate breakend pairing (joining two `BND` rows that are halves of one
  translocation) is deferred. Each `BND` row produces its own
  `StructuralVariant`; consumers can match `MATEID` themselves.
- `parse_structural_variants=False` is the default. Without the flag,
  symbolic ALTs are dropped with a warning that names the flag.

## Downstream consumers

`MutantTranscript` is the prediction-boundary type for downstream
neoantigen pipelines (topiary reads `mt.mutant_protein_sequence`;
vaxrank consumes the `EffectCollection` + protein pair to score
neoantigens). RNA-evidence callers (isovar, Exacto) plug in either as
registered annotators or via the `RNAEvidenceResolver` protocol —
see [Germline-aware annotation](germline.md) for the resolver pattern,
which the same evidence shape uses across germline / phase / RNA.

[i297]: https://github.com/openvax/varcode/issues/297
