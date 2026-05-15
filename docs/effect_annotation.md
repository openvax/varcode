# Effect annotation

How varcode turns a variant into one or more `MutationEffect` objects.

## How it composes

A single DNA event can produce one or more plausible mutant proteins.
varcode represents each concrete mutant as a `MutantTranscript`. The
annotator turns a variant into one or more of these; the classifier
turns each `MutantTranscript` into a typed `MutationEffect`. When the
DNA alone admits multiple plausible outcomes — splice ambiguity, SV
breakpoint resolution, unphased germline-overlapping codons — the
results are packaged in a `MultiOutcomeEffect` whose `candidates`
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
# SpliceOutcomeSet(...) replaces the splice effect.
# .candidates is a tuple[EffectCandidate, ...], in producer order.
# Each candidate's .effect is a SpliceMechanismEffect subclass:
#   EffectCandidate(effect=ExonSkipping(affected_exon=..., in_frame=True,
#                                       aa_ref="KGYK...", ...))
#   EffectCandidate(effect=IntronRetention(retained_intron_start=...,
#                                          side="donor", ...))
#   EffectCandidate(effect=CrypticDonor(affected_exon=..., ...))
#   EffectCandidate(effect=NormalSplicing(coding_effect=Substitution(...)))
```

`SpliceOutcomeSet` replaces the splice effect with a set of
candidate mechanisms. Class identity = mechanism — `NormalSplicing`,
`ExonSkipping`, `IntronRetention`, `CrypticDonor`, `CrypticAcceptor`.
Each is a `SpliceMechanismEffect` subclass that carries its own
protein vocab on the instance (`aa_ref`, `aa_alt`,
`mutant_protein_sequence`, `mutant_transcript`); these are `None`
when the protein math couldn't resolve (e.g. intron retention
without a genomic-sequence provider), populated otherwise. Each
mechanism also carries `splice_signal` — the underlying
`SpliceDonor` / `SpliceAcceptor` / `IntronicSpliceSite` /
`ExonicSpliceSite` effect describing *where* the disruption was.

Downstream consumers dispatch by class:

```python
for c in splice_set.candidates:
    if isinstance(c.effect, ExonSkipping):
        print(c.effect.affected_exon.exon_id, c.effect.in_frame)
    elif isinstance(c.effect, IntronRetention):
        print(c.effect.side, c.effect.retained_intron_start)
```

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
subclasses, so consumers iterate `.candidates` (a tuple of
`EffectCandidate` objects) uniformly without caring about which form
they're holding. `alternate_effect` works on both: on
`ExonicSpliceSite` it's the splicing-proceeds outcome directly; on
`SpliceOutcomeSet` it resolves to the inner effect of the
`NormalSplicing` candidate (or `None` when that candidate is just
a placeholder). `candidate.effect.short_description` is uniform
across both forms.

With RNA evidence, splice sets are reconciled rather than merely
extended. `SpliceOutcomeSet.with_rna_evidence(...)` returns a new set
whose `candidates` are the RNA-observed mechanisms, while
`dna_candidates`, `rna_evidence`, `excluded_candidates`,
`added_candidates`, and `candidate_rna_evidence` preserve the audit
trail. Use `splice_set.rna_evidence_for(candidate)` to inspect the
observations supporting one current candidate.

### Candidate provenance

There is no `plausibility` or `probability` field in the shared
candidate wrapper. The old splice-specific `plausibility` value was a
DNA-only ordering heuristic, not evidence. Varcode now keeps that
ordering only as producer order.

Producer-specific support belongs in `candidate.evidence` under
explicit names: `read_count`, `junction_id`, `psi`, `motif_score`,
`donor_score`, `acceptor_score`, and so on. Varcode stores evidence
as opaque provenance and does not normalize it into a probability.

### Picking a single candidate

When you need to collapse a multi-outcome effect to one Effect, two
notions of "best" are available — pick consciously:

| Accessor | Returns | Meaning |
|---|---|---|
| `.most_likely_candidate` | `EffectCandidate` | First candidate after producer ordering |
| `.most_likely_effect` | `MutationEffect` | Inner effect of the above |
| `.highest_priority_candidate` | `EffectCandidate` | Top by `effect_priority` (most protein-disruptive) |
| `.highest_priority_effect` | `MutationEffect` | Inner effect of the above |

The `_candidate` accessors keep the provenance wrapper (`.source`,
`.evidence`); the `_effect` accessors peel it off. The two "top by"
notions coincide whenever producer ordering and priority ranking
agree, which is common — but for clinical / functional filtering
("flag if any candidate is at least a frameshift") prefer
`highest_priority_*`: a disruptive candidate behind a less-disruptive
primary candidate should still light up.

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
subclasses — `e.candidates` exposes the candidate ORFs / cryptic-splice
outcomes as a tuple of `EffectCandidate` objects in producer order.
External evidence producers (RNA evidence, long-read assembly)
plug in via `apply_rna_evidence_to_effects` to append observed
candidates; see [Germline-aware annotation](germline.md)
for the same composition pattern applied to germline.

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
