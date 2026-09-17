# Effect annotation

How varcode turns a variant into one or more `MutationEffect` objects.

## Basic usage

```python
import varcode

variants = varcode.load_maf("my_variants.maf")

effects = variants.effects()
effects.top_priority_effect()
```

For ordinary use, this is all the selection you need. The default handles
point variants and structural variants through the same interface; you do
not need to choose an annotator. See [structural variant loading](#structural-variants)
if your input includes SV records.

The result is an `EffectCollection` with a prediction for each relevant
`(variant, transcript)` pair. Effects include `Substitution`, `Silent`,
`FrameShift`, structural consequences and unresolved predictions.

## How it composes

An effect describes a predicted consequence. It may also carry a
`MutantTranscript`: edits or reference-segment structure, provenance, and
optional mutant cDNA and protein sequences. This is **not necessarily a known
protein or complete RNA molecule**. It can describe a partial retained fragment,
or a structure whose sequence is still unknown. Not every effect has an attached
`MutantTranscript`; some effects expose a predicted protein directly.

When several outcomes are plausible, ordinary outcome sets expose alternatives
through `MultiOutcomeEffect.candidates`. Their order reflects the producer's
rules or evidence, not necessarily calibrated probabilities. An optional
RNA-evidence resolver can refine the alternatives or add observed outcomes.
Experimental result shapes are described in the [advanced section](#annotator-selection).

## The four primitives

| Primitive | What it represents | Module |
|---|---|---|
| `MutationEffect` (and subclasses) | A predicted consequence, unresolved result, or outcome set | `varcode.effects.effect_classes` |
| `MutantTranscript` | Edits or transcript structure with provenance; cDNA and protein sequences are optional | `varcode.mutant_transcript` |
| `MultiOutcomeEffect` | A set of candidate effects in producer-defined order | `varcode.effects.effect_classes` |
| `EffectAnnotator` | An implementation that predicts effects; ordinary callers use the default | `varcode.annotators` |

## Splice-disrupting variants

A single nucleotide change near an exon-intron boundary can hit
the splice signal *and* the coding sequence at the same time.
The splice surface captures both possibilities, gives every
splice-disrupting variant a uniform candidate-set shape, and
exposes accessors for the "what if splicing still proceeds?"
question.

### When splice disruption is in play

The classifier is **position-based**: it fires when a variant
lands in the canonical splice window around an exon-intron
boundary. The window is asymmetric — the donor consensus
(`MAG|GURAGU`) is wider on both sides than the acceptor
consensus (`YAG|R`):

- **exonic side**: the last 3 bases of an exon (donor side) or
  the first base of the next exon (acceptor side)
- **intronic side**: positions +1..+6 of the intron (donor side)
  and positions -3..-1 (acceptor side), including the canonical
  `GT` at +1/+2 and `AG` at -2/-1

Four classes record *where* in this window the variant landed:

| Class | Position |
|---|---|
| `ExonicSpliceSite` | Last 3 bases of an exon (donor side) or the first base of the next exon (acceptor side) |
| `SpliceDonor` | Canonical `GT` at intronic +1 / +2 |
| `SpliceAcceptor` | Canonical `AG` at intronic -2 / -1 |
| `IntronicSpliceSite` | Intronic +3..+6 (donor side) or -3 (acceptor side); also `+1/+2` or `-1/-2` when the reference base isn't the canonical `GT` / `AG` |

Variants outside this window are **not** flagged as
splice-disrupting, even when they may affect splicing
biologically — ESE/ESS motifs mid-exon, branch points ~20–50 bp
upstream of the acceptor, deep intronic cryptic activation.
Detecting those requires ML predictors or direct RNA evidence;
see [Limitations](#limitations).

### Splice and coding effects can co-occur

A variant in an exon sits on a coding base by definition — it
rewrites a codon. If that same exonic base is **also** in the
splice window (the exonic positions in the table above), the
same nucleotide change disrupts the splice signal *and* changes
the protein. varcode represents this duality as
**`ExonicSpliceSite`**:

- on the default 2-outcome shape, splice disruption is the
  primary effect; the coding consequence (a `Substitution`,
  `Silent`, etc.) hangs off `.alternate_effect`
- on the opt-in `SpliceOutcomeSet` shape, the same coding
  consequence is the `coding_effect` of the `NormalSplicing`
  candidate, reachable through
  `splice_set.effect_if_splicing_unchanged`

For purely **intronic** disruptions (`SpliceDonor`,
`SpliceAcceptor`, `IntronicSpliceSite`), there is no codon to
rewrite — the variant doesn't change a coding base. The default
shape doesn't expose `alternate_effect` on these classes; the
opt-in shape's `effect_if_splicing_unchanged` returns `None`.

For coding variants **outside** the splice window, varcode emits
a plain coding effect (`Substitution`, `Silent`, `FrameShift`,
…) with no splice annotation attached. The variant may still
disrupt splicing through a non-canonical mechanism, but varcode
won't flag it — see Limitations.

### The `SpliceOutcomeSet` shape

Every splice-disrupting variant emits a `SpliceOutcomeSet` — there
is no "bare splice class" path at the user-facing API as of
varcode 6.0.

```python
variant = Variant("17", 43082575 - 5, "C", "T", "GRCh38")
splice_set = variant.effect_on_transcript(transcript)
# SpliceOutcomeSet(disrupted_signal_class=ExonicSpliceSite, ...)
# .candidates is a tuple[EffectCandidate, ...] in producer order.
# Each candidate's .effect is a SpliceMechanismEffect subclass:
#   EffectCandidate(effect=NormalSplicing(coding_effect=Substitution(...)))
#   EffectCandidate(effect=ExonSkipping(affected_exon=..., in_frame=True,
#                                       aa_ref="KGYK...", ...))
#   EffectCandidate(effect=IntronRetention(retained_intron_start=...,
#                                          side="donor", ...))
#   EffectCandidate(effect=CrypticDonor(affected_exon=..., ...))
```

`SpliceOutcomeSet` carries:

- `disrupted_signal_class` — the `SpliceSite` subclass (`SpliceDonor`,
  `SpliceAcceptor`, `ExonicSpliceSite`, or `IntronicSpliceSite`)
  identifying where in the splice window the variant landed
- `candidates` — a tuple of `EffectCandidate` objects in producer
  order, one per plausible mechanism
- `effect_if_splicing_unchanged` — the coding consequence that
  applies if the spliceosome still splices normally (the
  `NormalSplicing` candidate's `coding_effect`), or `None` for
  purely intronic disruptions where the nucleotide change doesn't
  touch a coding base. Also exposed as `alternate_effect` for
  back-compat with code that read `ExonicSpliceSite.alternate_effect`

Each candidate's `.effect` is a `SpliceMechanismEffect` subclass
that carries its own protein vocab on the instance (`aa_ref`,
`aa_alt`, `mutant_protein_sequence`, `mutant_transcript`). Fields
are `None` when the protein math couldn't resolve (e.g. intron
retention without a `genomic_sequence` provider), populated
otherwise. Each mechanism also exposes `splice_signal` — the
underlying raw `SpliceDonor` / `SpliceAcceptor` /
`IntronicSpliceSite` / `ExonicSpliceSite` effect describing *where*
the disruption was.

**Lazy construction.** Only the cheap `NormalSplicing` candidate
is built eagerly when the set is constructed; `ExonSkipping`,
`IntronRetention`, and `CrypticDonor`/`CrypticAcceptor` materialise
on first `.candidates` access and are cached. Filter pipelines
that drop variants early via `modifies_protein_sequence` /
`effect_priority` never trigger the expensive candidates.

Downstream consumers dispatch by class:

```python
for c in splice_set.candidates:
    if isinstance(c.effect, ExonSkipping):
        print(c.effect.affected_exon.exon_id, c.effect.in_frame)
    elif isinstance(c.effect, IntronRetention):
        print(c.effect.side, c.effect.retained_intron_start)
```

### Common questions

A cheat sheet for the simple splice use cases. `splice_set` is a
`SpliceOutcomeSet` (every splice-disrupting variant produces one).

**Is this variant splice-disrupting?**

```python
from varcode import MultiOutcomeEffect, SpliceOutcomeSet

# Splice-specific check:
isinstance(effect, SpliceOutcomeSet)

# Or by disrupted signal class:
isinstance(effect, SpliceOutcomeSet) and effect.disrupted_signal_class is SpliceDonor

# Broader: any multi-outcome effect, including SV outcomes
# (LargeDeletion, GeneFusion, ...) — use when you want one
# uniform handler for splice + SV ambiguity.
isinstance(effect, MultiOutcomeEffect)
```

**What coding consequence applies if splicing still proceeds?**

```python
coding = splice_set.effect_if_splicing_unchanged   # canonical
coding = splice_set.alternate_effect               # back-compat alias

# Either returns the NormalSplicing candidate's coding_effect (a
# Substitution / Silent / PrematureStop / ...), or None for purely
# intronic disruptions where the variant doesn't change a coding base.
```

**What's the most likely splice mechanism?**

```python
splice_set.most_likely_effect                # SpliceMechanismEffect
splice_set.most_likely_candidate             # EffectCandidate (.effect + .source/.evidence)
```

**What are all candidate outcomes?**

```python
for candidate in splice_set.candidates:
    candidate.effect      # SpliceMechanismEffect (ExonSkipping, IntronRetention, ...)
    candidate.source      # producer name
    candidate.evidence    # opaque dict of provenance fields
```

**Which outcome is the most disruptive?**

```python
splice_set.highest_priority_effect           # most protein-disruptive
splice_set.highest_priority_candidate
```

Use this for clinical / functional filtering ("flag if any
candidate is at least a frameshift") — a disruptive candidate
ranked below a less-disruptive primary should still light up.
See [Picking a single candidate](#picking-a-single-candidate)
for the "most likely" vs "most disruptive" distinction.

**What protein sequences could result?**

```python
splice_set.candidate_proteins                # {ExonSkipping: "MA...", IntronRetention: "", ...}
splice_set.mutant_protein_sequences          # set[str] of distinct non-empty sequences
```

Empty string means the mechanism's protein math couldn't resolve
(typically: no `genomic_sequence` provider, so `IntronRetention`
and `CrypticDonor` stay predicted-only).

**Where on the transcript is the splice signal?**

```python
for candidate in splice_set.candidates:
    candidate.effect.splice_signal           # SpliceDonor / SpliceAcceptor / IntronicSpliceSite / ExonicSpliceSite
```

### RNA evidence reconciliation

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

Sequence-based splice signals are not flagged: exonic splicing
enhancer/silencer disruption mid-exon (~6-10nt SR-protein motifs),
branch points (~20-50nt upstream of the acceptor), deep intronic
cryptic sites. Detecting these needs ML predictors (SpliceAI,
Pangolin, MMSplice, SpliceTransformer) or direct RNA evidence;
tracked in [#297][i297].

<a id="annotator-selection"></a>

## Advanced: annotators and implementation limits

Skip this section unless you are comparing implementations, using an experiment,
or writing an integration. Ordinary callers should keep using `variants.effects()`.

### What the optional implementations change

| Selection | How it predicts | Supported scope and limits |
|---|---|---|
| Default (omit `annotator=`) | Established point-edit prediction, internal structural handling, and patient-baseline comparison for point edits | Point variants and SVs; germline/phase context for point edits. General SV-plus-germline composition is not supported. |
| `annotator="protein_diff"` | Constructs and translates point-edited transcripts, then compares proteins; reuses the default's splice/location and germline helpers | Experimental alternative for point edits, not broader biological coverage. SV inputs are unsupported. |
| `annotator="transcript_model"` | Enumerates phase/splice hypotheses, constructs transcript products, compares against the patient's baseline, then merges equivalent results | Experimental point edits and local DEL/DUP/INV layouts. Insertions require `alt_assembly`; CNVs are unsupported. BNDs use the existing structural helper, but BND-plus-haplotype composition is unresolved. Nonlocal DUP/INV handling remains limited; see #449 below. |

The default's historical registry name is `fast`; you may see it in provenance
headers. You do not need to pass that name. There is no separate selectable
structural annotator or outside router switching implementations for you.
**Getting protein sequences does not require selecting `protein_diff`.** Protein
comparison is also a shared helper used by the other prediction paths.

Given a `variant` and one of its `transcript` objects, selection is explicit:

```python
effect = variant.effect_on_transcript(transcript)  # ordinary use
comparison = variant.effect_on_transcript(transcript, annotator="protein_diff")
experimental = variant.effect_on_transcript(transcript, annotator="transcript_model")
```

The transcript model can use `germline=` and `phase_resolver=` context. Canonical
and exon-skip paths need only transcript annotation; a genome with reference
FASTA additionally supplies sequence for intron retention and cryptic splice
sites. Without a calibrated scorer, mechanism preference is an ordering rule,
not a probability. Selecting this experiment does not guarantee that every
structural or combined input is supported.

### Reading sequences and alternatives

For a single effect, distinguish an effect's predicted protein from its optional
transcript model:

```python
protein = effect.mutant_protein_sequence  # may be None
model = effect.mutant_transcript         # may be None even when protein is available
cdna = model.cdna_sequence if model is not None else None
evidence = model.evidence if model is not None else None
```

`None` does not mean an unchanged protein. A model containing reference segments
can still be partial: for example, a BND's retained fragment is labeled
`sequence_status="retained_reference_fragment"`, with full cDNA/protein unknown.
Do not concatenate that fragment and call it a complete allele. See
[observed RNA structures](structural_variants.md#importing-observed-rna-structures)
for sequence completeness and prediction-versus-observation details.

The ordinary and experimental candidate wrappers are not interchangeable:

| Candidate type | Shared access | Additional information |
|---|---|---|
| `EffectCandidate` (ordinary splice/SV/RNA outcome sets) | `candidate.effect` | `source`, `evidence`; sequences are on the effect or its optional `mutant_transcript` |
| `RealizedEffectCandidate` (transcript-model hypothesis pipeline) | `candidate.effect` | `outcomes`, `hypotheses`, `probability`, `ordinal_key`; each outcome has `baseline` and `mutant` products with `cdna_sequence`, `protein_sequence` and `evidence` |

The experiment returns an ordinary top effect with candidates attached, not
necessarily a `MultiOutcomeEffect`. Its BND delegation uses ordinary candidates;
early noncoding, incomplete or unresolved results may have no candidates at all.
Inspect the returned shape rather than assuming it from the annotator name:

```python
from varcode import EffectCandidate
from varcode.effect_hypotheses import RealizedEffectCandidate

print(experimental.short_description)
for candidate in getattr(experimental, "candidates", ()):
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

An absent candidate list does not mean no effect: the returned effect still
applies. Missing sequences and uncalibrated probabilities remain explicit unknowns.

### Current boundaries and legacy names

- **Combined haplotypes are not fully owned by the selected annotator.**
  `VariantCollection.effects()` annotates individual variants, then calls the
  shared haplotype builder when a `phase_resolver` is supplied. Selecting an
  experiment therefore does not control all combined predictions. Relocating
  this active behavior is tracked in [#437](https://github.com/openvax/varcode/issues/437).
  The separate `predict_transcript_model_effect(variants, transcript, ...)`
  function accepts explicitly known-cis variants for the experimental joint
  pipeline; it does not change the collection workflow.
- **The experiment can overstate what a nonlocal SV establishes.**
  [#449](https://github.com/openvax/varcode/issues/449) tracks a DUP spanning CFTR's
  coding sequence that is clipped to a local model and reported as `FivePrimeUTR`.
  This is an open correctness bug, not evidence that the event is harmless.
- **SV consequence/filtering semantics remain incomplete.** A false
  `modifies_protein_sequence` flag is not sufficient evidence of an unchanged
  SV protein; `drop_silent_and_noncoding()` can discard unresolved structural
  effects ([#418](https://github.com/openvax/varcode/issues/418)). Keep those
  results for separate review rather than treating unknown as noncoding.
- **Old names are aliases, not another engine.** `realized`,
  `RealizedEffectAnnotator` and `predict_realized_effect` still select the
  transcript model; new code should use `transcript_model` names. Its
  `RealizedEffectCandidate` product wrapper is still distinct from
  `EffectCandidate`, as shown above; the aliases do not unify those result shapes.

### Writing an annotator

Custom annotators register via the registry:

```python
varcode.register_annotator(my_annotator)
variant.effects(annotator=my_annotator.name)
```

The protocol requires only `name` and `annotate_on_transcript`; `version` is
optional provenance. Return a `MutationEffect` when a prediction is available
or Python's `NotImplemented` singleton when this particular input is outside
the implementation. No `supports` list is required:

```python
class MyExperiment:
    name = "my_experiment"

    def annotate_on_transcript(self, variant, transcript):
        prediction = my_model.predict(variant, transcript)
        if prediction is None:
            return NotImplemented
        return prediction  # a MutationEffect
```

Public APIs turn `NotImplemented` into
`Unresolved(mechanism="unsupported_annotation", reason=...)`, preserving the
selected annotator's provenance. They do not silently substitute the default.
For example, selecting `protein_diff` for an SV returns `Unresolved`, even
inside `use_annotator("protein_diff")`. `None` is a plugin error, and exceptions
retain the usual `raise_on_error` behavior. Unknown never means harmless.

An optional `annotate_with_context(variant, transcript, germline_ctx,
phase_resolver=None)` method has the same return contract. Without that method,
nonempty germline context produces `Unresolved` rather than being ignored or
sent to another implementation. The default retains the existing germline
point-edit path; SV plus germline composition remains experimental in
`transcript_model`. Empty context uses `annotate_on_transcript` as usual.

`variant.effect_on_transcript(transcript, annotator=..., germline=...)` and
`predict_variant_effect_on_transcript` now use the same selection rules as
`effects()`, including the current scoped default.

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
symbolic ALTs (`<DEL>`, `<INS:ME:ALU>`, `<CN0>`), breakends and single
breakends (`.ACGT` / `ACGT.`, loaded as a `BND` with no mate) as
`StructuralVariant` objects rather than dropping them. SVs get the
genome and contig-name settings passed to `load_vcf`, mates included.
Callers such as esvee and GRIDSS write deletions, duplications and
inversions as breakend pairs labeled `SVTYPE=DEL` / `DUP` / `INV`. Each
row loads as the breakend it is; `varcode.transforms.pair_breakends`
joins the pair into one SV of that type, spanning the event, when the
two halves agree on the label and their kept sides fit it.

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

Fusions follow breakend orientation and strand. [Structural variant
annotation](structural_variants.md) covers every case: breakends between
genes and intergenic space, strand combinations for deletions,
duplications and inversions, where a breakpoint lands, and which effect
class comes back.

Limitations:

- Each breakend row produces its own `StructuralVariant`;
  `varcode.transforms.pair_breakends` joins the two rows of a pair (see
  [Transforms](transforms.md)).
- `parse_structural_variants=False` is the default. Without the flag,
  symbolic ALTs are dropped with a warning that names the flag.
- The fusion partner is the first protein-coding transcript at the other
  breakpoint with the right orientation, not a ranked choice
  ([#406](https://github.com/openvax/varcode/issues/406)).

## Downstream consumers

`MutantTranscript` is the prediction-boundary type for downstream
neoantigen pipelines (topiary reads `mt.mutant_protein_sequence`;
vaxrank consumes the `EffectCollection` + protein pair to score
neoantigens). RNA-evidence callers (isovar, Exacto) plug in either as
registered annotators or via the `RNAEvidenceResolver` protocol —
see [Germline-aware annotation](germline.md) for the resolver pattern,
which the same evidence shape uses across germline / phase / RNA.

[i297]: https://github.com/openvax/varcode/issues/297
