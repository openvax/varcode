# Effect annotation

Use `effects()` to predict consequences, then inspect the transcript, protein
sequence, and any alternatives. Start with [Getting started](getting_started.md)
if you still need reference data or a first runnable example.

## Basic usage

```python
import varcode

variants = varcode.load_vcf("variants.vcf", genome=81)  # GRCh38
effects = variants.effects()
for variant, effect in effects.top_priority_effect_per_variant().items():
    print(variant.short_description, effect.short_description)
```

For ordinary use, this is all the selection you need. The default handles
point variants and structural variants through the same interface; you do
not need to choose an annotator. See [structural variant loading](#structural-variants)
if your input includes SV records.

The result is an `EffectCollection` with a prediction for each relevant
`(variant, transcript)` pair. Effects include `Substitution`, `Silent`,
`FrameShift`, structural consequences and unresolved predictions.

## Read an effect

An effect's `short_description` summarizes the change. Its `transcript` identifies
the prediction's transcript when one applies; intergenic effects have none.
`top_priority_effect_per_variant()` gives one summary per variant, while
`top_priority_effect()` selects one effect from the whole collection. Priority
is consequence ordering, not a probability or clinical classification.

Given an effect from the collection, read its predicted protein directly:

```python
protein = effect.mutant_protein_sequence  # may be None
model = effect.mutant_transcript         # may be None even when protein is available
cdna = model.cdna_sequence if model is not None else None
evidence = model.evidence if model is not None else None
```

`None` means unavailable, not unchanged. A transcript model may contain only
partial structure; see [how it composes](#how-it-composes). For a complete
single-variant example, see [reading the protein](getting_started.md#read-the-predicted-protein).

## Reading alternatives

Some predictions are outcome sets rather than a single known consequence:

```python
from varcode import MultiOutcomeEffect

if isinstance(effect, MultiOutcomeEffect):
    for candidate in effect.candidates:
        print(candidate.effect.short_description)
        print(candidate.effect.mutant_protein_sequence)
        print(candidate.source, candidate.evidence)
```

Candidate order follows the producer's rules or evidence. `most_likely_effect`
returns the first candidate; `highest_priority_effect` returns the most severe
by Varcode's ordering. Neither establishes that the outcome occurred. Keep the
candidate wrapper when you need its provenance.

## Choose a deeper topic

- [Splice outcomes](#splice-disrupting-variants): normal splicing, exon skips,
  intron retention, and cryptic sites.
- [Structural variants](structural_variants.md): loading SVs and interpreting
  fusion, partial-sequence, and RNA-supported results.
- [Germline and phasing](germline.md): patient-specific baselines and joint effects.
- [Experimental annotators and integrations](#annotator-selection): supported
  inputs, result-shape differences, and the extension contract.
- [Effect types](effect_types.md) and [API reference](api.md): individual fields
  and classes.

The rest of this page is reference detail; ordinary annotation does not require
an understanding of the internal representations.

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

An exonic variant can affect a codon and a splice signal at the same time.
Not every exonic base is coding: exons also include untranslated regions
([Ensembl glossary](https://www.ensembl.org/Help/Glossary?id=521)).

At the public API, splice disruptions are wrapped in `SpliceOutcomeSet`.
For a coding exonic disruption, the `NormalSplicing` candidate contains the
coding consequence if splicing proceeds unchanged, accessible through
`splice_set.effect_if_splicing_unchanged`. The underlying `ExonicSpliceSite`
signal also records it as `alternate_effect`; that signal is not a separate
default result shape.

For purely intronic disruptions (`SpliceDonor`, `SpliceAcceptor`,
`IntronicSpliceSite`), there is no directly changed coding base.
`effect_if_splicing_unchanged` returns `None`.

For coding variants **outside** the splice window, varcode emits
a plain coding effect (`Substitution`, `Silent`, `FrameShift`,
…) with no splice annotation attached. The variant may still
disrupt splicing through a non-canonical mechanism, but varcode
won't flag it — see Limitations.

### The `SpliceOutcomeSet` shape

Every splice-disrupting variant emits a `SpliceOutcomeSet`; no opt-in flag is
needed. The raw signal classes describe the disrupted site inside the set.

```python
from varcode import Variant

variant = Variant("7", 117_531_114, "G", "T", genome=81)
transcript = variant.genome.transcript_by_id("ENST00000003084")
splice_set = variant.effect_on_transcript(transcript)
print(type(splice_set).__name__)  # SpliceOutcomeSet
print(splice_set.disrupted_signal_class.__name__)  # ExonicSpliceSite
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
from varcode import ExonSkipping, IntronRetention

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
from varcode import MultiOutcomeEffect, SpliceDonor, SpliceOutcomeSet

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

Use this when reviewing the most disruptive predicted alternative, even if it
is not the producer's first choice. It is not evidence of clinical significance.
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
agree. Use `highest_priority_*` to review the most disruptive predicted
alternative, and inspect evidence before interpreting it as an observed outcome.

### Limitations

Sequence-based splice signals are not flagged: exonic splicing
enhancer/silencer disruption mid-exon (~6-10nt SR-protein motifs),
branch points (~20-50nt upstream of the acceptor), deep intronic
cryptic sites. Detecting these needs ML predictors (SpliceAI,
Pangolin, MMSplice, SpliceTransformer) or direct RNA evidence;
tracked in [#297][i297].

## Structural variants

Load symbolic alleles and breakends with `parse_structural_variants=True`.
The ordinary `effects()` interface handles them; there is no separate annotator
to select. Follow the [structural variant guide](structural_variants.md) for
loading, pairing, protein access, and RNA-supported models.

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
raises a warning on load. CSV loading re-annotates with the selected
implementation; selecting the recorded annotator does not restore historical
predictions or evidence. See [serialization limits](csv.md#csv-vs-json).

## The four primitives

| Primitive | What it represents | Module |
|---|---|---|
| `MutationEffect` (and subclasses) | A predicted consequence, unresolved result, or outcome set | `varcode.effects.effect_classes` |
| `MutantTranscript` | Edits or transcript structure with provenance; cDNA and protein sequences are optional | `varcode.mutant_transcript` |
| `MultiOutcomeEffect` | A set of candidate effects in producer-defined order | `varcode.effects.effect_classes` |
| `EffectAnnotator` | An implementation that predicts effects; ordinary callers use the default | `varcode.annotators` |

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

Use the [ordinary sequence accessors](#read-an-effect) for a single effect.

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

## Downstream consumers

Integration authors can use the [annotator contract](annotator_contract.md) and
[RNA evidence API](api.md#rna-evidence). Preserve effect candidates, transcript
identity, sequence completeness, and provenance when passing predictions to
other tools; a sequence string alone does not establish a complete expressed
protein.

[i297]: https://github.com/openvax/varcode/issues/297
