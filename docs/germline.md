# Germline-aware annotation

*New in varcode 4.19.0 ([#268](https://github.com/openvax/varcode/issues/268)).*

When you analyze somatic variants in a real patient, the biologically
correct "before" state is the **patient's germline** — not the
reference genome. The same somatic mutation in the same codon can
produce a different protein change depending on the patient's
germline at neighboring positions. varcode's `germline=` kwarg
threads patient germline through every step of effect prediction so
the output reflects the patient's actual baseline.

## Why this matters: a worked example

Reference codon at some position is `GCC` (Alanine).

Patient is heterozygous for a germline variant `V_g` changing the
middle base `G→T`:

- Haplotype A: `GCC` (Ala) — reference.
- Haplotype B: `GTC` (Val) — germline.

A somatic variant `V_s` changes the first base `G→A`. The somatic
effect depends on which haplotype `V_s` landed on:

| Haplotype | Codon before V_s | After V_s | Amino acid change |
|---|---|---|---|
| In trans with V_g (V_s on A) | `GCC` (Ala) | `ACC` (Thr) | **Ala→Thr** |
| In cis with V_g (V_s on B) | `GTC` (Val) | `ATC` (Ile) | **Val→Ile** |

Reference-relative annotation reports **Ala→Thr** for both cases —
correct for trans, wrong for cis. Without phase data, the honest
output is the **possibility set** `{Val→Ile, Ala→Thr}`. With phase
data, the set collapses to one. varcode produces both shapes
automatically.

## How it composes

Germline is a transcript modifier, applied before the somatic variant
is classified. Effect prediction with `germline=ctx` builds the
patient haplotype from any germline variants in the somatic's window,
then asks the annotator how the somatic changes the patient protein.
When phase between the somatic and one or more germline variants in
the same codon is unknown, the result is a `PhaseAmbiguousEffect` —
one classified `MutationEffect` per phase hypothesis, packaged
together. A `phase_resolver=` collapses the set; an `rna_resolver=`
collapses it further (or adds observed-isoform outcomes).

## Quickstart

```python
from varcode import load_vcf, GermlineContext

somatic = load_vcf("tumor.vcf", genome="GRCh38")
germline = GermlineContext.from_germline_vcf("normal.vcf")

effects = somatic.effects(germline=germline)

for effect in effects:
    print(effect.short_description)
    # Possibility-set effects expose .outcomes; iterate to see all
    # phase hypotheses with per-haplotype evidence.
    for outcome in getattr(effect, "outcomes", ()):
        hap = outcome.evidence.get("haplotype")
        print(f"  [{hap}] {outcome.short_description}")
```

`germline=None` (the default) is byte-identical to today's
reference-relative behavior — pinned by the test suite. The new
path engages only when you explicitly opt in.

## Input routes: how do I construct a `GermlineContext`?

In real cancer pipelines, "the germline VCF" can be one of several
shapes with different semantics for *absence* of a call. The
`GermlineContext` constructors mirror those shapes; see
`Completeness` for the absence-handling rules.

### Route 1: full germline call set (most common)

```python
ctx = GermlineContext.from_germline_vcf("normal.vcf")
```

Use this for output from a real germline caller run on the normal
BAM (DeepVariant, GATK HaplotypeCaller, Strelka2 germline). Defaults
to `Completeness.COMPLETE` — absence at a position implies ref/ref.

### Route 2: multi-sample VCF, extract a column

```python
from varcode import GermlineContext, Completeness

ctx = GermlineContext.from_multi_sample_vcf(
    "tumor_normal.vcf",
    sample="NORMAL",
    completeness=Completeness.SPARSE,  # required, no default
)
```

`completeness=` is a required keyword. Multi-sample VCFs from
somatic callers (Mutect2's `NORMAL` column, Strelka2 somatic) are
**sparse** — they only emit rows where the somatic caller saw
tumor-vs-normal signal. Pure-germline multi-sample VCFs (1000
Genomes batch) are complete. Forcing the caller to declare
prevents silently mis-treating a sparse column as if absence
implied ref/ref.

### Route 3: tumor-only with no germline

```python
ctx = GermlineContext.empty()  # explicit fallback
effects = somatic.effects(germline=ctx)
# Equivalent to germline=None: today's reference-relative output.
```

`empty()` is an explicit no-germline opt-in. Useful when your
pipeline always passes a `germline=` kwarg and you want to
document "we have no germline for this patient" rather than
silently omitting the kwarg.

Population-frequency germline (using gnomAD/ExAC as a probabilistic
substitute) is not in v1; tracked as a future direction in #268.

### Route 4: direct construction (tests, custom pipelines)

```python
from varcode import Variant

ctx = GermlineContext.from_variants(
    [Variant("7", 117_531_101, "T", "C", genome="GRCh38")],
    reference_name="GRCh38",
)
```

Useful when you already have variants in memory and don't need to
re-parse a VCF.

## Possibility sets: what does the output look like?

When a somatic variant and one or more germline variants share a
codon (or splice signal) and phase between them is unknown, the
honest output is a **possibility set** — one effect per phase
hypothesis. varcode emits a `PhaseAmbiguousEffect`:

```python
from varcode.effects.effect_classes import PhaseAmbiguousEffect

effects = somatic.effects(germline=germline)

for effect in effects:
    if isinstance(effect, PhaseAmbiguousEffect):
        # Most-likely candidate (highest priority class) is .most_likely.
        print("primary:", effect.most_likely.short_description)
        # Full set of hypotheses on .outcomes.
        for outcome in effect.outcomes:
            ev = outcome.evidence
            print(
                f"  haplotype={ev['haplotype']} "
                f"phase={ev['phase_state']} "
                f"germline_in_cis={[g.short_description for g in ev['germline_variants']]} "
                f"effect={outcome.effect.short_description}")
```

Each `Outcome.evidence` carries:

- `phase_state`: `"phased"`, `"implicit"` (hemizygous), `"unknown"`
  (enumerated), or `"too_many_hypotheses"` (cap exceeded).
- `haplotype`: an opaque tag (`"A"`, `"B"`, `"A_mixed_<n>"`, or
  `"unknown"`) used for cross-axis matching with RNA evidence.
- `germline_variants`: the cis germline variants on this
  hypothesis's haplotype (empty tuple for the all-trans case).

`PhaseAmbiguousEffect` is a `MultiOutcomeEffect` subclass — the
same machinery the splice-outcomes (#262) and SV (#264) work uses.
Consumers iterating `.outcomes` for splice variants already iterate
the same shape for germline-aware variants.

`effect.short_description` for an ambiguous effect is prefixed with
`?` (e.g. `?p.L159M`) to signal "this is the most-likely candidate
of a possibility set."

## Cross-VCF phase is structurally unknown by default

`PS` tags in a VCF describe phase **within** that VCF. If germline
is in one VCF and somatic is in another, the cross-VCF phase
relationship is unknown even when both VCFs have rich `PS` tags
internally — the somatic caller didn't see the germline calls when
it emitted its rows.

In practice, for a tumor/normal short-read pipeline, possibility
sets are the **guaranteed** case when somatic + germline share a
codon, not just the common case. To collapse a possibility set,
you need explicit phase data spanning both: typically running
WhatsHap or HapCUT2 on the merged VCF and feeding the resulting
phase into a `PhaseResolver`.

```python
from varcode import VCFPhaseResolver

phaser = VCFPhaseResolver(merged_vcf_path="phased.vcf")
effects = somatic.effects(germline=germline, phase_resolver=phaser)
```

When the resolver answers cis/trans for each (somatic, germline)
pair, the possibility set collapses to a single hypothesis and
the output is a regular `MutationEffect` (not a
`PhaseAmbiguousEffect`).

Long-read data (PacBio HiFi, ONT) typically resolves phase across
whole genes — so possibility sets become rare. Assembly-based
pipelines (hifiasm, Verkko, Exacto) provide full haplotype
sequences, which flow in via `MutantTranscript` directly.

## Loss of heterozygosity (LOH) detection

When a "somatic" variant call shares position+alt with a germline
het call, the apparent somatic is really a zygosity change:
patient was germline het, tumor lost the reference allele, now
appears alt-only in tumor. Real biology, not a true new mutation.

```python
effect = somatic.effect_on_transcript(transcript)
if getattr(effect, "is_loh", False):
    # This 'somatic' call is LOH at a germline het site.
    print("LOH event:", effect.short_description)
```

LOH detection runs whenever you pass `germline=`. The flag is
attached to the resulting effect; the underlying classification
proceeds normally (so you can still see the predicted protein
change at this site).

## Sparseness: how absence-of-a-call is interpreted

| Completeness | Typical pipeline | Absence at a position |
|---|---|---|
| `COMPLETE` | Germline caller (DeepVariant, HaplotypeCaller, Strelka2 germline) on normal BAM | ⇒ ref/ref |
| `SPARSE` | `NORMAL` column of a somatic VCF (Mutect2, Strelka2 somatic) | ⇒ unknown — somatic caller didn't query this position |
| `HOTSPOTS_ONLY` | Panel-of-normals filter list, ClinVar pathogenic list | ⇒ definitely unknown |
| `EMPTY` | Explicit fallback ("no germline data") | n/a |

When the context is `SPARSE` or `HOTSPOTS_ONLY` and a somatic
variant lands in a window with no germline calls, varcode flags
the resulting effect with `effect.germline_unknown = True` rather
than silently treating absence as ref/ref. Consumers can filter or
surface the uncertainty.

```python
ctx = GermlineContext.from_multi_sample_vcf(
    "tumor_normal.vcf", sample="NORMAL",
    completeness=Completeness.SPARSE,
)
for effect in somatic.effects(germline=ctx):
    if getattr(effect, "germline_unknown", False):
        # Patient germline state at this site is unknown.
        # Surface this in your downstream report.
        ...
```

## Composing germline + phase + RNA

The three axes layer through `Outcome.evidence`. Each resolver adds
information without overriding earlier layers:

```python
from varcode import (
    load_vcf,
    GermlineContext,
    VCFPhaseResolver,
)

somatic = load_vcf("tumor.vcf")
germline = GermlineContext.from_germline_vcf("normal.vcf")
phaser = VCFPhaseResolver(merged_vcf_path="phased_merged.vcf")
rna = MyIsovarResolver(reads=tumor_rna_reads)  # your implementation

effects = somatic.effects(
    germline=germline,
    phase_resolver=phaser,
    rna_resolver=rna,
)
```

Order matters:

1. **Germline first** — modifies the patient's transcript.
2. **Phase next** — collapses possibility sets when cis/trans is
   known.
3. **RNA last** — collapses sets to observed isoforms, or appends
   observed-only outcomes via the resolver protocol from #259.

The cross-axis key is `Outcome.evidence["haplotype"]`. An RNA
observation tagged `evidence={"haplotype": "B"}` aligns with the
haplotype-B germline-aware outcome. Consumers filtering by
`haplotype` get a self-consistent picture across all three axes.

## Cross-VCF build mismatch

Hard error by default when the germline and somatic VCFs were
called against different reference builds. Coordinates from the
two cannot be composed without lift-over.

```python
from varcode import GenomeBuildMismatchError

try:
    effects = somatic.effects(germline=germline)
except GenomeBuildMismatchError as e:
    print(f"Mismatch: somatic={e.somatic_reference} germline={e.germline_reference}")
    # Lift one VCF over and retry.
```

If you've explicitly lifted over and know the builds agree, opt
out:

```python
effects = somatic.effects(germline=germline, validate_reference=False)
```

## Lower-level helpers

The high-level `effects(germline=...)` path delegates to
`predict_germline_aware_effect`. Both that function and
`apply_germline_to_transcript` (which returns a `MutantTranscript`
with germline edits applied) are public — call them directly if you
have a custom annotator or want to translate the patient protein
without going through full effect prediction. See `varcode.germline`.

## Limitations

- **Germline-disrupted splice sites get no explicit downgrade.** When
  germline already breaks a splice signal that the somatic also
  targets, the patient transcript carries the germline edit, so
  classification runs against the patient signal — but there's no
  separate "germline already broke this; downgrade severity" path.
- **Subclonal somatic and CNV dosage are not modeled.** Every somatic
  variant is treated as present in 100% of cells; copy-number changes
  don't affect per-haplotype protein prediction.
- **Hypothesis cap of 8** by default when phase is unknown across
  multiple germline variants in a window. Raise via `max_hypotheses=`
  if your pipeline tolerates more.
- **Normalization mismatch between germline and somatic VCFs** (left-
  alignment, MNV split) causes apparent position mismatches.
  Normalize both with the same tool before loading.
- **Population-frequency germline (gnomAD/ExAC) as a substitute for
  patient germline** is not in v1.

## See also

- [#268](https://github.com/openvax/varcode/issues/268) — umbrella issue with the v1 acceptance criteria.
- [Effect annotation](effect_annotation.md) — pipeline overview.
- [Genotypes & sample queries](genotype.md) — sample-aware filtering.
