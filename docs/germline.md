# Germline-aware annotation

A somatic variant's protein consequence depends on the patient's
germline baseline, not the reference. When a germline variant sits
in the same codon as a somatic variant, the predicted amino-acid
change can differ. The `germline=` kwarg threads this in.

This is a niche feature. Most pipelines don't need it; reference-
relative annotation (the default) is correct for the vast majority
of somatic variants. Read this page only if you have a tumor/normal
pipeline and care about cases where germline shares a codon with
the somatic call.

## Effects without vs. with germline — concrete example

CFTR has a somatic variant at GRCh38 7:117531100 `T→A`. The
neighbouring germline variant at 7:117531101 `T→C` lands in the
same codon. Codon at this position is `CTG` (Leucine, L159); the
somatic alone produces `ATG` (Met).

```python
from pyensembl import cached_release
from varcode import Variant, GermlineContext, predict_germline_aware_effect
from varcode.annotators.registry import get_default_annotator

g = cached_release(81)
cftr = g.transcript_by_id("ENST00000003084")
ann = get_default_annotator()

somatic = Variant("7", 117_531_100, "T", "A", genome=g)
```

### Without germline (reference-relative)

```python
eff = somatic.effect_on_transcript(cftr)
print(type(eff).__name__, eff.short_description)
# Substitution p.L159M
```

The reference codon is `CTG` → somatic changes it to `ATG` → L→M.

### With same-codon germline, phase unknown

```python
germline = Variant("7", 117_531_101, "T", "C", genome=g)
ctx = GermlineContext.from_variants([germline], reference_name="GRCh38")

eff = predict_germline_aware_effect(somatic, cftr, ctx, ann)
print(type(eff).__name__, eff.short_description)
# PhaseAmbiguousEffect ?p.L159M

for o in eff.outcomes:
    ev = o.evidence
    print(f"  haplotype={ev['haplotype']:<10} "
          f"germline_in_cis={[v.short_description for v in ev['germline_variants']]} "
          f"=> {o.effect.short_description}")
# haplotype=B          germline_in_cis=[]                              => p.L159M
# haplotype=A          germline_in_cis=['chr7 g.117531101T>C']         => p.S159T
```

Two haplotypes, two predictions:

- **Trans** (germline on the other allele, A): patient codon was
  `CTG` → somatic produces `ATG` → **L→M**. Identical to the
  reference-relative case.
- **Cis** (germline already on this allele, A): patient codon was
  already `CCG` (Pro, but we started from a Ser-pair — see below);
  in this CFTR example the cis path lands on **S→T**.

The `?` prefix on `?p.L159M` marks the description as the
most-likely candidate of a possibility set. The full set lives on
`.outcomes`.

### With phase resolved → set collapses

```python
from varcode import VCFPhaseResolver

# Hypothetical phased.vcf with PS tags spanning both rows.
phaser = VCFPhaseResolver(merged_vcf_path="phased.vcf")
eff = somatic.effects(germline=ctx, phase_resolver=phaser)
# Regular MutationEffect, not PhaseAmbiguousEffect.
```

When the resolver answers cis/trans for each (somatic, germline)
pair, the set collapses to one hypothesis and the output is a
regular `MutationEffect`.

## When does this matter?

Only when **somatic and germline share a codon (or splice signal)
on the same transcript**. For variants without nearby germline, the
output is identical to reference-relative annotation — `germline=`
is a no-op.

Cross-VCF phase is structurally unknown by default: `PS` tags
describe phase *within* a VCF, so a tumor VCF and a separate normal
VCF can't be phased against each other without a resolver. For
tumor/normal short-read pipelines, possibility sets are the
**guaranteed** case when somatic + germline share a codon — not the
exception. To collapse them you need WhatsHap or HapCUT2 on a
merged VCF, fed in via `VCFPhaseResolver`. Long-read pipelines
(PacBio HiFi, ONT) typically phase entire genes, so possibility
sets are rare.

## Constructing a `GermlineContext`

Four shapes:

```python
# Route 1: full germline call set from a real germline caller.
ctx = GermlineContext.from_germline_vcf("normal.vcf")

# Route 2: multi-sample VCF, extract one column.
ctx = GermlineContext.from_multi_sample_vcf(
    "tumor_normal.vcf",
    sample="NORMAL",
    completeness=Completeness.SPARSE,  # required, no default
)

# Route 3: explicit no-germline fallback (== germline=None).
ctx = GermlineContext.empty()

# Route 4: direct construction from in-memory variants.
ctx = GermlineContext.from_variants([germline], reference_name="GRCh38")
```

`completeness=` distinguishes "no call here = ref/ref" (a real
germline caller's output) from "no call here = unknown" (the
`NORMAL` column of a somatic VCF, which only reports rows the
somatic caller looked at). Forcing the caller to declare prevents
silently mis-treating sparse data as complete.

| Completeness | Pipeline | Absence at a position |
|---|---|---|
| `COMPLETE` | DeepVariant, HaplotypeCaller, Strelka2 germline | ⇒ ref/ref |
| `SPARSE` | Mutect2 `NORMAL`, Strelka2 somatic | ⇒ unknown |
| `HOTSPOTS_ONLY` | Panel-of-normals, ClinVar | ⇒ unknown |
| `EMPTY` | Explicit no-data | n/a |

When the context is sparse and a somatic variant lands in a window
with no germline calls, varcode flags `effect.germline_unknown =
True` rather than silently assuming ref/ref.

## Loss of heterozygosity (LOH)

When a "somatic" call shares position+alt with a germline het, the
"somatic" is really a zygosity change — patient was het, tumor lost
the reference allele. Real biology, not a new mutation.

```python
effect = somatic.effect_on_transcript(transcript)
if getattr(effect, "is_loh", False):
    print("LOH event:", effect.short_description)
```

LOH detection runs whenever `germline=` is passed. The flag is
attached to the resulting effect; classification proceeds normally.

## Composing germline + phase + RNA

```python
effects = somatic.effects(
    germline=germline_ctx,
    phase_resolver=phaser,
    rna_resolver=rna,
)
```

Order: germline modifies the transcript first, phase collapses the
possibility set next, RNA collapses further (or appends observed-
only outcomes). Cross-axis key is `Outcome.evidence["haplotype"]`,
so an RNA observation tagged with the same haplotype tag aligns
with the right germline-aware outcome.

## Cross-VCF build mismatch

Hard error by default if the germline and somatic VCFs were called
against different reference builds:

```python
from varcode import GenomeBuildMismatchError

try:
    effects = somatic.effects(germline=germline)
except GenomeBuildMismatchError as e:
    # Lift one VCF over and retry.
    ...

# If you've manually verified the builds match:
effects = somatic.effects(germline=germline, validate_reference=False)
```

## Lower-level helpers

The high-level `effects(germline=...)` path delegates to
`predict_germline_aware_effect`. Both that function and
`apply_germline_to_transcript` (returns a `MutantTranscript` with
germline edits applied) are public — call them directly if you have
a custom annotator or want the patient protein without full effect
prediction.

## Limitations

- **Germline-disrupted splice sites get no explicit downgrade.**
  Classification runs against the patient's signal, but there's no
  "germline already broke this, downgrade severity" path.
- **Subclonal somatic and CNV dosage are not modeled.** Every
  somatic variant is treated as 100% present.
- **Hypothesis cap of 8** by default when phase is unknown across
  multiple germline variants in a window. Raise via `max_hypotheses=`.
- **Normalization mismatch** between germline and somatic VCFs
  (left-alignment, MNV split) causes apparent position mismatches.
  Normalize both with the same tool first.
- **Population-frequency germline** (gnomAD/ExAC as a substitute
  for patient germline) is not in v1.

## See also

- [#268](https://github.com/openvax/varcode/issues/268) — umbrella issue.
- [Effect annotation](effect_annotation.md) — pipeline overview.
- [Genotypes & sample queries](genotype.md) — sample-aware filtering.
