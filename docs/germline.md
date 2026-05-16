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

## How varcode handles unknown phase

Varcode does no phasing itself. When somatic and germline share a
codon and the relative phase is unknown, varcode **enumerates the
possibility set** — one classified effect per haplotype hypothesis
— and returns a `PhaseCandidateSet`. To collapse the set, pass a
`phase_resolver=` that knows the answer:

- `VCFPhaseResolver(merged_phased.vcf)` — reads `PS` tags from a
  WhatsHap- or HapCUT2-phased merged VCF.
- `MolecularPhaseResolver(source)` — wraps any RNA-phasing source
  (typically an Isovar adapter shipped by `openvax/isovar`) to check
  which haplotype the somatic was observed on in RNA reads.
  For RNA-seq BAMs without assembly, use
  `MolecularPhaseResolver(RNAReadPhasingSource("tumor.rna.bam"))` to phase
  by direct read/fragment co-occurrence. Raw BAM phasing does not
  provide observed `MutantTranscript`s; assembly-backed sources can.
  `ReadPhaseResolver` remains as the varcode 5.0 compatibility name.
- Anything implementing `in_cis(v1, v2, transcript) -> bool | None`.

Phase known → single `MutationEffect`. Phase unknown →
`PhaseCandidateSet` with `.candidates` for the full set.

## Concrete example: same codon, three scenarios

CFTR has a somatic at GRCh38 7:117531100 `T→A`. A neighbouring
germline at 7:117531101 `T→C` lands in the same codon.

```python
from pyensembl import cached_release
from varcode import Variant, GermlineContext, predict_germline_aware_effect
from varcode.annotators.registry import get_default_annotator

g = cached_release(81)
cftr = g.transcript_by_id("ENST00000003084")
ann = get_default_annotator()

somatic = Variant("7", 117_531_100, "T", "A", genome=g)
germline = Variant("7", 117_531_101, "T", "C", genome=g)
ctx = GermlineContext.from_variants([germline], reference_name="GRCh38")
```

### Scenario 1 — no germline (reference-relative, the default)

```python
eff = somatic.effect_on_transcript(cftr)
print(type(eff).__name__, eff.short_description)
# Substitution p.L159M
```

### Scenario 2 — germline passed, phase unknown → possibility set

```python
eff = predict_germline_aware_effect(somatic, cftr, ctx, ann)
print(type(eff).__name__, eff.short_description)
# PhaseCandidateSet ?p.L159M

for c in eff.candidates:
    ev = c.evidence
    print(f"  haplotype={ev['haplotype']:<2} "
          f"germline_in_cis={[v.short_description for v in ev['germline_variants']]} "
          f"=> {c.effect.short_description}")
# haplotype=B  germline_in_cis=[]                            => p.L159M
# haplotype=A  germline_in_cis=['chr7 g.117531101T>C']       => p.S159T
```

The `?` prefix on `?p.L159M` flags the description as the
most-likely candidate of a `PhaseCandidateSet`. Real consumers
should `isinstance(eff, PhaseCandidateSet)` and iterate
`eff.candidates` (a tuple of `EffectCandidate` objects carrying
per-hypothesis evidence keys); use `eff.effects` if only the
inner classified `MutationEffect`s are needed.

### Scenario 3 — force phasing → single effect

A real pipeline gets the cis/trans answer from a phased VCF or an
RNA assembly. For this demo, a hand-rolled resolver makes the
collapse visible:

These stubs only implement `in_cis(...)`, which is all the
codon-collapse path consults. Richer pipelines implement more of
the `varcode.phasing.PhaseResolver` protocol (`has_contig`,
`mutant_transcript`, `phased_partners`, ...).

```python
class ForceCis:
    source = "demo"
    def in_cis(self, v1, v2, transcript=None): return True

class ForceTrans:
    source = "demo"
    def in_cis(self, v1, v2, transcript=None): return False

eff_cis = predict_germline_aware_effect(somatic, cftr, ctx, ann,
                                         phase_resolver=ForceCis())
print("cis   ->", type(eff_cis).__name__, eff_cis.short_description)
# cis   -> Substitution p.S159T

eff_trans = predict_germline_aware_effect(somatic, cftr, ctx, ann,
                                           phase_resolver=ForceTrans())
print("trans ->", type(eff_trans).__name__, eff_trans.short_description)
# trans -> Substitution p.L159M
```

In real code use `VCFPhaseResolver(merged_phased.vcf)` or
`MolecularPhaseResolver(source)` — same `phase_resolver=` slot.

## When does this matter?

Only when **somatic and germline share a codon (or splice signal)
on the same transcript**. For variants without nearby germline, the
output is identical to reference-relative annotation — `germline=`
is a no-op.

Cross-VCF phase is structurally unknown by default: `PS` tags
describe phase *within* a VCF, so a tumor VCF and a separate normal
VCF can't be phased against each other without a resolver. For
tumor/normal short-read pipelines, `PhaseCandidateSet` is the
**guaranteed** output when somatic + germline share a codon — not
the exception. To collapse to a single effect you need WhatsHap or
HapCUT2 on a merged VCF, fed in via `VCFPhaseResolver`. Long-read
pipelines (PacBio HiFi, ONT) typically phase entire genes, so
candidate sets are rare there.

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
candidate set next, and RNA refines multi-candidate effects. Splice
mechanism sets are reconciled to observed mechanisms; other
multi-outcome effects append observed-only candidates. Cross-axis key is
`EffectCandidate.evidence["haplotype"]`, so an
RNA observation tagged with the same haplotype tag aligns with the
right germline-aware outcome.

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
