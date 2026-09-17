# Germline-aware annotation

Use patient germline calls when you want a somatic variant classified against
the patient's baseline rather than only the reference. Nearby germline changes
can alter the predicted amino-acid consequence.

## Basic usage

```python
from varcode import GermlineContext, load_vcf

somatic_variants = load_vcf("tumor.vcf", genome=81)
germline_ctx = GermlineContext.from_germline_vcf("normal.vcf", genome=81)
effects = somatic_variants.effects(germline=germline_ctx)
```

Use matching genome builds. `normal.vcf` should be a germline call set, not a
sparse normal column from a somatic caller; see
[constructing a context](#constructing-a-germlinecontext) for that distinction.

Unknown relative phase can produce several candidate effects. The worked example
below shows how to read them. Continue to [phase evidence](#how-varcode-handles-unknown-phase)
for phased VCF/RNA inputs, or [known deletion haplotypes](#known-deletion-haplotypes-in-rna-alignments)
for the specialized RNA-alignment workflow. No annotator selection is needed.

## Concrete example: same codon, three scenarios

CFTR has a somatic at GRCh38 7:117531100 `T→A`. A neighbouring
germline at 7:117531101 `T→C` lands in the same codon.

```python
from pyensembl import cached_release
from varcode import Variant, VariantCollection, GermlineContext

g = cached_release(81)
cftr = g.transcript_by_id("ENST00000003084")

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
eff = somatic.effect_on_transcript(cftr, germline=ctx)
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
the duck-typed resolver interface, including `mutant_transcript` and
`phased_partners`. The [phasing API](api.md#phasing) documents the built-in
resolvers and source protocols; no public `PhaseResolver` class is required.

```python
class ForceCis:
    source = "demo"
    def in_cis(self, v1, v2, transcript=None): return True

class ForceTrans:
    source = "demo"
    def in_cis(self, v1, v2, transcript=None): return False

eff_cis = somatic.effect_on_transcript(
    cftr, germline=ctx, phase_resolver=ForceCis())
print("cis   ->", type(eff_cis).__name__, eff_cis.short_description)
# cis   -> Substitution p.S159T

eff_trans = somatic.effect_on_transcript(
    cftr, germline=ctx, phase_resolver=ForceTrans())
print("trans ->", type(eff_trans).__name__, eff_trans.short_description)
# trans -> Substitution p.L159M
```

In real code use `VCFPhaseResolver("merged_phased.vcf")` or
`MolecularPhaseResolver(source)` — same `phase_resolver=` slot.

## How varcode handles unknown phase

When somatic and germline share a
codon and the relative phase is unknown, varcode **enumerates the
possibility set** — one classified effect per haplotype hypothesis
— and returns a `PhaseCandidateSet`. To collapse the set, pass a
`phase_resolver=` that knows the answer:

- `VCFPhaseResolver("merged_phased.vcf")` — reads `PS` tags from a
  WhatsHap- or HapCUT2-phased merged VCF.
- `MolecularPhaseResolver(source)` — wraps any RNA-phasing source
  (typically an Isovar adapter shipped by `openvax/isovar`) to check
  which haplotype the somatic was observed on in RNA reads.
  For RNA-seq BAMs without assembly, use
  `MolecularPhaseResolver(RNAReadPhasingSource("tumor.rna.bam"))` to phase
  by direct read/fragment co-occurrence. Raw BAM phasing does not
  provide observed `MutantTranscript`s; assembly-backed sources can.
  `ReadPhaseResolver` remains a compatibility name.
- Anything implementing `in_cis(v1, v2, transcript) -> bool | None`.

Phase known → single `MutationEffect`. Phase unknown →
`PhaseCandidateSet` with `.candidates` for the full set.

## When does this matter?

Nearby germline variants can change the baseline codon or splice signal used
to classify a somatic variant. The current default uses a local window, not a
complete patient-genome reconstruction; see [limitations](#limitations).

Phase-set tags describe phase within their source VCF. Do not equate tags
from independently produced tumor and normal files. Use a resolver backed by a
jointly phased VCF or molecular evidence. Short- and long-read data can both
leave phase unknown; coverage and linked alleles determine whether a particular
pair is resolved.

## Constructing a `GermlineContext`

Choose the source that matches your input:

```python
from varcode import Completeness, GermlineContext

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

When a somatic call matches a germline allele, Varcode attaches an `is_loh`
flag. This is an allele-overlap heuristic, not evidence on its own that a tumor
lost the other allele; interpreting LOH requires additional data. This naming
limitation is tracked in [#454](https://github.com/openvax/varcode/issues/454).

```python
# Reuse the example germline allele as an overlapping tumor call.
overlap = VariantCollection([germline])
for effect in overlap.effects(germline=ctx):
    if getattr(effect, "is_loh", False):
        print("Germline-overlap flag:", effect.short_description)
```

The check requires `germline=`; annotation without that context cannot attach
the flag.

## Composing germline + phase + RNA

Given a somatic `VariantCollection`, patient context, phase resolver, and RNA
resolver:

```python
effects = somatic_variants.effects(
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

Use a `VariantCollection` for collection-level build validation. It raises
`GenomeBuildMismatchError` when its somatic variants and germline context have
different reference builds:

```python
from varcode import GenomeBuildMismatchError, VariantCollection

somatic_variants = VariantCollection([somatic])
try:
    effects = somatic_variants.effects(germline=ctx)
except GenomeBuildMismatchError as error:
    print(error.somatic_reference, error.germline_reference)
    raise  # correct the input builds before retrying
```

The one-variant API does not provide the same collection-level preflight.
`validate_reference=False` is an option on `VariantCollection.effects()`,
not `Variant.effects()`; it suppresses that check and does not convert
coordinates. Prefer correcting mismatched inputs.

## Known deletion haplotypes in RNA alignments

Splice-aware aligners can encode a known deletion's RNA sequence with `N`
rather than `D`, or split the gap among mismatches and smaller deletions.
Opt in to a local **sequence hypothesis** instead of interpreting the gap as
proof of a DNA deletion:

```python
from varcode import MolecularPhaseResolver, RNAReadPhasingSource

source = RNAReadPhasingSource("tumor.rna.bam", min_alt_reads=1)
source.register_haplotype([known_deletion, adjacent_snv_1, adjacent_snv_2])
resolver = MolecularPhaseResolver(source)
effects = variants.effects(phase_resolver=resolver)
```

The variants must share one genome dataset and contig. The genome must provide
reference sequence across the interval, through a FASTA or annotated transcript.
Every retained base between five-base reference flanks must match one alignment
and pass the configured quality/read-edge filters. No mates or separate partial
haplotypes are stitched together. Registration invalidates cached counts.

Registration does not assert cis phase: it supplies a candidate to test. For
registered variants, only full-context matches count as alternate support;
nonmatches remain unknown, not trans evidence. The current local mode supports
nonoverlapping substitutions and deletions. Register competing combinations
separately. Unregistered variants keep the ordinary CIGAR-based behavior.

For the audited GRCh38 MAP2 example, the combination of `209694769 C>A`,
`209694770 T>G`, and the 28-base deletion at `209694773` gives the same local
sequence under `28D`, `28N`, and `22D2M6D` alignments. An unrelated 19,449-base
exon skip lacks the required local anchors and is not supporting evidence.

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
  for patient germline) is not supported.

## See also

- [#268](https://github.com/openvax/varcode/issues/268) — umbrella issue.
- [Effect annotation](effect_annotation.md) — pipeline overview.
- [Genotypes & sample queries](genotype.md) — sample-aware filtering.
