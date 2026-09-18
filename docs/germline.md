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

Unknown relative phase can produce several candidate effects. See
[two variants in one codon](phasing.md#two-variants-in-one-codon) for a worked
example, or [phased VCF input](phasing.md#phased-vcf) to supply evidence.
No annotator selection is needed.

## Constructing a `GermlineContext`

Choose the source that matches your input:

```python
from varcode import Completeness, GermlineContext

# Route 1: full germline call set from a real germline caller.
ctx = GermlineContext.from_germline_vcf("normal.vcf", genome=81)

# Route 2: multi-sample VCF, extract one column.
ctx = GermlineContext.from_multi_sample_vcf(
    "tumor_normal.vcf",
    sample="NORMAL", genome=81,
    completeness=Completeness.SPARSE,  # required, no default
)

# Route 3: explicit no-germline fallback (== germline=None).
ctx = GermlineContext.empty()

# Route 4: direct construction from an existing collection of germline variants.
ctx = GermlineContext.from_variants(germline_variants, reference_name="GRCh38")
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

## Cross-VCF build mismatch

Use a `VariantCollection` for collection-level build validation. It raises
`GenomeBuildMismatchError` when its somatic variants and germline context have
different reference builds:

```python
from varcode import GenomeBuildMismatchError

try:
    effects = somatic_variants.effects(germline=germline_ctx)
except GenomeBuildMismatchError as error:
    print(error.somatic_reference, error.germline_reference)
    raise  # correct the input builds before retrying
```

The one-variant API does not provide the same collection-level preflight.
`validate_reference=False` is an option on `VariantCollection.effects()`,
not `Variant.effects()`; it suppresses that check and does not convert
coordinates. Prefer correcting mismatched inputs.

<a id="when-does-this-matter"></a>

## How germline changes predictions

Nearby germline variants can change the baseline codon or splice signal used
to classify a somatic variant. The current default uses a local window, not a
complete patient-genome reconstruction; see [limitations](#limitations).

Phase-set tags describe phase within their source VCF. Do not equate tags
from independently produced tumor and normal files. Use a resolver backed by a
jointly phased VCF or molecular evidence. Short- and long-read data can both
leave phase unknown; coverage and linked alleles determine whether a particular
pair is resolved.

## Loss of heterozygosity (LOH)

When a somatic call matches a germline allele, Varcode attaches an `is_loh`
flag. This is an allele-overlap heuristic, not evidence on its own that a tumor
lost the other allele; interpreting LOH requires additional data. This naming
limitation is tracked in [#454](https://github.com/openvax/varcode/issues/454).

```python
for effect in effects:
    if getattr(effect, "is_loh", False):
        print("Germline-overlap flag:", effect.short_description)
```

Here `effects` is the collection annotated with germline context above.
The check requires `germline=`; annotation without that context cannot attach
the flag.

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

## Lower-level helpers

The high-level `effects(germline=...)` path delegates to
`predict_germline_aware_effect`. Both that function and
`apply_germline_to_transcript` (returns a `MutantTranscript` with
germline edits applied) are public — call them directly if you have
a custom annotator or want the patient protein without full effect
prediction.

## Related guides

- <a id="concrete-example-same-codon-three-scenarios"></a><a id="scenario-1-no-germline-reference-relative-the-default"></a><a id="scenario-2-germline-passed-phase-unknown-possibility-set"></a><a id="scenario-3-force-phasing-single-effect"></a>[Two variants in one codon](phasing.md#two-variants-in-one-codon).
- <a id="how-varcode-handles-unknown-phase"></a>[Phasing from VCF or RNA](phasing.md).
- <a id="composing-germline-phase-rna"></a>[Combining germline, phase, and RNA evidence](phasing.md#combining-evidence).
- <a id="known-deletion-haplotypes-in-rna-alignments"></a>[Known deletion haplotypes in RNA](phasing.md#known-deletion-haplotypes-in-rna-alignments).
- <a id="see-also"></a>[Genotypes and sample queries](genotype.md), [effect annotation](effect_annotation.md), and the [germline API](api_phasing.md#germline-aware-annotation).
- [Germline development roadmap (#268)](https://github.com/openvax/varcode/issues/268).
