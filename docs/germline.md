# Germline-aware annotation

By default, Varcode compares each variant with the reference genome. Pass the
patient's germline calls to classify somatic variants against the patient's
own sequence instead. This matters when an inherited variant sits nearby: if
it already changed a codon, a somatic mutation in that codon can produce a
different amino acid than the reference codon would suggest.

<a id="phase-enumeration-limit"></a>

!!! warning "Phase enumeration limit"
    When the phase-hypothesis cap is exceeded, the current implementation
    classifies one all-cis assignment and marks
    `germline_phase_state="too_many_hypotheses"`. That result does not resolve
    phase and should not be treated as a definitive consequence. The default
    cap is eight hypotheses. A correction that preserves uncertainty is
    tracked in [#503](https://github.com/openvax/varcode/issues/503).

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

When a reported somatic allele matches the supplied germline, Varcode returns
`GermlineAlleleOverlap` with `is_germline_overlap=True`, `is_loh=None`, and
`loh_status="not_assessed"`. It does not apply the inherited allele a second
time or claim a new somatic protein sequence. The allele may still have a
functional effect relative to the genome reference.

LOH requires evidence about loss of an allele previously present in the normal
sample ([NCI definition](https://www.cancer.gov/publications/dictionaries/genetics-dictionary/def/loss-of-heterozygosity)).
Varcode's sequence annotation does not perform that tumor/normal allelic-state
analysis. An allele match, including unchanged heterozygous or inherited
homozygous calls, cannot establish LOH. Integrate a dedicated LOH caller's
evidence separately.

```python
for effect in effects:
    if getattr(effect, "is_germline_overlap", False):
        print("Germline-overlap flag:", effect.short_description)
```

Here `effects` is the collection annotated with germline context above.
The check requires `germline=`; annotation without that context cannot attach
the flag. Use `detect_germline_overlap` for an allele-match query. The deprecated
`detect_loh` now returns `None` (not assessed) and warns, rather than making a
zygosity claim. Mixed inherited/somatic groups in the experimental transcript
model remain explicitly unresolved until their allele-aware baseline can be
represented.

## Limitations

- **Germline-disrupted splice sites get no explicit downgrade.**
  Classification runs against the patient's signal, but there's no
  "germline already broke this, downgrade severity" path.
- **Subclonal somatic and CNV dosage are not modeled.** Every
  somatic variant is treated as 100% present.
- **Hypothesis cap of 8** by default when phase is unknown across
  multiple germline variants in a window. Raise via `max_hypotheses=`.
  Above the cap, the result is not a definitive consequence; see the
  [phase enumeration limit](#phase-enumeration-limit).
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
