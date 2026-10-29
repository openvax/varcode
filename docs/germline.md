# Germline-aware annotation

By default, Varcode compares each variant with the reference genome. Pass the
patient's germline calls to classify somatic variants against the patient's
own sequence instead. This matters when an inherited variant sits nearby: if
it already changed a codon, a somatic mutation in that codon can produce a
different amino acid than the reference codon would suggest.

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
With many unphased germline variants nearby, the result is unresolved instead;
see the [phase enumeration limit](#phase-enumeration-limit).
No annotator selection is needed.

The default lookup covers the full transcript span, including introns that a
splice hypothesis may retain. `effect.baseline_transcript` and
`effect.mutant_transcript` retain the actual products compared. The baseline
includes inherited alleles assigned in cis by each hypothesis. The same
`germline=`, `phase_resolver=` and `rna_resolver=` inputs compose through the
[shared transcript pipeline](experimental_annotators.md#comparing-annotators).

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

When the context is sparse, Varcode flags `effect.germline_unknown =
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

Supplied germline variants can change baseline codons or splice signals used
to classify a somatic variant. This is a transcript comparison, not a complete
patient-genome reconstruction; see [limitations](#limitations).

Phase-set tags describe phase within their source VCF. Do not equate tags
from independently produced tumor and normal files. Use a resolver backed by a
jointly phased VCF or molecular evidence. Short- and long-read data can both
leave phase unknown; coverage and linked alleles determine whether a particular
pair is resolved.

## Phase enumeration limit

Each germline variant of unknown phase doubles the number of cis/trans
combinations. Homozygous germline variants are always cis, and variants phased
to each other (for example in one VCF phase set) count once, so both reduce the
count.

When the combinations exceed the cap (default eight, so more than three
independent unphased germline variants across the transcript), Varcode does not guess a phase. It
returns a `HypothesisLimit` effect, a kind of `Unresolved`:

- Its `phase` records which germline variants are known cis, known trans, or
  unphased.
- Its `reference_effect` is the effect without germline context. The
  `HypothesisLimit` is ranked like that effect.

To raise the cap, set it on the context:

```python
import dataclasses

ctx = dataclasses.replace(ctx, max_phase_hypotheses=16)
```

## When a patient comparison cannot be built

Each candidate retains its baseline and mutant products. Missing intronic
sequence, conflicting edits and unsupported rearrangements remain explicit
unresolved results rather than reference-relative substitutes. Unknown phase
and splice choices stay in the same `candidate.outcomes` records; inspect
`outcome.hypothesis.phase`, `baseline_splice_plan`, and `splice_plan`.

A germline-disrupted site retains all compatible baseline mechanisms. An
unchanged active site uses the same mechanism on both sides of the comparison;
a site deleted or changed in the tumor does not force one guessed baseline.
RNA support can select a matching predicted product without inventing phase.

When all inherited edits are in trans, the affected haplotype has the reference
baseline. Missing germline information cannot establish that the patient has
that baseline. For unavailable fields, see
[sequence availability](transcript_models.md#sequence-access).

The explicit `annotator="fast"` compatibility path retains local germline
lookup, `PhaseCandidateSet` results, and `Unresolved` with
`mechanism="germline_composition"` for combinations it cannot compose. Its
`reference_effect` is a diagnostic without germline context, not a patient
consequence. Existing serialized results remain readable.

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
the flag.

Use `detect_germline_overlap` for an allele-match query. The deprecated
`detect_loh` now returns `None` (not assessed) and warns, rather than making a
zygosity claim.

The default transcript model composes known-cis groups that mix inherited
and somatic alleles, applying each inherited allele once; see
[mixed inherited and somatic haplotypes](experimental_annotators.md#mixed-inherited-and-somatic-haplotypes).

## Limitations

- **Splice predictions remain hypotheses.** Sequence changes do not establish
  which transcript is expressed or automatically reduce splice severity.
- **Subclonal somatic and CNV dosage are not modeled.** Every
  somatic variant is treated as 100% present.
- **Hypothesis cap of 8** by default when phase is unknown across
  multiple germline variants across the transcript; above it the effect is
  [unresolved](#phase-enumeration-limit). Change it with
  `GermlineContext.max_phase_hypotheses`.
- **Normalization mismatch** between germline and somatic VCFs
  (left-alignment, MNV split) causes apparent position mismatches.
  Normalize both with the same tool first.
- **Population-frequency germline** (gnomAD/ExAC as a substitute
  for patient germline) is not supported.

## Lower-level helpers

The high-level `effects(germline=...)` path uses the shared transcript model.
`predict_transcript_model_effect` accepts an explicit known-cis group and
supplied germline variants. The existing `predict_germline_aware_effect` and
`apply_germline_to_transcript` helpers remain available for compatibility and
custom small-edit annotators; they retain their local/cDNA-edit contracts.

## Related guides

- <a id="concrete-example-same-codon-three-scenarios"></a><a id="scenario-1-no-germline-reference-relative-the-default"></a><a id="scenario-2-germline-passed-phase-unknown-possibility-set"></a><a id="scenario-3-force-phasing-single-effect"></a>[Two variants in one codon](phasing.md#two-variants-in-one-codon).
- <a id="how-varcode-handles-unknown-phase"></a>[Phasing from VCF or RNA](phasing.md).
- <a id="composing-germline-phase-rna"></a>[Combining germline, phase, and RNA evidence](phasing.md#combining-evidence).
- <a id="known-deletion-haplotypes-in-rna-alignments"></a>[Known deletion haplotypes in RNA](phasing.md#known-deletion-haplotypes-in-rna-alignments).
- <a id="see-also"></a>[Genotypes and sample queries](genotype.md), [effect annotation](effect_annotation.md), and the [germline API](api_phasing.md#germline-aware-annotation).
- [Germline development roadmap (#268)](https://github.com/openvax/varcode/issues/268).
