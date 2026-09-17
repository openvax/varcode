# Troubleshooting

Start by checking the input build, reference alleles, and sample names.
Do not suppress an error before understanding which predictions would be lost.

| Problem | First check |
|---|---|
| Reference data missing | Install the release used by `genome=`; see [setup](getting_started.md#reference-data) |
| `ReferenceMismatchError` | Input assembly, REF allele, and forward-strand convention |
| `GenomeBuildMismatchError` | Somatic and germline inputs must use the same assembly |
| `SampleNotFoundError` | Inspect `variants.samples` for the available names |
| Missing SV results | Load with `parse_structural_variants=True`; see [SV loading](structural_variants.md#basic-usage) |
| No protein sequence | May be unresolved or noncoding, not an exception; see [result access](effect_annotation.md#read-an-effect) |

The exception details below support programmatic handling. The domain-specific
exceptions retain standard `ValueError` or `KeyError` base classes.

## `ReferenceMismatchError`

Raised when a variant's reported `ref` allele doesn't match the
reference genome at the variant's position:

```python
import varcode

v = varcode.Variant("7", 117531114, "T", "A", "GRCh38")
# The real + strand ref at chr7:117531114 is G, not T.
v.effects()
```

Produces:

```
varcode.errors.ReferenceMismatchError:
Reference allele mismatch for Variant(contig='7', start=117531114, ref='T', alt='A', reference_name='GRCh38')
on Transcript(...) at transcript offset 620 (chromosome positions 117531114:117531114):
variant reports ref='T' but the reference genome has 'G' at this position.
This usually means the variant was called against a different genome
build, the ref field was filled in with the patient's germline allele
rather than the reference, or the variant is on the wrong strand.
Pass raise_on_error=False to .effects() to receive a Failure effect
instead of raising.
```

Subclasses `ValueError`, so `except ValueError` keeps working. For
programmatic handling, the structured fields are:

```python
try:
    v.effects()
except varcode.ReferenceMismatchError as e:
    e.variant           # the Variant
    e.transcript        # the Transcript being compared against
    e.expected_ref      # what the genome has
    e.observed_ref      # what the variant claims
    e.transcript_offset # position in the transcript
```

### Three common causes

1. **Wrong genome build.** A VCF called against GRCh37 annotated with
   GRCh38 (or vice versa) produces these errors at positions where
   the builds differ.
2. **Germline allele in the `ref` field.** VCF requires `ref` to match
   the reference genome. Patient-specific germline variants at the
   position should be encoded as separate variants, not by changing
   `ref`.
3. **Strand confusion.** The variant is specified on the negative
   strand but varcode expects positive-strand coordinates.

### The `raise_on_error=False` escape hatch

If you'd rather continue past these errors instead of surfacing them,
pass `raise_on_error=False` to `.effects()`. Each mismatched variant
produces a `Failure` effect instead of raising:

```python
from varcode.effects import Failure

effects = v.effects(raise_on_error=False)
assert any(isinstance(e, Failure) for e in effects)
```

Use this in batch pipelines only if you record and review the `Failure` results;
it does not correct the input or make those rows successfully annotated.

## `GenomeBuildMismatchError`

Raised by `VariantCollection.effects(germline=...)` when the somatic
collection and the germline context were called against different
reference genome builds (e.g. GRCh37 vs GRCh38). Distinct from
`ReferenceMismatchError`, which is per-variant: this is a top-level
pre-flight that fails fast on the whole pair so you don't see N
per-variant errors caused by a build mismatch.

```python
try:
    effects = somatic.effects(germline=germline)
except varcode.GenomeBuildMismatchError as e:
    e.somatic_reference   # the somatic collection's reference (e.g. 'GRCh38')
    e.germline_reference  # the germline context's reference (e.g. 'GRCh37')
```

Subclasses `ValueError`. `validate_reference=False` suppresses this check on
`VariantCollection.effects()`; it does not lift over coordinates or verify
that the data match. Prefer correcting the input reference metadata and builds.

## `SampleNotFoundError`

Raised by `VariantCollection` genotype/filter methods when the sample
name isn't in the collection's source VCF(s):

```python
vc = varcode.load_vcf("tumor_normal.vcf", genome="GRCh38")
vc.samples
# ['normal', 'tumor']

vc.for_sample("patient_01")
# SampleNotFoundError: Sample 'patient_01' not found.
# Available samples: ['normal', 'tumor']
```

Subclasses `KeyError` for back-compat. The early-fail-on-typo behavior
is intentional: silently returning an empty collection when a sample
name is misspelled would hide real bugs in analysis scripts.

## Debugging tips

Both errors include the specific variant and the transcript that
triggered them. When investigating:

1. Check the genome build. `v.reference_name` (e.g. `"GRCh38"`) should
   match what the VCF was called against.
2. Check the strand. If you expect a reverse-strand gene, the cDNA is
   the reverse complement of the + strand — a common confusion.
3. For sample errors, the exception message lists the available
   samples, so misspellings are easy to spot.
