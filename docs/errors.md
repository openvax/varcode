# Error handling

Varcode raises domain-specific exceptions for the two most common ways
analysis can fail on real data: the variant doesn't match the
reference genome, and a sample name doesn't exist in the loaded VCF.
Both subclass standard Python exception types so callers that already
catch `ValueError` or `KeyError` keep working.

## `ReferenceMismatchError`

*New in varcode 2.2.1
([#215](https://github.com/openvax/varcode/issues/215),
[#246](https://github.com/openvax/varcode/issues/246)).*

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

This is the right choice for batch pipelines that would rather log and
skip a bad row than halt.

## `SampleNotFoundError`

*New in varcode 2.3.0
([#267](https://github.com/openvax/varcode/issues/267)).*

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
