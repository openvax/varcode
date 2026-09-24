# Troubleshooting

Most errors come from a mismatch between your input and the reference data:
the wrong genome build, a REF allele that doesn't match the reference, or a
misspelled sample name. Understand an error before suppressing it, because
suppressing it loses those predictions.

| Problem | First check |
|---|---|
| Reference data missing | Install the release used by `genome=`; see [setup](getting_started.md#reference-data) |
| `ReferenceMismatchError` | Input genome build, REF allele, and forward-strand convention; see [below](#referencemismatcherror) |
| `GenomeBuildMismatchError` | Somatic and germline inputs must use the same genome build |
| `SampleNotFoundError` | Inspect `variants.samples` for the available names |
| Missing SV results | Python: pass `parse_structural_variants=True`. The CLIs already load SVs; inspect skip warnings. See [SV loading](structural_variants.md#basic-usage) |
| No protein sequence | Not an error: the effect may be unresolved or noncoding; see [protein sequences](effect_annotation.md#protein-sequences) |
| A few bad records stop a batch run | See [continuing past errors](#continuing-past-errors) |

The domain-specific exceptions below subclass the standard `ValueError` or
`KeyError`, so existing `except` clauses keep working.

## `ReferenceMismatchError`

Raised when a variant's `ref` allele doesn't match the reference transcript
sequence at the variant's position:

```python
import varcode

v = varcode.Variant("7", 117531114, "T", "A", genome=81)
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

<a id="three-common-causes"></a>

### Common causes

1. **Wrong genome build.** A VCF called against GRCh37 but annotated with
   GRCh38 (or vice versa) produces these errors wherever the builds differ.
   Check that `v.reference_name` (e.g. `"GRCh38"`) matches the build your
   variants were called against.
2. **Germline allele in the `ref` field.** VCF requires `ref` to match the
   reference genome. Encode a patient's germline variant at that position as
   a separate variant rather than changing `ref`.
3. **Strand confusion.** Varcode expects alleles on the reference genome's
   forward (+) strand. For a gene on the reverse strand, the cDNA is the
   reverse complement of what you should enter.

### Handling it in code

The exception subclasses `ValueError` and has structured fields:

```python
try:
    v.effects()
except varcode.ReferenceMismatchError as e:
    e.variant           # the Variant
    e.transcript        # the Transcript being compared against
    e.expected_ref      # reference bases in transcript orientation
    e.observed_ref      # variant REF bases in transcript orientation
    e.transcript_offset # position in the transcript
```

On minus-strand transcripts, `expected_ref` and `observed_ref` are
reverse-complemented, so they won't match the forward-strand REF text in your
VCF. The current error message calls them genome bases; that misleading
wording is tracked in [#434](https://github.com/openvax/varcode/issues/434).

## `GenomeBuildMismatchError`

Raised by `VariantCollection.effects(germline=...)` when the somatic
variants and the germline context use different genome builds (e.g. GRCh37
vs GRCh38). Unlike `ReferenceMismatchError`, which is per variant, this check
runs once for the whole collection, so a build mismatch fails immediately
instead of producing one error per variant.

```python
try:
    effects = somatic.effects(germline=germline)
except varcode.GenomeBuildMismatchError as e:
    e.somatic_reference   # the somatic collection's reference (e.g. 'GRCh38')
    e.germline_reference  # the germline context's reference (e.g. 'GRCh37')
```

It subclasses `ValueError`. `validate_reference=False` on
`VariantCollection.effects()` skips this check, but it does not convert
coordinates or verify that the data match. Prefer fixing the input builds.

## `SampleNotFoundError`

Raised by `VariantCollection` genotype and filter methods when a sample name
isn't in the collection's source VCFs. The message lists the available names:

```python
vc = varcode.load_vcf("tumor_normal.vcf", genome="GRCh38")
vc.samples
# ['normal', 'tumor']

vc.for_sample("patient_01")
# SampleNotFoundError: Sample 'patient_01' not found.
# Available samples: ['normal', 'tumor']
```

It subclasses `KeyError` for backward compatibility. Failing on a typo is
intentional: silently returning an empty collection for a misspelled sample
name would hide bugs in analysis scripts.

<a id="the-raise_on_errorfalse-escape-hatch"></a>
<a id="debugging-tips"></a>

## Continuing past errors

In Python, pass `raise_on_error=False` to `.effects()`. Each failed
annotation becomes a `Failure` effect instead of an exception, and other
transcripts can still annotate successfully:

```python
from varcode.effects import Failure

effects = v.effects(raise_on_error=False)
failures = [e for e in effects if isinstance(e, Failure)]
assert failures
print(failures[0].error)  # retained diagnostic
```

On the command line, `--skip-errors` does the same:

- Effects CSVs contain `Failure` rows and an `annotation_error` column.
- Gene CSVs add `annotation_status` and `annotation_error` columns, with the
  gene fields left empty for failures.
- Failures stay visible even with `--only-coding` and `--one-per-variant`.
  The latter selects one successful effect per variant and also keeps each
  failed transcript.
- Input parsing errors still stop the command.

Either way, the failed records are not corrected or annotated; record and
review them. Filtered records are a separate matter: the CLI skips non-`PASS`
records unless you pass `--include-filtered`, which you should use only when
you intend to annotate them.
