# Troubleshooting

Start by checking the input build, reference alleles, and sample names.
Do not suppress an error before understanding which predictions would be lost.

| Problem | First check |
|---|---|
| Reference data missing | Install the release used by `genome=`; see [setup](getting_started.md#reference-data) |
| `ReferenceMismatchError` | Input assembly, REF allele, and forward-strand convention |
| `GenomeBuildMismatchError` | Somatic and germline inputs must use the same assembly |
| `SampleNotFoundError` | Inspect `variants.samples` for the available names |
| Missing SV results | Python: pass `parse_structural_variants=True`. The CLIs already load SVs; inspect skip warnings. See [SV loading](structural_variants.md#basic-usage) |
| No protein sequence | May be unresolved or noncoding, not an exception; see [protein sequences](effect_annotation.md#protein-sequences) |

The exception details below support programmatic handling. The domain-specific
exceptions retain standard `ValueError` or `KeyError` base classes.

## `ReferenceMismatchError`

Raised when a variant's reported `ref` allele doesn't match the
reference transcript sequence at the variant's position:

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

Subclasses `ValueError`, so `except ValueError` keeps working. For
programmatic handling, the structured fields are:

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

On minus-strand transcripts these fields contain reverse-complemented bases,
so they do not directly match the forward-strand REF text in a VCF. The
current error message calls them genome bases; that misleading wording is
tracked in [#434](https://github.com/openvax/varcode/issues/434).

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
produces a `Failure` for the failed annotation instead of raising. Other
transcripts can still annotate successfully:

```python
from varcode.effects import Failure

effects = v.effects(raise_on_error=False)
failures = [e for e in effects if isinstance(e, Failure)]
assert failures
print(failures[0].error)  # retained diagnostic
```

Use this in batch pipelines only if you record and review the `Failure` results;
it does not correct the input or make those rows successfully annotated.
The CLI equivalent is `--skip-errors`; its output retains Failure rows and
their error text. FILTER exclusions are separate: use `--include-filtered`
only when you intend to annotate those records.

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

`ReferenceMismatchError` identifies the variant and transcript. Build and
sample errors instead identify the incompatible references or sample names.
When investigating:

1. Check the genome build. `v.reference_name` (e.g. `"GRCh38"`) should
   match what the VCF was called against.
2. Check the strand. If you expect a reverse-strand gene, the cDNA is
   the reverse complement of the + strand — a common confusion.
3. For sample errors, the exception message lists the available
   samples, so misspellings are easy to spot.
