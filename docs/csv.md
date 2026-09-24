# Saving and reloading results

Save results as CSV for a table you can inspect or share, as
[JSON](#csv-vs-json) to archive supported objects, or write variants back to
[VCF](#vcf-export). Whichever you choose, keep the original inputs and
evidence too: a results table does not preserve every alternative outcome,
sequence, or annotation context.

## Writing

```python
from varcode import load_vcf

vc = load_vcf("variants.vcf", genome="GRCh38")
vc.to_csv("variants.csv")

effects = vc.effects()
effects.to_csv("effects.csv")
```

By default, both writers prepend `#`-prefixed provenance lines so the
file is self-describing:

```
# varcode_version=10.2.1
# reference_name=GRCh38
chr,start,ref,alt,gene_name,gene_id
17,43082575,C,T,BRCA1,ENSG00000012048
...
```

Pass `include_header=False` for fast consumers that can't tolerate
comment lines:

```python
vc.to_csv("plain.csv", include_header=False)
```

## Reading

The examples here apply to ordinary point-variant tables. Structural-variant
CSVs are summaries only: `from_csv` rejects them rather than reconstructing
incomplete SVs. See [structural exports](sv_reference.md#alleles-coordinates-and-exports).

When the header is present, `from_csv` recovers the reference genome
automatically:

```python
from varcode import VariantCollection, EffectCollection

vc = VariantCollection.from_csv("variants.csv")
effects = EffectCollection.from_csv("effects.csv")
```

`EffectCollection.from_csv` re-runs annotation. It does not recover the original
germline/phase/RNA evidence, candidate sets, or historical predictions. A
reference-name header also does not pin the original annotation dataset; pass
the same explicit `genome` when that matters.

When the CSV has no header (written by an older varcode or with
`include_header=False`), pass `genome` explicitly:

```python
vc = VariantCollection.from_csv("plain.csv", genome="GRCh38")
```

Missing both the header *and* an explicit `genome` produces a clear error:

```
ValueError: from_csv needs a reference genome: pass the `genome`
argument explicitly, or write the CSV with
`to_csv(include_header=True)` so `# reference_name=...` is recorded
in the header. Neither was found at plain.csv.
```

## Column-name flexibility

`VariantCollection` historically writes a `chr` column while
`EffectCollection` writes `contig`. Both readers accept either spelling for
that column. This does not make their full schemas interchangeable:
effect tables also need effect/transcript fields.

## Version drift

Because `EffectCollection.from_csv` re-runs annotation on read, a
collection serialized by one major varcode version and loaded under
another can produce different effects. `from_csv` emits
a `UserWarning` when the header's `varcode_version` differs in major
version from the currently-installed version. Minor and patch drift
is silent, but bug fixes or a different annotation dataset can still change
predictions; API compatibility is not a guarantee of identical scientific output.

## CSV vs JSON

| | CSV | JSON |
|---|---|---|
| Intended use | Inspect or share a table | Serialize supported objects |
| Reloading effects | Re-annotates on read | Restores serialized state where supported |
| Structural variants | Export summaries; import rejected | Restore effect/candidate graphs and mutant transcript models |
| Carries annotator version header | Yes | (via the serialized object) |
| Preserves all effect-specific fields | No | Depends on the effect type and serialization support |

JSON (`to_json` / `from_json`) avoids CSV's re-annotation for supported objects.
Structural results use a versioned effect graph: self candidates and shared
effect links survive, as do primary, cryptic, splice, and external candidates,
their source/evidence fields, fusion partners, affected exons, and attached
mutant transcript sequences/segments. Collection annotator provenance is also
preserved. Deserialization restores recorded predictions without recomputing
them. Older Varcode versions cannot read this graph format.

This remains an object archive with Python class and reference-dataset
dependencies. Arbitrary extra attributes and third-party sequence providers
need their own serialization support. Test the round-trip for the actual result
types you use, and retain source variants, reference release, and original
RNA/phase/germline evidence independently.

<a id="vcf-export-limitations"></a>

## VCF export

`varcode.vcf_output.variants_to_vcf` writes small variants back to VCF:

```python
from varcode import load_vcf
from varcode.vcf_parsing import VCFHeader
from varcode.vcf_output import variants_to_vcf

variants = load_vcf("input.vcf", genome=81)
header = VCFHeader.from_path("input.vcf")
with open("export.vcf", "w") as out:
    variants_to_vcf(variants, variants.metadata, out=out, header=header)
```

What is preserved:

- **Sample identities.** Each value is looked up by sample name and FORMAT key,
  so values can't shift between sample columns. Sample order is deterministic,
  GT comes first, missing values are written as `.`, and samples present in
  only some records are kept.
- **Multi-allelic records.** `load_vcf` remembers each record's original ALT
  list and indexes, so export rebuilds the complete record in its original
  order, with GT and allele-indexed values intact.
- **INFO/FORMAT declarations,** when you pass the original `VCFHeader` as
  above. Without it, standard FORMAT definitions and the observed value types
  determine new declarations; custom descriptions and cardinalities can't be
  recovered from values alone.
- **Ordering.** Positions sort numerically within each contig.

What is refused or not kept:

- A subset that drops one of a record's original ALT alleles is rejected
  before any output is written. Keep all alleles, or explicitly remap the
  genotype and allele-indexed fields first.
- Records are not merged just because their IDs match.
- Hand-constructed records without an ALT index are treated as biallelic;
  their GT and standard allele-depth/likelihood fields must agree with that.
- Metadata saved by older Varcode versions lacks the full ALT list; reload it
  from the source VCF.
- Header lines other than INFO/FORMAT declarations are not archived,
  conflicting source records are not resolved, and `StructuralVariant`
  objects are not exported.

Keep the original VCF as the source of record.

## Annotation provenance

Annotated effect collections record `annotator`, `annotator_version`, and
`annotated_at` (an ISO-8601 UTC timestamp). Filtering and grouping preserve
these fields. CSV export writes them into the metadata header, and loading
recovers the original values.

A mismatch between the recorded annotator and the current default raises a
warning on load. Selecting the recorded implementation does not restore the
original predictions: CSV loading still re-annotates. Keep the original reference
release and evidence separately; see [CSV vs JSON](#csv-vs-json).

## Custom header fields

`read_metadata_header` is available in `varcode.csv_helpers` for
tools that want to add their own metadata lines:

```python
from varcode.csv_helpers import read_metadata_header, write_metadata_header

meta = read_metadata_header("variants.csv")
# OrderedDict([('varcode_version', '10.2.1'), ('reference_name', 'GRCh38')])
```

Annotator provenance fields (`annotator`, `annotator_version`) use
the same `# key=value` convention.
