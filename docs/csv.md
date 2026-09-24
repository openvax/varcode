# Saving and reloading results

Use CSV for a table you can inspect or share. Keep original variants and
evidence as well: a result table does not preserve every alternative, sequence,
or annotation context. Start with writing below; see [CSV vs JSON](#csv-vs-json)
before choosing an archive format.

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
# varcode_version=9.2.5
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

`varcode.vcf_output.variants_to_vcf` writes small variants with a deterministic
sample order. Each sample value is looked up by sample name and FORMAT key;
GT comes first, missing fields are written as `.`, and samples present in only
some records remain in the output. Positions sort numerically within contigs.

```python
from varcode import load_vcf
from varcode.vcf_parsing import VCFHeader
from varcode.vcf_output import variants_to_vcf

variants = load_vcf("input.vcf", genome=81)
header = VCFHeader.from_path("input.vcf")
with open("export.vcf", "w") as out:
    variants_to_vcf(variants, variants.metadata, out=out, header=header)
```

Pass the original `VCFHeader` to retain INFO/FORMAT declarations. Otherwise,
standard FORMAT definitions and observed value types determine new declarations;
original custom cardinality and descriptions cannot be recovered from values
alone. Other source header lines are not archived by this helper.

`load_vcf` retains original ALT lists and indexes. Export reconstructs complete
multi-allelic records in that order, preserving GT and allele-indexed values.
A subset missing an original ALT is rejected before output is written: retain
all alleles or explicitly remap genotype/allele-indexed metadata first. Records
are not merged solely because their IDs agree. Indexed metadata from older versions must be reloaded from its source VCF
to recover the full ALT list. Hand-constructed records without an ALT index
are treated as explicitly biallelic; GT and standard allele-depth/likelihood
fields must agree with that representation.

These sample/FORMAT and ALT-order protections resolve
[#502](https://github.com/openvax/varcode/issues/502) in 10.1.3. Retain original
VCFs as the source of record; the helper does not archive all header metadata,
resolve conflicting source records, or export StructuralVariant objects.

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
# OrderedDict([('varcode_version', '9.2.5'), ('reference_name', 'GRCh38')])
```

Annotator provenance fields (`annotator`, `annotator_version`) use
the same `# key=value` convention.
