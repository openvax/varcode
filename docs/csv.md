# CSV round-trip and metadata headers

*New in varcode 2.1.0, refined in 2.2.0.*

`VariantCollection` and `EffectCollection` can both round-trip through
CSV. The CSV format is human-readable and easy to inspect in a
spreadsheet; the JSON format (`to_json` / `from_json`) is the
byte-for-byte exact round-trip and is considerably faster for large
collections (it skips re-annotation on read).

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
# varcode_version=2.3.0
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

When the header is present, `from_csv` recovers the reference genome
automatically:

```python
from varcode import VariantCollection, EffectCollection

vc = VariantCollection.from_csv("variants.csv")
effects = EffectCollection.from_csv("effects.csv")
```

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
`EffectCollection` writes `contig`. As of 2.2.0 both readers accept
either alias, so CSVs from the two sides are interchangeable. Writers
are unchanged.

## Version drift

Because `EffectCollection.from_csv` re-runs annotation on read, a
collection serialized by one major varcode version and loaded under
another can produce different effects. As of 2.2.0, `from_csv` emits
a `UserWarning` when the header's `varcode_version` differs in major
version from the currently-installed version. Minor and patch drift
is silent (semver guarantees compatibility).

## CSV vs JSON

| | CSV | JSON |
|---|---|---|
| Human-readable | Yes | No |
| Byte-for-byte round-trip | **No** (re-annotates on read) | Yes |
| Speed on large collections | Slow (per-row Variant construction + annotation) | Fast |
| Carries annotator version header | Yes | (via the serialized object) |
| Preserves all effect-specific fields | No | Yes |

Use CSV when you want to inspect or edit the output manually. Use
JSON (`to_json` / `from_json`) for exact round-trip and for large
collections.

## Custom header fields

`read_metadata_header` is available in `varcode.csv_helpers` for
tools that want to add their own metadata lines:

```python
from varcode.csv_helpers import read_metadata_header, write_metadata_header

meta = read_metadata_header("variants.csv")
# OrderedDict([('varcode_version', '2.3.0'), ('reference_name', 'GRCh38')])
```

Annotator provenance fields (`annotator`, `annotator_version`) use
the same `# key=value` convention (tracked in
[#271](https://github.com/openvax/varcode/issues/271)).
