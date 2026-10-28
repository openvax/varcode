# API conventions and return values

Use `Variant.effects()` for all overlapping transcripts,
`Variant.effect_on_transcript()` for one selected transcript, and
`VariantCollection.effects()` for a batch. They share annotator selection
rules. The default is `fast`; experimental implementations do not silently
fall back to it for unsupported inputs.

## Return shapes

| Operation | Result | Empty or unavailable case |
|---|---|---|
| `Variant.effects()` or `VariantCollection.effects()` | EffectCollection | May be empty after selection; handled errors remain Failure effects |
| `Variant.effect_on_transcript()` | MutationEffect | Unsupported annotation becomes Unresolved; annotation errors propagate |
| `EffectCollection.top_priority_effect()` | MutationEffect | Raises ValueError on an empty collection |
| `predict_orf_effects(...)` | Tuple of ORFAnnotation objects | Empty tuple if no catalogue model is selected |
| `annotate_orf_variants(...)` | ORFAnnotation | Inspect unresolved reason and candidate outcomes |
| `annotate_orf_rna_product(...)` | ORFProductAnnotation | Retains mapped observations and unresolved comparisons |
| `annotate_orf_product(...)` | Dictionary | Compatibility JSON export of the observed-product annotation |
| `resolve_contig_name(...)` | ContigResolution | Inspect status; unresolved results have `resolved_name=None` |

An effect can itself contain multiple candidate outcomes. Selecting one
highest-priority effect does not choose the expressed outcome within it.
An unavailable protein sequence (`None`) is different from an empty sequence
or an unchanged protein. See [reading effects](effect_annotation.md) and
[ORF results](api_orfs.md).

Check whether a collection is nonempty before asking for its top effect:

```python
from varcode import EffectCollection

effects = EffectCollection([])
summary = effects.top_priority_effect() if effects else None
assert summary is None
```

## Error defaults

| Entry point | Default | Explicit alternatives |
|---|---|---|
| `Variant.effects()` | `raise_on_error=True`: annotation errors raise | False retains per-transcript failures |
| `VariantCollection.effects()` | `raise_on_error=None`: retain missing-contig failures and raise other annotation errors | True is strict; False retains annotation failures |
| `Variant.effect_on_transcript()` | Errors propagate | Use the collection-returning workflow when failures must be collected |

These policies apply to annotation; invalid configuration and file parsing can
still raise before annotation begins. Missing reference data does not establish
an intergenic or reference-identical result. See [troubleshooting](errors.md).

## Names, coordinates, and compatibility

Use documented root imports for common operations. Topic modules expose
additional specialized APIs. Names beginning with an underscore are internal;
their signatures are not extension contracts.

`genome=` is the canonical reference argument. The older `ensembl=` argument is
retained for compatibility. `load_vcf_fast` is deprecated in favor of
`load_vcf`. The old `realized` annotator names refer to the experimental
transcript-model implementation. The [uORF wrappers](api_orfs.md#upstream-orf-compatibility-api)
retain their upstream-only defaults, so changing names alone can change which
catalogue is used.

Variant positions and reference lookup helpers use one-based genomic
coordinates. cDNA/protein offsets used for Python indexing are zero-based;
formatted protein descriptions use one-based residue numbers. Each API's
parameter documentation specifies its interval convention. Chromosome-name
resolution changes spelling within a supplied vocabulary; it does not convert
coordinates between assemblies.

## Selecting an annotator

Pass `annotator=` explicitly when different calls need different implementations.
`set_default_annotator()` changes the process-wide default. `use_annotator()`
temporarily changes that same global state and restores it on exit, including
when an exception occurs. It is not isolated between threads or asynchronous
tasks. Its `with ... as ...` value is None.

An extension returns a MutationEffect or Python's `NotImplemented` singleton.
Returning None is an error. Public prediction methods convert an explicit
refusal into Unresolved. See the [annotator contract](annotator_contract.md).

## Type checking

Varcode ships a `py.typed` marker and adopts annotations gradually. File-loading
inputs, chromosome resolution, core annotation return paths, base effect
sequence/change fields, and the
annotator protocol/registry are checked with static consumer examples.
`register_annotator()` preserves the concrete type of the registered object.

Coverage is not yet complete. Inherited collection iteration/indexing, many
effect subclass properties and optional integration arguments can still expose Any.
An accepted type-check is not proof that every step in a pipeline was checked.
Runtime-checkable protocols check attribute presence, not complete signatures
or biological validity.

Contributors can run `./lint.sh` to check implementation and consumer contracts.
The [maintenance guide](dev/api-quality.md) records remaining gaps and how to
extend coverage without masking them with Any.
