# VariantCollection transforms

A **transform** is a pure function: given a `VariantCollection` and
optional auxiliary inputs (reference, phase resolver, ...), it returns
a new `VariantCollection`. Cardinality may be preserved, reduced, or
increased. Composition is by application — no registry, no protocol.

The pattern lives in `varcode.transforms`. As of 4.20.0 it ships one
transform, `pair_breakends`; the module is designed to grow.

## The contract

Every transform owes three things, documented in its docstring:

| Field | Meaning |
|---|---|
| **Cardinality** | `preserves`, `reduces`, or `increases`. |
| **Provenance** | Every output variant carries `source_variants: tuple[Variant, ...]`. Empty tuple for pass-through; one element for derived-from-one; two or more for combined. Not part of hash/equality. |
| **Metadata behavior** | Explicit rule for how `source_to_metadata_dict` entries flow through (which fields are inherited from which source, which require agreement, what happens on disagreement). |

Transforms are **idempotent on inputs they don't recognize**. Running
`pair_breakends` twice produces the same VC; the second pass finds no
unpaired BNDs to combine because every combined row's `source_variants`
is already populated.

Composition is just function application:

```python
import varcode
from varcode.transforms import pair_breakends

vc = varcode.load_vcf("tumor.vcf", genome="GRCh38",
                     parse_structural_variants=True)
vc = pair_breakends(vc)
effects = vc.effects()
```

## `pair_breakends`

Merges MATEID-paired BND rows into a single combined
`StructuralVariant`. **Reduces.**

A VCF can represent the same translocation event two ways:

| Caller pattern | Input rows | After `pair_breakends` |
|---|---|---|
| Manta / DELLY / SVABA / newer GRIDSS | Two BND rows linked by `MATEID` | One combined row, `source_variants=(a, b)` |
| Older GRIDSS | Two BND rows linked by `PARID` | One combined row (PARID treated as MATEID alias) |
| BreakDancer / CREST / older DELLY | One row, `SVTYPE=TRA`, `CHR2`/`END` in INFO | Pass-through, `source_variants=()` |
| GRIDSS unresolved single-end | One BND row, no `MATEID` | Pass-through |

The point is **caller-uniformity after the transform**: regardless of
which caller produced the VCF, effect prediction sees one variant per
rearrangement event.

### Usage

```python
from varcode import load_vcf
from varcode.transforms import pair_breakends

vc = load_vcf("manta.vcf", genome="GRCh38",
              parse_structural_variants=True)
vc = pair_breakends(vc)

# Each combined variant carries both endpoints + provenance.
for v in vc:
    if v.sv_type == "BND" and v.source_variants:
        bnd_a, bnd_b = v.source_variants
        print(f"{v.contig}:{v.start} <-> {v.mate_contig}:{v.mate_start} "
              f"from rows {bnd_a.info.get('paired_with')} + "
              f"{v.info.get('paired_with')}")

effects = vc.effects()
```

### Pairing rules

1. **Primary key**: `MATEID` on the variant's `info` dict matched
   against VCF row IDs captured at load time.
2. **Alias**: `PARID` (older GRIDSS) is treated as `MATEID`.
3. **Symmetric**: A.mateid must equal B.id **and** B.mateid must equal
   A.id. Asymmetric references log a warning and pass through.
4. **In-degree 1**: each ID must be referenced by exactly one other
   row's MATEID. If any ID has incoming degree > 1, the whole
   connected component is left unpaired with a warning.
5. **Mate present**: if a MATEID points to an ID not in this
   collection (filtered out, chunked load), the half passes through
   with a warning.

### Metadata merge

The combined variant's per-source metadata entry is built fresh:

| Field | Rule |
|---|---|
| `id` | A's ID (lex-earlier of the two). |
| `qual` | `min(A.qual, B.qual)` when both present. |
| `filter` | Union of FILTER tokens; `PASS` drops out if any non-PASS label appears. |
| `info` | A's values; `MATEID`/`PARID` removed (no longer meaningful); `paired_with` added pointing at B's ID. |
| `sample_info` | Per-sample: GT must match across A and B (raises on disagreement); other FORMAT fields taken from A. |
| `alt_assembly` | One populated -> use it. Both equal -> use shared. Both differ -> A's wins with warning. |

The strict-GT rule is intentional: both halves of a legitimate paired
BND describe the same biological event, so disagreement indicates
either a caller bug (asymmetric filtering, separate re-genotyping per
half) or a real analytical concern. The transform raises with both
row IDs and both GT values so the problem surfaces.

### Trade-off: single fusion direction post-collapse

A reciprocal translocation produces two derivative chromosomes
(`der(15)t(15;19)` and `der(19)t(15;19)` for BRD4-NUTM1). Pre-pair,
varcode emits one `GeneFusion` effect per half × overlapping
transcript — so both fusion directions are represented. Post-pair,
the combined variant is anchored at the lex-earlier endpoint, so
effects represent that single direction.

The other direction is reachable: `combined.source_variants` returns
both originals, and you can annotate the other half directly:

```python
combined = next(v for v in vc if v.source_variants)
other_direction = combined.source_variants[1].effects()
```

If you want both directions in the same effect collection without
running `pair_breakends`, just don't run it — the parser already
emits both halves.

## `left_align_indels`

Shifts indels to their canonical leftmost equivalent position.
**Preserves cardinality** (1:1, value may change).

Two pipelines that exchange variants need a canonical representation
per indel. `CTT→T` at position 10 inside a `CT`-repeat means the
same biological event as `CT→_` at position 8 — but tools that
compare variants by `(contig, start, ref, alt)` see those
representations as distinct. Left-alignment normalizes to the
leftmost equivalent position.

### Usage

```python
from varcode import load_vcf, Genome
from varcode.transforms import left_align_indels

# Default (transcript-cDNA coverage only) — exonic indels normalize,
# intronic/intergenic pass through unchanged.
vc = load_vcf("tumor.vcf", genome="GRCh38")
vc = left_align_indels(vc)

# Full coverage — every indel normalizes.
g = Genome(81, fasta="/path/to/GRCh38.fa")
vc = load_vcf("tumor.vcf", genome=g)
vc = left_align_indels(vc)
```

No `reference` parameter — `left_align_indels` reads bases via the
genome the variants already carry (see
[varcode.Genome](api.md#varcodegenome)). Coverage depends on which
genome shape was passed.

### Behavior

| Variant kind | Behavior |
|---|---|
| Pure SNV / MNV / complex (`ATG→GCC`) | Pass-through, `source_variants=()` |
| Already-canonical indel (no equivalent leftward position) | Pass-through, `source_variants=()` |
| Indel inside an exon, no FASTA | Shifts via transcript cDNA |
| Indel in homopolymer / STR repeat | Shifts to leftmost equivalent position |
| Indel in intron / intergenic, no FASTA | Pass-through (no reference coverage) |
| Indel in intron / intergenic, FASTA attached | Shifts via FASTA |
| Indel that starts exonic but would shift past the exon boundary, no FASTA | **Partial shift** — moves left until the boundary, sets `info["left_align_partial"] = True` |
| Indel at chromosome start | Shift bounded by `start > 1` |

### Metadata

Shifted variants carry:

- `source_variants = (original,)` — the pre-shift variant.
- Source-path metadata re-keyed under the shifted variant; the
  original key is removed.
- `info["original_start"]` — the pre-shift start position, for
  round-trip diagnostics.
- `info["left_align_partial"] = True` — set only when the shift was
  bounded by reference coverage (not by the reference disagreeing).
  Tells downstream tools that a chromosome FASTA might have shifted
  the variant further.

### Idempotence

Running `left_align_indels` twice produces the same result on the
second call — every variant is already at its canonical leftmost
position after the first call. Composes cleanly with `pair_breakends`
in either order; SVs and BNDs are not indels and pass through.

## Roadmap

Transforms planned for future PRs. Each lands as its own ticket; all
follow the contract above.

| Transform | Cardinality | Brief |
|---|---|---|
| `combine_cis_snvs(vc, phase_resolver)` | reduces | Adjacent in-codon SNVs sharing a phase set merge into MNVs. |

See the [API reference](api.md) for
`varcode.transforms.pair_breakends` and
`varcode.transforms.left_align_indels`.
