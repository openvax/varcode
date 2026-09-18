# Structural variant annotation

Load structural variants (SVs) with an explicit VCF option, then use the same
`effects()` interface as for small variants.

## Basic usage

```python
from varcode import load_vcf

variants = load_vcf(
    "structural_variants.vcf",
    genome=81,  # GRCh38; use the annotation chosen for your input
    parse_structural_variants=True,
)
effects = variants.effects()
for effect in effects:
    if effect.variant.is_structural:
        print(effect.variant.sv_type, effect.short_description)
        print(effect.mutant_protein_sequence)  # may be None
```

Without `parse_structural_variants=True`, symbolic alleles and breakends are
skipped with a warning. No separate annotator selection is needed.
[Reference setup](getting_started.md#reference-data) explains the genome argument.

## Reading SV results

Effects include `LargeDeletion`, `LargeDuplication`, `Inversion`, `GeneFusion`,
and `TranslocationToIntergenic`. They may contain alternatives in `.candidates`;
see [alternative outcomes](effect_annotation.md#alternative-outcomes).

A DNA rearrangement does not by itself establish a complete expressed fusion
protein. Sequence may be unknown or partial; `None` is not an unchanged protein.
Do not rely on `drop_silent_and_noncoding()` to retain unresolved SV effects:
their protein-change flags remain incomplete
([#418](https://github.com/openvax/varcode/issues/418)).

## Which transcripts get annotated

`effects()` produces one effect per overlapping transcript, as for any
variant. What "overlapping" means depends on the SV:

| Variant | Transcripts annotated |
|---|---|
| Breakend record (`BND`) | Those at its own breakpoint. The mate's transcripts are annotated from the mate's record. |
| `DEL`, `DUP`, `INV`, `INS`, `CNV` | Every transcript overlapping `start`..`end`. |

Two outcomes apply before any SV logic:

- No gene at the variant's position(s) → a single `Intergenic` effect.
- A non-coding transcript → `NoncodingTranscript`.

## Pairing breakend records

`pair_breakends` joins the two records of a breakend pair:

- If both carry the same `SVTYPE` of `DEL`, `DUP` or `INV` and their kept
  sides fit it, the result is that typed event, annotated over its whole
  span. See [pairing rules and examples](transforms.md#pair_breakends).
- Otherwise the result is one `BND` anchored at the record whose VCF ID
  sorts first, and only that end's transcripts are annotated. If that end
  is intergenic the combined variant reports just `Intergenic`; the other
  end is still reachable through `combined.source_variants`.

## Effect classes

| Class | Meaning |
|---|---|
| `GeneFusion` | The SV joins this transcript sense-to-sense with a coding transcript in another gene. |
| `TranslocationToIntergenic` | A breakend in this transcript that doesn't form a gene fusion. |
| `LargeDeletion` | A deletion (or `<CN0>`) removing one or more exons. |
| `LargeDuplication` | A duplication, insertion or CNV overlapping exons. |
| `Inversion` | An inversion overlapping exons. |
| `Intronic` | A span inside the transcript that overlaps no exon. |
| `Intergenic` | No gene at the variant's position. |
| `NoncodingTranscript` | The transcript isn't protein-coding. |

## Limitations

- The partner isoform isn't ranked
  ([#406](https://github.com/openvax/varcode/issues/406)).
- Chains of several SVs aren't assembled into one allele, and regulatory
  effects (promoter or enhancer hijacking) aren't modeled.
- Annotating multi-megabase spans can be slow
  ([#407](https://github.com/openvax/varcode/issues/407)).

<a id="contents"></a>

## Related guides

- <a id="alleles-coordinates-and-exports"></a>[SV coordinates and exports](sv_reference.md#alleles-coordinates-and-exports).
- <a id="junctions-which-side-of-each-breakpoint-is-kept"></a><a id="5-and-3-roles"></a>[Junction orientation and fusion partners](sv_reference.md#junction-orientation).
- <a id="breakends-every-combination"></a><a id="deletions-duplications-and-inversions"></a><a id="where-in-the-gene-a-breakpoint-lands"></a>[Breakend, span, and protein outcomes](sv_reference.md#breakend-outcomes).
- <a id="importing-observed-rna-structures"></a>[Import observed RNA structures](rna_structures.md).
- <a id="import-exactos-protein-predictions"></a>[Import Exacto protein predictions](rna_structures.md#import-exactos-protein-predictions).
