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

Local deletion and tandem-duplication models are translated and classified as
`StartLoss`, `Deletion`, `Insertion`, `FrameShift`, or another protein consequence.
The annotated start is mapped through the retained transcript segments; a
remaining downstream ATG does not rescue a deleted start. UTR-only edits retain
`FivePrimeUTR` / `ThreePrimeUTR` labels when their protein is unchanged.

A single resolved model returns its consequence directly. When a breakpoint cuts
an exon, or motif/splice alternatives exist, a `StructuralVariantEffect` contains
the consequences in `.candidates`. Its description and priority reflect those
consequences. Candidate evidence records `sv_type`, `splice_ambiguous`, and the
`reference_splicing` assumption; DUP models also record `tandem_duplication`.
`mutant_transcript` retains the modeled cDNA, translated protein, and evidence.
See [alternative outcomes](effect_annotation.md#alternative-outcomes).

`GeneFusion` and `TranslocationToIntergenic` remain distinct consequences.
Unmaterialized inversions, unspecified insertions/CNVs, and assemblies without a
mapped CDS return `Unresolved` candidates. A deletion on an incompletely annotated
coding transcript falls back to `ExonLoss`.

**Migrating to 10.0:** `LargeDeletion`, `LargeDuplication`, and `Inversion` remain
importable for legacy objects but are no longer emitted. Read the DNA event from
`effect.variant.sv_type`; use the consequence or its candidates for protein impact.
Not every SV effect has `.candidates`: check `isinstance(effect, MultiOutcomeEffect)`.

For CFTR (`ENST00000003084`, Ensembl 81), deleting exons 1–3 loses the annotated
start, deleting exon 5 yields an in-frame `Deletion` (1480 → 1450 aa), and deleting
exons 5–6 yields `FrameShift` (171 aa). Reading-frame classification uses the
spliced edit, consistent with [Ensembl's consequence definitions](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html).

A DNA rearrangement does not by itself establish a complete expressed fusion
protein. Sequence may be unknown or partial; `None` is not an unchanged protein.

## Filtering by protein change

```python
retained = effects.drop_silent_and_noncoding()
# Require a positive prediction of protein change, excluding unknowns:
resolved_changes = effects.drop_silent_and_noncoding(keep_unresolved=False)
```

SV effects report `modifies_protein_sequence` and `modifies_coding_sequence` as
`True` (changed), `False` (unchanged), or `None` (unresolved). A changed CDS can
still encode the same protein. Comparisons use the transcript being annotated;
for a fusion, the 5′ partner supplies the initiation site.

The flags cover the whole candidate set: any changed alternative makes the flag
`True`; otherwise an unresolved alternative makes it `None`. The filter retains
the original set, including provenance, rather than selecting one candidate.
Unresolved effects are kept by default. To test for a *known* unchanged SV, use
`effect.modifies_protein_sequence is False`, not `not effect.modifies_protein_sequence`.

These flags describe existing predictions; they do not establish expression.
Partial BND fragments and assemblies without a mapped ORF stay unresolved.
So do partial observations (`protein_completeness` other than `start_to_stop`)
unless a mapped, in-frame observed codon differs from the reference.
The filter does not construct missing proteins or change the effect class.

Selenocysteine (`U` in the Ensembl reference protein) is encoded by UGA,
which is decoded as Sec only with a SECIS element in the mRNA's 3′ UTR.
Ensembl doesn't annotate SECIS positions. A fusion's stored protein reads Sec
unless no selenoprotein 3′ UTR remains ([transcript models](transcript_models.md#selenocysteine)).
For the flags, UGA is read as Sec where a model keeps
the transcript intact from that codon through its 3′ end, and as a stop where
no selenoprotein 3′ UTR remains. Otherwise the protein flag is reported only if
both readings agree, else `None`; for example, a partial 3′ UTR deletion leaves
the protein unresolved. The coding flag compares CDS bases, which don't depend
on how Sec is decoded.
A start codon other than ATG counts as the initiator methionine, even though
Ensembl writes CTG and TTG starts as `L`.

## Fusion protein candidates

For an annotated effect, inspect every compatible partner isoform rather than
only the protein on the first result:

```python
from varcode.effects import GeneFusion

for candidate in effect.candidates:
    fusion = candidate.effect
    if isinstance(fusion, GeneFusion):
        print(fusion.five_prime_transcript.id, fusion.three_prime_transcript.id)
        print(fusion.mutant_protein_sequence)  # None when unresolved
```

There is no candidate-count cap. Distinct transcript pairs remain separate even
when their predicted proteins match. The first candidate follows annotation
order, not measured likelihood. These are reference-isoform predictions, not
every possible splice/phase combination or evidence that a fusion is expressed.
See [fusion rules](sv_reference.md#fusion-partners) and [RNA imports](rna_structures.md).

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
| `StartLoss` | The annotated initiation codon is not retained. |
| `Deletion` / `Insertion` | An in-frame protein deletion / insertion. |
| `FrameShift` | A spliced coding edit changes the reading frame. |
| `StructuralVariantEffect` | A set of conditional consequences or unresolved models. |
| `Unresolved` | Available sequence/coordinates cannot establish the consequence. |
| `ExonLoss` | Deleted exons on a transcript whose CDS annotation is incomplete. |
| `Intronic` | A span inside the transcript that overlaps no exon. |
| `Intergenic` | No gene at the variant's position. |
| `NoncodingTranscript` | The transcript isn't protein-coding. |

## Limitations

- Partner isoforms are enumerated but not ranked by RNA support or likelihood.
- Local DEL/DUP models with inserted or unreadable paired-breakend alleles
  remain unresolved until their junction bases can be incorporated
  ([#491](https://github.com/openvax/varcode/issues/491)). Fusion models already
  retain resolved exonic junction inserts.
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

For exhaustive callset reconciliation, see [Compare SV calls](sv_comparison.md).
