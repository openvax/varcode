# Effect annotation

Call `effects()` to predict the consequences of your variants. If you have not
installed reference data yet, start with [Getting started](getting_started.md).

<a id="basic-usage"></a>

## Annotate variants

```python
import varcode

variants = varcode.load_vcf("variants.vcf", genome=81)  # GRCh38
effects = variants.effects()
for variant, effect in effects.top_priority_effect_per_variant().items():
    print(variant.short_description, effect.short_description)
```

The collection holds one prediction per variant per overlapping transcript,
and a variant's effect can differ between transcripts. The loop above prints
one summary per variant, picking the most severe effect by Varcode's
consequence ordering, not by likelihood or clinical significance.
`effects.top_priority_effect()` instead selects one effect from the entire
collection. See [how to read results](getting_started.md#how-to-read-results).

The same interface handles structural variants. To include them when loading
a VCF, pass `parse_structural_variants=True`; see [SV loading](structural_variants.md#basic-usage).

## Choose transcripts

By default, `effects()` includes every overlapping transcript. Use
`transcript_filter` to select transcripts before prediction:

```python
# Restrict annotation to a known set of reference transcript IDs.
selected_ids = {"ENST00000003084", "ENST00000357654"}
effects = variants.effects(transcript_filter=lambda t: t.id in selected_ids)

# Or select by the reference annotation's biotype.
effects = variants.effects(
    transcript_filter=lambda t: t.biotype == "protein_coding")
```

The same argument works on one `Variant` and on
`varcode.effects.predict_variant_effects()`. It receives the original reference
transcript. Use a deterministic predicate; it may be called again for another
variant or ORF. The predicate chooses the primary transcripts sent to the
annotator; an SV fusion prediction can still refer to a partner transcript
discovered by that annotator.

When supplying `extra_orfs`, selection also applies to each ORF's original
transcript placement. Excluded placements do not generate coding predictions
or joint haplotype effects. Unplaced ORF failures remain visible because there
is no resolved transcript to select. See [additional ORFs](uorf_annotations.md).

Selection does not change how a transcript is interpreted. With the built-in
annotators, ordinary point-variant prediction still reports
`NoncodingTranscript` outside the `protein_coding` biotype, and
`IncompleteTranscript` for incomplete coding annotations. A stop codon that
polyadenylation completes does not make a transcript incomplete: see
[mitochondrial stop codons](#mitochondrial-stop-codons). Selecting an NMD,
IG/TR or pseudogene transcript does not by itself request a different coding
model. Additional supported coding frames are supplied through ORFdata's
`extra_orfs` API. [Ensembl](https://www.ensembl.org/info/genome/genebuild/biotypes.html)
and [GENCODE](https://www.gencodegenes.org/pages/biotypes.html) define the
reference biotypes separately from the sequence consequences calculated here.

If all overlapping transcripts are rejected, their effects are absent and
the collection can be empty. Filtering does not turn an overlap into an
`Intergenic` or `Intragenic` result; genuine transcript-free overlaps remain.
Errors inside the predicate follow `raise_on_error`, like annotation errors.
A non-callable `transcript_filter` always raises `TypeError`.

<a id="read-an-effect"></a>

## Protein sequences

Given an effect from the collection:

```python
protein = effect.mutant_protein_sequence  # may be None
```

`None` means the sequence could not be determined, not that the protein is
unchanged. The prediction belongs to `effect.transcript`; intergenic effects
have no transcript. For cDNA, partial structures, and sequence evidence, see
[Transcript models](transcript_models.md).

To keep only effects that may change the protein, use
`effects.drop_silent_and_noncoding()`. It removes known silent and noncoding
predictions but keeps unresolved ones by default; see
[filtering by protein change](structural_variants.md#filtering-by-protein-change).

<a id="mitochondrial-stop-codons"></a>

### Mitochondrial stop codons

MT-ND1, MT-ND2, MT-CO3, MT-ND3, MT-ND4 and MT-CYB end their reading frame in
`T` or `TA`, and the poly(A) tail adds the rest of a `TAA` stop codon
([NC_012920.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1)). Ensembl
annotates no stop codon for them. Varcode translates each as a complete
transcript, with the tail finishing the last partial codon, so m.3460G>A in
MT-ND1 is `p.A52T` and m.11778G>A in MT-ND4 is `p.R340H`.

A change to the `T` or `TA` that the genome does encode can lose the stop. The
last partial codon is read with the `A`s that finish it, and the reading frame
ends there, since the tail past it is not in the genome:

| variant | codon read | effect |
|---|---|---|
| MT-ND1 m.4261T>C (`TA`) | `CAA` | `StopLoss`, `p.*319Q` |
| MT-ND1 m.4262A>G (`TA`) | `TGA`, Trp in the mitochondrial code | `StopLoss`, `p.*319W` |
| MT-ND4 m.12137T>A (`T`) | `AAA` | `StopLoss`, `p.*460K` |

This needs PyEnsembl 2.34.1 or later. With an earlier PyEnsembl these genes are
still `IncompleteTranscript`.

`MT-CO1` and `MT-ND6` instead carry annotated `AGA` and `AGG` `stop_codon`
features in Ensembl. Varcode decodes them with NCBI table 2's static stop set
(`TAA`, `TAG`, `AGA`, `AGG`) rather than modeling the `-1` mitoribosomal
frameshift onto a preceding `U` reported by Temperley et al. (2010), so
`MT-CO1` m.7445A>G (`AGA` to `AGG`) and m.7442T>C (`TCT` to `TCC`, `p.S513=`)
are `Silent`, while m.7444G>A (`AGA` to `AAA`, `p.*514K`) is `StopLoss`.

## Filter and group a collection

```python
from varcode.effects import Substitution

print(variants.gene_counts())  # number of variants overlapping each gene name
cftr_effects = variants.groupby_gene_name()["CFTR"].effects()
most_severe = effects.drop_silent_and_noncoding().top_priority_effect_per_gene_id()
substitution_or_worse = effects.filter_by_effect_priority(Substitution)
```

`groupby_gene_name()` returns a dict from gene name to the variants that
overlap that gene; a variant overlapping several genes appears in each group.
`top_priority_effect_per_gene_id()` returns a dict from gene ID to that gene's
highest-priority effect. `filter_by_effect_priority(cls)` keeps the effects at
least as severe as `cls` in Varcode's
[priority order](api_effects.md#varcode.effect_priority).

<a id="reading-alternatives"></a>

## Alternative outcomes

Every `MutationEffect` exposes the same candidate interface — `.candidates`,
`.effects`, `.most_likely_candidate`, `.most_likely_effect`,
`.highest_priority_candidate`, `.highest_priority_effect`, and
`.add_candidates(...)`. Bare `"fast"` effects return `()` / `None`:

```python
for candidate in effect.candidates:
    comparison = candidate.effect
    print(comparison.short_description)
    print(comparison.original_protein_sequence)
    print(comparison.mutant_protein_sequence)
```

`effect.most_likely_effect` is the preferred comparison and
`effect.highest_priority_effect` is the most severe one. Preference can
reflect a splice rule or supplied RNA support; it is not a calibrated
probability. Each `EffectCandidate` carries `source`, `mechanisms`, `phase`,
`uncertainty`, `evidence`, and (when produced by the transcript-model
composition path) paired baseline/mutant `outcomes` that survive JSON
round-trip. See [the result contract](experimental_annotators.md#transcript-model-results).

Collections filter and rank across all candidates while retaining the complete
result. Use `effect.add_candidates(...)` to attach external observations (for
example from RNA or custom splice scorers) without wrapping or mutating private
attributes.

<a id="choose-a-deeper-topic"></a>

## Related guides

- <a id="splice-disrupting-variants"></a><a id="when-splice-disruption-is-in-play"></a><a id="splice-and-coding-effects-can-co-occur"></a><a id="the-spliceoutcomeset-shape"></a><a id="common-questions"></a><a id="rna-evidence-reconciliation"></a><a id="candidate-provenance"></a><a id="picking-a-single-candidate"></a><a id="limitations"></a>[Splice variants](splice_variants.md): affected signals, candidate mechanisms, and RNA evidence.
- <a id="structural-variants"></a>[Structural variants](structural_variants.md): deletions, duplications, inversions, and fusions.
- [Germline](germline.md) and [phasing](phasing.md): patient-specific baselines and linked variants.
- <a id="how-it-composes"></a><a id="the-four-primitives"></a>[Transcript models](transcript_models.md): cDNA, partial structures, and missing sequences.
- <a id="provenance"></a>[Saving results](csv.md#annotation-provenance): tables, provenance, and round-trip limits.
- <a id="annotator-selection"></a><a id="advanced-annotators-and-implementation-limits"></a><a id="what-the-optional-implementations-change"></a><a id="reading-sequences-and-alternatives"></a><a id="current-boundaries-and-legacy-names"></a>[Annotation pipeline and compatibility](experimental_annotators.md): default results, migration, and implementation limits.
- <a id="writing-an-annotator"></a><a id="downstream-consumers"></a>[Writing an annotator](annotator_contract.md): custom implementations and unsupported inputs.
- [Effect types](effect_types.md) and [API reference](api.md): class definitions and parameters.
