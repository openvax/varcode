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
`IncompleteTranscript` for incomplete coding annotations. Selecting an NMD,
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

<a id="reading-alternatives"></a>

## Alternative outcomes

The default pipeline attaches alternatives to an ordinary classified effect.
Iterate `candidates` without requiring a wrapper class:

```python
for candidate in getattr(effect, "candidates", ()):
    comparison = candidate.effect
    print(comparison.short_description)
    print(comparison.original_protein_sequence)
    print(comparison.mutant_protein_sequence)
```

For shared transcript results, `most_likely_candidate.effect` is the preferred
comparison and `highest_priority_candidate.effect` is the most severe one.
Preference can reflect a splice rule or supplied RNA support; it is not a
calibrated probability. Each candidate retains its `outcomes`, including the
baseline/mutant products and phase/splice hypotheses. See
[the result contract](experimental_annotators.md#transcript-model-results).

Existing structural builders and explicit `annotator="fast"` can return
`MultiOutcomeEffect` wrappers. Those expose `most_likely_effect`,
`highest_priority_effect`, and candidates with `source` and `evidence` fields.
Keep this provenance when reporting a candidate's sequence. Collections filter
and rank across all candidates while retaining the complete result.

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
