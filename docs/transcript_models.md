# Transcript models

An effect may carry a `MutantTranscript` describing edits, retained reference
segments, or an observed RNA structure. It can contain cDNA and protein sequence,
but it may also be partial or have no resolved sequence.

You usually want the effect's `mutant_protein_sequence`; a transcript model is
useful when you also need the RNA structure or its provenance. Not every effect
has a model, even when a predicted protein is available.

## Sequence access

Given an annotated `effect`:

```python
protein = effect.mutant_protein_sequence
model = effect.mutant_transcript
cdna = model.cdna_sequence if model is not None else None
evidence = model.evidence if model is not None else None
baseline = effect.baseline_transcript
baseline_cdna = baseline.cdna_sequence if baseline is not None else None
baseline_protein = effect.original_protein_sequence
```

`None` means unavailable. It does not establish an unchanged protein, an absent
transcript, or a harmless variant.

For a composed comparison, `baseline_transcript` and `mutant_transcript` are
the actual products used to classify the effect. `original_protein_sequence`
is then the patient baseline, including inherited changes, and amino-acid
offsets refer to that baseline. The annotation reference remains available as
`effect.reference_transcript`. These models and sequences survive JSON and
pickle round trips.

The transcript-model annotator assembles selected exon, retained-intron and
cryptic-splice intervals through one cDNA/translation path. Unavailable paths
have `None` sequences; they never substitute the canonical transcript. Each
candidate exposes its own products through the same effect accessors. Products
with different distant sequences remain distinct candidates, even if their
local amino-acid changes match.

The default considers supplied germline calls across the whole transcript span.
Unknown phase produces alternatives or a hypothesis limit; missing germline
calls cannot be recovered. The explicit compatibility selections `fast` and
`protein_diff` retain their local germline lookup. See
[annotation choices](experimental_annotators.md).

## Selenocysteine

Ensembl marks selenocysteine (Sec) as `U` in the reference protein, encoded by
an in-frame UGA. UGA is decoded as Sec only with a SECIS element in the same
mRNA's 3′ UTR, whose position isn't annotated. Predicted proteins read an
annotated Sec codon as `U` when the model still maps it (no edit touches it)
and keeps some of a selenoprotein's 3′ UTR. Where no selenoprotein 3′ UTR
remains, for example a fusion downstream of Sec, UGA is read as a stop.
A sequence without reference coordinates, such as an imported RNA assembly,
is translated literally. A partly kept 3′ UTR still gets the Sec reading; the
[SV change flags](structural_variants.md#filtering-by-protein-change) report
such cases as unresolved.

## Partial structures

A model with reference segments is not necessarily a complete mutant RNA.
For example, a breakend can retain only a transcript's local 5′ prefix or 3′
suffix. That model has `evidence["sequence_status"] = "retained_reference_fragment"`,
while full cDNA and protein remain unknown. Concatenating its segments does
not reconstruct the missing partner or establish the full allele.

Likewise, a sequence supplied by an RNA assembler may be only a junction
fragment. Keep its completeness and source information with the sequence;
neither observed RNA nor sequence-predicted ORFs alone prove translation.

## Splice junctions and open rearrangement templates

For downstream RNA verification (such as Isovar), `MutantTranscript` exposes
junction coordinates and open rearrangement templates directly:

- `model.is_complete`: `True` when a full spliced `cdna_sequence` is resolved;
  `False` for open rearrangements (large insertions `>= 50` bp or translocations
  without a supplied `alt_assembly`).
- `model.novel_junctions`, `model.back_junctions`, `model.wildtype_junctions`:
  1-based inclusive ascending genomic tuples `(contig, start, end, strand)` with
  `start <= end` on both `+` and `-` strands, directly matching BAM CIGAR `N`
  intron spans (or chimeric tandem-duplication back-junctions).
- `model.mate_region` and `model.exit_frame`: oriented 1-based mate search
  window `(contig, start, end, strand)` and exiting 5′ CDS reading frame
  (`0`, `1`, `2`, or `None`) for open rearrangements.
- `model.template_sequence(sequence_provider)`: materializes the oriented
  5′-to-3′ mutant pre-mRNA template sequence on demand from a genomic sequence
  callable `(contig, start_1based, end_1based) -> str`.

## Related objects

Effects produced by `annotator="dna"` (the default) and
`annotator="transcript_model"` retain their candidates, paired `outcomes`,
baseline and mutant products, evidence, and ranking through
`EffectCollection.to_json()` / `from_json()` and individual effect
serialization (`effect_to_json()` / `effect_from_json()`). Extra-ORF effects
also round-trip natively through the versioned `orf_effect_schema` graph. Older
effect JSON archives remain readable. Transcript references still use the
existing PyEnsembl serialization contract; the archive is not a replacement for
reference annotation data.

| Object | Contents |
|---|---|
| `MutationEffect` | A predicted consequence or unresolved result, exposing `.candidates`, `.most_likely_candidate`, `.highest_priority_candidate`, and `.add_candidates()` |
| `EffectCandidate` | One candidate outcome with `effect`, `source`, `evidence`, `mechanisms`, `phase`, `uncertainty`, and paired `outcomes` |
| `MutantTranscript` | Edits or transcript structure, provenance, optional sequences, splice junctions, and open rearrangement templates |

The [API reference](api_rna.md#mutant-transcripts) lists model fields and
construction helpers; see [result access](experimental_annotators.md#transcript-model-results)
for comparing `dna` and `transcript_model` outputs.


Ambiguous reference codons retain their uniquely determined amino acid, or `X`
when no unique residue is established. `X` is an unknown residue, including a
possible stop for mixed sense/stop codons; it is not a concrete amino acid.
The original cDNA remains available. `translate_sequence` in
`varcode.effects.codon_tables` applies this rule under every supported NCBI
genetic code.
