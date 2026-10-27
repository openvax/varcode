# Effect types

Use this page to look up a class after reading the
[effect annotation guide](effect_annotation.md). Definitions below are generated
from the source docstrings; this is reference material, not a required first read.

<a id="find-an-effect-family"></a>

## Effect families

| Family | Examples |
|---|---|
| [Coding changes](#coding-effects) | `Substitution`, `Insertion`, `Deletion`, `FrameShift`, `PrematureStop`, `Silent` |
| [Noncoding and untranslated regions](#noncoding-effects) | `Intronic`, `FivePrimeUTR`, `ThreePrimeUTR`, `NoncodingTranscript`, `Intergenic` |
| [Splice effects](#splice-effects) | `SpliceOutcomeSet`, with `NormalSplicing`, `ExonSkipping`, `IntronRetention`, or cryptic-site candidates |
| [Structural effects](#structural-effects) | `StructuralVariantEffect`, `GeneFusion`, `TranslocationToIntergenic` |
| [Linked variants and phase](#haplotype-effects) | `HaplotypeEffect`, `PhaseCandidateSet` |
| [Inherited allele overlap](germline.md#loss-of-heterozygosity-loh) | `GermlineAlleleOverlap`: no new allele sequence; LOH unassessed |
| [Unknown results and failures](#unknown-results) | `Unresolved`, `HypothesisLimit`, `IncompleteTranscript`, `Failure` |

These labels describe predictions and their limits, not clinical significance.

Some effects describe a change to the protein; others only say where the
variant landed. Coding effects (`Substitution`, `FrameShift`, `PrematureStop`,
...) and splice mechanism effects (`NormalSplicing`, `ExonSkipping`, ...)
describe a protein change. Splice-signal disruption effects (`SpliceDonor`,
`SpliceAcceptor`, `IntronicSpliceSite`, `ExonicSpliceSite`, all subclasses of
`SpliceSite`) and region effects (`Intronic`, `FivePrimeUTR`, ...) carry no
protein consequence on their own.

Most effects are a single answer. Subclasses of `MultiOutcomeEffect`
(`SpliceOutcomeSet`, `ExonicSpliceSite`, the structural-variant effects,
`HaplotypeEffect`, `PhaseCandidateSet`) hold several candidate effects when the
protein-level outcome isn't determined. Each provides `.candidates`,
`.most_likely_effect`, and `.highest_priority_effect`.

Severity ordering across types is set by
[`effect_priority`](api_effects.md#varcode.effect_priority).

<a id="class-reference"></a>
<a id="varcode.effects.effect_classes"></a>

## Coding effects

### Protein deletion labels

By default, `Deletion.short_description` includes the complete deleted sequence
and its one-based starting position. `p.RSP3del` means deleting `RSP` starting at
residue 3. To display the first and last deleted residues instead:

```python
effect.short_description                                      # "p.RSP3del"
effect.format_short_description(deletion_notation="range")     # "p.R3_P5del"
effect.format_short_description(deletion_notation="sequence")  # "p.RSP3del"

effects.to_dataframe(deletion_notation="range")
effects.to_csv("effects.csv", deletion_notation="range")
```

Single-residue deletions use the same label in both forms, e.g. `p.K405del`.
The option also follows nested deletion effects in normal-splicing, splice-set,
phase and structural summaries. Other consequence and splice-mechanism labels
retain their existing notation. The endpoint form follows the
[HGVS deletion range convention](https://hgvs-nomenclature.org/stable/recommendations/protein/deletion/);
these short labels are not a complete HGVS representation.

This is a per-call display option. It leaves `aa_ref`, zero-based effect offsets,
repeat normalization and predicted proteins unchanged. JSON archives and CSV
reloads retain the biological effect, whose `short_description` uses the default
sequence form. Use the formatter again when a range label is wanted.

::: varcode.effects.effect_classes.Silent

::: varcode.effects.effect_classes.AlternateStartCodon

::: varcode.effects.effect_classes.StartLoss

### Start-loss nucleotide details

`StartLoss` retains the known sequence change while leaving
`mutant_protein_sequence` as `None` and the summary as `p.M1? (start-loss)`.
For example, with Ensembl 81 / GRCh38:

```python
from pyensembl import cached_release
from varcode import Variant

genome = cached_release(81)
transcript = genome.transcript_by_id("ENST00000003084")
effect = Variant("7", 117480096, "T", "C", genome).effect_on_transcript(transcript)
assert (effect.ref_codon, effect.alt_codon) == ("ATG", "ACG")
assert effect.alt_codon_amino_acid == "T"
assert effect.mutant_protein_sequence is None
edit, = effect.mutant_transcript.edits
```

Codons and edits use transcript orientation, including reverse-strand
transcripts. `ref_codon_offset` and `alt_codon_offset` are zero-based offsets
in the baseline and mutant cDNA. Each resolved `TranscriptEdit` carries its
reference cDNA interval, alternate bases and source variant. Realized compound
products retain their complete cDNA and mapping evidence; they carry a single
edit only when that edit explains the entire product.

`codon_change_status` describes whether the alternate triplet can be reported:

| Status | Meaning |
| --- | --- |
| `mapped` | The start position maps to a complete triplet in the mutant cDNA. |
| `start_removed` | The mapped original start position was removed. |
| `start_ambiguous` | Multiple retained copies prevent a unique mapping. |
| `edit_unresolved` | The edit or start mapping could not be established. |
| `sequence_unavailable` | The required sequence or complete triplet is unavailable. |

Unavailable codons and alternate offsets are `None`; available edits and cDNA
remain accessible through `mutant_transcript`. `alt_codon_amino_acid` gives the
triplet's **internal decoding**, using the transcript's genetic code (`"*"`
for a stop, `None` for an unknown codon). This does not establish translation
initiation or an alternate protein. See [NCBI's genetic-code tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
for the distinction between internal decoding and initiator methionine.
Recognized alternative initiators retain `AlternateStartCodon` classification.

Native `to_json()` / `from_json()` retain these details. Collections containing
start-loss effects expose codon fields in `to_dataframe()` and include the
native `effect_json` archive in `to_csv()`; `EffectCollection.from_csv()` restores
these saved results without reannotation.

::: varcode.effects.effect_classes.Substitution

::: varcode.effects.effect_classes.ComplexSubstitution

::: varcode.effects.effect_classes.Insertion

::: varcode.effects.effect_classes.Deletion

::: varcode.effects.effect_classes.PrematureStop

::: varcode.effects.effect_classes.StopLoss

::: varcode.effects.effect_classes.FrameShift

::: varcode.effects.effect_classes.FrameShiftTruncation

::: varcode.effects.effect_classes.ExonLoss

## Noncoding effects

::: varcode.effects.effect_classes.Intergenic

::: varcode.effects.effect_classes.NoncodingTranscript

::: varcode.effects.effect_classes.FivePrimeUTR

::: varcode.effects.effect_classes.ThreePrimeUTR

::: varcode.effects.effect_classes.Intronic

## Splice effects

<a id="splice-outcome-container"></a>

See the [splice guide](splice_variants.md) for candidate access.

::: varcode.SpliceOutcomeSet

::: varcode.effects.effect_classes.SpliceSite

::: varcode.effects.effect_classes.IntronicSpliceSite

::: varcode.effects.effect_classes.SpliceDonor

::: varcode.effects.effect_classes.SpliceAcceptor

::: varcode.effects.effect_classes.ExonicSpliceSite

::: varcode.effects.effect_classes.SpliceMechanismEffect

::: varcode.effects.effect_classes.NormalSplicing

::: varcode.effects.effect_classes.ExonSkipping

::: varcode.effects.effect_classes.IntronRetention

::: varcode.effects.effect_classes.CrypticSpliceSiteEffect

::: varcode.effects.effect_classes.CrypticDonor

::: varcode.effects.effect_classes.CrypticAcceptor

::: varcode.effects.effect_classes.CrypticExonCandidate

## Structural effects

::: varcode.effects.effect_classes.StructuralVariantEffect

Legacy event-type classes remain importable; the default annotator emits
[transcript consequences](structural_variants.md#reading-sv-results).

::: varcode.effects.effect_classes.LargeDeletion

::: varcode.effects.effect_classes.LargeDuplication

::: varcode.effects.effect_classes.Inversion

::: varcode.effects.effect_classes.GeneFusion

::: varcode.effects.effect_classes.TranslocationToIntergenic

## Haplotype effects

::: varcode.effects.effect_classes.HaplotypeEffect

::: varcode.effects.effect_classes.PhaseCandidateSet

::: varcode.effects.effect_classes.GermlineAlleleOverlap

## Unknown results

::: varcode.effects.effect_classes.Unresolved

::: varcode.effects.effect_classes.HypothesisLimit

::: varcode.effects.effect_classes.IncompleteTranscript

::: varcode.effects.effect_classes.Failure

## Base classes

::: varcode.effects.effect_classes.MutationEffect

::: varcode.effects.effect_classes.MultiOutcomeEffect

::: varcode.effects.effect_classes.Intragenic

::: varcode.effects.effect_classes.TranscriptMutationEffect

::: varcode.effects.effect_classes.Exonic

::: varcode.effects.effect_classes.CodingMutation

::: varcode.effects.effect_classes.NonsilentCodingMutation

::: varcode.effects.effect_classes.KnownAminoAcidChange
