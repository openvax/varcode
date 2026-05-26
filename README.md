[![Tests](https://github.com/openvax/varcode/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/varcode/actions/workflows/tests.yml)
[![Coverage Status](https://coveralls.io/repos/openvax/varcode/badge.svg?branch=main&service=github)](https://coveralls.io/github/openvax/varcode?branch=main)
[![PyPI](https://img.shields.io/pypi/v/varcode.svg?maxAge=1000)](https://pypi.python.org/pypi/varcode/)
[![PyPI downloads](https://img.shields.io/pypi/dm/varcode.svg)](https://pypistats.org/packages/varcode)

# Varcode

Varcode is a library for working with genomic variant data in Python and predicting the impact of those variants on protein sequences.

## Installation

You can install varcode using [pip](https://pip.pypa.io/en/latest/quickstart.html):

```bash
pip install varcode
```

You can install required reference genome data through [PyEnsembl](https://github.com/openvax/pyensembl) as follows:

```bash
# Downloads and installs the Ensembl releases (75 and 76)
pyensembl install --release 75 76
```

## Example

```python
import varcode

# Load TCGA MAF containing variants from their
variants = varcode.load_maf("tcga-ovarian-cancer-variants.maf")

print(variants)
### <VariantCollection from 'tcga-ovarian-cancer-variants.maf' with 6428 elements>
###  -- Variant(contig=1, start=69538, ref=G, alt=A, genome=GRCh37)
###  -- Variant(contig=1, start=881892, ref=T, alt=G, genome=GRCh37)
###  -- Variant(contig=1, start=3389714, ref=G, alt=A, genome=GRCh37)
###  -- Variant(contig=1, start=3624325, ref=G, alt=T, genome=GRCh37)
###  ...

# you can index into a VariantCollection and get back a Variant object
variant = variants[0]

# groupby_gene_name returns a dictionary whose keys are gene names
# and whose values are themselves VariantCollections
gene_groups = variants.groupby_gene_name()

# get variants which affect the TP53 gene
TP53_variants = gene_groups["TP53"]

# predict protein coding effect of every TP53 variant on
# each transcript of the TP53 gene
TP53_effects = TP53_variants.effects()

print(TP53_effects)
### <EffectCollection with 789 elements>
### -- PrematureStop(variant=chr17 g.7574003G>A, transcript_name=TP53-001, transcript_id=ENST00000269305, effect_description=p.R342*)
### -- ThreePrimeUTR(variant=chr17 g.7574003G>A, transcript_name=TP53-005, transcript_id=ENST00000420246)
### -- PrematureStop(variant=chr17 g.7574003G>A, transcript_name=TP53-002, transcript_id=ENST00000445888, effect_description=p.R342*)
### -- FrameShift(variant=chr17 g.7574030_7574030delG, transcript_name=TP53-001, transcript_id=ENST00000269305, effect_description=p.R333fs)
### ...

premature_stop_effect = TP53_effects[0]

print(str(premature_stop_effect.mutant_protein_sequence))
### 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMF'

print(premature_stop_effect.aa_mutation_start_offset)
### 341

print(premature_stop_effect.transcript)
### Transcript(id=ENST00000269305, name=TP53-001, gene_name=TP53, biotype=protein_coding, location=17:7571720-7590856)

print(premature_stop_effect.gene.name)
### 'TP53'
```

If you are looking for a quick start guide, you can check out [this iPython book](https://github.com/openvax/varcode/blob/main/examples/varcode-quick_start.ipynb) that demonstrates simple use cases of Varcode.

## Further reading

Feature guides live in [`docs/`](https://github.com/openvax/varcode/tree/main/docs) and on the [docs site](https://openvax.github.io/varcode/):

- [**Effect annotation**](https://openvax.github.io/varcode/effect_annotation/) — how variants
  become effects, splice outcome representations, pluggable annotators,
  and structural variants.
- [**Genotypes and sample-aware queries**](https://openvax.github.io/varcode/genotype/) —
  per-sample zygosity on multi-sample VCFs.
- [**VariantCollection transforms**](https://openvax.github.io/varcode/transforms/) —
  `pair_breakends`, `left_align_indels`, and the transform contract.
- [**CSV round-trip and metadata headers**](https://openvax.github.io/varcode/csv/) —
  `to_csv` / `from_csv` with `#`-prefixed provenance headers.
- [**Error handling**](https://openvax.github.io/varcode/errors/) — `ReferenceMismatchError`,
  `GenomeBuildMismatchError`, `SampleNotFoundError`, and
  `raise_on_error=False`.
- [**Germline-aware annotation**](https://openvax.github.io/varcode/germline/) — niche.
  Classify somatic variants against the patient's germline-applied
  transcript when somatic and germline share a codon.

See [`CHANGELOG.md`](https://github.com/openvax/varcode/blob/main/CHANGELOG.md) for the release history.

## Effect Types

Every concrete `MutationEffect` subclass that varcode may emit, grouped
by biological context. Each row links to the class definition in
[`varcode/effects/effect_classes.py`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py)
via a browser text-fragment URL — links survive line-number drift as
the source file evolves. Severity ordering across types is set by
[`effect_priority()`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_ordering.py);
the abstract bases — `MutationEffect`, `TranscriptMutationEffect`,
`CodingMutation`, `NonsilentCodingMutation`, `SpliceMechanismEffect`,
`StructuralVariantEffect` — define the shared interface;
`MutationEffect` and `NonsilentCodingMutation` have dedicated entries
in the [API reference](https://openvax.github.io/varcode/api/), and
`MultiOutcomeEffect` is described in its own section below.

### Effects that carry multiple possibilities

Several effects don't have a single deterministic protein-level
outcome — splice-signal disruption can resolve as normal splicing,
exon skipping, intron retention, or cryptic-site use; an exon-edge
variant might be a routine coding effect *or* a splice disruption;
a structural variant might affect multiple transcripts or have
multiple plausible breakpoint resolutions; two or more cis variants
on one transcript can compose into a joint mutant protein; and when
phase is unknown between somatic and germline variants sharing a
window, the somatic effect depends on which haplotype it landed on.
Varcode wraps these as
[`MultiOutcomeEffect`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20MultiOutcomeEffect%28)
instances. Every multi-outcome effect exposes the same surface:

- `.candidates` → `tuple[EffectCandidate, ...]`
- Each [`EffectCandidate`](https://github.com/openvax/varcode/blob/main/varcode/effect_candidates.py#:~:text=class%20EffectCandidate%28)
  wraps an inner effect (`.effect`), a producer tag (`.source`,
  e.g. `"varcode"`, `"rna_evidence"`), and free-form `.evidence`.
- `.effects` → `tuple[MutationEffect, ...]` — convenience that unwraps `.candidates` to inner effects when provenance isn't needed.
- `.most_likely_candidate` / `.most_likely_effect` — producer-ordered top pick.
- `.highest_priority_candidate` / `.highest_priority_effect` — most severe by `effect_priority()`.

The `MultiOutcomeEffect` containers appear in the sub-tables below
where they're emitted:
[`SpliceOutcomeSet`](https://github.com/openvax/varcode/blob/main/varcode/splice_outcomes.py#:~:text=class%20SpliceOutcomeSet%28),
[`ExonicSpliceSite`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20ExonicSpliceSite%28),
the [`StructuralVariantEffect`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20StructuralVariantEffect%28)
sub-hierarchy,
[`HaplotypeEffect`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20HaplotypeEffect%28),
and
[`PhaseCandidateSet`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20PhaseCandidateSet%28).

### Coding region — in-frame changes

| Effect type | Description |
| --- | --- |
| [`Substitution`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20Substitution%28) | Coding mutation which causes simple substitution of one amino acid for another. |
| [`Insertion`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20Insertion%28) | Coding mutation which causes insertion of amino acid(s). |
| [`Deletion`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20Deletion%28) | Coding mutation which causes deletion of amino acid(s). |
| [`ComplexSubstitution`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20ComplexSubstitution%28) | Insertion and deletion of multiple amino acids. |
| [`Silent`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20Silent%28) | Mutation in coding sequence which does not change the amino acid sequence of the translated protein. |
| [`AlternateStartCodon`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20AlternateStartCodon%28) | Replace annotated start codon with alternative start codon (_e.g._ "ATG>CAG"); a `Silent` subclass since the initiator tRNA still loads Met. |

### Coding region — frame-disrupting / truncating

| Effect type | Description |
| --- | --- |
| [`FrameShift`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20FrameShift%28) | Out-of-frame insertion or deletion of nucleotides, causes novel protein sequence and often premature stop codon. |
| [`FrameShiftTruncation`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20FrameShiftTruncation%28) | A frameshift which leads immediately to a stop codon (no novel amino acids created). |
| [`PrematureStop`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20PrematureStop%28) | Insertion of stop codon, truncates protein. |
| [`StartLoss`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20StartLoss%28) | Mutation causes loss of start codon, likely result is that an alternate start codon will be used down-stream (possibly in a different frame). |
| [`StopLoss`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20StopLoss%28) | Loss of stop codon, causes extension of protein by translation of nucleotides from 3' UTR. |

### Splice-site disruption — *where* the signal was hit

DNA-level locations: these classes identify *where* a variant landed
in the canonical splice window. They no longer appear as top-level
effects — every splice-disrupting variant is wrapped in a
`SpliceOutcomeSet` (varcode 6.0+), and these classes survive as the
wrapper's `disrupted_signal_class` (a type) and as the
`splice_signal` reference on each candidate mechanism (an instance).
All four share the
[`SpliceSite`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20SpliceSite%28)
base.

| Effect type | Description |
| --- | --- |
| [`SpliceDonor`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20SpliceDonor%28) | Mutation at canonical donor `GT` (intronic +1/+2). |
| [`SpliceAcceptor`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20SpliceAcceptor%28) | Mutation at canonical acceptor `AG` (intronic -2/-1). |
| [`IntronicSpliceSite`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20IntronicSpliceSite%28) | Other intronic positions in the splice window (+3..+6 donor, -3 acceptor; also +1/+2 or -1/-2 when the reference signal isn't canonical). |
| [`ExonicSpliceSite`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20ExonicSpliceSite%28) | Last 3 bases of an exon (donor side) or first base of an exon (acceptor side); changes a codon *and* disrupts the splice signal. Carries `alternate_effect` (the coding consequence if splicing proceeds). |

### Splice mechanism — *what the spliceosome does* in response

**These are the splice effects that carry a protein consequence.** The
protein-level outcome of a splice-signal hit is not deterministic from
DNA alone, so varcode wraps every splice-signal disruption in a
[`SpliceOutcomeSet`](https://github.com/openvax/varcode/blob/main/varcode/splice_outcomes.py#:~:text=class%20SpliceOutcomeSet%28)
(a `MultiOutcomeEffect`) carrying these mechanisms as candidates.
Wrapping is always-on as of varcode 6.0 and lazy — only the cheap
`NormalSplicing` candidate is built eagerly; the rest materialise on
`.candidates` access. Each mechanism carries the originating
disruption on its `.splice_signal` attribute (a `SpliceSite`
*instance*), so you can always recover *where* the hit was off any
mechanism. The set also records the disruption's *class* on
`.disrupted_signal_class` (the `SpliceSite` subclass, e.g.
`SpliceDonor` — a type, not an instance) for priority lookup.

| Effect type | Description |
| --- | --- |
| [`NormalSplicing`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20NormalSplicing%28) | Splice signal hit but splicing proceeds normally; protein consequence (if any) is whatever the underlying nucleotide change would produce. |
| [`ExonSkipping`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20ExonSkipping%28) | Affected exon excluded from the mature transcript; in-frame skip deletes amino acids, out-of-frame skip propagates a frameshift. |
| [`IntronRetention`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20IntronRetention%28) | Intron stays in the mature transcript; translation usually hits a premature stop inside the retained intron. |
| [`CrypticDonor`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20CrypticDonor%28) | Disrupted canonical donor replaced by a nearby cryptic GT donor; exon extended or truncated. |
| [`CrypticAcceptor`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20CrypticAcceptor%28) | Disrupted canonical acceptor replaced by a nearby cryptic AG acceptor; exon extended or truncated. |

### Non-coding regions and unclassifiable contexts

| Effect type | Description |
| --- | --- |
| [`FivePrimeUTR`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20FivePrimeUTR%28) | Variant affects 5' untranslated region before start codon. |
| [`ThreePrimeUTR`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20ThreePrimeUTR%28) | Variant affects 3' untranslated region after stop codon of mRNA. |
| [`Intronic`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20Intronic%28) | Variant occurs between exons and is unlikely to affect splicing. |
| [`NoncodingTranscript`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20NoncodingTranscript%28) | Transcript doesn't code for a protein. |
| [`IncompleteTranscript`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20IncompleteTranscript%28) | Can't determine effect since transcript annotation is incomplete (often missing either the start or stop codon). |
| [`Intergenic`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20Intergenic%28) | Occurs outside of any annotated gene. |
| [`Intragenic`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20Intragenic%28) | Within the annotated boundaries of a gene but not in a region that's transcribed into pre-mRNA. |
| [`Failure`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20Failure%28) | Placeholder effect emitted when annotation failed but a non-empty effect list is required (`raise_on_error=False`). |

### Exon-level and structural-variant effects

`ExonLoss` is a plain `Exonic` effect. The structural-variant
effects below (`LargeDeletion` through `TranslocationToIntergenic`)
are `MultiOutcomeEffect`s — their `.candidates` may include
cryptic-exon outcomes, RNA-evidence-ranked alternatives, and so on.
`CrypticExonCandidate` typically appears as a candidate inside those
SV effects rather than standalone.

| Effect type | Description |
| --- | --- |
| [`ExonLoss`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20ExonLoss%28) | Deletion of an entire exon, significantly disrupts protein. |
| [`LargeDeletion`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20LargeDeletion%28) | Structural deletion (`<DEL>` / `<CN0>`) removing one or more exons or an entire gene. |
| [`LargeDuplication`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20LargeDuplication%28) | Tandem duplication (`<DUP>`) overlapping exons; may yield copy-number increase or a fused reading frame. |
| [`Inversion`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20Inversion%28) | Inversion (`<INV>`) flipping a stretch of a transcript; consequence depends on whether breakpoints fall in exons or introns. |
| [`GeneFusion`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20GeneFusion%28) | Breakend (`<BND>`) whose mate lies in another protein-coding gene — the canonical fusion shape. |
| [`TranslocationToIntergenic`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20TranslocationToIntergenic%28) | Breakend whose mate lies in intergenic space; consequence depends on cryptic splice / ORF signals downstream. |
| [`CrypticExonCandidate`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20CrypticExonCandidate%28) | An SV brings novel sequence into range of a transcript and motif scoring flags a plausible new splice acceptor / donor pair; attached as additional candidates on SV effects. |

### Multi-variant / phase-dependent effects

Both are `MultiOutcomeEffect`s, emitted alongside per-variant effects
(additive, not a replacement) when a phase resolver groups cis
variants together or when phase between somatic and germline variants
is unknown.

| Effect type | Description |
| --- | --- |
| [`HaplotypeEffect`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20HaplotypeEffect%28) | Joint effect of two or more cis variants on the same transcript; the combined mutant cDNA is built and translated as one unit. |
| [`PhaseCandidateSet`](https://github.com/openvax/varcode/blob/main/varcode/effects/effect_classes.py#:~:text=class%20PhaseCandidateSet%28) | Possibility set across phase hypotheses when a somatic variant and one or more germline variants share a window on a transcript and phase between them is unknown. |

## Coordinate System

Varcode currently uses a "base counted, one start" genomic coordinate system, to match the Ensembl annotation database. We are planning to switch over to "space counted, zero start" (interbase) coordinates, since that system allows for more uniform logic (no special cases for insertions). To learn more about genomic coordinate systems, read this [blog post](http://alternateallele.blogspot.com/2012/03/genome-coordinate-conventions.html).
