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
see [reading alternatives](effect_annotation.md#reading-alternatives).

A DNA rearrangement does not by itself establish a complete expressed fusion
protein. Sequence may be unknown or partial; `None` is not an unchanged protein.
Do not rely on `drop_silent_and_noncoding()` to retain unresolved SV effects:
their protein-change flags remain incomplete
([#418](https://github.com/openvax/varcode/issues/418)).

## Find the detail you need

- [Which transcripts get annotated](#which-transcripts-get-annotated) and
  [effect classes](#effect-classes).
- [Pairing breakend records](#pairing-breakend-records), including when pairing
  changes which partner transcripts appear.
- [Alleles, coordinates, and exports](#alleles-coordinates-and-exports).
- [Junction orientation](#junctions-which-side-of-each-breakpoint-is-kept),
  [fusion partners](#5-and-3-roles), and [breakpoint position](#where-in-the-gene-a-breakpoint-lands).
- [Observed RNA import](#importing-observed-rna-structures) and
  [Exacto protein predictions](#import-exactos-protein-predictions).

The detailed biological examples below use Ensembl 95.

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

## Alleles, coordinates, and exports

Structural `ref`/`alt` and `original_ref`/`original_alt` expose the actual record
alleles, including symbolic or breakend ALT. Small-edit flags (`is_snv`,
`is_indel`, `is_insertion`, `is_deletion`, `is_transition`, `is_transversion`)
are false; use `is_structural` and `sv_type` for event classification.
SV-containing tables add `sv_type`, `end`, `mate_contig`, `mate_start`,
`affected_start`, and `affected_end`. These CSVs are summaries, not lossless
structural archives: `from_csv` rejects them. Retain original VCF/JSON variants
and their RNA evidence; point-only and empty table columns are unchanged.

Symbolic VCF span records retain `POS` as the padding/junction coordinate;
`affected_start..affected_end` is the inclusive changed interval `POS+1..END`.
This applies to DEL, DUP, INV and CNV (including CN0/CN3). This parser requires
a nonempty span with `END > POS`. Direct `StructuralVariant(...)` construction
keeps its existing explicit-coordinate defaults: pass `affected_start`
separately when `start` includes padding. Paired breakends already supply it.

## Junctions: which side of each breakpoint is kept

An SV creates one or more *junctions*, each joining two breakpoints. At
each breakpoint the junction keeps either the bases to the left (up to
and including the position) or those to the right.
`StructuralVariant.junctions` lists them.

For a breakend record, the ALT says which sides are kept:

| ALT | This record keeps | Mate keeps |
|---|---|---|
| `t[p[` | left | right |
| `t]p]` | left | left |
| `]p]t` | right | left |
| `[p[t` | right | right |

For the other SV types the junctions follow from the type (`start` is the
base before the event, as in VCF):

| Type | Junctions |
|---|---|
| `DEL` | `start` (keeps left) joined to `end + 1` (keeps right) |
| `DUP` | `end` (keeps left) joined to `start + 1` (keeps right) |
| `INV` (symbolic) | two: `start`–`end` keeping left, and `start + 1`–`end + 1` keeping right |
| `INV` built by `pair_breakends` | only the junction its breakend pair observed |
| `INS`, `CNV`, single breakend | none |

A breakend with no readable ALT has unknown sides, except that
`mate_orientation` `"[["` / `"]]"` still gives the mate's side (right /
left).

## 5' and 3' roles

Keeping a side of a transcript keeps either its 5' end or its 3' end,
depending on its strand:

| Strand | Keeps left | Keeps right |
|---|---|---|
| forward (`+`) | its 5' end → **5' partner** | its 3' end → **3' partner** |
| reverse (`−`) | its 3' end → **3' partner** | its 5' end → **5' partner** |

A junction makes a `GeneFusion` for a transcript when all of these hold:

1. The transcript contains one end of the junction.
2. Its gene doesn't also contain the other end (that would be intragenic).
3. At the other end there's a protein-coding transcript in a different
   gene, whose gene doesn't contain this end either.
4. The two roles are opposite: one 5' partner, one 3' partner. Two 5'
   ends (head to head) or two 3' ends (tail to tail) can't be read
   through, so they don't fuse.

The partner is the first transcript at the other end that satisfies
these, in pyensembl's order; varcode doesn't rank isoforms
([#406](https://github.com/openvax/varcode/issues/406)).

`GeneFusion.transcript` is the transcript being annotated, which can be
either partner. `five_prime_transcript` and `three_prime_transcript` say
which is which, and `partner_transcript` is the other one.

## Breakends: every combination

Effects on the transcripts at the record's own breakpoint:

| This breakpoint | Mate | Effect |
|---|---|---|
| Intergenic | anything | `Intergenic` |
| In a coding gene | intergenic | `TranslocationToIntergenic` |
| In a coding gene | in a non-coding gene only (e.g. MALAT1) | `TranslocationToIntergenic` |
| In a coding gene | in a coding gene, opposite roles | `GeneFusion` |
| In a coding gene | in a coding gene, head to head or tail to tail | `TranslocationToIntergenic` |
| In a coding gene | in the same gene | `TranslocationToIntergenic` |
| In a coding gene | single breakend (no mate) | `TranslocationToIntergenic` |
| In a coding gene | ALT isn't a breakend | `GeneFusion` treating this transcript as the 5' partner, with a warning |

`TranslocationToIntergenic` therefore means "a breakend that doesn't form
a gene fusion", not only "the mate is intergenic".

Its reference-derived `mutant_transcript` contains only the retained local
5' prefix or 3' suffix, including the base at an exonic breakpoint. Intronic
breakpoints retain the corresponding annotated exons. `A.` keeps genomic left;
`.A` keeps genomic right, even without a mate. Transcript cDNA is already
strand-oriented, so the segment itself has `strand="+"` on either gene strand.
The fragment is labeled `evidence["sequence_status"] = "retained_reference_fragment"`;
full `cdna_sequence` and `mutant_protein_sequence` remain None. Concatenating this
one fragment does not establish the full allele or prove a mature RNA/protein.
Unknown local orientation (including mate-orientation metadata alone) leaves
`mutant_transcript` None. A supplied `alt_assembly` still takes precedence.

**Intergenic ↔ gene.** The record at the intergenic end reports
`Intergenic`; the gene is annotated only from its own record, which
reports `TranslocationToIntergenic`. varcode doesn't model an intergenic
promoter or enhancer driving a gene.

**Same gene.** A breakend pair with both ends in one gene is intragenic,
so it isn't a fusion. To get `LargeDeletion` / `LargeDuplication` /
`Inversion` instead, pair the records with `pair_breakends` so the
caller's `SVTYPE` types the event.

**Direction.** Both records of a junction report the same fusion
direction, each on its own gene. For CFTR (chr7:117,485,000, `+`, keeping
left) joined to BRCA1 (chr17:43,120,000, `−`, keeping left):

| Record | Annotated on | Result |
|---|---|---|
| `N]17:43120000]` at chr7:117,485,000 | CFTR | `GeneFusion`, 5' = CFTR, 3' = BRCA1 |
| `N]7:117485000]` at chr17:43,120,000 | BRCA1 | `GeneFusion`, 5' = CFTR, 3' = BRCA1 |
| `N[17:43120000[` at chr7:117,485,000 | CFTR | `TranslocationToIntergenic` (head to head) |
| `]17:43120000]N` at chr7:117,485,000 | CFTR | `TranslocationToIntergenic` (tail to tail) |

The two records' fused proteins can differ (712 aa vs 1,854 aa here),
because each pairs its own transcript with the first matching isoform of
the other gene.

## Deletions, duplications and inversions

Whether a span fuses the genes at its two ends depends on their strands:

| Type | Genes on the same strand | Genes on opposite strands |
|---|---|---|
| `DEL` | fusion: the gene upstream in transcription is 5' | no fusion |
| `DUP` | fusion: the gene downstream in transcription is 5' | no fusion |
| `INV` (symbolic) | no fusion | two reciprocal fusions; each gene reports the one driven by its own promoter |

Examples from chr7, with CFTR (`+`), LSM8 (`+`, downstream) and CTTNBP2
(`−`, downstream):

| Event | On CFTR | On the other gene |
|---|---|---|
| `DEL` CFTR..LSM8 | `GeneFusion` CFTR → LSM8 | `GeneFusion` CFTR → LSM8 |
| `DUP` CFTR..LSM8 | `GeneFusion` LSM8 → CFTR | `GeneFusion` LSM8 → CFTR |
| `INV` CFTR..LSM8 | `Inversion` | `Inversion` |
| `DEL` CFTR..CTTNBP2 | `LargeDeletion` | `LargeDeletion` |
| `DUP` CFTR..CTTNBP2 | `LargeDuplication` | `LargeDuplication` |
| `INV` CFTR..CTTNBP2 | `GeneFusion` CFTR → CTTNBP2 | `GeneFusion` CTTNBP2 → CFTR |

On the reverse strand the same rules give TMPRSS2 → ERG for the chr21
deletion behind that prostate-cancer fusion.

Per transcript:

| Transcript | Effect |
|---|---|
| Wholly inside the span | `LargeDeletion` / `LargeDuplication` / `Inversion` |
| Holds both ends | the same, or `Intronic` if no exon overlaps |
| Holds one end, fusion conditions met | `GeneFusion`, with the span's effect as a further candidate |
| Holds one end, no partner | the span's effect (e.g. a truncating `LargeDeletion`) |

When a boundary transcript gets a `GeneFusion`, what the span does to its
exons isn't lost:

```python
from varcode.effects import LargeDeletion

fusion = deletion.effect_on_transcript(tmprss2)
(span_effect,) = [c.effect for c in fusion.candidates
                  if isinstance(c.effect, LargeDeletion)]
span_effect.affected_exons
```

`INS` and `CNV` never fuse; they report `LargeDuplication` when they
overlap exons.

Local DUP/INV transcript models require every junction end to lie in the
selected transcript. An event extending outside it (including a span enclosing
the entire transcript) does not establish duplicated or inverted cDNA within
that transcript. Without a fusion partner or supplied `alt_assembly`, the
DNA-level `LargeDuplication`/`Inversion` remains but `mutant_transcript` is None.
This means unresolved sequence, not an unchanged protein or a proven truncation.
The same restriction applies to a span candidate attached to a fusion; it does
not remove the fusion's own transcript model. A longer isoform of the same gene
does not establish the shorter isoform's structure.

## Where in the gene a breakpoint lands

The fused cDNA is the 5' partner's cDNA up to its breakpoint followed by
the 3' partner's from its breakpoint on.

| Breakpoint | What's kept |
|---|---|
| In an intron | The 5' side ends with the last complete exon before it; the 3' side starts with the first exon after it. For CFTR intron 1 that's exon 1's 185 bases. |
| In an exon | The base at the breakpoint stays on the side that keeps it, so the cut is mid-exon: 51 bases for a breakpoint 50 bases into CFTR exon 1. |

The protein is translated from the 5' partner's start codon, through the
junction, to the first stop. That gives these cases, all still reported
as `GeneFusion`:

| 5' partner breakpoint | `mutant_transcript.mutant_protein_sequence` |
|---|---|
| Before its start codon (5' UTR or an early intron) | `None`. BRCA1's start codon is in exon 2, so a break in intron 1 loses it. |
| Within its coding sequence, in frame with the 3' partner | 5' partner's N-terminus plus the 3' partner's C-terminus (CPEB2::FAM193A, 1,500 aa) |
| Within its coding sequence, out of frame | 5' partner's N-terminus, then the 3' partner read in a shifted frame to the first stop (OTX1::KIF3C, 35 aa) |
| After its stop codon (3' UTR) | the 5' partner's unchanged protein (CFTR, 1,480 aa) |

There is no separate in-frame flag; compare the protein with the partners'
reference proteins. A 5' partner transcript without a complete coding
sequence also gives `None`.

Breakpoints near exon boundaries add candidates to any SV effect:

- **Splice outcomes**, when the SV's `start` or `end` is within 6 bp of an
  exon edge (`source="varcode_splice"`).
- **Cryptic exons**, scanned in ±500 bp windows around every breakpoint
  (`source="varcode_motif"`), from a chromosome FASTA when one is attached.

## Pairing breakend records

`pair_breakends` joins the two records of a breakend pair:

- If both carry the same `SVTYPE` of `DEL`, `DUP` or `INV` and their kept
  sides fit it, the result is that typed event, annotated over its whole
  span as above.
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

## Importing observed RNA structures

Use an existing DNA `StructuralVariant` as the anchor, so an RNA splice junction
is not mistaken for a genomic breakpoint. Multiple observed RNA models can be
attached to that same variant:

```python
from varcode import load_exacto_fusions

rna = load_exacto_fusions(
    "transcript_structures.tsv", "integrated_variants.tsv",
    variants_by_id={"42": structural_variant},
)
effects = variants.effects(rna_resolver=rna)
for candidate in rna.candidates:
    sequence = candidate.effect.mutant_transcript.cdna_sequence
    model_id = candidate.evidence["transcript_model_id"]
```

`variants_by_id` maps **Exacto DNA call IDs**, not assumed VCF record IDs, to
already loaded variants. Only these IDs are selected from the integration file.
The model key includes both `transcript_model_id` and the reference transcript
ID set. Different models/reference interpretations and their DNA links are kept;
duplicate integration rows do not duplicate a candidate.

The adapter follows Exacto's `index` ordering and inclusive `read_start` /
`read_end` base-row coordinates. Sequences already run 5' to 3' even for reverse
strand genes: they are never reverse-complemented a second time. Event rows
preserve splice structure but do not invent missing bases. The complete original
structure and integration rows are in `candidate.evidence` and the mutant
transcript's evidence. Row counts are **not** supporting-read counts.

This is a deliberately limited adapter for linear, two-locus, SV-linked models
with an annotated sense 5' anchor. A missing or antisense 3' annotation gives
`TranslocationToIntergenic`, whose existing meaning includes non-sense joins.
Unknown transcript IDs and incomplete/unsupported structures raise errors.
This does not implement all Exacto DNA/RNA variant formats, circular RNA,
or multi-gene paths. Those remain in the broader
[#260](https://github.com/openvax/varcode/issues/260) roadmap.

For other producers, use the same existing effect classes directly:

```python
from varcode import RNAEvidence, make_fusion_outcome

candidate = make_fusion_outcome(
    structural_variant, five_prime_transcript,
    partner_transcript=three_prime_transcript,  # only when sense-oriented
    sequence=observed_sequence,               # already 5' to 3'
    transcript_model_id="assembled-model-7",
    source="rna_assembler", read_count=6,
    extra_evidence={"sequence_status": "junction_fragment"},
)
effects = variants.effects(rna_resolver=RNAEvidence([candidate]))
```

Do not supply a sense partner merely because a gene name appears near a
breakpoint. Leave `partner_transcript=None` for unresolved/intergenic/antisense
partners and record the observed loci/orientation in `extra_evidence`.
The importer preserves the supplied sequence, including insertions, without
appending reference exons. Read counts do not become probabilities or change
default candidate ranking. DNA predictions and RNA observations coexist.

Protein sequence is `None` by default. A caller can pass `cds_start=<0-based
offset>` to `make_fusion_outcome`, or a `cds_starts` mapping keyed by
`("model-id", ("sorted-reference-id-1", "sorted-reference-id-2"))` to the Exacto
loader. This requires a valid unambiguous start-to-stop ORF in the supplied
sequence. It produces a **sequence-predicted protein**, not evidence of
translation, full-length RNA, tumor-cell identity, or a mature two-gene fusion.

### Import Exacto's protein predictions

When Exacto has already selected reading frames, supply its native
primary-structures table instead of choosing `cds_starts` yourself:

```python
rna = load_exacto_fusions(
    "transcript_structures.tsv", "integrated_variants.tsv",
    variants_by_id=variants_by_id,
    primary_structures_path="primary_structures.tsv",
)
for candidate in rna.candidates:
    print(candidate.evidence.get("exacto_peptide_id"))
    print(candidate.effect.mutant_protein_sequence)
    print(candidate.evidence.get("protein_completeness"))
```

Each `(model, reference-transcript group, peptide_id)` remains separate. The
loader validates native `primary_structure_index`, `amino_acid_index`, and
`codon_index`, checks each nucleotide against its `transcript_structure_index`
and inclusive read coordinates, and checks the amino acid against its codon.
The amino acid repeated on all three base rows is emitted once. Terminal `*`
is retained in provenance but omitted from `mutant_protein_sequence`.

Partial peptides are available, with `protein_completeness` set to
`partial_start`, `partial_end`, or `partial_both`; start-to-stop predictions
are labeled `start_to_stop`. Trailing incomplete codons cannot claim an amino
acid. Missing peptide rows leave that model's protein `None`. Per-base variant
IDs, frameshift state and all original fields remain in `exacto_primary_structure`.
Exacto's standard genetic-code choice is recorded as `protein_translation_table=1`.

These sequences are **predictions from RNA**, not evidence of translation.
A complete start-to-stop ORF need not cross the rearrangement or establish a
full-length fusion transcript. Downstream users must check completeness and
source coordinates rather than treating every protein string as a complete
expressed fusion protein.

The small osteosarc regression fixtures deliberately exercise the negative
case: GABBR1 joins sequence upstream of SLC29A1, OTUD7A joins an antisense FMN1
intron, and the KLF15-side reads are intronic. Their RNA junctions are retained
without manufacturing coding fusions or protein sequences.

Format reference: [Exacto's structure translation source, pinned revision
307c086](https://github.com/pirl-unc/exacto/blob/307c08670d5e706734bddf393bcebc84db497f9f/exacto/exacto-translator/src/algorithms/translation.rs).

## Limitations

- The partner isoform isn't ranked
  ([#406](https://github.com/openvax/varcode/issues/406)).
- Chains of several SVs aren't assembled into one allele, and regulatory
  effects (promoter or enhancer hijacking) aren't modeled.
- Annotating multi-megabase spans can be slow
  ([#407](https://github.com/openvax/varcode/issues/407)).
