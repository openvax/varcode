# SV coordinates and fusion rules

For loading and result access, start with [Structural variants](structural_variants.md).
This reference describes how the current implementation handles junctions and
transcript boundaries. The biological examples use Ensembl 95.

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

## Junction orientation

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

## Fusion partners

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

## Breakend outcomes

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

## Breakpoint position and protein sequence

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
