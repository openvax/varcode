# Structural variant annotation

This page describes what varcode reports for each kind of structural
variant (SV): which transcripts get annotated, how a fusion's direction is
decided, what depends on strand and on where a breakpoint lands, and which
effect class comes back in each case. Examples use Ensembl 95.

Load SVs with `load_vcf(..., parse_structural_variants=True)`; see
[Effect annotation](effect_annotation.md#structural-variants) for loading
and [Transforms](transforms.md) for `pair_breakends`.

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
fusion = annotator.annotate_on_transcript(deletion, tmprss2)
(span_effect,) = [c.effect for c in fusion.candidates
                  if isinstance(c.effect, LargeDeletion)]
span_effect.affected_exons
```

`INS` and `CNV` never fuse; they report `LargeDuplication` when they
overlap exons.

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

## Limitations

- The partner isoform isn't ranked
  ([#406](https://github.com/openvax/varcode/issues/406)).
- A `DUP` or `INV` with one end in a transcript and no partner builds its
  mutant transcript as if the whole event were inside
  ([#405](https://github.com/openvax/varcode/issues/405)).
- Symbolic `<DEL>` / `<DUP>` count the VCF padding base as part of the
  event ([#404](https://github.com/openvax/varcode/issues/404)).
- Chains of several SVs aren't assembled into one allele, and regulatory
  effects (promoter or enhancer hijacking) aren't modeled.
- Annotating multi-megabase spans can be slow
  ([#407](https://github.com/openvax/varcode/issues/407)).
