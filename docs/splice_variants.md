# Splice variants

Variants near an exon–intron boundary can change how a transcript is spliced.
Varcode predicts compatible mechanisms, realizes their transcript products,
and compares each with the patient baseline. DNA alone does not establish
which transcript is expressed.

## Example

This CFTR donor variant uses Ensembl 81 (GRCh38); see
[reference setup](getting_started.md#reference-data).

```python
from varcode import Variant

variant = Variant("7", 117_531_115, "G", "A", genome=81)
transcript = variant.genome.transcript_by_id("ENST00000003084")
effect = variant.effect_on_transcript(transcript)
print(type(effect).__name__)  # Deletion: the preferred exon-skip comparison
for candidate in effect.candidates:
    print(candidate.effect.short_description)
    print(candidate.effect.mutant_protein_sequence)
```

No opt-in flag is needed. The top effect is one classified comparison;
`candidates` retains the alternatives, including unresolved products. Read the
candidates when evaluating uncertainty; the top description is not an audit of
all hypotheses. JSON preserves the full products and evidence.

## Candidate outcomes

Mechanisms include `normal`, `exon_skip`, `intron_retention`, and `cryptic`.
The effect class describes the resulting protein change, not the mechanism:
exon skipping can produce a `Deletion` or a frameshift, for example.

```python
for candidate in effect.candidates:
    print(candidate.splice_mechanisms, candidate.novel_junctions)
    print(candidate.baseline_cdna, candidate.mutant_cdna)
```

Canonical and exon-skip paths can use transcript annotation. Retention and
cryptic paths need the relevant genomic sequence; absent sequence remains
`None` with an unresolved consequence. First and last exons cannot be skipped
by this internal-exon mechanism. Plans exclude incompatible combinations such
as retaining an intron while skipping an exon it requires.

## Coding effects when splicing is unchanged

An exonic variant can alter both a codon and a splice signal. To inspect the
conditional coding consequence, select candidates whose mutant plan uses only
`normal` splicing. With germline splice uncertainty, several patient baseline
plans may still give different comparisons.

```python
normal_effects = [
    candidate.effect
    for candidate in effect.candidates
    if candidate.splice_mechanisms == ("normal",)
]
```

## Ranking and protein sequences

`most_likely_candidate.effect` is the first ranked comparison;
`highest_priority_candidate.effect` is the most severe consequence. Mechanism
preference is ordinal, not a measured probability. Candidate `probability`
remains `None` when splice choices are uncalibrated. Iterate candidates to
retain each sequence's hypotheses and baseline.

Unchanged active splice sites share their mechanism between patient baseline
and mutant. If a site changes or disappears, all compatible baseline choices
remain. Neither a germline splice change nor an arbitrary preferred mechanism
establishes that the patient's normal transcript lacks an exon.

## Splice signals

The default uses positional windows, not a learned splice score:

| Region | Positions in transcript direction |
|---|---|
| Donor motif | Intronic +1/+2 |
| Acceptor motif | Intronic −2/−1 |
| Intronic splice region | Donor +3…+6 or acceptor −3 |
| Exonic splice region | Last three exon bases at a donor or first exon base at an acceptor |

Removed junction coordinates and changed motif bases produce stronger splice
hypotheses than peripheral changes. Unchanged annotated noncanonical sites,
including GC–AG and AT–AC introns, remain intact; these are established splice
classes ([Burset et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC113136/)).
An altered site still retains normal splicing as a possibility.

## RNA evidence

Supplied full RNA models that exactly match a predicted cDNA support every
compatible comparison and rank those candidates first. Other alternatives
remain. Partial or unmatched models retain their supplied sequence as an
unresolved comparison with an unknown baseline. Read counts do not become
probabilities; missing reads do not establish normal splicing.

Read collection, alignment and assembly belong upstream in Isovar or another
RNA tool. Varcode consumes models via `rna_resolver=` or a phase resolver; it
does not fill missing observed sequence with reference exons. See
[combining evidence](phasing.md#combining-evidence) and the
[RNA evidence API](api_rna.md#rna-evidence).

## Compatibility with Varcode 10

Explicit `annotator="fast"` retains `SpliceOutcomeSet`, its mechanism-specific
classes, `effect_if_splicing_unchanged`, and its legacy RNA reconciliation
contract. Existing imports and serialized sets remain supported. See the
[SpliceOutcomeSet reference](effect_types.md#splice-outcome-container) and
[migration guide](experimental_annotators.md#migrating-from-varcode-10).

## Bounded cryptic-site search

By default, `predict_transcript_model_effect` scans ±50 bp around a disrupted
canonical boundary and retains the single highest-scoring heuristic motif
(`CrypticSplicePolicy(scan_flank=50, max_sites_per_junction=1, min_score=0.0)`).
Experimentally observed cryptic sites outside that window or below a competing
heuristic score — such as the mouse *Ap2s1* `del17` donor 51 nt upstream of
exon 2 ([PMID 29479578](https://pmc.ncbi.nlm.nih.gov/articles/PMC5824975/)) —
are not enumerated under the default bounds.

Pass an explicit `cryptic_policy=CrypticSplicePolicy(scan_flank=..., max_sites_per_junction=...)`
to enumerate multiple distinct cryptic boundaries per disrupted site while
keeping the combined search bounded by `max_hypotheses` and preserving
baseline/mutant coupling on shared disrupted sites. All DNA-enumerated cryptic
candidates are tagged `discovery_status="heuristic_dna_motif"` and are never
labeled experimentally observed; RNA-supported products continue to work
independently of DNA hypothesis coverage.

## Limitations

Variants outside the windows can still affect splicing. Mid-exon regulatory
sites, branch points and deep-intronic activation require broader modeling,
tracked in [#297](https://github.com/openvax/varcode/issues/297). Mechanism and
protein predictions alone do not prove expression or translation.
