# Splice variants

Variants near an exon–intron boundary can change how a transcript is spliced,
for example by skipping an exon or keeping an intron. DNA alone doesn't show
which will happen, so Varcode returns a `SpliceOutcomeSet` listing the possible
outcomes, each with its own predicted protein where it can be determined.

## Example

This CFTR variant affects an exonic splice position in Ensembl 81 (GRCh38).
See [reference setup](getting_started.md#reference-data) if needed.

```python
from varcode import Variant

variant = Variant("7", 117_531_114, "G", "T", genome=81)
transcript = variant.genome.transcript_by_id("ENST00000003084")
splice_set = variant.effect_on_transcript(transcript)
print(type(splice_set).__name__)  # SpliceOutcomeSet
print(splice_set.disrupted_signal_class.__name__)  # ExonicSpliceSite
```

No opt-in flag is needed. `disrupted_signal_class` identifies the affected
signal; `candidates` describes what might happen to the transcript.

## Candidate outcomes

```python
for candidate in splice_set.candidates:
    print(type(candidate.effect).__name__)
    print(candidate.effect.mutant_protein_sequence)
    print(candidate.source, candidate.evidence)
```

Candidates can include `NormalSplicing`, `ExonSkipping`, `IntronRetention`,
`CrypticDonor`, and `CrypticAcceptor`. Their sequences and amino-acid changes
are available when the required sequence can be resolved; otherwise the
fields may be `None`. Intron retention and cryptic-site predictions can need
genomic sequence beyond the annotated transcript.

To inspect mechanism-specific fields:

```python
from varcode import ExonSkipping, IntronRetention

for candidate in splice_set.candidates:
    mechanism = candidate.effect
    if isinstance(mechanism, ExonSkipping):
        print(mechanism.affected_exon.exon_id, mechanism.in_frame)
    elif isinstance(mechanism, IntronRetention):
        print(mechanism.side, mechanism.retained_intron_start)
```

Each mechanism's `splice_signal` points to the underlying site-disruption
effect. Only `NormalSplicing` is built eagerly; the other candidates are
created and cached when `candidates` is first accessed.

## Coding effects when splicing is unchanged

An exonic variant can affect a codon and a splice signal at the same time.
Exons also contain untranslated sequence, so exonic does not always mean coding
([Ensembl glossary](https://www.ensembl.org/Help/Glossary?id=521)).

```python
coding = splice_set.effect_if_splicing_unchanged
if coding is not None:
    print(coding.short_description)
```

This returns the `NormalSplicing` candidate's coding effect when one applies.
It is `None` for a purely intronic disruption with no directly changed coding
base. `alternate_effect` is a compatibility alias.

## Ranking and protein sequences

| Accessor | Returns |
|---|---|
| `most_likely_effect` | First effect in the producer's candidate order |
| `highest_priority_effect` | Effect with the highest consequence priority |
| `most_likely_candidate`, `highest_priority_candidate` | The same choices with their provenance wrappers |

These answer different questions: the producer's preferred mechanism may not
be the most disruptive. Neither ordering establishes that a mechanism occurred,
and candidate order is not a calibrated probability.

```python
preferred = splice_set.most_likely_effect
most_disruptive = splice_set.highest_priority_effect
proteins = splice_set.mutant_protein_sequences
```

`mutant_protein_sequences` contains distinct nonempty sequences. The separate
`candidate_proteins` mapping is keyed by mechanism class and uses an empty
string for unresolved protein sequence. To retain the evidence for each
sequence, iterate the candidates instead of using either summary.

## Splice signals

The default classifier uses a positional window around each exon–intron
boundary, not a learned splice score:

| Signal class | Position |
|---|---|
| `ExonicSpliceSite` | Last 3 exonic bases at a donor, or the first exonic base at an acceptor |
| `SpliceDonor` | Canonical `GT` at intronic +1/+2 |
| `SpliceAcceptor` | Canonical `AG` at intronic −2/−1 |
| `IntronicSpliceSite` | Intronic +3…+6 at a donor or −3 at an acceptor; also the canonical positions when the supplied reference signal is not canonical |

These classes identify the site inside the outcome set; they are not separate
top-level predictions. Coding variants outside this window receive their coding
effect without a splice annotation.

## RNA evidence

`SpliceOutcomeSet.with_rna_evidence(...)` returns a new set whose current
candidates are the RNA-observed mechanisms. The original predictions remain
in `dna_candidates`; `rna_evidence`, `excluded_candidates`, `added_candidates`,
and `candidate_rna_evidence` retain the comparison. Use
`splice_set.rna_evidence_for(candidate)` to inspect support for one candidate.

Candidate evidence can include named measurements such as `read_count`,
`junction_id`, `psi`, or `motif_score`. Varcode preserves these as provenance;
it does not normalize them into a probability. The ordinary `EffectCandidate`
wrapper has no `probability` or `plausibility` field.

See [germline, phase, and RNA inputs](phasing.md#combining-evidence)
for collection-level use and the [RNA evidence API](api_rna.md#rna-evidence)
for integration details.

## Limitations

A variant outside the positional window can still affect splicing biologically.
The default does not flag mid-exon enhancer/silencer disruption, branch-point
changes, or deep-intronic cryptic-site activation. Broader sequence-based
prediction is tracked in [#297](https://github.com/openvax/varcode/issues/297).
RNA evidence can help distinguish predicted alternatives, but missing evidence
does not establish normal splicing.

For fields and methods, see the [SpliceOutcomeSet reference](effect_types.md#splice-outcome-container).
