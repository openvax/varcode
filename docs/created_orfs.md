# Discovering newly created ORF starts

`discover_created_orfs` searches **paired RNA products** for newly created start
codons. It is opt-in and leaves the ordinary transcript consequences unchanged.
Use `extra_orfs=` to annotate supplied reference ORFs, including ORFs on noncoding
host transcripts; discovery answers the separate question of a newly created
initiation hypothesis.

```python
from pyensembl import cached_release
from varcode import Variant, ORFDiscoveryPolicy, discover_created_orfs

genome = cached_release(95)
transcript = genome.transcript_by_id("ENST00000304494")
variant = Variant("9", 21974861, "C", "A", genome=genome)
effect = variant.effect_on_transcript(transcript)
report = discover_created_orfs(effect, policy=ORFDiscoveryPolicy(
    start_codons=("ATG",), max_scan_bases=100_000,
    max_protein_length=1_000, max_candidates=32, max_products=256,
))
for assessment in report.assessments:
    print(assessment.status, assessment.reasons)
for candidate in report.effects:
    print(candidate.transcript_id, candidate.discovered_orf.start,
          candidate.mutant_protein_sequence)
```

This is the published **CDKN2A c.-34G>T** allele (genomic forward-strand C>A on
GRCh38). It creates an upstream AUG. Its conditional peptide is
`MGRAGSGRRGAAWSRRRGAAWSLRLTGWPRPRPGVG`, independently reproduced from
**NM_000077.4**. The Ensembl and RefSeq RNAs have different 5' extents but the same
new coding sequence. [Liu et al., 1999](https://doi.org/10.1038/5082) established
the initiation defect experimentally; this function does not infer the measured
reduction in canonical translation or treat the entire predicted peptide as an
observed patient peptide.

## Common effect fields and the new model

Each candidate is a `varcode.effects.StartGain` with the normal host identity and
sequence accessors. The focal fields distinguish a new initiation hypothesis:

| Field | Meaning |
|---|---|
| `reference_transcript`, `transcript_id`, `gene_id` | Original host RNA and gene |
| `baseline_transcript`, `mutant_transcript` | The actual paired RNA sequences |
| `original_protein_sequence` | `None`: baseline initiation in this frame is unassessed |
| `mutant_protein_sequence` | Conditional peptide, possibly a marked prefix |
| `modifies_protein_sequence` | `None`: the candidate does not establish translation |
| `orf` | `None`: discovery has not supplied a reference ORF |
| `discovered_orf` | Mutant start/end, origins, genetic code, termination and context |
| `hypothesis`, `source_evidence`, `discovery_policy` | Phase, provenance and search assumptions |

Offsets and intervals on the mutant RNA are zero-based and half-open. An `end`
includes the stop codon when one is present. `termination` is `observed_stop`
(a stop in the sequence under the selected genetic code), `sequence_boundary`,
`protein_length_limit`, or `ambiguous_sequence`. The latter three retain a
conditional peptide prefix and do not claim a complete protein.

Kozak context records available flanking sequence. It is neither a probability
nor an exclusion threshold. Even a start at supplied offset zero has unknown
biological 5' completeness. Accepted non-AUG codons can be requested explicitly;
conditional initiation supplies methionine, while subsequent codons use the
reference genetic code. Known recoded transcripts require an explicit model.

## Scope, limits and mapped RNA

`report.complete` means every supplied product was searched under the declared
policy. It does not mean the transcriptome or every possible ORF was searched.
Inspect every assessment, including those with no candidates:

- `complete`: every triplet in this supplied RNA was assessed under this policy.
- `bounded`: a scan, candidate or product limit left part of the search unassessed.
- `partial`: sequence or ancestry ambiguity prevents a complete assessment.
- `unavailable`: required paired sequence, ancestry or host reference is missing.

Existing transcript-model phase and splice hypotheses are retained separately.
An unchanged start shifted by an insertion is excluded using its mapped ancestry.
A new junction or newly included RNA segment can create a start candidate.
Neither hypothesis ranking nor phase mass is a translation probability.

`discover_created_orfs_in_products(baseline, mutant, ...)` accepts two
`RealizedTranscriptProduct` objects, including mapped products of a noncoding
host. Both must identify the same host/reference. Their `evidence.base_origins`
contains one `(contig, one_based_position, strand)` or `None` per nucleotide.
The ordinary engine's `genomic_positions` are also accepted. `None` is unknown
unless its offset appears in `inserted_base_offsets`. If both products have
inserted/unmapped sequence, new adjacency cannot always be established; the
report retains that ambiguity. No sequence alignment guesses missing origins.

The initial adapter requires existing paired products. Legacy `fast` results,
ordinary noncoding results without materialized RNA, partial observed RNAs,
unmapped nonlocal structural products and hosts without a Genome return explicit
unavailable assessments. Supplied reference ORFs still use the existing ORF
annotation engine. This bounded addition leaves the broader partial-product and
discovery work in [#583](https://github.com/openvax/varcode/issues/583) open.

## Export

```python
from varcode import EffectCollection, ORFDiscovery, to_proteome

restored = ORFDiscovery.from_json(report.to_json())
effects = EffectCollection(report.effects)
effects.to_csv("created-orfs.csv")
proteome = to_proteome(effects)
proteome.save("created-orfs.json.gz")
```

Native JSON/pickle and CSV archives retain the candidate model, paired RNA,
policy, hypothesis and pinned reference. CSV requires its `effect_json` column;
canonical reannotation cannot reconstruct a discovered frame. Proteome entries
are labeled `nc_orf` / `variant_created` and carry start/stop completeness.
Preserve the full discovery report as well: exporting only its effects omits
empty and unavailable assessments.
