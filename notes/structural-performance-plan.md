# Structural annotation performance (#407)

## Scope

Reproduce the reported ESVEE breakend and large-span slowdowns on current
Varcode using the existing public Osteosarc callsets and Ensembl 115 cache.
Keep this investigation and any implementation changes within Varcode.

## Method

1. Record exact source records, versions, cold/warm timings, and Python/SQL
   profiles for the named GNAS/CNOT3 breakends and a large typed event.
2. Identify repeated or unnecessarily broad work. Preserve transcript and
   candidate coverage, biological interpretations, and uncertainty semantics.
3. If the bottleneck is in Varcode, implement the smallest measured correction
   with output-parity and behavioral regression tests. Report dependency
   bottlenecks separately for the owning agent.
4. Record before/after evidence. For a code fix, bump the version, update release
   notes, pass lint/full tests and CI, merge, and deploy from clean main.

## Acceptance

Explain the observed cost with a reproducible profile. Any performance fix must
preserve complete effect output on the benchmark inputs and avoid a new global
unbounded cache or silently reduced annotation scope.

## Findings and implementation

The original minute-long cases do not reproduce on current main (`cb49a63`,
10.17.2), PyEnsembl 2.28.1, Datacache 1.29.0, Python 3.12.6. Earlier commit
`07873df` already reduced the reference-range fallback to one locus query and
cached variant-level mate lookups. The remaining hotspot is repeated
`spliced_offset` calls inside that range: every base is mapped through every
overlapping transcript, repeatedly walking the same exons.

A coverage-only early rejection was insufficient: some breakpoint windows are
fully covered by complementary transcript exons. The final implementation reads
each intersecting exon as a contiguous slice, combining answers in the same
transcript order and retaining disagreement diagnostics. Uncovered unions return
empty immediately; missing/short sequences and adapters without usable exon
intervals retain the existing all-or-nothing behavior and fallback. Endpoint
offsets now reject intron-spanning slices before loading cDNA.

Coordinate and sequence contracts were checked against the installed PyEnsembl
`Transcript.exon_intervals`, `spliced_offset`, and `sequence` implementations and
the [primary API documentation](https://pyensembl.readthedocs.io/en/latest/pyensembl.html#module-pyensembl.transcript).
No biological prediction or candidate-selection rule changes.

### Measurements

Same local ESVEE files, Ensembl 115, no chromosome FASTA; one unprofiled run per
checkout. Cache state and machine load affect timing, particularly the first
case, which loads transcript/protein dictionaries. These are illustrative
measurements, not portable timing assertions.

| Input | Before (s) | After (s) | Transcript effects |
| --- | ---: | ---: | ---: |
| CNOT3 19:54153839 → 54155308 (first, dictionaries unloaded) | 3.030 | 1.609 | 24 |
| CNOT3 19:54155308 → 54153839 | 0.355 | 0.042 | 24 |
| 20:36148058 → 36148629 | 0.028 | 0.013 | 51 |
| GNAS 20:58909683 → 58909577 | 0.692 | 0.039 | 63 |
| INV 2:25955847–63053516 (37 Mb) | 1.148 | 1.131 | 6,373 |

All five complete serialized effect graphs compare equal after removing only
the annotation timestamp and Varcode version. This includes structural primary,
cryptic and splice candidates, model sequences, and evidence. The large inversion
continues to annotate every transcript in its span; its cost is proportional to
that output, not a reason to silently switch to breakpoint-only annotation.

The GNAS profile falls from 100,386 to 904 `spliced_offset` calls. Its reference
lookup falls from 2.135 s to 0.055 s under cProfile. Regression tests check work
counts instead of wall time: the uncovered-flank case drops from 40,060 offset
calls to 40 on each strand. The old implementation fails that test and the test
that forbids loading cDNA for an already-proven intron-spanning range.

### Reproduction

`benchmark_structural_407.py` accepts the public T1 organoid and T2 tumor ESVEE
somatic VCFs from [Osteosarc](https://osteosarc.com/data/), installed human Ensembl
115, and an output directory. It records exact source SHA256s, dependency
versions, repeated timings, and complete effect graphs. `--profile` also writes
cProfile dumps and summaries. Compare the generated JSON graphs across
checkouts; run benchmarks separately from the test suite to limit interference.

```sh
PYTHONPATH=/path/to/checkout python notes/benchmark_structural_407.py \
    IPISRC044_tumor_T1_organoid_ucla.esvee.somatic.vcf.gz \
    IPISRC044_tumor_T2_ucla.esvee.somatic.vcf.gz /tmp/varcode-407-results --profile
```
