# Compare SV calls across samples and callers

Different SV callers, and different samples, often report the same event in
different forms: a pair of reciprocal breakends or a symbolic `<DEL>`,
explicit or symbolic alleles, slightly different coordinates.
`varcode.sv_comparison` groups calls that report the same junctions, or nearby
ones, so you can see where callers agree and disagree. It compares what was
reported; it does not validate the calls.

```sh
python -m varcode.sv_comparison --calls calls.csv --output comparison --max-distance 100
```

`calls.csv` has one row per source ALT allele, with unique `call_id` and columns
`caller,sample,build,chrom,pos,ref,alt`. Coordinates are one-based VCF coordinates.
Use `GRCh37` or `GRCh38` for `build`; a trailing parenthesized note, such as
`GRCh38 (from caller header)`, is accepted and preserved. Optional `info` is the original VCF INFO text;
`end` supplies END when absent from INFO. Include `record_id`, `mate_id`,
`source_url`, `source_id`, filter status and `duplicate_export_of` to retain
source provenance. All additional columns survive unchanged in each member's
`raw` mapping. No PASS filter or record deduplication is applied.

Outputs, written to the `--output` directory:

- `groups.csv`: comparison group, exact representation IDs, samples, callers,
  original call IDs, per-breakpoint range/spread, insertion disagreement and
  sequence/imprecision status.
- `members.csv`: every input row, normalized oriented breakpoints and insertion
  sequence, exact ID, comparison group, or explicit unresolved reason.
- `comparison.json`: the same exhaustive data plus input SHA256, package version
  and distance threshold. The output directory must be new.

The Python API is `varcode.sv_comparison.compare_sv_calls(rows, max_distance=100)`.
It needs no reference sequence downloads. Synthetic unit fixtures are constructed
inline in `tests/test_sv_comparison.py`, including a complete CSV/CLI example;
there is no machine-specific fixture-generation prerequisite.

## What a group means

`exact_reported_allele` means identical reported oriented adjacencies and inserted
bases, with no reported positional uncertainty. Reciprocal BNDs and equivalent
explicit/symbolic deletions share that representation. Multi-junction events
must match all junctions: one BND cannot establish a complete inversion.
This is representation matching, not biological validation of a caller.

`same_reported_adjacency` retains uncertainty about inserted bases or coordinates.
`nearby_candidate` contains distinct exact IDs. Every corresponding breakpoint of
every pair must differ by at most the stated distance. Greedy complete-link
clustering, ordered by coordinates and stable identifiers, prevents an A–B–C
chain from grouping A with a distant C. There may be multiple valid partitions;
the deterministic partition is not a biological event-resolution claim. At equal
coordinates, different inserted sequences still have different exact IDs, and
`insertion_disagreement` is true. Original positions, confidence intervals and
alleles remain available for review.

Different assemblies, contigs, orientations or numbers of adjacencies never
share a group. Reference liftover, repeat-aware left-alignment and inferred
complex-event decomposition are not performed. Single breakends and unsupported
or malformed alleles remain `unresolved` singleton rows rather than disappearing.
They do not establish a complete adjacency even if their local positions agree.

Counts are source rows. Reciprocal records, repeat exports, caller dependencies,
and multiple samples do not establish independent molecular support. Use the
preserved source metadata when interpreting support. Assemblies come from the
input declarations; reference alleles are not checked against a genome FASTA.

The input contract and coordinate normalization follow the
[VCF specification](https://samtools.github.io/hts-specs/VCFv4.5.pdf).
