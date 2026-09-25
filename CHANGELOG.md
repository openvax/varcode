# Change Log

## [v10.4.2](https://github.com/openvax/varcode/tree/v10.4.2) (2026-09-25)

- Allow osteosarc `>=0.2.3,<0.3` in the `test-data` extra, the range Isovar and
  Vaxrank also require. The 0.7.0 pin in 10.3.1 made a development environment
  with Varcode's test data and Isovar impossible to install. The snapshot check
  accepts any version in the range, so osteosarc patch releases no longer
  break it.

## [v10.4.1](https://github.com/openvax/varcode/tree/v10.4.1) (2026-09-25)

- Ranking a structural candidate set that lists itself as a candidate, for
  example through attached RNA evidence, no longer recurses without end. The
  self-candidate ranks as the set's own class. Found while running Vaxrank on
  Varcode 10.

## [v10.4.0](https://github.com/openvax/varcode/tree/v10.4.0) (2026-09-25)

- `MolecularPhaseResolver` no longer reports trans from a source's partner
  lists (#517). Variants that are not co-observed are unknown (`None`) unless
  the source implements its own `in_cis` from fragments that show one alt
  allele with the other's reference allele. Previously a germline variant
  outside an RNA source's input, which the source never examined, was
  reported as trans and collapsed germline-aware annotation to the
  reference-relative result. Co-observation still establishes cis, so
  haplotype grouping is unchanged.

## [v10.3.1](https://github.com/openvax/varcode/tree/v10.3.1) (2026-09-25)

- Pin the optional `test-data` extra to osteosarc 0.7.0, the current release.
  The bundled snapshot checks pass unchanged on it; the 0.2.3 pin made them
  fail in any environment with a newer osteosarc.

## [v10.3.0](https://github.com/openvax/varcode/tree/v10.3.0) (2026-09-24)

- Unenumerated phase no longer becomes a precise all-cis consequence (#503).
  When germline phase needs more hypotheses than the cap, the default path and
  the experimental `transcript_model` return a new `HypothesisLimit` effect (a
  kind of `Unresolved`). It records the phase partition and the effect without
  germline context, ranks as that effect would in top-priority selection, and
  round-trips through JSON.
- New public helpers: `partition_germline_by_phase` returns a `PhasePartition`
  whose `hypotheses()` enumerates only what is unknown. Homozygous germline
  variants are always cis, phase follows through germline-to-germline answers,
  and variants phased to each other flip as one block, so partial evidence is
  kept and impossible combinations are not generated. `query_in_cis` gives
  every caller one policy for asking a resolver (NumPy booleans accepted;
  resolver errors logged and treated as unknown). Exceeding the cap raises a
  picklable `HypothesisLimitError`; `enumerate_phase_hypotheses` now raises it
  instead of returning an all-cis placeholder.
- The phase cap is `GermlineContext.max_phase_hypotheses` (default 8), so
  `effects(germline=...)` can raise it. Caps must be positive integers (NumPy
  integers accepted). `transcript_model` uses the context's cap through
  `effects()`, and its combined phase/splice overflow now returns
  `HypothesisLimit` instead of raising.
- `PhaseCandidateSet` now ranks as its most severe candidate instead of below
  every other effect. Noncoding and incompletely annotated transcripts skip
  phase enumeration, and haplotype labels depend only on the cis set.

## [v10.2.1](https://github.com/openvax/varcode/tree/v10.2.1) (2026-09-24)

- Documentation readability pass: the README and docs home now describe what
  Varcode does and how to read its results (per-transcript predictions,
  severity vs. likelihood, `None` as unknown). Getting started gains a
  glossary of the main objects.
- Reorganize troubleshooting, sample identity checks, VCF export, SV results,
  and RNA imports around plain-language summaries, with status/threshold tables
  and collapsible notes for edge cases.
- Remove stale limitations already fixed in #437 and #449, note that
  genotype set operations treat uncalled normals as absent, and move the dated
  quality audit out of the user reference section.
## [v10.2.0](https://github.com/openvax/varcode/tree/v10.2.0) (2026-09-24)

- Add `varcode check-samples` and a Python API for pairwise donor genotype
  concordance and somatic tumor overlap (#510), with explicit sample roles,
  matched-normal evidence, expected identity flags, JSON and TSV reports.
- Report sparse or uninformative evidence as inconclusive; never infer reference
  genotypes from missing records. Tolerate tumor LOH and partial somatic overlap
  with documented, configurable screening heuristics and provenance/QC counts.
- Include reproducible synthetic VCF builders and checks for swaps, mixed
  tumor/normal inputs, missing quality, allele ordering and precise SV matching.

## [v10.1.3](https://github.com/openvax/varcode/tree/v10.1.3) (2026-09-24)

- Preserve sample identities, per-sample FORMAT values, missing samples, and
  original multi-allelic ALT order during VCF export (#502).
- Write INFO/FORMAT declarations, accept original field definitions, and reject
  incomplete/ambiguous ALT reconstructions before writing output.
- Support generator/empty exports, numeric position ordering, missing list
  elements, and per-sample filters; add independent identity round-trip tests.

## [v10.1.2](https://github.com/openvax/varcode/tree/v10.1.2) (2026-09-24)

- Enforce VCF allele limits across simple and structural calls and retain
  structural ALT indexes; propagate parser failures instead of silently
  truncating input (#504).
- Isolate default collection metadata (#505), preserve transcript failure
  diagnostics (#506), and run pytest with the interpreter used to probe
  optional plugins (#490).
- Enable undefined-name lint checks, remove mutable parser defaults and
  deprecated logging calls, and preserve column order during MAF normalization.
- Update developer/release instructions, fixture provenance guidance,
  Osteosarc compatibility, API navigation, and documented export/phase limits.

## [v10.1.1](https://github.com/openvax/varcode/tree/v10.1.1) (2026-09-24)

- Do not infer protein truncation from an internal-Met reference suffix when
  an observed transcript's 5′ completeness is unknown (#467). Preserve the
  start-to-stop ORF, its evidence and mapped local changes, while retaining
  full comparisons for predictions with a mapped annotated initiator.

## [v10.1.0](https://github.com/openvax/varcode/tree/v10.1.0) (2026-09-23)

- Report identical inherited alleles as `GermlineAlleleOverlap`, without
  manufacturing a new somatic sequence change or claiming LOH (#454).
- Add `detect_germline_overlap`; deprecate the allele-only `detect_loh` query,
  which now returns `None` (not assessed). LOH needs independent allelic-state
  evidence. Mixed inherited/somatic transcript-model groups remain unresolved.

## [v10.0.4](https://github.com/openvax/varcode/tree/v10.0.4) (2026-09-23)

- Leave overlapping DUP/INV events unresolved in the experimental transcript
  model when their junctions exceed its finite layout, instead of classifying
  a clipped allele as a UTR-only or unchanged product (#449). Preserve supplied
  assemblies and fully represented local rearrangements.

## [v10.0.3](https://github.com/openvax/varcode/tree/v10.0.3) (2026-09-23)

- Archive structural effects as versioned graphs so self candidates, fusion
  partners, primary/cryptic/splice/external candidates, mutant transcript
  models, and annotation provenance survive JSON round trips (#438).

## [v10.0.2](https://github.com/openvax/varcode/tree/v10.0.2) (2026-09-23)

- Load structural variants from both CLIs, normalize explicit chr contigs,
  report filtered VCF record counts, and add `--include-filtered` (#433).
- Add `--skip-errors` with retained Failure rows and error details, including
  failed initial annotation lookups and failures alongside coding filters.

## [v10.0.1](https://github.com/openvax/varcode/tree/v10.0.1) (2026-09-23)

- Reconstruct inserted junction bases in local DEL/DUP transcripts on both
  strands, including insertions at exon anchors flanking a deleted intron (#491).
  Protein consequences now use the net spliced edit with the retained insert.
- Share insertion-retention rules with fusion models: retain exonic inserts,
  exclude intronic inserts under reference splicing, and preserve uncertainty
  for mixed retention or conflicting/unreadable reciprocal alleles.

## [v10.0.0](https://github.com/openvax/varcode/tree/v10.0.0) (2026-09-23)

**Changed**
- Classify local SV transcript models by their protein consequence (#420):
  start loss, in-frame deletion/insertion, frameshift, and other changes now
  replace event-type labels. Preserve conditional consequences in candidate
  sets with splice assumptions, event provenance, and mutant sequences.
- Return unresolved candidates for unspecified insertion/CNV structure,
  unmaterialized inversions, and unmapped assemblies. Keep fusion consequences
  and legacy class imports; callers should read DNA type from `variant.sv_type`.
- Map retained CDS starts on both strands and preserve selenocysteine/SECIS
  uncertainty and independent coding/protein change flags.

## [v9.5.0](https://github.com/openvax/varcode/tree/v9.5.0) (2026-09-23)

**Added**
- Compare SV call tables across samples and callers with a reusable API and
  command. Preserve every input record, distinguish exact reported alleles
  from nearby candidates, and export breakpoint/insertion disagreements,
  source IDs, unresolved calls, and reproducible run provenance.

## [v9.4.2](https://github.com/openvax/varcode/tree/v9.4.2) (2026-09-23)

**Fixed**
- Retain inserted junction bases when translating exonic breakend fusions
  (#485), with strand-correct sequence from either reciprocal record and
  paired structural events. Mixed exon/intron insertion retention and
  conflicting or unreadable inserted alleles remain explicitly unresolved.
- Record the reference-splicing assumption when an insertion between two
  intronic breakpoints is excluded from the predicted fusion transcript.

## [v9.4.1](https://github.com/openvax/varcode/tree/v9.4.1) (2026-09-22)

**Testing**
- Bundle the verified historical Osteosarc metadata snapshot and run the
  native 0.1.4 adapter checks in CI on Python 3.10+. A fresh checkout with
  `.[test-data]` now runs these tests offline without a manually prepared
  cache or local sibling repository (#483).

## [v9.4.0](https://github.com/openvax/varcode/tree/v9.4.0) (2026-09-22)

**Fixed**
- Delegate known-cis groups to the selected annotator's optional
  `annotate_haplotype` method (#437). Retain unsupported groups as unresolved
  effects instead of silently dropping them or predicting with another backend.
- Route experimental transcript-model joint effects through its existing
  multi-variant engine, including patient germline context from every member's
  window and separately retained observed RNA models. The default retains
  point-edit haplotypes and explicitly declines joint germline composition.

**Testing**
- Update the optional `test-data` extra to Osteosarc 0.1.4, retaining Python 3.10+
  compatibility and the historical fixture's package/snapshot provenance (#478).

## [v9.3.7](https://github.com/openvax/varcode/tree/v9.3.7) (2026-09-21)

**Testing**
- Collect 177 real site variants and five unresolved entries from the pinned
  osteosarc snapshot into a portable fixture. Ordinary offline tests now
  annotate every ready allele with `fast` and `protein_diff`, including the
  mitochondrial variant and corrected MAP2 complex allele (#464).
- Add an explicit offline regeneration command and native test-variant loader,
  preserving original alleles, reference identity, correction notes, and
  source hashes. Optional snapshot checks verify byte-for-byte regeneration.

## [v9.3.6](https://github.com/openvax/varcode/tree/v9.3.6) (2026-09-19)

**Fixed**
- The `fast` and `transcript_model` annotators now classify an insertion
  immediately before the retained CDS start on a reverse-strand transcript
  as 5′ UTR. The insertion no longer produces a false coding- or
  protein-sequence-change flag; all three annotators agree (#474).

## [v9.3.5](https://github.com/openvax/varcode/tree/v9.3.5) (2026-09-19)

**Fixed**
- Genomic-layout translation only recodes an annotated selenocysteine when
  all three codon origins remain contiguous. Insertions or deletions that
  create a new TGA from part of that codon now terminate translation (#473).
- `protein_diff` inspects the mutant start codon at the same mapped CDS
  offset used for translation. An insertion immediately before the retained
  start codon is classified as 5′ UTR, without a false alternate-start or
  coding-sequence-change call (#473).

**Testing**
- Optional offline corpus checks use `osteosarc==0.1.0` and a pinned public
  snapshot, preserving native alleles and provenance and distinguishing the
  corrected MAP2 complex allele from the older deletion (#464).

## [v9.3.4](https://github.com/openvax/varcode/tree/v9.3.4) (2026-09-18)

**Fixed**
- Predicted proteins no longer stop at selenocysteine (#470). Mutant
  transcripts, splice outcomes, fusions and the experimental `transcript_model`
  layouts now read an annotated Sec UGA (`U` in the Ensembl reference protein)
  as Sec when no edit touches it and some selenoprotein 3′ UTR (SECIS) remains.
  Where none remains, as in a fusion downstream of Sec, UGA still terminates.
  Across the 60 complete selenoprotein transcripts in Ensembl 81, `protein_diff`
  had called most coding variants `PrematureStop` (176 of 240 probes), and
  `transcript_model` had called Sec→Trp a stop loss and missense changes after
  Sec silent (112 of 240). Both now agree with the default annotator on all 240.
  Opaque RNA imports without reference coordinates are still translated
  literally. When the 3′ UTR is only partly kept, the SV change flags ignore a
  supplied protein and stay unresolved.
- Joint translation of several edits (germline context, phasing, splice
  outcomes) now finds the CDS start after a 5′ UTR indel instead of reading
  the wrong frame (#471).
- `translate_sequence` accepts `selenocysteine=` offsets of TGA codons to
  read as `U`.

## [v9.3.3](https://github.com/openvax/varcode/tree/v9.3.3) (2026-09-18)

**Fixed**
- SV sequence-change flags no longer call unchanged selenoproteins
  protein-changing (#468). A TGA mapped onto an annotated selenocysteine
  (`U` in the Ensembl reference protein) is read as Sec where the model keeps
  that transcript through its 3′ end, which holds the SECIS element. It is read
  as a stop where no selenoprotein 3′ UTR remains (for example a fusion
  downstream of Sec, or a deletion of the whole 3′ UTR, which truncates the
  protein). When the SECIS may be only partly lost (a 3′ UTR deletion or
  duplication, a fusion in the 3′ UTR, or an unmapped import), the protein flag
  is reported only if both readings agree; otherwise it is `None`. The coding
  flag compares CDS bases and does not depend on Sec decoding. A supplied
  protein that ends exactly at a Sec residue no longer counts as a truncation,
  and partial observations continue past a decoded Sec codon. All 60 complete
  selenoprotein transcripts in Ensembl 81 now read unchanged for a
  reference-identical model.
- Complete-ORF comparisons read any start codon as the initiator methionine
  on both sides. Ensembl writes CTG/TTG initiators as `L`, so the 68 complete
  non-ATG transcripts in Ensembl 81 were previously called protein-changing
  when unchanged.

## [v9.3.2](https://github.com/openvax/varcode/tree/v9.3.2) (2026-09-18)

**Fixed**
- SV sequence-change flags no longer call a partial observation
  (`protein_completeness` of `partial_start`, `partial_end`, `partial_both`, or
  an explicit unknown) protein-changing just because it is shorter than the
  full reference protein (#462). Missing sequence is neither unchanged nor a
  truncation: these flags are `True` only when ORF bounds and reference-transcript
  segments place an observed codon in frame on a differing reference CDS codon,
  including premature stops and stop loss, and otherwise stay `None`.
  Only `start_to_stop`, or no label, permits a whole-protein comparison or
  the complete-ORF fallback; any other label fails closed. Exacto partial
  peptides, which have no reference coordinates, remain unresolved. Complete
  predictions and 9.3.0 candidate aggregation/filtering are unchanged.

## [v9.3.1](https://github.com/openvax/varcode/tree/v9.3.1) (2026-09-18)

**Documentation**
- Explain the shared Varcode / Isovar / Vaxrank responsibility split and evidence
  handoffs, with short README links and a focused integration guide. Distinguish
  available functionality from planned RNA/SV reconciliation (Isovar #305) and
  remaining downstream candidate-selection limits. No annotation behavior changed.

## [v9.3.0](https://github.com/openvax/varcode/tree/v9.3.0) (2026-09-18)

**Fixed**
- Structural effects now report coding/protein changes from available sequence
  and retained reference ORFs, rather than inheriting `False` (#418). Compare
  every candidate; retain a set if any alternative predicts a protein change.
  Coding deletions and the seven audited CPEB2–FAM193A fusion proteins survive
  filtering, while unchanged and synonymous proteins do not.

**Changed**
- SV sequence-change flags can now be `None` for unknown, distinct from `False`
  for unchanged. `Unresolved` and unclassified cryptic-exon candidates also use
  `None`. `drop_silent_and_noncoding()` retains unknowns by default; pass
  `keep_unresolved=False` to require a positive protein-change prediction.
  Candidate ordering, evidence, effect classes, and stored proteins are unchanged.

## [v9.2.8](https://github.com/openvax/varcode/tree/v9.2.8) (2026-09-18)

**Fixed**
- Retain all compatible fusion partner isoforms and junctions as candidates
  instead of selecting only the first (#406). Each carries its predicted
  protein when available; missing transcript sequence stays unresolved.
  Preserve the primary ordering and existing span/splice/RNA candidates, with
  no arbitrary count cap and no merging of distinct isoforms by protein alone.
- Add the audited CPEB2–FAM193A 1,541/1,500-aa alternatives, a 257-isoform case,
  and strand, provenance, assembly, and missing-sequence regressions. These
  are annotated-isoform predictions, not exhaustive splice/phase hypotheses
  or evidence of expression.

## [v9.2.7](https://github.com/openvax/varcode/tree/v9.2.7) (2026-09-17)

**Documentation**
- Split structural annotation, SV reference rules, and observed RNA imports;
  separate germline setup from phasing workflows. Replace the long API page
  with an index and five topic references, without duplicate object headings.
- Group effect definitions by family, simplify headings, and move transform
  contributor notes out of the user guide. Preserve existing section links.
- Correct the phased-VCF resolver example and cover it with executable tests.
  No annotation behavior or public APIs changed (#457, #458).

## [v9.2.6](https://github.com/openvax/varcode/tree/v9.2.6) (2026-09-17)

**Documentation**
- Shorten the effect annotation guide and split splice variants, transcript
  models, and experimental annotators into separate pages. Consolidate plugin
  instructions and provenance, use specific headings, and update navigation
  and example tests. Existing section links still lead to the relevant guides.
  No annotation behavior or public APIs changed.

## [v9.2.5](https://github.com/openvax/varcode/tree/v9.2.5) (2026-09-17)

**Documentation**
- Make the README a short introduction and add a reproducible getting-started
  guide. Group navigation into everyday tasks, advanced workflows, and reference;
  move implementation details behind usage examples (#415, #453).
- Correct germline examples and resolver descriptions (#413), stale splice and
  serialization explanations, and the missing Genome API link (#451).
- Lead the annotation guide with ordinary `effects()` usage. Explain optional
  implementations, supported-input limits, partial transcript models, sequence
  and candidate access, legacy aliases, and current haplotype ownership in an
  advanced section. No annotation behavior or public APIs changed.

## [v9.2.4](https://github.com/openvax/varcode/tree/v9.2.4) (2026-09-17)

**Fixed**
- BND fallbacks retain only the strand/orientation-correct reference cDNA
  prefix or suffix, not the entire transcript (#447). Single breakends retain
  their known local side; unknown local orientation leaves the model unresolved.
  Fragments are explicitly labeled partial, full cDNA/protein remain unknown,
  and supplied assemblies and existing coding-fusion predictions are preserved.

## [v9.2.3](https://github.com/openvax/varcode/tree/v9.2.3) (2026-09-17)

**Fixed**
- DUP/INV events with junction ends outside the selected transcript no longer
  fabricate a local duplicated/inverted transcript (#405). The DNA event class
  remains, but its unresolved mutant transcript is None. Existing fusion
  predictions and explicitly supplied allele assemblies are preserved.

## [v9.2.2](https://github.com/openvax/varcode/tree/v9.2.2) (2026-09-17)

**Fixed**
- Structural variants expose their actual REF/ALT, not the constructor's
  temporary nucleotide placeholders; small-edit flags are always false (#417).
- Mixed/SV tables include type, endpoint, mate and affected-span columns;
  point-only and empty table schemas remain unchanged. Structural CSV import
  raises explicitly instead of reconstructing misleading point variants.

## [v9.2.1](https://github.com/openvax/varcode/tree/v9.2.1) (2026-09-17)

**Fixed**
- Parsed symbolic DEL/DUP/INV/CNV spans exclude the retained VCF padding
  base from exon and mutant-cDNA annotation (#404). POS and junctions are
  unchanged; direct constructors retain their explicit affected-span defaults.
- Symbolic span records with missing/non-increasing END raise instead of
  being interpreted as a one-base event at the padding position.

## [v9.2.0](https://github.com/openvax/varcode/tree/v9.2.0) (2026-09-16)

**Added**
- `RNAReadPhasingSource.register_haplotype` tests an explicitly supplied local
  allele combination against anchored, quality-filtered RNA sequence. Equivalent
  deletion, splice-gap and split-gap alignments can support a known sequence
  without treating arbitrary RNA skips as DNA deletions (#441).
- `load_exacto_fusions(..., primary_structures_path=...)` imports native Exacto
  peptide predictions, including separate ORFs and partial proteins. Codons and
  coordinates are checked against observed RNA; original per-base provenance and
  completeness remain explicit. Prediction is not evidence of translation.

**Fixed**
- RNA read support and pair grouping use `(RG, QNAME)`, preventing unrelated
  libraries with reused names from creating false cis evidence or being
  collapsed into one supporting fragment (#443).

## [v9.1.0](https://github.com/openvax/varcode/tree/v9.1.0) (2026-09-16)

**Added**
- `load_exacto_fusions` imports selected SV-linked transcript structures and
  DNA/RNA integration rows without guessing an ORF or choosing one isoform.
  Original oriented sequence, splice/path rows, model IDs, and variant links
  remain available on existing `MutantTranscript`/effect candidates (#259, #261).
- `RNAEvidence` is a small concrete resolver for imported candidates;
  `make_fusion_outcome` imports a caller-specified observed sequence using
  existing structural effects. Partial and antisense/intergenic observations
  need not be promoted to coding fusions. Translation requires an explicit,
  complete start-to-stop ORF and is labeled prediction, not protein evidence.
- `MutantTranscript.from_sequence` shares external-sequence construction with
  the existing assembled-SV path, removing the duplicate sequence wrapper.

**Fixed**
- RNA observations attached to single-outcome structural predictions such as
  `Intronic` are retained alongside the DNA prediction. Deterministic
  point-variant behavior and default annotation remain unchanged.

## [v9.0.0](https://github.com/openvax/varcode/tree/v9.0.0) (2026-09-15)

**Breaking changes**
- Removed `StructuralVariantAnnotator`, the `varcode.annotators.structural_variant`
  module, and `annotator="structural_variant"`. Use `effects()` or
  `annotator="fast"` for both point variants and SVs; direct callers can use
  `FastEffectAnnotator().annotate_on_transcript(variant, transcript)`.
- Removed `UnsupportedVariantError` and its package/registry exports. Partial
  annotators return `NotImplemented`; public prediction APIs expose an
  `Unresolved` effect without silently selecting another annotator.
- Pickles referring to classes in the removed structural-annotator module or
  to `UnsupportedVariantError` must be read with varcode 8 before migrating.

**Changed**
- The default annotator calls internal structural prediction helpers directly,
  without instantiating another annotator. The experimental transcript model
  shares those helpers for fusion prediction. Classification, sequence assembly,
  candidates, and existing provenance are unchanged.
- `fast` remains the default; `protein_diff` and `transcript_model` remain
  experimental. Joint haplotype construction is unchanged (tracked in #437).

## [v8.0.3](https://github.com/openvax/varcode/tree/v8.0.3) (2026-09-15)

**Fixed**
- Contig validation now uses the actual Genome annotation dataset instead of
  sharing a cache by assembly name. Subsets, custom annotations and releases
  named GRCh38 no longer accept or reject chromosomes based on loading order
  (#402).
- Suppressed gene/transcript lookup errors return an empty `EffectCollection`
  with annotator provenance, rather than a plain list, and log the error.

## [v8.0.2](https://github.com/openvax/varcode/tree/v8.0.2) (2026-09-15)

**Fixed**
- `varcode` and `varcode-genes` no longer fail at startup when `pkg_resources`
  is unavailable. Both commands load their packaged logging configuration with
  the standard library's `importlib.resources`, including from zipped packages.
- Added subprocess regression tests for both commands with `pkg_resources`
  imports blocked, so older setuptools in a test environment cannot mask the
  missing runtime dependency.

## [v8.0.1](https://github.com/openvax/varcode/tree/v8.0.1) (2026-09-15)

**Changed**
- Renamed the experimental `realized` annotator to `transcript_model`, with
  `TranscriptModelEffectAnnotator` and `predict_transcript_model_effect` in
  `varcode.transcript_model`. The old registry name, class, function and module
  imports remain compatibility aliases. New annotation provenance uses
  `transcript_model`, including when selected through the old alias.
- The transcript model remains experimental and opt-in. The `fast` default,
  biological model, candidate ordering and unsupported-input behavior are
  unchanged.

## [v8.0.0](https://github.com/openvax/varcode/tree/v8.0.0) (2026-09-15)

**Changed**
- One built-in default annotator (`fast`) now owns point-edit and structural
  routing. Explicit `annotator="fast"` handles SVs just like `effects()`.
  Single-transcript APIs accept annotator and germline context and honor the
  same scoped selection as collection annotation.
- Partial experimental annotators return `NotImplemented` for unsupported
  inputs, without declaring `supports`. Public APIs expose those results as
  `Unresolved` with a reason and retain the selected annotator's provenance;
  they never silently substitute another annotator. `None` is an invalid
  plugin result. `UnsupportedVariantError` remains a compatibility import.
- Nonempty germline context requires `annotate_with_context`. The default
  retains established germline point-edit prediction; SV/germline composition
  remains experimental in `realized`. `protein_diff` declines SVs before
  reading placeholder nucleotide alleles. Realized-layout limitations return
  explicit unknowns without swallowing other errors.
- `protein_diff` and `realized` remain optional experiments; shared protein
  comparison helpers and point-edit parity tests remain in use. This release
  changes the annotator contract, not the default biological model or ranking.

**Fixed**
- The registry-default test now restores the prior annotator rather than
  leaking `protein_diff` into later tests.

## [v7.3.0](https://github.com/openvax/varcode/tree/v7.3.0) (2026-09-15)

**Added**
- Added the opt-in `realized` effect annotator and
  `predict_realized_effect`. It composes unknown germline phase, splice-site
  choices, point variants and local structural variants on one genomic
  layout, realizes the patient baseline and mutant product, classifies their
  protein difference, merges identical comparisons, and returns an ordinary
  top effect with ordered alternatives in `effect.candidates`.
- Added graph-native exon projection and splicing on forward and reverse
  strands. Canonical and exon-skip products resolve without a FASTA; with
  genomic sequence, intron-retention and cryptic-site products are translated
  from the mutated haplotype. Sequence-dependent tier-0 branches are reported
  as `Unresolved`.
- Added BND support to the realized annotator through the independently
  validated fusion assembler, including translated cross-gene products.

**Fixed**
- Splice choices whose required exon was removed by another choice are no
  longer enumerated, including the impossible combination of retaining an
  intron while skipping the exon that terminates it.
- The 64-hypothesis safety limit now applies to the combined phase × splice
  product instead of separately to each phase branch.
- Exonic alleles excluded from a realized mRNA are reported as `Silent` with
  `excluded_from_mrna=True`, rather than being mislabeled as genomic
  `Intronic` variants.
- Structural effects now expose the protein on their realized
  `mutant_transcript` through `effect.mutant_protein_sequence`.

## [v7.2.0](https://github.com/openvax/varcode/tree/v7.2.0) (2026-09-15)

**Added**
- Added the internal coordinate-aware genomic-layout foundation for the
  realized effect engine. Layouts preserve base origins through point edits,
  deletions, duplications, inversions, reverse strands, and cross-contig
  breakend joins without expanding large introns into per-base objects.
- Added composable phase/splice hypothesis and realized-product types.
  Splice hypotheses are validated as graph rewrites, rule preferences remain
  ordinal unless a calibrated probability is supplied, and candidate merging
  now keys on both the patient baseline and mutant product.

## [v7.1.2](https://github.com/openvax/varcode/tree/v7.1.2) (2026-09-15)

**Fixed**
- Cryptic splice-site scans now score the mutated allele and use one
  interbase boundary convention across donors, acceptors, and both strands.
  A destroyed canonical donor can no longer nominate itself as a cryptic
  replacement.
- FASTA-backed cryptic-splice and intron-retention outcomes now realize the
  originating exonic or intronic allele on the same mutant transcript instead
  of constructing the splice mechanism from reference sequence alone.

## [v7.1.1](https://github.com/openvax/varcode/tree/v7.1.1) (2026-09-14)

**Fixed**
- Four regression tests no longer silently skip. Three splice-outcome tests
  now inspect the always-on `SpliceOutcomeSet` wrapper and assert its
  `disrupted_signal_class`; the explicit 50-base deletion test now derives
  its coordinates from the pinned CFTR transcript instead of using a
  position one base before the intended exon.

## [v7.1.0](https://github.com/openvax/varcode/tree/v7.1.0) (2026-09-11)

**Changed**
- Structural-variant fusion annotation follows breakend orientation. The
  side of each breakpoint that's kept, with transcript strand, decides
  which transcript is the 5' and which the 3' partner, and a partner has
  to be in another gene and join sense-to-sense. `GeneFusion` is reported
  on either partner
  (`GeneFusion.transcript` can be the 3' partner) and carries
  `five_prime_transcript` / `three_prime_transcript`. A breakend whose
  ALT keeps the wrong sides for a fusion now gives
  `TranslocationToIntergenic`.
- A `DEL` / `DUP` / `INV` with one end in a transcript and a
  sense-to-sense partner at the other end gives a `GeneFusion` (e.g. the
  TMPRSS2-ERG deletion). The deletion, duplication or inversion of that
  transcript's exons follows the fusion in its `candidates`.
- `pair_breakends` builds a typed `DEL` / `DUP` / `INV` from a breakend
  pair that carries that `SVTYPE` (esvee, GRIDSS) when both halves agree
  on the label and their kept sides fit it, so `effects()` covers the
  whole span. Each record still loads as the breakend it is, so
  annotating an unpaired collection doesn't report the event twice.
- `StructuralVariant.junctions` exposes the novel adjacencies a variant
  creates as pairs of `Breakend` ends (position plus the side kept), and
  `breakpoints` lists their positions. Fusion assembly, transcript
  containment and cryptic-exon scanning all read from them.
- The reverse-complement orientation warning is replaced by a warning
  when a breakend with a mate has no breakend ALT to read orientation
  from.

**Fixed**
- VCF single breakends (`.ACGT` / `ACGT.`) no longer crash `load_vcf`;
  they load as `BND`s with no mate.
- SVs from `load_vcf` use the genome and contig-name settings passed to
  it. Previously they fell back to the default GRCh38 and kept `chr`
  names, so `effects()` raised on UCSC-named VCFs.
- `StructuralVariant.mate_contig` is normalized like `contig`, and
  `Variant._convert_ucsc_contig_name_to_ensembl` converts its argument
  rather than the variant's own contig.
- `StructuralVariant` equality compares every field the record carries
  (type, span, ALT, mate, assembled allele and confidence intervals), so
  `load_vcf(distinct=True)` no longer merges different SV records at the
  same position.
- `StructuralVariant.to_dict` covers every SV field, so JSON and pickle
  round-trip SVs (pickling one previously raised `TypeError`).
- A fusion keeps the base at each breakpoint, so a breakpoint inside an
  exon no longer drops one base from the 5' partner (and frameshifts the
  predicted protein).
- A fusion partner must be in a gene that doesn't span both ends of the
  junction, so an event inside one gene isn't reported as a fusion with a
  nested or antisense gene at its far end.
- A `GeneFusion` reported on its 3' partner carries a mutant transcript
  whose `reference_transcript` is the transcript being annotated.
- `reference_range` reads a range with one locus query instead of one
  per position when no single transcript spans it, which is the path
  cryptic-exon scoring takes on a genome with no chromosome FASTA.

**Added**
- A structural variant annotation guide (`docs/structural_variants.md`)
  covering which transcripts are annotated, how strand and kept sides
  decide fusion direction, every breakend combination, deletions,
  duplications and inversions by strand, and where a breakpoint lands.
- Fusion regression tests from the public osteosarc.com osteosarcoma
  dataset, validated against LINX (`tests/test_osteosarc_fusions.py`,
  `tests/data/osteosarc_esvee_somatic.vcf`).

## [v7.0.0](https://github.com/openvax/varcode/tree/v7.0.0) (2026-07-08)

**Changed**
- The default effect annotator is now **`fast`** (the offset-based
  classifier varcode has shipped since 2.0.0) instead of `protein_diff`.
  `Variant.effects()` / `VariantCollection.effects()` with no explicit
  `annotator=` now route through `fast`. This is a behavior change for
  callers that relied on the default: for SNVs / indels / MNVs the two
  annotators are fully reconciled (see `tests/test_protein_diff_parity.py`,
  `tests/test_annotator_parity_adversarial.py`, and
  `tests/test_annotator_divergence_scenarios.py`), so nearly all output is
  unchanged, but any residual divergence now resolves to `fast`'s
  classification. `protein_diff` stays available via
  `annotator="protein_diff"` or `varcode.use_annotator("protein_diff")`, and
  remains the substrate the `MutantTranscript` / splice-outcome / germline
  machinery builds on. Rationale: `fast` is the more battle-tested path —
  during the `protein_diff` bring-up it was effectively the correctness oracle
  `protein_diff` was reconciled against
  ([#318](https://github.com/openvax/varcode/issues/318)–[#321](https://github.com/openvax/varcode/issues/321)),
  and it has the cleaner bug history
  ([#397](https://github.com/openvax/varcode/issues/397)). Major version bump
  because default behavior changes.

**Added**
- Regression pins in `tests/test_annotator_divergence_scenarios.py` for the
  frameshift coincidental-shared-suffix class
  ([#396](https://github.com/openvax/varcode/pull/396) /
  [#397](https://github.com/openvax/varcode/issues/397)): a frameshift whose
  novel C-terminus coincidentally ends with the reference protein's terminal
  residue(s) must retain the full tail under **both** annotators (CFTR
  `p.L127fs` on the + strand, BRCA1 `p.R71fs` on the − strand).

## [v6.0.2](https://github.com/openvax/varcode/tree/v6.0.2) (2026-07-08)

**Fixed**
- The default protein-diff annotator (`classify_from_protein_diff`) no
  longer drops the last 1–2 residues of a frameshift/stop-loss novel
  C-terminus when that tail coincidentally ends with the same residue(s)
  as the reference protein's own C-terminus. Reading-frame-changing
  branches now trim only the shared *prefix* (`trim_shared_prefix`)
  instead of `trim_shared_flanking_strings`, which also stripped a shared
  *suffix* — correct for an in-frame indel, but wrong for a novel tail
  running to a new stop codon. Example: ATM p.F61fs (GRCh38
  `chr11:g.108227882delT`, `ENST00000675843`) now reports `74 aa … FRKKQNV`
  from `Variant.effects()`, matching `annotator="fast"` and the true ORF,
  instead of the truncated `73 aa … FRKKQN`. `annotator="fast"` and
  `predict_variant_effect_on_transcript()` were never affected
  ([#396](https://github.com/openvax/varcode/pull/396),
  [#397](https://github.com/openvax/varcode/issues/397)).

## [v6.0.1](https://github.com/openvax/varcode/tree/v6.0.1) (2026-06-18)

**Fixed**
- In-frame deletions that remove the stop codon of a transcript with no
  3' UTR sequence (e.g. MAPK3-006 / `ENST00000395199`, whose
  `three_prime_utr_sequence` is `""`) no longer raise
  `ValueError: If no amino acids added by StopLoss then it should be Silent`.
  With no readthrough sequence to translate into, the effect is now
  classified as a C-terminal `Deletion` instead of an invalid `StopLoss`
  with an empty `aa_alt`. Both annotators agree on this: the in-frame
  predictor (`FastEffectAnnotator`) no longer constructs the invalid
  `StopLoss`, and the default protein-diff classifier
  (`classify_from_protein_diff`) no longer mislabels the truncated protein
  as a `PrematureStop` — there is no stop codon in the mutant CDS, so a
  premature stop is incorrect. The earlier
  [#246](https://github.com/openvax/varcode/issues/246) fix only covered
  transcripts with a non-empty 3' UTR
  ([#394](https://github.com/openvax/varcode/issues/394)).

## [v6.0.0](https://github.com/openvax/varcode/tree/v6.0.0) (2026-05-26)

**Fixed**
- `apply_variants_to_transcript` now refuses an insertion abutting
  another edit at the same cDNA offset. The previous overlap check
  only caught range overlap and the two-insertions-at-the-same-offset
  case, so an insertion paired with a substitution or deletion
  starting at the same offset slipped through and produced an
  order-dependent joint cDNA. User-visible impact: phased haplotype /
  germline pipelines (`varcode/phasing.py`, `varcode/germline.py`)
  that previously returned a silently order-dependent joint effect for
  these inputs now return independent per-variant effects computed
  against the reference instead.
- IUPAC ambiguity codes (R, Y, S, W, K, M, D, V, H, B) in a
  reverse-strand variant's alleles are now reverse-complemented via
  the full IUPAC translation table. The private helper in
  `varcode/mutant_transcript.py` only covered A/C/G/T/N and silently
  passed ambiguity codes through unchanged; `_resolve_variant_edit`
  now delegates to `varcode.nucleotides.reverse_complement`.

**Added**
- Added `RNAReadPhasingSource`, a BAM-backed `ReadPhasingSource`
  implementation for RNA read/fragment co-occurrence. It is consumed
  through `MolecularPhaseResolver(source)` and lives behind the optional
  `varcode[rna]` / `pysam` dependency. `ReadPhaseResolver` remains as
  the varcode 5.0 compatibility name.
- `SpliceOutcomeSet.effect_if_splicing_unchanged` — the canonical
  "alternative outcome" accessor: the coding consequence that applies
  if splicing proceeds normally (the `NormalSplicing` candidate's
  `coding_effect`). Unlike the legacy `ExonicSpliceSite.alternate_effect`,
  it works for intronic splice disruptions — returning `None` when
  the nucleotide change leaves the protein untouched (i.e. there is
  no coding consequence to attach). Sits alongside `most_likely_effect`
  and `candidates` as the three-accessor surface on `SpliceOutcomeSet`
  ([#391](https://github.com/openvax/varcode/issues/391)).

**Breaking**
- `SpliceOutcomeSet` is now always-on for splice-disrupting variants
  ([#391](https://github.com/openvax/varcode/issues/391)). Every variant
  that lands in the canonical splice window — `SpliceDonor`,
  `SpliceAcceptor`, `ExonicSpliceSite`, `IntronicSpliceSite` — is
  wrapped in a `SpliceOutcomeSet` carrying the candidate mechanisms.
  Specifically:
  - The `splice_outcomes=True` flag on `Variant.effects()` /
    `VariantCollection.effects()` / `predict_variant_effects()` is
    **removed**. Callers passing it explicitly get a `TypeError`.
    Migration: drop the keyword — wrapping is unconditional.
  - `Variant.effect_on_transcript(transcript)` and the
    `FastEffectAnnotator` / `ProteinDiffEffectAnnotator` per-transcript
    paths return a `SpliceOutcomeSet` for splice-disrupting variants
    instead of the raw `ExonicSpliceSite` / `SpliceDonor` / etc. class.
    Migration: replace `isinstance(effect, ExonicSpliceSite)` with
    `isinstance(effect, SpliceOutcomeSet) and effect.disrupted_signal_class is ExonicSpliceSite`.
    `effect.alternate_effect` still works as a back-compat alias for
    `effect.effect_if_splicing_unchanged`, so attribute access keeps
    working through the wrapper.
  - `SpliceOutcomeSet.modifies_protein_sequence` is hardcoded to
    `True` (a splice disruption is always *potentially* protein-
    modifying via a non-`NormalSplicing` candidate). Closes a long-
    standing filter bug where `drop_silent_and_noncoding()` silently
    dropped exonic-splice-site variants whose `NormalSplicing.coding_effect`
    happened to be `Silent`.
  - Candidate construction is lazy: only the cheap `NormalSplicing`
    candidate is built eagerly; `ExonSkipping`, `IntronRetention`,
    `CrypticDonor` / `CrypticAcceptor` materialise on first
    `.candidates` access. Pipelines that filter on
    `modifies_protein_sequence` / `effect_priority` and never read
    `.candidates` pay only the eager cost.
  - `SpliceOutcomeSet` is now a `TranscriptMutationEffect` subclass
    (alongside `MultiOutcomeEffect`), so it carries `gene` /
    `transcript` and matches the standard `isinstance(effect, TranscriptMutationEffect)`
    filter used by downstream consumers.
- Unified the multi-outcome machinery: `SpliceCandidate` deleted;
  `MultiOutcomeEffect.outcomes` accessor + `_with_extra_outcomes`
  helper + `_extra_outcomes` slot removed
  ([#382](https://github.com/openvax/varcode/issues/382)).
  - `SpliceOutcomeSet.candidates` now returns
    `tuple[EffectCandidate, ...]` — the same shape every other
    `MultiOutcomeEffect` subclass exposes. Each entry wraps an
    inner `SpliceMechanismEffect` (`NormalSplicing`,
    `ExonSkipping`, `IntronRetention`, `CrypticDonor`, or
    `CrypticAcceptor`). Mechanism identity now lives on
    `type(candidate.effect)`, and the mechanism object carries
    fields like `affected_exon`, `side`, `cryptic_genomic_position`,
    `aa_ref`, `aa_alt`, and `mutant_transcript`. The previous
    `candidate.outcome` / `candidate.plausibility` /
    `candidate.coding_effect` / `candidate.predicted_class_name` /
    `candidate.mutant_transcript` / `candidate.has_protein`
    fields are gone — read provenance off the `EffectCandidate`
    (`.source`, `.evidence`) and mechanism/protein state off the
    inner effect (`candidate.effect`,
    `candidate.effect.mutant_transcript`).
    `candidate.plausibility` has no one-for-one semantic replacement:
    it was the old splice-specific name for a DNA-only ordering
    heuristic, not evidence. There is no shared
    `EffectCandidate.probability`; producer-specific support belongs
    in `candidate.evidence` under explicit names.
  - `MultiOutcomeEffect.candidates` is the single accessor on every
    subclass (`SpliceOutcomeSet`, `StructuralVariantEffect`,
    `PhaseCandidateSet`, `ExonicSpliceSite`, `HaplotypeEffect`).
    A new `MultiOutcomeEffect.effects` convenience property unwraps
    to `tuple(c.effect for c in self.candidates)` for callers that
    don't need provenance.
  - The post-hoc attachment slot renamed `_extra_outcomes` →
    `_extra_candidates`; the merge helper renamed
    `_with_extra_outcomes` → `_combine_with_extra_candidates`.
    `apply_rna_evidence_to_effects` still uses `_extra_candidates`
    for non-splice multi-outcome effects; splice mechanism sets now
    reconcile RNA evidence into a replacement set with
    `rna_evidence`, `added_candidates`, `excluded_candidates`, and
    `candidate_rna_evidence` audit fields. External integrations
    that touched these private names must rename.
  - `StructuralVariantEffect.__init__` parameter renamed
    `candidates=` → `primary_effects=` (carries the inner
    `MutationEffect` tuple; the `candidates` accessor now lifts
    to `EffectCandidate` automatically). Same on `LargeDeletion`,
    `LargeDuplication`, `GeneFusion`.
  - `MultiOutcomeEffect.most_likely` is **removed**. Replaced by
    four explicit accessors so callers never confuse "wrapped vs
    unwrapped" or "likeliest vs most-disruptive":
    - `.most_likely_candidate` → `EffectCandidate` (same as
      `candidates[0]`)
    - `.most_likely_effect` → inner `MutationEffect` of the above
    - `.highest_priority_candidate` → `EffectCandidate` with the
      highest `effect_priority` among candidates (worst-case
      classification, independent of producer order)
    - `.highest_priority_effect` → inner `MutationEffect` of the above
    Callers doing `effect.most_likely.mutant_protein_sequence` or
    `effect.most_likely.aa_ref` should switch to
    `effect.most_likely_effect.mutant_protein_sequence` etc.;
    callers that want the wrapper (with `.source`, `.evidence`) use
    `effect.most_likely_candidate`.
  - Splice mechanisms promoted to first-class `MutationEffect`
    classes, each carrying its own protein vocab (`aa_ref` /
    `aa_alt` / `mutant_protein_sequence` / `mutant_transcript`) on
    the instance — no more wrapping a separate coding effect.
    New hierarchy:
    - `SpliceMechanismEffect(TranscriptMutationEffect)` — base,
      carries `splice_signal` referencing the underlying
      `SpliceDonor` / `SpliceAcceptor` / `IntronicSpliceSite` /
      `ExonicSpliceSite` so each mechanism knows *where* the
      disruption was.
    - `NormalSplicing` — splice signal hit but splicing proceeds;
      carries `coding_effect` for the underlying nucleotide-level
      change (or `None` for purely intronic disruption).
    - `ExonSkipping` — affected exon excluded; carries
      `affected_exon` and `in_frame`.
    - `IntronRetention` — intron retained; carries
      `retained_intron_start`, `retained_intron_end`, `side`.
    - `CrypticDonor` / `CrypticAcceptor` — cryptic site replaces
      canonical; carry `affected_exon`,
      `cryptic_genomic_position`, `motif_score`,
      `exon_length_delta`.
    Unresolved state is "protein fields are `None`" — no parallel
    placeholder class hierarchy. Class identity = mechanism;
    consumers dispatch on `isinstance(candidate.effect,
    ExonSkipping)` instead of evidence-key checks. Resolved
    mechanisms retain their classified protein consequence as
    `protein_effect` so severity queries (`effect_priority`,
    `modifies_protein_sequence`) behave like ordinary coding effects
    while preserving mechanism identity.
  - **Deleted**: `SpliceOutcome` enum (replaced by class identity),
    `PredictedIntronRetention` (subsumed by `IntronRetention`),
    `PredictedCrypticSpliceSite` (split into `CrypticDonor` +
    `CrypticAcceptor`), `SpliceOutcomeSet.to_dict` /
    `.from_dict` overrides (no enum left to stringify),
    `evidence["splice_outcome"]` / `evidence["placeholder"]` /
    `evidence["description"]` keys (info now lives on the
    mechanism Effect — `type(candidate.effect)`,
    `candidate.effect.resolved`, `candidate.effect.short_description`),
    `_placeholder_effect_for_outcome` and `_make_splice_candidate`
    helpers.
  - **`SpliceOutcomeSet.candidate_proteins`** now keys by mechanism
    class (e.g. `proteins[ExonSkipping]`) instead of `SpliceOutcome`
    enum value.
- `varcode.Outcome` renamed to `varcode.EffectCandidate`. The helper
  `outcomes_from_candidates` renamed to `candidates_from_effects`.
  The module `varcode.outcomes` renamed to
  `varcode.effect_candidates`. The test file `tests/test_outcomes.py`
  renamed to `tests/test_effect_candidates.py`. The `description`
  field on the wrapper class is removed — use
  `candidate.effect.short_description` (it was always a passthrough).
  Also removed from `make_rna_outcome(description=...)`. No
  deprecation alias; update imports. The class ships with a module
  docstring explaining *why* the wrapper exists: the same
  `MutationEffect` instance can appear in multiple multi-outcome
  contexts with different per-context provenance (e.g. a splice
  candidate re-surfaced by the SV annotator with a different
  `source` tag and `sv_type` in evidence); putting metadata on the
  wrapper instead of the Effect lets the Effect stay shared while
  the labels diverge. (Per #382, `.outcomes` has since been removed
  and `MultiOutcomeEffect.candidates` is the single accessor across
  every subclass.) Aspirational
  `"isovar"` / `"exacto"` / `"longread_assembly"` example tags
  scrubbed from varcode docstrings.
- Phasing API generalized; Isovar-named identifiers removed from the
  varcode public surface ([#378](https://github.com/openvax/varcode/issues/378)).
  Varcode no longer imports or names any upstream tool — implementations
  of the new generic Protocols live in their respective packages
  (e.g. `isovar.IsovarReadPhasing`, openvax/isovar#183).
  - `IsovarAssemblyProvider` (Protocol) **removed**, split into:
    - `ReadPhasingSource` with `has_evidence(variant) -> bool` and
      `partners_in_cis(variant) -> Sequence[Variant]`.
    - `MutantTranscriptSource` with
      `mutant_transcript(variant, transcript) -> Optional[MutantTranscript]`.
  - `IsovarPhaseResolver` renamed to `ReadPhaseResolver`. Constructor
    accepts any `ReadPhasingSource`; routes `mutant_transcript(...)`
    to the wrapped source when it also satisfies `MutantTranscriptSource`.
    Returns `None` otherwise instead of raising.
  - Resolver `source` tag changed from `"isovar"` to `"read_phasing"`
    on `ReadPhaseResolver`. Consumers filtering effects by phase
    source need to update their filter values.
- `varcode.effects.effect_classes.PhaseAmbiguousEffect` renamed to
  `PhaseCandidateSet`. No deprecation alias — update imports
  ([#376](https://github.com/openvax/varcode/pull/376)).
  Per #382, the public surface is `.candidates` (a
  `tuple[EffectCandidate, ...]` with per-hypothesis evidence keys),
  `.most_likely_candidate` / `.most_likely_effect` /
  `.highest_priority_candidate` / `.highest_priority_effect`, and
  `.short_description`.

**Changed**
- `SpliceOutcomeSet` serialization migrated onto `serializable>=1.1.0`'s
  standard introspection (the parallel `SpliceCandidate` dataclass that
  also lived on this path has since been deleted per #382). The
  `__effect_class__` tagging and hand-rolled class registries
  (`_CODING_EFFECT_CLASS_REGISTRY`, `_SPLICE_SIGNAL_CLASS_REGISTRY`,
  `_rehydrate_coding_effect`) are gone; JSON round-trip now flows
  through `serializable.helpers`' standard `__class__` / `__module__`
  stamping ([#343](https://github.com/openvax/varcode/issues/343)).
  `SpliceOutcomeSet.to_dict` / `from_dict` are overridden to stringify
  the `SpliceOutcome` enum stored under
  `candidate.evidence["splice_outcome"]` (and rehydrate it on the way
  back) without mutating `self` mid-call, and to emit a single
  `candidates` key on the wire (no parallel `_candidates`).
  The JSON wire format is unchanged, but **the pre-#305 migration
  shim for the internal `_ExonSkipFrameshiftEffect` class has been
  removed**. JSON produced by varcode 2.4.x or earlier (which could
  contain `"__effect_class__": "_ExonSkipFrameshiftEffect"` tags) will
  no longer rehydrate; re-emit from the current annotator or hand-patch
  those tags to `FrameShift` before loading. Anyone still reading such
  payloads can pin `varcode<4.7` or resurrect the shim in user code.

**Changed**
- README "Effect Types" section rewritten: 14 missing concrete
  effect classes added (`Failure`, the splice mechanisms, the
  structural-variant effects, `CrypticExonCandidate`,
  `HaplotypeEffect`, `PhaseCandidateSet`), grouped into 7 sub-tables
  by biological context, with a new intro explaining
  `MultiOutcomeEffect` and how multi-possibility effects are
  represented. Each class name links to its source definition via
  GitHub text-fragment URLs that survive line-number drift.
- Documented and streamlined the splice-effect model: clarified that
  `SpliceSite` effects (`SpliceDonor` / `SpliceAcceptor` /
  `IntronicSpliceSite` / `ExonicSpliceSite`) describe *where* a splice
  signal was hit and carry no protein consequence on their own, while
  `SpliceMechanismEffect` subclasses describe *what* the spliceosome
  does and carry the protein change. `SpliceSite` (previously an
  undocumented marker base) is now the documented, load-bearing type:
  `enumerate_splice_outcomes` gates on `isinstance(effect, SpliceSite)`
  rather than enumerating the four subclasses. Documented the
  `SpliceOutcomeSet.disrupted_signal_class` (a `SpliceSite` *subclass*,
  i.e. a type, for priority lookup) vs each candidate's
  `effect.splice_signal` (a `SpliceSite` *instance*) distinction.
  No behavior change.

**Added**
- Docs site now has an "Effect types" page that auto-renders every
  class in `varcode.effects.effect_classes` via mkdocstrings, so the
  documented catalog stays in sync with the code (previously only the
  abstract bases were in the API reference, while the full list lived
  only in the hand-maintained README table).
- Exported the splice-signal disruption effects at the package root
  for parity with the already-public splice mechanism effects:
  `SpliceSite` (the shared base), `SpliceDonor`, `SpliceAcceptor`,
  `IntronicSpliceSite`, and `ExonicSpliceSite`. `from varcode import
  SpliceSite` now works, so `isinstance(effect, SpliceSite)` can be
  used to catch any splice-site disruption without reaching into
  `varcode.effects.effect_classes`.

## [v2.3.0](https://github.com/openvax/varcode/tree/v2.3.0) (2026-04-13)

**Added**
- Per-sample genotype / zygosity access ([#267](https://github.com/openvax/varcode/issues/267)).
  `Genotype` frozen dataclass, `Zygosity` enum (`ABSENT`/`HETEROZYGOUS`/`HOMOZYGOUS`/`MISSING`),
  new `VariantCollection` methods `.samples`, `.genotype(variant, sample)`,
  `.zygosity(variant, sample)`, `.for_sample(name)`, `.heterozygous_in(name)`,
  `.homozygous_alt_in(name)`. Multi-allelic aware: each split Variant
  reports zygosity relative to its own alt.
- `varcode.SampleNotFoundError(KeyError)` raised on typoed sample names.

## [v2.2.1](https://github.com/openvax/varcode/tree/v2.2.1) (2026-04-13)

**Fixed**
- Ref-vs-genome mismatches now raise a dedicated
  `varcode.ReferenceMismatchError` (subclass of `ValueError`) with an
  actionable message naming the likely causes and pointing at
  `raise_on_error=False` ([#215](https://github.com/openvax/varcode/issues/215),
  [#246](https://github.com/openvax/varcode/issues/246)).

## [v2.2.0](https://github.com/openvax/varcode/tree/v2.2.0) (2026-04-12)

**Added**
- `from_csv` now accepts either `chr` or `contig` as the contig column
  name on both `VariantCollection` and `EffectCollection`, so CSVs are
  interchangeable between the two types
  ([#274](https://github.com/openvax/varcode/issues/274)).
- `from_csv` warns on major `varcode_version` drift recorded in the
  CSV header ([#275](https://github.com/openvax/varcode/issues/275)).

**Changed**
- `from_csv` docstrings now point users at `from_json` for byte-for-byte
  round-trip or larger collections
  ([#276](https://github.com/openvax/varcode/issues/276)).

## [v2.1.0](https://github.com/openvax/varcode/tree/v2.1.0) (2026-04-12)

**Added**
- `VariantCollection.from_csv` and `EffectCollection.from_csv` for
  round-trip deserialization ([#273](https://github.com/openvax/varcode/pull/273)).
- `to_csv` prepends a `# key=value` metadata header by default
  (`varcode_version`, `reference_name`); `from_csv` reads it so the
  `genome` argument becomes optional. Pass `include_header=False` for
  legacy consumers.

**Fixed**
- `VariantCollection.variants` and `EffectCollection.effects` now
  match the collection's iteration order instead of holding the raw
  pre-sort / pre-dedup input list
  ([#220](https://github.com/openvax/varcode/issues/220)).

## [v2.0.0](https://github.com/openvax/varcode/tree/v2.0.0) (2026-04-11)

Major release — several backward-incompatible fixes. See
[#263](https://github.com/openvax/varcode/pull/263) and
[#265](https://github.com/openvax/varcode/pull/265) for full details.

**Breaking**
- `Silent.short_description` returns HGVS `p.{ref}{pos}=` (e.g.
  `p.R6=`) instead of the literal `"silent"`
  ([#217](https://github.com/openvax/varcode/issues/217)).
- `Silent.aa_pos` no longer includes the shared-prefix offset; it now
  points at the actual synonymous codon
  ([#208](https://github.com/openvax/varcode/issues/208)).
- `PrematureStop.short_description` returns `p.{pos}ins{alt}*` when
  `aa_ref` is empty instead of the ambiguous `p.{pos}{alt}*`
  ([#216](https://github.com/openvax/varcode/issues/216)).
- `EffectCollection` is sorted by effect priority (most severe first)
  by default. Pass `sort_key=False` to disable or a custom callable
  to override ([#227](https://github.com/openvax/varcode/issues/227)).
- Intronic splice classification is sequence-aware: variants at
  `+1`/`+2` or `-1`/`-2` with a non-canonical reference base are
  classified as `IntronicSpliceSite` rather than
  `SpliceDonor`/`SpliceAcceptor`
  ([#262](https://github.com/openvax/varcode/issues/262)).

**Fixed**
- SNV in the stop codon with a stop-prefixed 3' UTR is correctly
  classified as `StopLoss` instead of `Insertion`
  ([#250](https://github.com/openvax/varcode/issues/250),
  [#205](https://github.com/openvax/varcode/issues/205)).
- Insertion before the stop codon that produces an identical protein
  is correctly classified as `Silent`
  ([#201](https://github.com/openvax/varcode/issues/201)).
- `changes_exonic_splice_site` now applies the mutation before
  checking the splice pattern ([#262](https://github.com/openvax/varcode/issues/262)).
- VCF loader skips symbolic alleles (`<DEL>`, `<CN0>`, `<INS:ME:ALU>`,
  ...) and breakend notation with a visible warning instead of
  crashing ([#88](https://github.com/openvax/varcode/issues/88)).
  Full SV support is tracked in [#264](https://github.com/openvax/varcode/issues/264).

## [v0.5.15](https://github.com/hammerlab/varcode/tree/v0.5.15) (2017-04-28)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.14...v0.5.15)

**Closed issues:**

- Allow contig to be empty to allow parsing of MAF with faulty mutation\(s\) [\#210](https://github.com/hammerlab/varcode/issues/210)

**Merged pull requests:**

- Fixes to load\_maf [\#223](https://github.com/hammerlab/varcode/pull/223) ([tavinathanson](https://github.com/tavinathanson))
- added raise\_on\_error option to load\_maf and load\_maf\_dataframe [\#221](https://github.com/hammerlab/varcode/pull/221) ([iskandr](https://github.com/iskandr))
- Optionally allow duplicated mutations when using load\_vcf or load\_maf. Fixes \#211  [\#212](https://github.com/hammerlab/varcode/pull/212) ([tuomastik](https://github.com/tuomastik))

## [v0.5.14](https://github.com/hammerlab/varcode/tree/v0.5.14) (2017-04-05)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.13...v0.5.14)

**Merged pull requests:**

- Adding 'distinct' as a parameter to load\_vcf. [\#222](https://github.com/hammerlab/varcode/pull/222) ([julia326](https://github.com/julia326))

## [v0.5.13](https://github.com/hammerlab/varcode/tree/v0.5.13) (2017-04-01)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.12...v0.5.13)

**Closed issues:**

- Effect prediction throws error \(even if raise\_on\_error=False\) [\#213](https://github.com/hammerlab/varcode/issues/213)
- Optionally allow duplicated mutations when using load\_vcf or load\_maf [\#211](https://github.com/hammerlab/varcode/issues/211)

**Merged pull requests:**

- install ensembl 87 on travis [\#219](https://github.com/hammerlab/varcode/pull/219) ([iskandr](https://github.com/iskandr))
- Allow user to affect the sorting of variants when loading a VCF or MAF. [\#218](https://github.com/hammerlab/varcode/pull/218) ([tuomastik](https://github.com/tuomastik))

## [v0.5.12](https://github.com/hammerlab/varcode/tree/v0.5.12) (2017-01-18)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.11...v0.5.12)

**Closed issues:**

- Make Varcode correctly infer genome for b37-decoy string [\#207](https://github.com/hammerlab/varcode/issues/207)
- Longer indels in random variants [\#47](https://github.com/hammerlab/varcode/issues/47)
- Predict coding sequence of StartLoss mutations [\#4](https://github.com/hammerlab/varcode/issues/4)

**Merged pull requests:**

- Add optional\_cols list to load\_maf [\#209](https://github.com/hammerlab/varcode/pull/209) ([tavinathanson](https://github.com/tavinathanson))

## [v0.5.11](https://github.com/hammerlab/varcode/tree/v0.5.11) (2016-12-05)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.10...v0.5.11)

**Fixed bugs:**

- Varcode noncoding variant in a drop\_silent\_and\_noncoding\(\) list [\#200](https://github.com/hammerlab/varcode/issues/200)

**Merged pull requests:**

- Adding aa\_ref argument to StopLoss for variants which delete codons before stop [\#203](https://github.com/hammerlab/varcode/pull/203) ([iskandr](https://github.com/iskandr))

## [v0.5.10](https://github.com/hammerlab/varcode/tree/v0.5.10) (2016-10-19)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.9...v0.5.10)

**Fixed bugs:**

- Variant pickling won't work for not-human and non-EnsemblRelease Genomes [\#147](https://github.com/hammerlab/varcode/issues/147)

**Closed issues:**

- Link on PyPI badge broken [\#191](https://github.com/hammerlab/varcode/issues/191)
- Reference incorrectly inferred when "b36" in reference file path [\#181](https://github.com/hammerlab/varcode/issues/181)
- Premature stop codon error [\#166](https://github.com/hammerlab/varcode/issues/166)

**Merged pull requests:**

- explicit args to \_\_init\_\_ of Intronic splice effects fixes serialization [\#199](https://github.com/hammerlab/varcode/pull/199) ([iskandr](https://github.com/iskandr))
- Update RELEASING.md, fixing tagging instructions [\#198](https://github.com/hammerlab/varcode/pull/198) ([julia326](https://github.com/julia326))

## [v0.5.9](https://github.com/hammerlab/varcode/tree/v0.5.9) (2016-10-11)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.8...v0.5.9)

**Fixed bugs:**

- StopLoss pickling is broken [\#188](https://github.com/hammerlab/varcode/issues/188)

**Closed issues:**

- One logger per module [\#196](https://github.com/hammerlab/varcode/issues/196)
- SNV results in deletion [\#193](https://github.com/hammerlab/varcode/issues/193)

**Merged pull requests:**

- One logger per module. [\#197](https://github.com/hammerlab/varcode/pull/197) ([julia326](https://github.com/julia326))
- Fix edge case where PrematureStop in last amino acid got interpreted as a Deletion  [\#194](https://github.com/hammerlab/varcode/pull/194) ([iskandr](https://github.com/iskandr))
- Fix inferred-reference-bug [\#182](https://github.com/hammerlab/varcode/pull/182) ([jburos](https://github.com/jburos))

## [v0.5.8](https://github.com/hammerlab/varcode/tree/v0.5.8) (2016-09-28)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.7...v0.5.8)

**Merged pull requests:**

- changed Markdown image links to HTML [\#192](https://github.com/hammerlab/varcode/pull/192) ([iskandr](https://github.com/iskandr))

## [v0.5.7](https://github.com/hammerlab/varcode/tree/v0.5.7) (2016-09-28)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.3...v0.5.7)

## [v0.5.3](https://github.com/hammerlab/varcode/tree/v0.5.3) (2016-09-28)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.2...v0.5.3)

**Merged pull requests:**

- Use conda to install pypandoc [\#189](https://github.com/hammerlab/varcode/pull/189) ([arahuja](https://github.com/arahuja))
- Ensure README.md is packaged [\#186](https://github.com/hammerlab/varcode/pull/186) ([arahuja](https://github.com/arahuja))
- Upgrade serializable dependency with tests [\#185](https://github.com/hammerlab/varcode/pull/185) ([arahuja](https://github.com/arahuja))
- Add pypi badge [\#184](https://github.com/hammerlab/varcode/pull/184) ([arahuja](https://github.com/arahuja))

## [v0.5.2](https://github.com/hammerlab/varcode/tree/v0.5.2) (2016-09-28)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.1...v0.5.2)

**Closed issues:**

- Make sure `Variant` works with any Genome \(not just a human EnsemblRelease\) [\#127](https://github.com/hammerlab/varcode/issues/127)

**Merged pull requests:**

- Move extraneous variables to properties for normalization [\#190](https://github.com/hammerlab/varcode/pull/190) ([tavinathanson](https://github.com/tavinathanson))
- Use is\_protein\_coding property of pyensembl.Transcript and pyensembl.Gene [\#180](https://github.com/hammerlab/varcode/pull/180) ([iskandr](https://github.com/iskandr))

## [v0.5.1](https://github.com/hammerlab/varcode/tree/v0.5.1) (2016-09-16)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.5.0...v0.5.1)

**Merged pull requests:**

- Add MutationEffect to \_\_init\_\_.py [\#178](https://github.com/hammerlab/varcode/pull/178) ([timodonnell](https://github.com/timodonnell))

## [v0.5.0](https://github.com/hammerlab/varcode/tree/v0.5.0) (2016-09-13)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.19...v0.5.0)

**Implemented enhancements:**

- Support collection.as\_dataframe\(\) [\#128](https://github.com/hammerlab/varcode/issues/128)

**Closed issues:**

- Substitution mis-annotated as stop-loss [\#176](https://github.com/hammerlab/varcode/issues/176)
- Wrong aa\_mutation\_end\_offset for insertion of stop codon [\#175](https://github.com/hammerlab/varcode/issues/175)
- Wrong aa\_ref for insertion of stop codon [\#174](https://github.com/hammerlab/varcode/issues/174)
- Insertions after the stop codon annotated as plain Insertions [\#172](https://github.com/hammerlab/varcode/issues/172)
- Mutations before the stop codon confused as StopLosses [\#171](https://github.com/hammerlab/varcode/issues/171)
- StopLosses do not translate into 3' UTR [\#170](https://github.com/hammerlab/varcode/issues/170)
- Insertion of stop codon is annotated as simple Insertion and not PrematureStop [\#169](https://github.com/hammerlab/varcode/issues/169)
- Synonimous FrameShift over stop codon not annotated as silent  [\#168](https://github.com/hammerlab/varcode/issues/168)
- Wrong offset for insertion of StopCodon [\#167](https://github.com/hammerlab/varcode/issues/167)
- Document release process [\#154](https://github.com/hammerlab/varcode/issues/154)
- compare variants that use different references [\#83](https://github.com/hammerlab/varcode/issues/83)
- Annotate with predicted pathogenicity [\#46](https://github.com/hammerlab/varcode/issues/46)

**Merged pull requests:**

- Reorganize effect prediction code, fixed annotation bugs/issues [\#173](https://github.com/hammerlab/varcode/pull/173) ([iskandr](https://github.com/iskandr))

## [v0.4.19](https://github.com/hammerlab/varcode/tree/v0.4.19) (2016-09-12)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.18...v0.4.19)

**Fixed bugs:**

- original\_start doesn't get pickled [\#141](https://github.com/hammerlab/varcode/issues/141)

**Closed issues:**

- replace `load\_vcf` with `load\_vcf\_fast` ? [\#144](https://github.com/hammerlab/varcode/issues/144)
- Add `annotate\_random\_variants` commandline script [\#49](https://github.com/hammerlab/varcode/issues/49)
- support filtering a variant collection to variants overlapping specified gene names [\#32](https://github.com/hammerlab/varcode/issues/32)
- Use SPANR to identify splicing misregulation [\#2](https://github.com/hammerlab/varcode/issues/2)

## [v0.4.18](https://github.com/hammerlab/varcode/tree/v0.4.18) (2016-08-08)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.17...v0.4.18)

**Closed issues:**

- vcf unit tests broken in python 3 \(?\) [\#164](https://github.com/hammerlab/varcode/issues/164)
- maximum recursion depth exceeded when loading a vcf from a URL [\#163](https://github.com/hammerlab/varcode/issues/163)

**Merged pull requests:**

- In load\_vcf, when passed a URL download it first to a local file then… [\#165](https://github.com/hammerlab/varcode/pull/165) ([timodonnell](https://github.com/timodonnell))
- Removed Collection from varcode, moved to separate 'sercol' repo instead [\#162](https://github.com/hammerlab/varcode/pull/162) ([iskandr](https://github.com/iskandr))

## [v0.4.17](https://github.com/hammerlab/varcode/tree/v0.4.17) (2016-08-05)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.16...v0.4.17)

**Merged pull requests:**

- Commandline interface, simplified serialization, merging VariantCollections [\#161](https://github.com/hammerlab/varcode/pull/161) ([iskandr](https://github.com/iskandr))

## [v0.4.16](https://github.com/hammerlab/varcode/tree/v0.4.16) (2016-07-30)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.15...v0.4.16)

## [v0.4.15](https://github.com/hammerlab/varcode/tree/v0.4.15) (2016-07-15)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.14...v0.4.15)

**Fixed bugs:**

- Fix clone\_with\_new\_elements for VariantCollection [\#159](https://github.com/hammerlab/varcode/pull/159) ([tavinathanson](https://github.com/tavinathanson))

**Closed issues:**

- load\_vcf\_fast fails when sample names contain spaces [\#158](https://github.com/hammerlab/varcode/issues/158)

**Merged pull requests:**

- Fix load\_vcf\_fast for sample names containing a space character [\#160](https://github.com/hammerlab/varcode/pull/160) ([timodonnell](https://github.com/timodonnell))

## [v0.4.14](https://github.com/hammerlab/varcode/tree/v0.4.14) (2016-06-07)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.12...v0.4.14)

**Merged pull requests:**

- Don't memoize so much [\#157](https://github.com/hammerlab/varcode/pull/157) ([iskandr](https://github.com/iskandr))

## [v0.4.12](https://github.com/hammerlab/varcode/tree/v0.4.12) (2016-05-28)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.13...v0.4.12)

## [v0.4.13](https://github.com/hammerlab/varcode/tree/v0.4.13) (2016-05-28)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.11...v0.4.13)

**Merged pull requests:**

- Fix versioneer prefix format [\#155](https://github.com/hammerlab/varcode/pull/155) ([armish](https://github.com/armish))

## [v0.4.11](https://github.com/hammerlab/varcode/tree/v0.4.11) (2016-05-27)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.10...v0.4.11)

## [v0.4.10](https://github.com/hammerlab/varcode/tree/v0.4.10) (2016-05-27)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.9...v0.4.10)

## [v0.4.9](https://github.com/hammerlab/varcode/tree/v0.4.9) (2016-05-27)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.8...v0.4.9)

**Closed issues:**

- Add serialization for EffectCollection and VariantCollection [\#71](https://github.com/hammerlab/varcode/issues/71)

**Merged pull requests:**

- Reorganized coding effects to use KnownAminoAcidChange base class [\#153](https://github.com/hammerlab/varcode/pull/153) ([iskandr](https://github.com/iskandr))

## [v0.4.8](https://github.com/hammerlab/varcode/tree/v0.4.8) (2016-05-27)
[Full Changelog](https://github.com/hammerlab/varcode/compare/v0.4.2...v0.4.8)

**Fixed bugs:**

- Potentially wrong translated sequence from frameshift on mm10 [\#151](https://github.com/hammerlab/varcode/issues/151)

**Closed issues:**

- ExonicSpliceSite mutations are classified as Noncoding [\#136](https://github.com/hammerlab/varcode/issues/136)
- Filter field is not saved after loading a VCF [\#89](https://github.com/hammerlab/varcode/issues/89)
- investigate porting read evidence module to use impala [\#69](https://github.com/hammerlab/varcode/issues/69)
- Attach genotypes and other sample information to Variants [\#30](https://github.com/hammerlab/varcode/issues/30)
- support determining the evidence for a variant in a bam [\#26](https://github.com/hammerlab/varcode/issues/26)

**Merged pull requests:**

- Added unit tests for Klf6 frameshift, fix bug in frameshift translation [\#152](https://github.com/hammerlab/varcode/pull/152) ([iskandr](https://github.com/iskandr))
- Add as\_dataframe to EffectCollection [\#150](https://github.com/hammerlab/varcode/pull/150) ([arahuja](https://github.com/arahuja))
- Use versioneer to manage version number [\#149](https://github.com/hammerlab/varcode/pull/149) ([arahuja](https://github.com/arahuja))
- Fix pyvcf error from passing \_parse\_samples a tuple instead of a list [\#148](https://github.com/hammerlab/varcode/pull/148) ([timodonnell](https://github.com/timodonnell))
- Fix variant pickling [\#146](https://github.com/hammerlab/varcode/pull/146) ([tavinathanson](https://github.com/tavinathanson))
- Parse and expose sample info, including for multisample VCFs [\#145](https://github.com/hammerlab/varcode/pull/145) ([timodonnell](https://github.com/timodonnell))
- Preserve contig name [\#140](https://github.com/hammerlab/varcode/pull/140) ([iskandr](https://github.com/iskandr))
- Quotes around nucleotides in Variant representation [\#139](https://github.com/hammerlab/varcode/pull/139) ([iskandr](https://github.com/iskandr))
- added is\_deletion, is\_insertion, and is\_indel properties to variants [\#138](https://github.com/hammerlab/varcode/pull/138) ([iskandr](https://github.com/iskandr))

## [v0.4.2](https://github.com/hammerlab/varcode/tree/v0.4.2) (2016-02-25)
**Implemented enhancements:**

- VariantCollection.high\_priority\_effect != Variant.top\_effect [\#58](https://github.com/hammerlab/varcode/issues/58)
- Improves the documentation for varcode [\#110](https://github.com/hammerlab/varcode/pull/110) ([armish](https://github.com/armish))
- Convert effect-type section into a sorted table [\#104](https://github.com/hammerlab/varcode/pull/104) ([armish](https://github.com/armish))
- Start highlighting Python syntax in README [\#103](https://github.com/hammerlab/varcode/pull/103) ([armish](https://github.com/armish))

**Fixed bugs:**

- Varcode requires pandas \>= 0.13.1, however it uses 0.15 functionality \#12 [\#92](https://github.com/hammerlab/varcode/issues/92)
- Varcode version 0.3.10 cannot be imported when installed through pip [\#90](https://github.com/hammerlab/varcode/issues/90)
- pip installing Varcode doesn't seem to work lately [\#84](https://github.com/hammerlab/varcode/issues/84)
- AttributeError: 'FrameShiftTruncation' object has no attribute 'aa\_alt' [\#70](https://github.com/hammerlab/varcode/issues/70)
- Use find\_packages correctly [\#85](https://github.com/hammerlab/varcode/pull/85) ([tavinathanson](https://github.com/tavinathanson))

**Closed issues:**

- memoize a bit less [\#131](https://github.com/hammerlab/varcode/issues/131)
- Intragenic variants do not have a short\_description field [\#129](https://github.com/hammerlab/varcode/issues/129)
- move read\_evidence module and Locus class to varlens [\#124](https://github.com/hammerlab/varcode/issues/124)
- Support Structural Variants [\#122](https://github.com/hammerlab/varcode/issues/122)
- PrematureStop called as Silent [\#116](https://github.com/hammerlab/varcode/issues/116)
- PrematureStop called as a Deletion [\#111](https://github.com/hammerlab/varcode/issues/111)
- UnboundLocalError in in\_frame\_coding\_effect.py [\#107](https://github.com/hammerlab/varcode/issues/107)
- Double mutations in a MAF file cause error [\#105](https://github.com/hammerlab/varcode/issues/105)
- varcode.load\_vcf\_fast used 0.16.1 Pandas options [\#101](https://github.com/hammerlab/varcode/issues/101)
- Configuring datacache default cache directory [\#98](https://github.com/hammerlab/varcode/issues/98)
- Improve the README to include some examples of working with Varcode in IPython [\#95](https://github.com/hammerlab/varcode/issues/95)
- support loading VCFs over HTTP [\#91](https://github.com/hammerlab/varcode/issues/91)
- Travis should include setup.py testing [\#86](https://github.com/hammerlab/varcode/issues/86)
- Make Variants pickle-able [\#77](https://github.com/hammerlab/varcode/issues/77)
- modifies\_coding\_sequence is always false [\#64](https://github.com/hammerlab/varcode/issues/64)
- AssertionError: aa\_ref and aa\_alt can't both be empty string [\#63](https://github.com/hammerlab/varcode/issues/63)
- Too many open files on error on getting top effect [\#62](https://github.com/hammerlab/varcode/issues/62)
- KeyError: 'reference' in load\_vcf [\#60](https://github.com/hammerlab/varcode/issues/60)
- Issue with n\_skip? [\#56](https://github.com/hammerlab/varcode/issues/56)
- Optional random seed argument for generating random variants [\#48](https://github.com/hammerlab/varcode/issues/48)
- An argument for using == and not \>= for requirements? [\#43](https://github.com/hammerlab/varcode/issues/43)
- deploy a test coverage tool [\#38](https://github.com/hammerlab/varcode/issues/38)
- Replace raise\_on\_error parameter to property of VariantCollection [\#36](https://github.com/hammerlab/varcode/issues/36)
- assertion error in infer\_coding\_effect [\#33](https://github.com/hammerlab/varcode/issues/33)
- add a memoized "highest\_priority\_effect" property to Variant [\#31](https://github.com/hammerlab/varcode/issues/31)
- support deep reloading varcode module [\#25](https://github.com/hammerlab/varcode/issues/25)
- handle multiallelic variants [\#22](https://github.com/hammerlab/varcode/issues/22)
- vcf.load\_vcf should provide an option to load all variants, regardless of whether filter is PASS [\#21](https://github.com/hammerlab/varcode/issues/21)
- empty variant collection when loading strelka vcf [\#16](https://github.com/hammerlab/varcode/issues/16)
- Incorrect handling of variants which run past the beginning/end of an exon's boundary [\#14](https://github.com/hammerlab/varcode/issues/14)
- Reference amino acid sequence sometimes empty for coding variants [\#12](https://github.com/hammerlab/varcode/issues/12)
- handle single-sample VCFs with INFO fields containing list values of size \> 1 [\#9](https://github.com/hammerlab/varcode/issues/9)
- Do FrameShift \(or StopGain\) mutations affect splicing? [\#6](https://github.com/hammerlab/varcode/issues/6)
- What to do with mutations that span the 5' UTR / CDS boundary? [\#5](https://github.com/hammerlab/varcode/issues/5)
- Annotate essential splice site mutations [\#1](https://github.com/hammerlab/varcode/issues/1)

**Merged pull requests:**

- Modest change to filtering of coding mutations include ExonicSpliceSite  [\#137](https://github.com/hammerlab/varcode/pull/137) ([iskandr](https://github.com/iskandr))
- Version bump [\#135](https://github.com/hammerlab/varcode/pull/135) ([tavinathanson](https://github.com/tavinathanson))
- Fix conda install on Travis [\#134](https://github.com/hammerlab/varcode/pull/134) ([iskandr](https://github.com/iskandr))
- Don't memoize EffectCollection.top\_priority\_effect\(\) [\#132](https://github.com/hammerlab/varcode/pull/132) ([timodonnell](https://github.com/timodonnell))
- All effects should have a default `short\_description` field [\#130](https://github.com/hammerlab/varcode/pull/130) ([armish](https://github.com/armish))
- Remove read\_evidence and locus modules [\#125](https://github.com/hammerlab/varcode/pull/125) ([timodonnell](https://github.com/timodonnell))
- Include a link to the iPython notebook in README.md [\#121](https://github.com/hammerlab/varcode/pull/121) ([armish](https://github.com/armish))
- Add varcode to Travis [\#120](https://github.com/hammerlab/varcode/pull/120) ([tavinathanson](https://github.com/tavinathanson))
- Minor problem in Variant.\_\_init\_\_ [\#119](https://github.com/hammerlab/varcode/pull/119) ([iskandr](https://github.com/iskandr))
- Update Varcode to work with new multi-species PyEnsembl  [\#118](https://github.com/hammerlab/varcode/pull/118) ([iskandr](https://github.com/iskandr))
- Fix \#116 and call PrematureStop when stop codon is added in the middle of an insertion [\#117](https://github.com/hammerlab/varcode/pull/117) ([leekaiinthesky](https://github.com/leekaiinthesky))
- Warn when variants in MAF file have wrong end position [\#115](https://github.com/hammerlab/varcode/pull/115) ([iskandr](https://github.com/iskandr))
- Bump pyensembl/varcode version [\#114](https://github.com/hammerlab/varcode/pull/114) ([tavinathanson](https://github.com/tavinathanson))
- fix logic for determining whether the protein length decreases [\#112](https://github.com/hammerlab/varcode/pull/112) ([leekaiinthesky](https://github.com/leekaiinthesky))
- decreasing 3' splice site to distance 3 from boundary [\#109](https://github.com/hammerlab/varcode/pull/109) ([iskandr](https://github.com/iskandr))
- fixed typo in effect inference, added breaking variant to unit tests [\#108](https://github.com/hammerlab/varcode/pull/108) ([iskandr](https://github.com/iskandr))
- Allow Varcode to work with mouse data via Genome [\#106](https://github.com/hammerlab/varcode/pull/106) ([tavinathanson](https://github.com/tavinathanson))
- Manually set compression in read\_vcf\_into\_dataframe [\#102](https://github.com/hammerlab/varcode/pull/102) ([timodonnell](https://github.com/timodonnell))
- Added examples to README [\#100](https://github.com/hammerlab/varcode/pull/100) ([iskandr](https://github.com/iskandr))
- depend on pandas \>= 0.15 [\#99](https://github.com/hammerlab/varcode/pull/99) ([iskandr](https://github.com/iskandr))
-  Faster VCFs loading, support HTTP, and refactored variant metadata [\#94](https://github.com/hammerlab/varcode/pull/94) ([timodonnell](https://github.com/timodonnell))
- Support for regular varcode variant instances in read evidence module [\#87](https://github.com/hammerlab/varcode/pull/87) ([timodonnell](https://github.com/timodonnell))
- Read and write json files [\#82](https://github.com/hammerlab/varcode/pull/82) ([iskandr](https://github.com/iskandr))
- JSON serialization for VariantCollection. [\#81](https://github.com/hammerlab/varcode/pull/81) ([timodonnell](https://github.com/timodonnell))
- Add short\_description field to intergenic variants [\#80](https://github.com/hammerlab/varcode/pull/80) ([timodonnell](https://github.com/timodonnell))
- Speed up PileupCollection.group\_by\_allele [\#79](https://github.com/hammerlab/varcode/pull/79) ([timodonnell](https://github.com/timodonnell))
- Variant serialization [\#78](https://github.com/hammerlab/varcode/pull/78) ([timodonnell](https://github.com/timodonnell))
- added option for genome name in load\_vcf [\#76](https://github.com/hammerlab/varcode/pull/76) ([iskandr](https://github.com/iskandr))
- Fix variant.effects\(\) to always return an EffectCollection [\#75](https://github.com/hammerlab/varcode/pull/75) ([timodonnell](https://github.com/timodonnell))
- Bump pysam dependency [\#74](https://github.com/hammerlab/varcode/pull/74) ([timodonnell](https://github.com/timodonnell))
- Cufflinks RNA filtering [\#73](https://github.com/hammerlab/varcode/pull/73) ([iskandr](https://github.com/iskandr))
- Read evidence tweaks [\#72](https://github.com/hammerlab/varcode/pull/72) ([timodonnell](https://github.com/timodonnell))
- Filter effect collection by expression [\#67](https://github.com/hammerlab/varcode/pull/67) ([iskandr](https://github.com/iskandr))
- Created EpitopeCollection, refactored effects, fix assertion failure while annotating silent stop codon [\#66](https://github.com/hammerlab/varcode/pull/66) ([iskandr](https://github.com/iskandr))
- Created EpitopeCollection, refactored effects [\#65](https://github.com/hammerlab/varcode/pull/65) ([iskandr](https://github.com/iskandr))
- include substitution in high priority effects [\#61](https://github.com/hammerlab/varcode/pull/61) ([arahuja](https://github.com/arahuja))
- don't annotate StopLoss variants that are immediately followed by another stop codon [\#57](https://github.com/hammerlab/varcode/pull/57) ([iskandr](https://github.com/iskandr))
- Refactor coding effect [\#55](https://github.com/hammerlab/varcode/pull/55) ([iskandr](https://github.com/iskandr))
- Add read\_evidence module [\#53](https://github.com/hammerlab/varcode/pull/53) ([timodonnell](https://github.com/timodonnell))
- Use transcript protein sequence [\#45](https://github.com/hammerlab/varcode/pull/45) ([iskandr](https://github.com/iskandr))
- Add contributing md [\#41](https://github.com/hammerlab/varcode/pull/41) ([iskandr](https://github.com/iskandr))
- Small coding effect refactoring and fixes [\#39](https://github.com/hammerlab/varcode/pull/39) ([iskandr](https://github.com/iskandr))
- Test problematic variants [\#37](https://github.com/hammerlab/varcode/pull/37) ([iskandr](https://github.com/iskandr))
- Typechecks and test fixes [\#35](https://github.com/hammerlab/varcode/pull/35) ([timodonnell](https://github.com/timodonnell))
- Fix maf parsing [\#34](https://github.com/hammerlab/varcode/pull/34) ([iskandr](https://github.com/iskandr))
- parse multiple alleles into distinct Variant records [\#29](https://github.com/hammerlab/varcode/pull/29) ([iskandr](https://github.com/iskandr))
- PEP8 & pyflakes fixes [\#28](https://github.com/hammerlab/varcode/pull/28) ([iskandr](https://github.com/iskandr))
- Remove pyfaidx [\#27](https://github.com/hammerlab/varcode/pull/27) ([iskandr](https://github.com/iskandr))
- Variant collection tweaks [\#24](https://github.com/hammerlab/varcode/pull/24) ([timodonnell](https://github.com/timodonnell))
- Improved vcf parsing [\#23](https://github.com/hammerlab/varcode/pull/23) ([timodonnell](https://github.com/timodonnell))
- Associate EnsemblRelease with each Variant object [\#20](https://github.com/hammerlab/varcode/pull/20) ([iskandr](https://github.com/iskandr))
- Variant collection filtering [\#19](https://github.com/hammerlab/varcode/pull/19) ([iskandr](https://github.com/iskandr))
- added IntronicSpliceSite, SpliceDonor, SpliceAcceptor effects [\#17](https://github.com/hammerlab/varcode/pull/17) ([iskandr](https://github.com/iskandr))
- collect effect annotation errors in dictionary, only look up overlapping... [\#13](https://github.com/hammerlab/varcode/pull/13) ([iskandr](https://github.com/iskandr))
- don't flatten INFO dictionary of VCF, lists are part of the field format [\#11](https://github.com/hammerlab/varcode/pull/11) ([iskandr](https://github.com/iskandr))
- Small fixes [\#10](https://github.com/hammerlab/varcode/pull/10) ([timodonnell](https://github.com/timodonnell))
- Add support for Python 3 [\#8](https://github.com/hammerlab/varcode/pull/8) ([timodonnell](https://github.com/timodonnell))
- Refactor core logic [\#7](https://github.com/hammerlab/varcode/pull/7) ([iskandr](https://github.com/iskandr))
- Classes for protein/transcript variant effects [\#3](https://github.com/hammerlab/varcode/pull/3) ([iskandr](https://github.com/iskandr))



\* *This Change Log was automatically generated by [github_changelog_generator](https://github.com/skywinder/Github-Changelog-Generator)*
