# Documentation and code quality audit — 2026-09-24

This audit started at `525499d` (Varcode 10.1.1). It combines a broad static
review with targeted execution of parsing, metadata, failure handling, and
phase-enumeration paths. The accompanying 10.1.2 maintenance change contains
bounded corrections; larger export and scientific-model problems are tracked
separately below.

## Coverage and method

The repository contained 67 library Python modules, 106 test/helper modules,
and 26 documentation pages before this work. The audit checked:

- README, contributor/release instructions, test-data provenance, navigation,
  troubleshooting, API links, and workflow examples against code and CI.
- VCF/MAF loading and VCF export: mixed allele types, limits, multi-allelic
  indexes, sample/FORMAT ordering, normalization, and round-trip assertions.
- Collection construction, mutable defaults, annotation error reporting,
  lint configuration, test interpreter selection, and fixture construction.
- Germline phase enumeration and its consumer, including an end-to-end CFTR
  reproduction where changing the cap changes the reported consequence.
- Existing issue and PR coverage, to avoid duplicating known problems.

The static pass checked undefined names and nearby disabled lint rules.
Reproductions used inline synthetic files or explicitly pinned Ensembl 81.
The audit does not establish biological validity for every annotator or
exhaustively review every source line. Large consequence modules still need
domain-specific reviews and independent scientific oracles.

## Major findings requiring dedicated work

| Priority | Finding and evidence | Disposition |
|---|---|---|
| High | VCF export derives header sample order from a set but writes row values in dictionary order. A two-row reproduction assigns tumor GT/DP to normal and emits `44:0/1` under `GT:DP`. Differing sample sets also discard sample data; headers and allele-indexed metadata lack a complete preservation contract. | Filed [#502](https://github.com/openvax/varcode/issues/502); documented in [saving results](csv.md#vcf-export-limitations). |
| High | The phase cap returns an all-cis placeholder that the consumer classifies as a precise effect. The existing CFTR pair yields cis/trans alternatives at cap 8 but only `p.S159T` at cap 1, with no molecular evidence resolving phase. Four unphased germline alleles exceed the default cap. | Filed [#503](https://github.com/openvax/varcode/issues/503); documented in [germline annotation](germline.md). |
| High | `deploy.sh` lacks the branch/clean-tree guards, version handling, tagging, and pushing described in AGENTS.md. It can upload from the wrong checkout if invoked without independent checks. | Existing [#414](https://github.com/openvax/varcode/issues/414); release instructions now describe the actual script and required manual checks. |
| Medium | The test cache probe accesses `EnsemblRelease.db`, which can download missing reference data while supposedly checking local installation. | Existing [#493](https://github.com/openvax/varcode/issues/493); contributor instructions distinguish reference provisioning from offline snapshot tests. |
| Medium | Some tests cannot detect failures, including the fast-path self-comparison and protocol/attribute surveys. VCF round-trip checks omitted ALT and complete record counts. | Existing [#479](https://github.com/openvax/varcode/issues/479); this change strengthens the VCF checks. The broader test review remains open. |
| Medium | Reference mismatch fields contain transcript-oriented bases but the message labels them as genome bases; minus-strand errors can be read backwards. | Existing [#434](https://github.com/openvax/varcode/issues/434); troubleshooting now explains the current orientation. |
| Medium | Historical changelog coverage is incomplete. | Existing [#411](https://github.com/openvax/varcode/issues/411); no invented reconstruction of missing release history. |

The export review used the primary
[VCF specification](https://samtools.github.io/hts-specs/VCFv4.5.pdf), particularly
header/sample ordering and FORMAT field interpretation. The phase finding is
an observed inconsistency between unresolved evidence and the implementation's
precise output, not a proposed replacement biological model.

## Small corrections in the maintenance change

| Area | Correction | Verification |
|---|---|---|
| VCF limits | Count successfully loaded alleles consistently across small/SV calls, ALT lists and chunks; zero means none; reject invalid caps; stop without consuming a later chunk. | Synthetic mixed VCF, multiple chunk sizes and limits, explicit deduplication boundary ([#504](https://github.com/openvax/varcode/issues/504)). |
| VCF errors | Remove StopIteration as an internal early-exit mechanism so unexpected parser failures propagate; validate dataframe columns with ValueError even under optimized Python. | Deliberately failing parser and malformed dataframe regressions. |
| Structural metadata | Retain the original ALT index for structural alleles using the same metadata path as small variants. | A `C,<DEL>` / `GT=0/2` fixture selects only the deletion for the tumor sample. |
| Collection state | Allocate independent default metadata dictionaries and source sets; remove nearby mutable parser defaults. | Mutating one collection cannot change another ([#505](https://github.com/openvax/varcode/issues/505)). |
| Error diagnostics | Keep the original exception text in the transcript failure helper, matching the collection annotation path. | Real pinned-reference mismatch checks variant, transcript and error ([#506](https://github.com/openvax/varcode/issues/506)). |
| MAF normalization | Rename case-only columns in place rather than moving them to the end; remove the Python 2 workaround. | Synthetic partly lowercased header retains canonical column order ([#507](https://github.com/openvax/varcode/issues/507)). |
| Test execution | Probe and invoke pytest with the same Python interpreter; preserve arguments and process exit status. | A deliberately conflicting pytest executable with and without xdist ([#490](https://github.com/openvax/varcode/issues/490)). |
| Logging and lint | Replace deprecated `warn` calls and enable undefined-name checks that the code already passes. | Ruff on library and tests. |
| Existing test quality | Check every round-tripped ALT and record count, clean temporary exports, and materialize parametrization iterables for pytest compatibility. | Existing VCF/MAF round-trip and parsing tests. |
| Developer workflow | Add environment, reference-data, test, parity and docs-build commands; make editable installation use the active interpreter. | Compare with CI and execute documented checks. |
| Release documentation | Use `main`, the actual version-file path and real script capabilities; document clean-main publishing and artifact verification. | Inspect script and release workflow. |
| Test data documentation | Correct current Osteosarc pin/support to 0.2.3 and Python 3.9+, while preserving historical 0.1.0 collection provenance and snapshot hashes. Require reproducible fixture recipes alongside tests. | Compare optional extra, CI matrix and integration fixture. |
| User documentation | Add missing inherited-overlap API links and SV comparison navigation; clarify CLI SV defaults, Failure handling, transcript-oriented mismatch bases, and large unresolved limitations. | Strict MkDocs build and local link checks. |

All new synthetic fixture construction is checked in beside its tests. No
machine-local data bundle or bespoke external generator is required.

## Existing maintenance PRs

[#409](https://github.com/openvax/varcode/pull/409) covers broader staging-note
cleanup and a fast-path test correction. [#410](https://github.com/openvax/varcode/pull/410)
covers unused developer configuration. [#431](https://github.com/openvax/varcode/pull/431)
covers CLI error presentation. This audit does not supersede those PRs; shared
files will need rebasing. In particular, this PR does not remove legacy config
files or rewrite all experimental-annotator documentation.

## Follow-up sequence

Fix sample identity/FORMAT export integrity (#502) and phase-cap uncertainty
(#503) first. Then make release guards executable (#414) and make reference
cache checks truly local (#493). Broader uncertainty composition (#421/#423),
shared RNA fixture adoption (#464), and test-oracle cleanup (#479) remain
foundational work; they should be reviewed as their own changes.

## Validation for the maintenance change

Local verification on Python 3.12 passed `./lint.sh`, all 2,735 tests through
`./test.sh`, and the 123-test focused parser/collection/error/script set.
The strict MkDocs build passed; checking the generated HTML found no broken
files or anchors among 8,055 local links. The PR records CI and publication
verification separately, including the supported Python 3.9–3.11 matrix.
