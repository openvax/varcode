# Public somatic benchmarks

These are complete, pinned benchmark call sets for the two requested human
cell lines, plus the complete published canine consequence table. They are
not an exhaustive union of every historical caller output. Sets overlap;
their record counts must not be added and described as unique mutations.

| Fixture | Records | Reference used in tests |
|---|---:|---|
| `hcc1395-seqc2-snv.vcf.gz` | 39,560 | GRCh38, Ensembl 93 |
| `hcc1395-seqc2-indel.vcf.gz` | 1,922 | GRCh38, Ensembl 93 |
| `hcc1395-lancet2-two-tech.vcf.gz` | 75,825 | GRCh38, Ensembl 93 |
| `colo829-nygc-v6.vcf.gz` | 43,192 | GRCh38, Ensembl 93 |
| `colo829-lancet2-two-tech.vcf.gz` | 38,850 | GRCh38, Ensembl 93 |
| `hcc1395-consensus-sv.tsv` | 1,788 | GRCh38, Ensembl 93 |
| `colo829-sv-truth.vcf` | 134 / 68 events | GRCh37, Ensembl 75 |
| `canine-osa-vep99.tsv` | 2,141 / 1,936 unique alleles | CanFam3.1, Ensembl 99 |

`manifest.json` records exact source URLs, content hashes, transformations,
assemblies, releases and row counts. `../SHA256SUMS` additionally covers every
fixture and supporting file. The source distribution includes this directory.
No reference data or benchmark downloads occur in the test code.

## Sources and selection

**SEQC2 HCC1395:** [Fang et al., 2021](https://doi.org/10.1038/s41587-021-00993-6),
from the NCBI Somatic Mutation Working Group release `v1.2.1`, whose superset
VCFs are named `v1.2`. The fixtures contain every record with the exact
`PASS` filter token and either `HighConf` or `MedConf`. `AllPASS` alone is not
`PASS`: the superset also contains normal-sample LOH calls that are not truth
somatic calls. We retain calls outside the high-confidence benchmark intervals,
rather than trimming the truth set to the smaller evaluation regions.
The first seven VCF fields are unchanged. Large per-run metrics are omitted;
`SOURCE_ROW` gives the one-based record in the checksum-pinned original.

**NYGC / Lancet2:** [Musunuri et al., 2026](https://doi.org/10.1093/nargab/lqag036),
[versioned benchmark archive](https://doi.org/10.5281/zenodo.19191074).
The complete NYGC v6 COLO829 VCF and final two-technology VCFs are unchanged,
including records with `FILTER=.`. GCS object generations and SHA256 hashes
pin the downloaded files, which also match the versioned Zenodo archive
byte for byte (archive MD5 and member names are recorded). The COLO829 VCF's original consequence annotations
identify VEP 93 / GRCh38.p12. Tests use that release for the GRCh38 call sets;
they do not assert human protein concordance with a different annotation model.

**HCC1395 structural variants:** [Talsania et al., 2022](https://doi.org/10.1186/s13059-022-02816-6),
Supplementary Table S4. Every data row is retained: 717 DEL, 551 INS, 230 DUP,
133 INV, 146 TRA and 11 complex BND events. The final prose footnote is removed
and header whitespace stripped. The source uses GRCh38 and describes clustered
consensus breakpoints, not necessarily exact allele boundaries. Its reported
length can differ from the endpoint distance. The test adapter annotates the
reported interval, sorts reversed interval endpoints, maps TRA to BND and
preserves the entire source row. It does not derive orientation from the
undocumented ID suffix or invent insertion sequence. The 11 complex BND rows
are tested as their two reported endpoints, not as fully reconstructed events.

**COLO829 structural variants:** [Espejo Valle-Inclan et al., 2022](https://doi.org/10.1016/j.xgen.2022.100139),
[original repository](https://github.com/UMCUGenetics/COLO829_somaticSV),
commit `f9caa273af7a09e9a199e570ca3fdc037156d67d`. The unchanged GRCh37 VCF has
66 reciprocal breakend pairs and two symbolic insertions: 134 records represent
68 events. Tests annotate every record, pair the breakends, check complete
source-record accounting, and annotate every resulting event.

**Canine osteosarcoma:** [Das et al., 2021](https://doi.org/10.1038/s42003-021-02683-0),
Supplementary Data 3 in `42003_2021_2683_MOESM4_ESM.xlsx`. All 2,141 sample-level
rows are retained, including two synonymous variants. Repeated alleles in
different samples remain separate source observations. The paper specifies
CanFam3.1 and VEP 99; the tests use the exact published Ensembl transcript IDs
with dog release 99. MAF-style insertion coordinates and `-` alleles are
interpreted as in Varcode's MAF loader.

## What the tests establish

Run after the reference installation in [CONTRIBUTING](../../../CONTRIBUTING.md):

```sh
python -m pytest -q tests/test_public_somatic.py
```

Every human source allele must load with its original coordinates and alleles
and produce nonempty annotation with errors enabled. Lazy SV alternatives and
their protein properties are evaluated too. Use `only_passing=False` for these
already selected fixtures: filtering again would lose the SEQC2 compound
filters and the two-technology `.` records.

The canine comparison independently applies the **published** protein HGVS
to the reference peptide, rather than using Varcode to generate expectations:

| Comparison | Sample-level records |
|---|---:|
| Protein changes agree (including frameshift prefix / known stop) | 1,954 |
| Stop-loss agrees; the source leaves extension termination unknown | 1 |
| Splice-site disruption agrees; the source supplies no exact protein | 1 |
| Initiation outcome uncertain (`p.M1?`; Varcode recognizes alternate TTG) | 1 |
| Protein prediction unavailable: partial CDS / IG_V_gene | 184 |

For frameshifts, HGVS specifies only the unchanged prefix, first new residue
and sometimes the stop distance; it cannot validate an unspecified tail.
For exonic splice sets, VEP's coding comparison uses the explicitly conditional
`effect_if_splicing_unchanged`, preserving the original alternative outcomes.

The 184 unresolved records (181 partial coding transcripts and 3 IG_V_gene
records) are listed by original spreadsheet row in `canine-unresolved.json`.
Tests verify their transcript identity, biotype and missing start/stop
annotations. They are **not counted as concordant**, skipped or assigned fake
protein predictions. This is the existing [#195 limitation](https://github.com/openvax/varcode/issues/195).
The initiation record is also not counted as protein concordance.

Successful SV annotation does not establish an exact mutant protein. Missing
insertion sequence, unresolved orientation, consensus breakpoint uncertainty,
incomplete transcripts and absent intronic DNA can leave products unresolved.
These public-data tests use transcript references; configured chromosome-FASTA
behavior is covered separately by the genome-DNA/structural regression tests.

## Regeneration and attribution

Download each manifest `source_url` to its `source_file` in a separate directory,
then run (requires `pandas` and `openpyxl` for the two spreadsheets):

```sh
python tests/data/public_somatic/build_fixtures.py /path/to/downloads /tmp/rebuilt
```

The offline builder verifies original SHA256 hashes before reading, verifies
row counts, and checks each regenerated fixture's exact SHA256 afterward.
Unchanged VCFs are copied; the SEQC2 gzip output is deterministic. Review source
updates and expectations before replacing the manifest or checksum inventory.
`canine-unresolved.json` is a reviewed limitation ledger, not a regenerated
scientific truth set.

Data remain attributed to their original authors and repositories:

- SEQC2: public NCBI benchmark; associated article
  [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
- Canine and HCC1395 SV supplementary data: articles
  [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); the HCC1395 article
  also dedicates its data to [CC0](https://creativecommons.org/publicdomain/zero/1.0/).
- NYGC/Lancet2 archive: BSD 3-Clause in the Zenodo metadata; the authors'
  license is retained in `LANCET2-LICENSE.txt`.
- Unmodified COLO829 SV VCF: original repository GPL-3.0, retained in
  `COLO829-SV-LICENSE.txt`. This data attribution does not change Varcode's
  Apache-2.0 source-code license.

The transformations above are Varcode's; the cited authors have not endorsed
these tests or their predictions.
