# Public somatic benchmarks

These include complete, pinned benchmark call sets for the two requested human
cell lines, the complete published CanFam3.1 canine consequence table, and
selected native ROS_Cfam_1.0 canine calls. They are not an exhaustive union
of every historical caller output. Sets overlap;
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
| `canine-ros-native.json` | 10 alleles / 11 observations, selected | ROS_Cfam_1.0, Ensembl 2020_10 |

`manifest.json` records the complete historical sets' source URLs, content
hashes, transformations, assemblies, releases and row counts. The curated
ROS_Cfam fixture embeds its own source and reference metadata.
`../SHA256SUMS` additionally covers every
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

## Native ROS_Cfam_1.0 canine variants

`canine-ros-native.json` adds **10 distinct alleles / 11 tumor observations**
on ROS_Cfam_1.0 (GCA_014441545.1). These are selected published observations,
not a complete cohort or a somatic truth set. They were called on this
assembly directly; no CanFam3.1 liftover is involved.

| Source | Coverage | Evidence |
|---|---|---|
| [Matsumura et al., 2025](https://doi.org/10.3389/fvets.2025.1535446) | All 4 ZDBF2 alleles in Table 4, in 5 tumor observations | Targeted sequencing; filtered against available blood controls |
| [Tsuruta et al., 2024](https://doi.org/10.1292/jvms.24-0096) | All 5 STAT3/JAK1 alleles in Tables 1–2 | Three pass the main VAF threshold; two remain explicitly below-threshold candidates. STAT3 Y640F was also detected by the hotspot assay |
| [Tagawa et al., 2026](https://doi.org/10.1292/jvms.25-0492) | One selected SPTBN1 deletion from the 15-entry Table 1 | Matched normal plus ddPCR confirmation in tumor DNA; published p.Asn795fs |

The melanoma table omits chromosomes and REF bases. For SPTBN1, the pinned
annotation identifies chromosome 10 and Ensembl genomic DNA supplies the
deleted `AC` at 56,527,053–56,527,054 (1-based inclusive). This derivation is
explicit in the fixture. The other validated marker, ARID1A, has insufficiently
specified insertion anchoring in the table and is not imported. This benchmark
does not claim coverage of the other melanoma calls.

The current Ensembl dog annotation is **2020_10**, checked against its
[published annotation dates](https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/014/441/545/1/ensembl/)
on 2026-10-20. This date records when the annotation was built.
The fixture pins GTF, cDNA and peptide SHA256s. Small genomic windows from
Ensembl REST preserve the assembly, interval, sequence and sequence SHA256;
tests require no network or full chromosome FASTA. These windows do not
provide full-genome SV coverage.

The papers do not specify transcript IDs for these records. Tests select
complete current models compatible with each published coding position,
check the published consequence, and independently translate the edited RNA
using Biopython. For SNVs, the expected edit comes from the published cDNA
position; for SPTBN1, it comes from the stated genomic deletion. Only STAT3
Y640F and SPTBN1 Asn795fs have explicit protein descriptions in these selected
source statements. Other peptide changes and the full frameshift tail are
computed predictions, not independently published peptide sequences.

The original article XML URLs and SHA256s are recorded alongside each source.
The fixture contains curated factual variant observations; articles retain
their original licenses (Matsumura: CC BY 4.0; Tsuruta and Tagawa: CC BY-NC-ND
4.0). Reacquire and check those statements and genomic windows with:

```sh
python tests/data/public_somatic/check_ros_sources.py /tmp/ros-sources --download
python -m pytest -q tests/test_canine_ros.py
```

Omit `--download` to verify an existing source directory offline. An upstream
source checksum change needs review; do not refresh expectations automatically.

## Regeneration and attribution: complete historical benchmarks

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
