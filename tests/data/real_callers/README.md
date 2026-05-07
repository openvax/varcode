# Real-world somatic-caller VCF fixtures

Tiny fixtures exercising what production cancer-pipeline callers actually
emit. Tests in `tests/test_real_world_callers.py` pin our parser against
hand-verified expected values; the parametrized oracle tests in
`tests/test_vcf_parsing.py` additionally compare every row against PyVCF3
element-by-element.

| File | Provenance | What it exercises |
|---|---|---|
| `mutect2_example.vcf` | **Verbatim copy** from `broadinstitute/gatk` | Per-alt arrays (`TLOD`/`POPAF` Number=A), per-allele arrays (`AD`/`F1R2`/`F2R1`/`MBQ`/`MFRL` Number=R), `RPA` Number=., `STR` Flag, quad-allelic site with four-way GT `0/1/2/3`, mitochondrial chrM contig, `##Mutect Version=2.1` metadata key with embedded whitespace, declared-but-unused FORMAT fields (`PGT`/`PID`/`SAAF`/`SAPP`) |
| `strelka2_somatic_snvs.vcf` | **Format-faithful reconstruction** following Strelka2's documented output | Per-base tier counts (`AU:CU:GU:TU` Number=2), per-sample format with **no GT field**, `SGT=GG->AG` notation, `##cmdline=` metadata with paths and CLI flags |
| `vep_annotated_csq.vcf` | **Format-faithful reconstruction** following Ensembl VEP's documented `--vcf` output | `CSQ` Number=. with multi-transcript comma-split, per-record pipe-delimited subfields (23 in this fixture's format spec), HGVS notation with embedded `>`/`:`/`/` and percent-encoded `%3D`, Description containing the format spec |

## Provenance details

### `mutect2_example.vcf` — verbatim from GATK

- **Source**: <https://github.com/broadinstitute/gatk>
- **Path**: `src/test/resources/org/broadinstitute/hellbender/tools/mutect/mito/unfiltered.vcf`
- **Pinned commit**: `3021e6924aeb84d9f3b333e5298abb1ec27d350a` (April 2020)
- **License**: BSD 3-Clause (compatible with varcode's Apache 2.0)
- **Direct URL**: <https://raw.githubusercontent.com/broadinstitute/gatk/3021e6924aeb84d9f3b333e5298abb1ec27d350a/src/test/resources/org/broadinstitute/hellbender/tools/mutect/mito/unfiltered.vcf>

7 records, single sample (NA12878), GATK4 MuTect2 mitochondrial DNA test
output. Reproduced byte-for-byte; the `##Mutect Version=2.1` metadata line
with a *whitespace inside the key* is a real example of why our regex uses
`.+?` (non-greedy) for the key match.

### `strelka2_somatic_snvs.vcf` — format-faithful reconstruction

Strelka2 ships under GPLv3, which is incompatible with varcode's Apache 2.0
license — we cannot bundle a verbatim excerpt. The fixture is a
reconstruction following Strelka2's documented output schema in:

- **Spec**: <https://github.com/Illumina/strelka/blob/master/docs/userGuide/README.md>
- **Sample data shape**: documented in
  <https://github.com/Illumina/strelka/tree/master/docs/userGuide#somatic>

Header fields and per-sample format match the schema verbatim; record
contents (positions, depths, SGT values) are synthetic but plausible. The
load-bearing format quirk this fixture pins — Strelka2 SNV records have **no
GT field** in their per-sample format — is a published characteristic of
the tool, not an invention.

### `vep_annotated_csq.vcf` — format-faithful reconstruction

Despite an extensive search across `Ensembl/ensembl-vep` test data,
`Ensembl/VEP_plugins`, `nf-core/test-datasets` (Sarek branch), `sigven/pcgr`,
`ga4gh/vrs-python`, `Illumina/hap.py`, `broadinstitute/gatk-sv`, and other
public bioinformatics corpora, no small public VEP-annotated VCF with a
`CSQ` field could be located: VEP-annotated VCFs are pipeline *outputs*, not
artifacts that anyone curates as test fixtures.

This fixture is therefore a reconstruction following Ensembl VEP's published
output format:

- **Format spec**: <https://www.ensembl.org/info/docs/tools/vep/vep_formats.html#vcfout>
- The 23 pipe-delimited subfield order
  (`Allele|Consequence|IMPACT|SYMBOL|Gene|...|HGNC_ID`) matches what VEP
  emits in default `--vcf` mode.
- The two real anchors in the records are biologically accurate: rs28934578
  is a real TP53 R282\* hotspot ClinVar variant, and ENST00000269305 is the
  canonical TP53 transcript.
- Position/depth/CSQ field-content values are otherwise synthetic but
  format-faithful.

If a small, public, real VEP-annotated VCF surfaces in a permissively-licensed
corpus, this fixture should be replaced with a verbatim excerpt at that point.
