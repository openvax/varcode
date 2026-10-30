# Native reference consequence fixtures

These are small, unrenumbered slices of **GRCm39 / Ensembl 115**,
**ROS_Cfam_1.0 / Ensembl 2020_10**, and **T2T-CHM13v2.0 / Ensembl rapid
2022_07**. The three JSON manifests retain assembly accession, species,
translation metadata, and SHA-256 checksums of the original GTF/cDNA/protein
files. GTF rows and FASTA records retain their native identifiers and sequences.
These slices are test references, not complete genome installations.

`cases.json` has 100 variant/transcript cases. Expected full coding sequences
and proteins were calculated from native RNA with a separate Biopython edit and
translation calculation, never Varcode predictions. Published protein positions
are additional constraints where available. A full expected sequence calculated
from a public variant is not a measured mutant protein.

## Public evidence

| Cases | Source and interpretation |
|---|---|
| Mouse Atp7a brindled | [MGI:1856098](https://www.alliancegenome.org/allele/MGI%3A1856098), [PMID 9147645](https://pubmed.ncbi.nlm.nih.gov/9147645/): curated GRCm39 placements of A514T and L801_A802del on the same allele. The joint test uses that allele-level cis evidence. |
| Mouse mt-Nd6 m1Dwa | [MGI:3777241](https://www.informatics.jax.org/allele/MGI%3A3777241), [Fan et al. 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC3049809/): minus-strand mitochondrial insertion. The original LA9/C3H description says stop at codon 79; translation of the pinned GRCm39 sequence gives 79 residues followed by a stop at codon 80. The discrepancy is retained, not silently reconciled. This is not proof that the source strain and reference mitochondrial haplotypes are identical. |
| Mouse Ap2s1 del17 | [PMID 29479578](https://pmc.ncbi.nlm.nih.gov/articles/PMC5824975/): RT-PCR and sequencing establish a 51-nt exon-2 deletion (amino acids 35–51). The test completes the rest from reference RNA and marks that completion; it does not assert observed full-length RNA. DNA-only cryptic-site discovery misses this junction: [#675](https://github.com/openvax/varcode/issues/675). |
| Dog GRIK2, VPS54, CCDC85A, BRINP3 | [Animals 13:810](https://doi.org/10.3390/ani13050810), Table 2: GRIK2 has native ROS_Cfam Sanger comparison. The other calls were originally obtained on CanFam3.1 and have published ROS_Cfam placements. Their origins remain explicit. |
| Dog FOXI3 duplication | [OMIA 000323](https://omia.org/OMIA000323/9615/), [PMID 18787161](https://pubmed.ncbi.nlm.nih.gov/18787161/): curated ROS_Cfam placement of the 7-bp duplication. Ensembl's protein has a 152-aa N-terminal extension relative to RefSeq NP_001129118.1; the entire shared sequence is checked. RefSeq p.A23Rfs*219 is therefore not copied verbatim onto the Ensembl isoform. |
| Dog tumour calls | [PLOS ONE 2024, doi:10.1371/journal.pone.0307792](https://doi.org/10.1371/journal.pone.0307792), supplementary S3A/S3B: native ROS_Cfam WGS. All **101 rows / 96 distinct variants** remain in `dog-tumour-source-audit.json`, including filters, genotypes and source row numbers. **74 rows** reproduce the published amino-acid change on a complete native transcript; four more provide splice predictions. Other rows have explicit model/position/format limitations. Source filters are retained, including low-confidence splice calls. These are published predictions, not validated mutant protein measurements. `0|1` without shared phase-set evidence is not treated as cis across rows. |
| CHM13 SweGen | [Genome Research, doi:10.1101/gr.279320.124](https://doi.org/10.1101/gr.279320.124), supplemental Table S5: nine selected native CHM13 calls marked PacBio CLR validated. Validation establishes the allele, not the protein consequence. The paper used SnpEff/RefSeq v5.1, while these tests use pinned Ensembl. ERAP1's insertion is at an Ensembl exon–intron boundary and DHRS4L2's deleted base is intronic; their unchanged proteins are conditional on normal splicing. TRIM51G and DHRS4L2 overlap the consortium segmental-duplication track. |
| CHM13 newly resolved sequence | [T2T consortium native 1KGP calls](https://github.com/marbl/CHM13): two variants at chr22:21798486 and chr22:21801591, each on RIMBP3B and RIMBP3C. Original VCF rows and hashes are retained. Both loci overlap the consortium `chm13v2-unique_to_hg38.bed` track. Distinct transcript identities must survive even though their coding sequences coincide. These population calls are not individually experimentally validated here. |

`structural-cases.json` adds two public curated deletion placements: mouse
[Pax3 Sp2H](https://archive.har.mrc.ac.uk/stock?id=239) (PMID 1682057), whose 32-nt
deletion gives the independently established 237-aa truncation, and dog
[COL9A2 Samoyed](https://omia.org/OMIA001523/9615/) (PMID 20686772), which removes
the annotated start. The COL9A2 endpoints span 1271 bp while the original report
says 1267 bp; the test verifies start removal without claiming exact patient
breakpoint sequence. The CHM13 structural deletion is a **controlled edit of
real reference sequence**, not a patient rearrangement. Inversion/fusion truth
is not supplied by these local deletion cases.

`dog-mito.json` adds four public mitochondrial variants from
[Kowal et al. 2024, Table 5](https://doi.org/10.3390/ijms252413414), independently
translated against the exact NCBI **NC_002008.4** sequence and CDS annotations.
ROS_Cfam Ensembl 2020_10 contains no mitochondrial transcripts, so this is a
separately identified reference, never relabeled ROS_Cfam. COX2 M187T, ND5 M74T
and ATP6 P136S agree with the paper. The paper's ND6 I106V label differs from
I107V on its stated NC_002008.4 reference; that discrepancy is explicit. ND6
tests the minus strand; ND5 tests the vertebrate mitochondrial ATA initiator.

## DNA and provenance

`dna-windows.json.gz` stores actual chromosome windows with their original
coordinates. The sparse reader returns unavailable sequence outside those
windows; it never pads with invented reference bases.

- Mouse: NCBI `NC_000073.7:16472335-16483219`, fetched with EFetch.
- Dog: NCBI ROS_Cfam accessions `NC_051824.1` (20), `NC_051805.1` (1),
  `NC_051836.1` (32); individual EFetch URLs accompany the windows.
- CHM13: consortium `chm13v2.0.fa.gz`; its `chr` names are explicitly recorded
  as `source_contig` and mapped to the pinned annotation's names without changing
  coordinates. The documented full installation uses the matching Ensembl FASTA.

Source papers and supplementary data require attribution. PLOS and the cited
Genome Research article are CC BY; these are derived frozen regression fixtures,
not updated copies of the articles. Consortium data are released as CC0; NCBI
records retain their accessions. Ensembl data and source metadata retain their
original attribution.

## Rebuilding

`notes/build_native_reference_fixtures.py` requires Biopython, openpyxl, PyEnsembl
2.32.3+, ORFdata 0.43+, installed mouse 115 and dog 2020_10 references, and the
pinned CHM13 installation. It never imports Varcode. Supply downloaded originals:

- `dog-carcinoma.xlsx` and `dog-sarcoma.xlsx`: the PLOS article's `.s005` and
  `.s006` supplementary files.
- `swedish-s5.json`: rows from `Supplemental Tables S5` of
  [Supplemental_Tables.xlsx](https://pmc-oa-opendata.s3.amazonaws.com/PMC12581843.1/Supplemental_Tables.xlsx),
  retaining native position/alleles/gene/sample/validation fields.
- `chm13-new-calls.vcf`: tabix query of the consortium 1KGP CHM13 chromosome-22
  allele-frequency VCF for `chr22:21797405-21802321` (not a lifted GRCh38 VCF).
- `chm13v2.0_SD.bed`, `chm13v2-unique_to_hg38.bed`: original consortium tracks.
- `mouse-ap2s1-dna.fa`, `dog-splice-dna.json`, `dog-FOXI3-NP_001129118.1.fa`:
  NCBI reference sequences, with original genomic coordinates as above.
- `dog-mito-NC_002008.4.gb`: NCBI EFetch `rettype=gbwithparts` for NC_002008.4.

```sh
python notes/build_native_reference_fixtures.py \
  --sources /path/to/public-inputs --chm13-cache /path/to/orfdata-cache \
  --chm13-dna /path/to/chm13v2.0.fa.gz --dna-contig-prefix chr
```

Tests consume only the committed slices and do not download these sources.
