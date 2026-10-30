# Public start-gain regression

`CDKN2A.json` pins the known CDKN2A c.-34G>T allele, GRCh38 chromosome 9
21974861 C>A, with source links and the SHA256 of NCBI's NM_000077.4 GenBank
record. Liu et al. (1999), doi:10.1038/5082, PMID 9916806, experimentally showed
creation of an upstream initiation codon and reduced canonical translation.
ClinVar RCV000493169.3 supplies the current coordinate/HGVS cross-reference.

The expected full peptide is an independent **sequence prediction**, not a claim
of observed patient translation: replace G with T at CDS-start minus 34 in the
frozen RefSeq RNA, begin at the preceding A (CDS-start minus 35), and translate
to the first stop with Biopython's standard code. The test repeats this oracle
and checks that the pinned Ensembl 95 ENST00000304494 coding path is identical.
The full RNAs differ in their 5' extent, so their start offsets are kept separate.
No liftover or Varcode-derived expectation is used.

The RefSeq nucleotide record is NCBI public sequence data. Source metadata and
checksum permit retrieval of the exact accession independently of this fixture.
