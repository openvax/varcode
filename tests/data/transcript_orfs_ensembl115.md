The adjacent JSON contains native Ensembl 115 GRCh38 transcript, CDS, exon, cDNA
and protein slices, with source URLs, sizes and SHA-256 hashes retained from the
ORFdata reference contract fixtures. RPL18 and KAZN additionally retain exact
GTF source lines, including partial-CDS tags and phases. Other GTF wrappers in
`tests/transcript_orf_fixture.py` are reconstructed from the native coordinates,
with exon numbers in RNA order. Protein sequences come independently from the
provider FASTA, not translation by the tested implementation.

The selected records cover POMC's complete CDS, RPL18's 5-prime incomplete CDS,
KAZN's 3-prime incomplete CDS, EIF4G2's GTG start, FGF2's alternative CTG starts,
SNURF/SNRPN and MIEF1/MIURF shared coding placements, and unsupported recoding.
FGF2 and EIF4G2 protein definitions also retain the pinned UniProt 2026_03 source
release, evidence text, source hashes and mappings. No patient evidence is added.

Ensembl data terms: https://www.ensembl.org/info/about/legal/disclaimer.html.
Partial locations and codon_start: https://www.insdc.org/submitting-standards/feature-table/.
