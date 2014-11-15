import argparse

import vcf
import pyensembl

parser = argparse.ArgumentParser()

parser.add_argument("input",
    help="VCF input file")

parser.add_argument("--output-vcf",
    help="VCF output file (keeping only header and coding variants)")

parser.add_argument("--output-csv",
    help="CSV output with chr/pos/ref/alt from VCF and annotation")


if __name__ == "__main__":
    # TODO: determine ensembl release from VCF metadata
    ensembl = pyensembl.EnsemblRelease(75)

    args = parser.parse_args()
    with open(args.input, 'r') as f:
        vcf_reader = vcf.Reader(f)
        for record in vcf_reader:
            chrom, pos = record.CHROM, record.POS
            ref, alt = record.REF, record.ALT
            gene_names =  ensembl.gene_names_at_locus(chrom, pos)
            print chrom, "@", pos, " :: ", ref, "->", alt
            for gene in gene_names:
                print "  ", gene
                # TODO: implement transcript_names_of_gene_name in pyensembl
                # and transcript_name_of_transcript_id
                transcript_ids = ensembl.transcript_ids_of_gene_name(gene)
                for transcript_id in transcript_ids:
                    print "    ", transcript_id
