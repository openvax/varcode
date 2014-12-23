import argparse
from collections import namedtuple

import maf
from variant_annotator import VariantAnnotator
from variant_collection import VariantCollection

import numpy as np
import pyensembl


parser = argparse.ArgumentParser()

parser.add_argument("input",
    help="VCF input file")

parser.add_argument("--output-csv",
    help="Comma-separated output with chr/pos/ref/alt from VCF and annotation")

parser.add_argument("--output-tsv",
    help="Tab-separated output with chr/pos/ref/alt from VCF and annotation")



if __name__ == "__main__":
    # TODO: determine ensembl release from VCF metadata
    annot = VariantAnnotator(ensembl_release=75)

    args = parser.parse_args()

    for record in VariantCollection(args.input):
        chrom, pos = record.contig, record.pos
        ref, alt = record.ref, record.alt
        print chrom, pos, ref, alt
        print "--", annot.describe_variant(chrom, pos, ref, alt)
        """
        genes_to_transcripts = {}
        for transcript_name in transcript_names:
            gene = ensembl.gene_name_of_transcript_name(transcript_name)
            if gene in genes_to_transcripts:
                genes_to_transcripts[gene].append(transcript_name)
            else:
                genes_to_transcripts[gene] = [transcript_name]
        for (gene, transcript_names) in genes_to_transcripts.iteritems():
            print "\t", gene
            for transcript_name in transcript_names:
                print "\t\t", transcript_name
        """