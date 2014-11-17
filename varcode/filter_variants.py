import argparse
from collections import namedtuple
import vcf
import pyensembl

parser = argparse.ArgumentParser()

parser.add_argument("input",
    help="VCF input file")

parser.add_argument("--output-vcf",
    help="VCF output file (keeping only header and coding variants)")

parser.add_argument("--output-csv",
    help="CSV output with chr/pos/ref/alt from VCF and annotation")


Gene = namedtuple("Gene",
    [
        'gene_name',
        'gene_id',
        'chrom',
        'start',
        'end',
        'strand',
        'biotype',
    ]
)

Transcript = namedtuple("Transcript",
    [
        'transcript_name',
        'transcript_id',
        'chrom',
        'start_position',
        'end_position',
        'start_codon_position',
        'stop_codon_position',
        'strand',
        'biotype',
        'exons',
    ]
)

Exon = namedtuple("Exon",
    [
        'exon_id',
        'start',
        'end',
        'strand',
        'frame',
    ]
)

AnnotatedVariant = namedtuple("AnnotatedVariant",
    [
        'chrom',
        'pos',
        'ref',
        'alt',
        'variant_type',
        # dictionary mapping each gene name to its Gene info struct
        'genes',
        # dictionary mapping each transcript name to its Transcript info struct
        'transcripts'
        # dictionary mapping each exon id to its Exon info struct
        'exons',
        # dict from transcript name to description of protein variant
        'protein_variants',
        'coding',
    ]
)

class VariantAnnotator(object):
    def __init__(self, ensembl_release):
        self.ensembl = pyensembl.EnsemblRelease(ensembl_release)


    def variant_genes(self, chrom, pos, number_modified_bases=1):
        """
        Parameters
        ----------

        chrom : str
            Chromosome or contig name

        pos : int
            Position in the chromosome

        number_modified_bases : int
            How many reference bases were changed or deleted?
        """
        gene_names = ensembl.gene_names_at_loci(
            chrom, pos, pos + number_modified_bases)
        genes = {}

        for gene_name in gene_names:
            self.ensembl.query(

            )
    def variant_transcripts(self, chrom, pos, number_modified_bases=1):
        """
        Parameters
        ----------

        chrom : str
            Chromosome or contig name

        pos : int
            Position in the chromosome

        number_modified_bases : int
            How many reference bases were changed or deleted?
        """
        return ensembl.transcript_names_at_loci(
            chrom, pos, pos + number_modified_bases)

    def describe_variant(self, chrom, pos, ref, alt):
        genes = self.variant_genes(chrom, pos, len(ref))

        if len(genes) == 0:
            return AnnotatedVariant(
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                variant_type='intergenic',
                genes={},
                transcripts={},
                exons={},
                protein_variants={},
                coding=False)

        transcripts = self.variant_transcripts(chrom, pos, len(ref))

        if len(transcripts) == 0:
            return AnnotatedVariant(
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                variant_type='intronic',
                genes=genes,
                transcripts={},
                exons={},
                protein_variants={},
                coding=False)


if __name__ == "__main__":
    # TODO: determine ensembl release from VCF metadata
    annot = VariantAnnotator(ensembl_release=75)

    args = parser.parse_args()
    with open(args.input, 'r') as f:
        vcf_reader = vcf.Reader(f)

        for record in vcf_reader:
            chrom, pos = record.CHROM, record.POS
            ref, alt = record.REF, record.ALT
            print chrom, pos, ref, alt
            print "--", annot.variant_gene_names(chrom, pos, ref, alt)
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