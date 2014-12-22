import argparse
from collections import namedtuple

import maf

import numpy as np
import pyensembl
from pyensembl.locus import normalize_chromosome
import vcf

parser = argparse.ArgumentParser()

parser.add_argument("input",
    help="VCF input file")

parser.add_argument("--output-csv",
    help="Comma-separated output with chr/pos/ref/alt from VCF and annotation")

parser.add_argument("--output-tsv",
    help="Tab-separated output with chr/pos/ref/alt from VCF and annotation")



# include all pseudonucleotides encoding repeats and uncertain bases
VALID_NUCLEOTIDES = {'A', 'C', 'T', 'G'}

EXTENDED_NUCLEOTIDES = {
    'A', 'C', 'T', 'G',
    'Y', # Pyrimidine (C or T)
    'R', # Purine (A or G)
    'W', # weak (A or T)
    'S', # strong (G or C)
    'K', # keto (T or G)
    'M', # amino (C or A)
    'D', # A, G, T (not C)
    'V', # A, C, G (not T)
    'H', # A, C, T (not G)
    'B', # C, G, T (not A)
    'X', # any base
    'N', # any base
}

def normalize_nucleotide_string(nucleotides, allow_extended_nucleotides=False):
    """
    Normalizes a nucleotide string by converting various ways of encoding empty
    strings into "", making all letters upper case, and checking to make sure
    all letters in the string are actually nucleotides.

    Parameters
    ----------
    nucleotides : str
        Sequence of nucleotides, e.g. "ACCTG"

    extended_nucleotides : bool
        Allow non-canonical nucleotide characters like 'X' for unknown base
    """
    # some MAF files represent deletions/insertions with NaN ref/alt values
    if isinstance(nucleotides, float) and np.isnan(nucleotides):
        return ""
    # VCFs sometimes have '.' ref/alt for insertions and deletions
    elif nucleotides == '.':
        return ""

    if not isinstance(nucleotides, (str, unicode)):
        raise TypeError(
                "Expected nucleotide string, got %s : %s" % (
                    nucleotides, type(nucleotides)))

    nucleotides = nucleotides.upper()

    for letter in set(nucleotides):
        if allow_extended_nucleotides:
            valid_nucleotides = EXTENDED_NUCLEOTIDES
        else:
            valid_nucleotides = VALID_NUCLEOTIDES
        if letter not in valid_nucleotides:
            raise ValueError(
                "Invalid character in nucleotide string: %s" % letter)

    return nucleotides

class Variant(object):
    def __init__(self, contig, pos, ref, alt):
        self.contig = normalize_chromosome(contig)
        self.pos = int(pos)
        self.ref = normalize_nucleotide_string(ref)
        self.alt = normalize_nucleotide_string(alt)

    @property
    def end_pos(self):
        return self.pos + len(self.ref) + 1

    def __str__(self):
        return "Variant(contig=%s, pos=%d, ref=%s, alt=%s)" % (
            self.contig, self.pos, self.ref, self.alt)

CodingEffect = namedtuple(
    "CodingEffect",
    [
        "OriginalProteinSequence",

        "MutantProteinSequence",

        # unless a frameshift results in a stop codon within the same
        # exon as the mutation, we leave the sequence incomplete
        # which is marked here
        "MutantSequenceComplete",

        "CodingVariantDescription",
    ]
)

class AnnotatedVariant(object):
    def __init__(
            self,
            variant,
            variant_type,
            genes,
            transcripts,
            coding_effects):
        """
        variant : Variant

        variant_type : str
            One of the following:
             - intergenic: not mapped to any gene in Ensembl
             - intronic: on a gene but not on or near an exon
             - exonic: on a gene inside an exon

        genes : list
            List of Gene objects

        transcripts : dict
            Dictionary mapping gene IDs to list of Transcript objects

        coding_effects : dict
            Dictionary from transcript ID to description of protein variant
        """
        self.variant = variant
        self.variant_type = variant_type
        self.genes = genes
        self.transcripts = transcripts
        self.coding_effects = coding_effects


class VariantAnnotator(object):
    def __init__(self, ensembl_release):
        self.ensembl = pyensembl.EnsemblRelease(ensembl_release)


    def variant_gene_ids(self, contig, pos, number_modified_bases=1):
        """
        Parameters
        ----------

        contig : str
            Chromosome or contig name

        pos : int
            Position in the chromosome

        number_modified_bases : int
            How many reference bases were changed or deleted?
        """
        return self.ensembl.gene_ids_at_locus(
            contig, pos, pos + number_modified_bases)


    def variant_transcript_ids(self, contig, pos, number_modified_bases=1):
        """
        Parameters
        ----------

        contig : str
            Chromosome or contig name

        pos : int
            Position in the chromosome

        number_modified_bases : int
            How many reference bases were changed or deleted?
        """
        return self.ensembl.transcript_ids_at_locus(
            contig, pos, pos + number_modified_bases)

    def make_intergenic(self, variant):
        return AnnotatedVariant(
            variant=variant,
            variant_type='intergenic',
            genes=[],
            transcripts={},
            coding_effects={})


    def make_intronic(self, variant, genes, transcripts):
        return AnnotatedVariant(
            variant=variant,
            variant_type='intronic',
            genes=genes,
            transcripts=transcripts,
            coding_effects={})

    def overlaps_any_exon(self, transcript, contig, start, end):
        return any(
            exon.overlaps(contig=contig, start=start, end=end)
            for exon in transcript.exons)

    def group_by(self, records, field_name):
        groups = {}
        for record in records:
            value = getattr(record, field_name)
            if value in groups:
                groups[value].append(record)
            else:
                groups[value] = [record]
        return groups

    def describe_variant(self, contig, pos, ref, alt):
        variant = Variant(contig=contig, pos=pos, ref=ref, alt=alt)

        overlapping_genes = self.ensembl.genes_at_locus(
            contig, pos, variant.end_pos)

        if len(overlapping_genes) == 0:
            return self.make_intergenic(variant)

        overlapping_transcripts = self.ensembl.transcripts_at_locus(
                contig, pos, variant.end_pos)

        assert len(overlapping_transcripts) > 0, \
            "No transcripts found for mutation %s:%d %s>%s" % (
                contig, pos, ref, alt)

        exonic = any(
            self.overlaps_any_exon(
                transcript, contig, start=pos, end=variant.end_pos)
            for transcript in
            overlapping_transcripts
        )

        # group transcripts by their gene ID
        overlapping_transcript_groups = self.group_by(
            overlapping_transcripts, field_name='gene_id')

        if not exonic:
            return self.make_intronic(
                contig, pos, ref, alt,
                genes=overlapping_genes,
                transcripts=overlapping_transcript_groups)

        protein_variants = {}
        for transcript in overlapping_transcripts:
            if not transcript.complete:
                protein_variants[transcript.id] = "incomplete"
            else:
                seq = transcript.coding_sequence
                original_aa = "V"
                aa_position = 600
                new_aa = "E"
                variant_string = "%s%d%s" % (original_aa, aa_position, new_aa)
                protein_variants[transcript.id] = variant_string

        if len(protein_variants) > 0:
            variant_type = "coding"
        else:
            variant_type = "coding-without-complete-transcripts"

        return AnnotatedVariant(
            variant=variant,
            variant_type=variant_type,
            genes=overlapping_genes,
            transcripts=overlapping_transcripts,
            coding_effects=protein_variants)


def infer_reference_name(path):
    # NCBI builds and hg releases aren't identical
    # but the differences are all on chrM and unplaced contigs
    candidates = {
        'GRCh36' : ['hg18', 'b36', 'B36', 'GRCh36'],
        'GRCh37' : ['hg19', 'b37', 'B37', 'GRCh37'],
        'GRCh38' : ['b38', 'B38', 'GRCh38'],
    }

    for name in sorted(candidates.keys(), reverse=True):
        aliases = candidates[name]
        for alias in aliases:
            if alias in path:
                return name

    assert False, "Failed to infer human genome assembly name for %s" % path



class VariantCollection(object):

    def __init__(self, filename):
        assert isinstance(filename, (str,unicode)), \
            "Expected filename to be str, got %s : %s" % (
                filename, type(filename)
            )

        if filename.endswith(".vcf"):
            vcf_reader = vcf.Reader(filename=filename)
            print vcf_reader.metadata
            self.reference = infer_reference_name(
                vcf_reader.metadata['reference'])
            print vcf_reader
            self.records = [
                Variant(x.CHROM, x.POS, x.REF, x.ALT[0].sequence)
                for x in vcf_reader
            ]
        elif filename.endswith(".maf"):
            maf_df = maf.load_maf_dataframe(filename)
            assert len(maf_df) > 0, "Empty MAF file %s" % filename
            ncbi_builds = maf_df.NCBI_Build.unique()
            assert len(ncbi_builds) == 1, \
                "Multiple NCBI builds: %s" % ncbi_builds
            ncbi_build = ncbi_builds[0]
            self.reference = infer_reference_name(ncbi_build)
            self.records = []
            for _, x in maf_df.iterrows():
                start_pos = x.Start_Position
                end_pos = x.End_Position
                contig = x.Chromosome
                ref = normalize_nucleotide_string(x.Reference_Allele)

                if x.Tumor_Seq_Allele1 != ref:
                    alt = x.Tumor_Seq_Allele1
                else:
                    assert x.Tumor_Seq_Allele2 != ref, \
                        "Both tumor alleles agree with reference: %s" % (x,)
                    alt = x.Tumor_Seq_Allele2

                alt = normalize_nucleotide_string(alt)

                record = Variant(contig, start_pos, ref, alt)
                self.records.append(record)

        else:
            assert False, "Unrecognized file type: %s" % filename

    def __iter__(self):
        return iter(self.records)

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