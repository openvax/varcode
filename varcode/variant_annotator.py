from common import reverse_complement
from variant import Variant
from variant_effect import VariantEffect
from apply_variant import apply_variant


import pyensembl
from pyensembl.biotypes import is_coding_biotype


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


    def describe_variant(self, contig, pos, ref, alt):
        variant = Variant(contig=contig, pos=pos, ref=ref, alt=alt)
        overlapping_genes, transcript_effects_groups = apply_variant(
            self.ensembl, variant)

        return VariantEffect(
            variant=variant,
            genes=overlapping_genes,
            transcript_effects=transcript_effects_groups)
