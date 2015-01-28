from variant import Variant
from variant_effect import VariantEffect
from infer_effect import infer_effects
from transcript_mutation_effects import top_priority_variant_effect

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

    def describe_variant(self, variant):
        overlapping_genes, transcript_effects_groups = infer_effects(
            self.ensembl, variant)

        # if our variant overlaps any genes, then choose the highest
        # priority transcript variant, otherwise call the variant "Intergenic"
        if len(overlapping_genes) == 0:
            variant_type = "Intergenic"
        else:
            all_variant_effects = []
            for _, variant_effects in transcript_effects_groups.iteritems():
                all_variant_effects.extend(variant_effects)
            summary_effect = top_priority_variant_effect(all_variant_effects)
            variant_type = summary_effect.__class__.__name__

        return VariantEffect(
            variant=variant,
            variant_type=variant_type,
            genes=overlapping_genes,
            gene_transcript_effects=transcript_effects_groups)
