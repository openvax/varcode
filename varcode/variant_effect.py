from transcript_effects import (
    Coding,
    Exonic,
    IncompleteTranscript,
)

from pyensembl.biotypes import is_coding_biotype


class VariantEffect(object):
    """
    A VariantEffect object is a container for all the TranscriptEffects of
    a particular mutation, as well as some properties that attempt to
    summarize the impact of the mutation across all genes/transcripts.
    """

    def __init__(
            self,
            variant,
            genes,
            transcripts,
            transcript_effects):
        """
        variant : Variant

        genes : list
            List of Gene objects

        transcripts : dict
            Dictionary mapping transcript IDs to Transcript objects

        transcript_effects : dict
            Dictionary from transcript ID to description of protein variant
        """
        self.variant = variant
        self.variant_type = variant_type
        self.genes = genes
        self.transcripts = transcripts
        self.transcript_effects = transcript_effects

    @property
    def coding_genes(self):
        """
        The `genes` property includes all sorts of esoteric entities
        like anti-sense annotations and pseudogenes.

        This property gives you just protein coding genes and functional RNAs
        """
        return [
            gene for gene in self.genes
            if is_coding_biotype(gene.biotype)
        ]

    @property
    def variant_type(self):
         """
         Returns one of the following strings:
             - intergenic: not mapped to any gene in Ensembl
             - intronic: on a gene but not on or near an exon
             - exonic: on a gene inside an exon
        """
        if len(self.genes) == 0:
            return "intergenic"

        for effect in self.transcript_effects.values():
            if isinstance(effect, Exonic):
                return "exonic"
        return "intronic"

    def __str__(self):
        fields = [
            ("variant", self.variant),
            ("n_genes", len(self.genes)),
            ("n_coding_genes", len(self.coding_genes)),
            ("coding_genes", self.coding_genes),
            ("transcript_effects", self.transcript_effects)
        ]
        return "VariantEffect(%s)" % (
            ", ".join(["%s=%s" % (k,v) for (k,v) in fields]))
