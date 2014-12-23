from pyensembl.biotypes import is_coding_biotype

class Annot(object):
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

    def __str__(self):
        fields = [
            ("variant", self.variant),
            ("variant_type", self.variant_type),
            ("n_genes", len(self.genes)),
            ("n_coding_genes", len(self.coding_genes)),
            ("coding_genes", self.coding_genes),
            ("coding_effects", self.coding_effects)
        ]
        return "Annot(%s)" % ", ".join(["%s=%s" % (k,v) for (k,v) in fields])
