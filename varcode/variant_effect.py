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
            variant_type,
            genes,
            gene_transcript_effects):
        """
        variant : Variant

        variant_type : str
            Summary of variant effect across all transcripts

        genes : list
            List of Gene objects

        gene_transcript_effects : dict
            Dictionary from gene ID to list of transcript variant effects
        """
        self.variant = variant
        # add "highest priority" this name
        self.variant_type = variant_type
        self.genes = genes
        self.gene_transcript_effects = gene_transcript_effects

        # dictionary mapping from transcript IDs to transcript mutation effects
        self.transcript_effects = {}
        for (_, transcript_effects) in self.gene_transcript_effects.iteritems():
            for effect in transcript_effects:
                self.transcript_effects[effect.transcript.id] = effect

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
            ("genes", [gene.name for gene in self.genes]),
            ("n_coding_genes", len(self.coding_genes)),
            ("transcript_effects", self.transcript_effects)
        ]
        return "VariantEffect(%s)" % (
            ", ".join(["%s=%s" % (k,v) for (k,v) in fields]))

    def __repr__(self):
        return str(self)
