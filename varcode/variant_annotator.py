from __future__ import print_function, division, absolute_import

from .common import group_by
from .core_logic import infer_transcript_effect
from .variant import Variant
from .variant_effect import VariantEffect

import pyensembl

class VariantAnnotator(object):
    def __init__(self, ensembl_release):
        self.ensembl = pyensembl.EnsemblRelease(ensembl_release)

    def describe_variant(self, variant):
        """
        Determine the effects of a variant on any transcripts it overlaps.
        Returns a VariantEffect object.
        """
        contig = variant.contig
        pos = variant.pos
        end_pos = variant.end_pos
        ref = variant.ref
        alt = variant.alt

        overlapping_genes = self.ensembl.genes_at_locus(contig, pos, end_pos)

        if len(overlapping_genes) == 0:
            return [], {}
        else:
            overlapping_transcripts = self.ensembl.transcripts_at_locus(
                    contig, pos, end_pos)

            assert len(overlapping_transcripts) > 0, \
                "No transcripts found for mutation %s:%d %s>%s" % (
                    contig, pos, ref, alt)

        # group transcripts by their gene ID
        overlapping_transcript_groups = group_by(
            overlapping_transcripts, field_name='gene_id')

        # dictionary from gene ID to list of transcript effects
        gene_transcript_effects_groups = {}
        for (gene_id, transcripts) in overlapping_transcript_groups.items():
            effects = []
            for transcript in transcripts:
                effect = infer_transcript_effect(variant, transcript)
                effects.append(effect)
            gene_transcript_effects_groups[gene_id] = effects

        return VariantEffect(
            variant=variant,
            genes=overlapping_genes,
            gene_transcript_effects=gene_transcript_effects_groups)
