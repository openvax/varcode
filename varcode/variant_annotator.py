from __future__ import print_function, division, absolute_import
import logging

from .common import group_by
from .core_logic import infer_transcript_effect
from .variant import Variant
from .variant_effect import VariantEffect

import pyensembl

class VariantAnnotator(object):
    def __init__(self, ensembl_release):
        self.ensembl = pyensembl.EnsemblRelease(ensembl_release)

    def describe_variant(self, variant, raise_on_error=True):
        """
        Determine the effects of a variant on any transcripts it overlaps.
        Returns a VariantEffect object.
        """
        contig = variant.contig
        pos = variant.pos
        end_pos = variant.end_pos
        ref = variant.ref
        alt = variant.alt

        overlapping_transcripts = self.ensembl.transcripts_at_locus(
            contig, pos, end_pos)

        if len(overlapping_transcripts) == 0:
            # intergenic variant
            return VariantEffect(
                variant=variant,
                genes=[],
                gene_transcript_effects={},
                errors={})

        # group transcripts by their gene ID
        overlapping_transcript_groups = group_by(
            overlapping_transcripts, field_name='gene_id')

        # dictionary from gene ID to list of transcript effects
        gene_transcript_effects_groups = {}

        # mapping from Transcript objects to errors encountered
        # while trying to annotated them
        errors = {}

        for (gene_id, transcripts) in overlapping_transcript_groups.items():
            effects = []
            for transcript in transcripts:
                try:
                    effects.append(infer_transcript_effect(variant, transcript))
                except (AssertionError, ValueError) as error:
                    if raise_on_error:
                        raise
                    else:
                        logging.warn(
                            "Encountered error annotating %s for %s: %s",
                            variant,
                            transcript,
                            error)
                    errors[transcript] = error
            gene_transcript_effects_groups[gene_id] = effects

        # construct Gene objects for all the genes containing
        # all the transcripts this variant overlaps
        overlapping_genes = [
            self.ensembl.gene_by_id(gene_id)
            for gene_id
            in overlapping_transcript_groups.keys()
        ]

        return VariantEffect(
            variant=variant,
            genes=overlapping_genes,
            gene_transcript_effects=gene_transcript_effects_groups,
            errors=errors)
