# Copyright (c) 2015. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function, division, absolute_import

from .effect import CodingSequenceMutation

def aggregate_gene_expression_levels(
        gene_fpkm_df,
        map_novel_genes_onto_ensembl=False,
        ensembl=None):
    """
    Create a dictionary mapping gene IDs to expression levels.

    Parameters
    ----------
    gene_fpkm_df : DataFrame
        DataFrame with columns:
            - 'id'
            - 'novel'
            - 'fpkm'
            - 'chr'
            - 'start'
            - 'end'
            - 'gene_names'
        IDs can be either from Ensembl or Cufflinks novel genes

    ensembl : pyensembl.EnsemblRelease
        Ensembl annotation database used to look up genes overlapping
        chromosomal locations.

    We can't just use the values from a Cufflinks gene.fpkm_tracking file and
    instead have to add up contributions from possibly multiple novel gene IDs
    all of which vouch for expression of annotated genes. For unclear reasons,
    Cufflinks sometimes assigns a novel gene ID (e.g. "CUFF.1") to regions of
    the genome which overlap an annotated gene.

    Example entry from a Cufflinks gene tracking file:

        CUFF.29167  -   -   CUFF.29167  BRAF    -   chr7:140424942-140624564

    This region overlaps BRAF entirely and thus the FPKM value of this
    novel gene should be assigned to BRAF's gene ID (ENSG00000157764). In the
    case that a novel gene entry overlaps multiple Ensembl genes, divide that
    entry's contribution by the number of overlapped genes.
    """

    # first drop any entries with FPKM == 0, since they don't
    # contribute to the sum of any gene

    gene_fpkm_df = gene_fpkm_df[gene_fpkm_df.fpkm > 0]

    # start with expression levels of known Ensembl genes
    mask = gene_fpkm_df.known
    known_df = gene_fpkm_df[mask]
    result_dict = dict(zip(known_df.id, known_df.fpkm))

    novel_df = gene_fpkm_df[~mask]
    for _, novel_row in novel_df.iterrows():
        # find genes overlapping the chromosomal positions spanned
        # by the novel gene constructed by Cufflinks
        overlapping_genes = ensembl.genes_at_locus(
            novel_row.chr,
            novel_row.start,
            novel_row.end)
        # some genes may be on the wrong strand so only use those
        # whose name was in the Cufflinks file
        matching_genes = [
            gene
            for gene in overlapping_genes
            if gene.name in novel_row.gene_names
        ]
        n_matching_genes = len(matching_genes)
        for gene in matching_genes:
            old_fpkm = result_dict.get(gene.id, 0.0)
            # split FPKM expression value across all matching genes
            # overlapping this locus, since we don't know which one to assign
            # expression and don't want to inflate loci that span lots of
            # genes
            result_dict[gene.id] = old_fpkm + novel_row.fpkm / n_matching_genes
    return result_dict

def expressed_gene_ids(
        gene_fpkm_df,
        fpkm_threshold,
        map_novel_genes_onto_ensembl=False,
        ensembl=None):
    """Return all gene IDs with expression level above the given threshold
    in the given Cufflinks file.
    """
    gene_fpkm_dict = aggregate_gene_expression_levels(
        gene_fpkm_df,
        map_novel_genes_onto_ensembl=map_novel_genes_onto_ensembl,
        ensembl=ensembl)
    return {
        gene_id
        for (gene_id, fpkm) in gene_fpkm_dict.iteritems()
        if fpkm >= fpkm_threshold
    }

def choose_principal_transcripts(
        variant_collection,
        gene_fpkm_df,
        gene_expression_threshold,
        transcript_fpkm_df,
        transcript_expression_threshold=0.0):

    # dictionary whose keys are Ensembl gene IDs and values are FPKM values
    gene_fpkm_dict = aggregate_gene_expression_levels(
        gene_fpkm_df, ensembl=variant_collection.ensembl)

    transcript_fpkm_dict = dict(
        zip(transcript_fpkm_df.id, transcript_fpkm_df.fpkm))

    principal_transcript_effects = []
    for variant_effect in variant_collection.variant_effects():

        # mapping from transcript ID to pair (gene fpkm, transcript fpkm)
        # we use this to first look at genes of high expression and then
        # choose their most highly expressed transcript
        transcript_expression_levels = {}
        transcript_effect_list = []
        for gene_id, transcript_effects in \
                variant_effect.gene_transcript_effects.iteritems():
            gene_fpkm = gene_fpkm_dict.get(gene_id)

            if gene_fpkm <= gene_expression_threshold:
                continue

            for transcript_effect in transcript_effects:
                if not isinstance(transcript_effect, CodingSequenceMutation):
                    continue

                transcript_id = transcript_effect.transcript.id
                transcript_fpkm = transcript_fpkm_dict.get(
                    transcript_id, 0.0)

                if transcript_fpkm <= transcript_expression_threshold:
                    continue

                fpkm_pair = (gene_fpkm, transcript_fpkm)
                transcript_expression_levels[transcript_id] = fpkm_pair
                transcript_effect_list.append(transcript_effect)

        def key(transcript_effect):
            return transcript_expression_levels[transcript_effect.transcript.id]

        transcript_effect_list.sort(key=key, reverse=True)

        if len(transcript_effect_list) > 0:
            best_transcript_effect = transcript_effect_list[0]
            principal_transcript_effects.append(best_transcript_effect)

    return principal_transcript_effects