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

def remap_novel_gene_expression_onto_ensembl_ids(gene_fpkm_df, ensembl):
    """Create a dictionary mapping annotated Ensembl gene IDs to expression
    levels. Rolls up novel detected genes into annotated genes which fully
    contain them.

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
