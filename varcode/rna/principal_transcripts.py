from .combine_genes import combine_gene_expression_levels

def choose_principal_transcripts(
        variant_collection,
        gene_fpkm_df,
        gene_expression_threshold,
        transcript_fpkm_df,
        transcript_expression_threshold=0.0):

    # dictionary whose keys are Ensembl gene IDs and values are FPKM values
    gene_fpkm_dict = combine_gene_expression_levels(
        gene_fpkm_df,
        ensembl=variant_collection.ensembl)

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
                transcript_id = transcript_effect.transcript.id
                if not transcript_id:
                    continue

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
