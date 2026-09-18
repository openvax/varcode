# API reference

Start with [Getting started](getting_started.md) for examples, or choose a topic
below for signatures, parameters, and return types. Most annotation uses the
same `effects()` interface; [annotator selection](experimental_annotators.md)
is optional.

## Variants and files

<a id="reference-genomes"></a><a id="varcodegenome"></a><a id="variants"></a><a id="varcodevariant"></a><a id="varcodevariantcollection"></a><a id="varcodestructuralvariant"></a><a id="varcodesv_types"></a><a id="genotypes"></a><a id="varcodegenotype"></a><a id="varcodezygosity"></a><a id="variantcollection-transforms"></a><a id="file-loading"></a><a id="varcodeload_vcf"></a><a id="exceptions"></a><a id="varcodereferencemismatcherror"></a><a id="varcodesamplenotfounderror"></a><a id="varcodegenomebuildmismatcherror"></a>[Variants and files reference](api_variants.md)

- <a id="varcode.Genome"></a><a id="varcode.Genome.__getattr__"></a><a id="varcode.Genome.__dir__"></a><a id="varcode.Genome.sequence"></a><a id="varcode.Genome.reference_base"></a><a id="varcode.Genome.reference_range"></a>[Genome](api_variants.md#varcode.Genome)
- <a id="varcode.Variant"></a><a id="varcode.Variant.ensembl"></a><a id="varcode.Variant.trimmed_ref"></a><a id="varcode.Variant.trimmed_alt"></a><a id="varcode.Variant.trimmed_base1_start"></a><a id="varcode.Variant.trimmed_base1_end"></a><a id="varcode.Variant.short_description"></a><a id="varcode.Variant.coding_transcripts"></a><a id="varcode.Variant.genes"></a><a id="varcode.Variant.gene_ids"></a><a id="varcode.Variant.gene_names"></a><a id="varcode.Variant.coding_genes"></a><a id="varcode.Variant.is_insertion"></a><a id="varcode.Variant.is_deletion"></a><a id="varcode.Variant.is_indel"></a><a id="varcode.Variant.is_snv"></a><a id="varcode.Variant.is_transition"></a><a id="varcode.Variant.is_transversion"></a><a id="varcode.Variant.__lt__"></a><a id="varcode.Variant.to_dict"></a><a id="varcode.Variant.effects"></a><a id="varcode.Variant.effect_on_transcript"></a><a id="varcode.Variant.clone_without_ucsc_data"></a>[Variant](api_variants.md#varcode.Variant)
- <a id="varcode.VariantCollection"></a><a id="varcode.VariantCollection.metadata"></a><a id="varcode.VariantCollection.samples"></a><a id="varcode.VariantCollection.to_dict"></a><a id="varcode.VariantCollection.clone_with_new_elements"></a><a id="varcode.VariantCollection.effects"></a><a id="varcode.VariantCollection.reference_names"></a><a id="varcode.VariantCollection.original_reference_names"></a><a id="varcode.VariantCollection.groupby_gene_name"></a><a id="varcode.VariantCollection.gene_counts"></a><a id="varcode.VariantCollection.filter_by_transcript_expression"></a><a id="varcode.VariantCollection.filter_by_gene_expression"></a><a id="varcode.VariantCollection.exactly_equal"></a><a id="varcode.VariantCollection.union"></a><a id="varcode.VariantCollection.intersection"></a><a id="varcode.VariantCollection.difference"></a><a id="varcode.VariantCollection.to_dataframe"></a><a id="varcode.VariantCollection.to_csv"></a><a id="varcode.VariantCollection.from_csv"></a><a id="varcode.VariantCollection.has_sample_data"></a><a id="varcode.VariantCollection.genotype"></a><a id="varcode.VariantCollection.zygosity"></a><a id="varcode.VariantCollection.for_sample"></a><a id="varcode.VariantCollection.heterozygous_in"></a><a id="varcode.VariantCollection.homozygous_alt_in"></a>[VariantCollection](api_variants.md#varcode.VariantCollection)
- <a id="varcode.StructuralVariant"></a><a id="varcode.StructuralVariant.symbolic_alt"></a><a id="varcode.StructuralVariant.junctions"></a><a id="varcode.StructuralVariant.breakpoints"></a><a id="varcode.StructuralVariant.length"></a><a id="varcode.StructuralVariant.to_dict"></a>[StructuralVariant](api_variants.md#varcode.StructuralVariant)
- <a id="varcodeparse_symbolic_alt"></a><a id="varcode.parse_symbolic_alt"></a>[parse_symbolic_alt](api_variants.md#varcode.parse_symbolic_alt)
- <a id="varcode.SV_TYPES"></a>[SV_TYPES](api_variants.md#varcode.SV_TYPES)
- <a id="varcode.Genotype"></a><a id="varcode.Genotype.is_called"></a><a id="varcode.Genotype.ploidy"></a><a id="varcode.Genotype.from_sample_info"></a><a id="varcode.Genotype.carries_alt"></a><a id="varcode.Genotype.copies_of_alt"></a><a id="varcode.Genotype.zygosity_for_alt"></a><a id="varcode.Genotype.depth_for_alt"></a>[Genotype](api_variants.md#varcode.Genotype)
- <a id="varcode.Zygosity"></a>[Zygosity](api_variants.md#varcode.Zygosity)
- <a id="varcodetransformspair_breakends"></a><a id="varcode.transforms.pair_breakends"></a>[transforms.pair_breakends](api_variants.md#varcode.transforms.pair_breakends)
- <a id="varcodetransformsleft_align_indels"></a><a id="varcode.transforms.left_align_indels"></a>[transforms.left_align_indels](api_variants.md#varcode.transforms.left_align_indels)
- <a id="varcode.vcf.load_vcf"></a>[load_vcf](api_variants.md#varcode.vcf.load_vcf)
- <a id="varcodeload_maf"></a><a id="varcode.load_maf"></a>[load_maf](api_variants.md#varcode.load_maf)
- <a id="varcode.load_maf_dataframe"></a>[load_maf_dataframe](api_variants.md#varcode.load_maf_dataframe)
- <a id="varcode.ReferenceMismatchError"></a>[ReferenceMismatchError](api_variants.md#varcode.ReferenceMismatchError)
- <a id="varcode.SampleNotFoundError"></a>[SampleNotFoundError](api_variants.md#varcode.SampleNotFoundError)
- <a id="varcode.GenomeBuildMismatchError"></a>[GenomeBuildMismatchError](api_variants.md#varcode.GenomeBuildMismatchError)

## Effects

<a id="varcodemutationeffect"></a><a id="varcodenonsilentcodingmutation"></a><a id="varcodemultioutcomeeffect"></a><a id="varcodeeffectcollection"></a><a id="varcodeeffectcandidate"></a><a id="priority-ordering"></a>[Effects reference](api_effects.md)

- <a id="varcode.MutationEffect"></a><a id="varcode.MutationEffect.short_description"></a><a id="varcode.MutationEffect.original_protein_sequence"></a><a id="varcode.MutationEffect.__lt__"></a>[MutationEffect](api_effects.md#varcode.MutationEffect)
- <a id="varcode.NonsilentCodingMutation"></a>[NonsilentCodingMutation](api_effects.md#varcode.NonsilentCodingMutation)
- <a id="varcode.MultiOutcomeEffect"></a><a id="varcode.MultiOutcomeEffect.effects"></a><a id="varcode.MultiOutcomeEffect.most_likely_candidate"></a><a id="varcode.MultiOutcomeEffect.most_likely_effect"></a><a id="varcode.MultiOutcomeEffect.highest_priority_candidate"></a><a id="varcode.MultiOutcomeEffect.highest_priority_effect"></a>[MultiOutcomeEffect](api_effects.md#varcode.MultiOutcomeEffect)
- <a id="varcode.EffectCollection"></a><a id="varcode.EffectCollection.gene_counts"></a><a id="varcode.EffectCollection.filter_by_transcript_expression"></a><a id="varcode.EffectCollection.filter_by_gene_expression"></a><a id="varcode.EffectCollection.filter_by_effect_priority"></a><a id="varcode.EffectCollection.drop_silent_and_noncoding"></a><a id="varcode.EffectCollection.detailed_string"></a><a id="varcode.EffectCollection.top_priority_effect"></a><a id="varcode.EffectCollection.top_priority_effect_per_variant"></a><a id="varcode.EffectCollection.top_priority_effect_per_transcript_id"></a><a id="varcode.EffectCollection.top_priority_effect_per_gene_id"></a><a id="varcode.EffectCollection.effect_expression"></a><a id="varcode.EffectCollection.top_expression_effect"></a><a id="varcode.EffectCollection.to_dataframe"></a><a id="varcode.EffectCollection.to_csv"></a><a id="varcode.EffectCollection.from_csv"></a>[EffectCollection](api_effects.md#varcode.EffectCollection)
- <a id="varcode.EffectCandidate"></a><a id="varcode.EffectCandidate.short_description"></a>[EffectCandidate](api_effects.md#varcode.EffectCandidate)
- <a id="varcode.effect_priority"></a>[effect_priority](api_effects.md#varcode.effect_priority)
- <a id="varcode.top_priority_effect"></a>[top_priority_effect](api_effects.md#varcode.top_priority_effect)

## Phasing and germline

<a id="phasing"></a><a id="varcodereadphasingsource"></a><a id="varcodemutanttranscriptsource"></a><a id="varcodemolecularphaseresolver"></a><a id="varcodereadphaseresolver"></a><a id="varcodernareadphasingsource"></a><a id="varcodevcfphaseresolver"></a><a id="germline-aware-annotation"></a><a id="varcodegermlinecontext"></a><a id="varcodecompleteness"></a>[Phasing and germline reference](api_phasing.md)

- <a id="varcode.ReadPhasingSource"></a><a id="varcode.ReadPhasingSource.has_evidence"></a><a id="varcode.ReadPhasingSource.partners_in_cis"></a>[ReadPhasingSource](api_phasing.md#varcode.ReadPhasingSource)
- <a id="varcode.MutantTranscriptSource"></a><a id="varcode.MutantTranscriptSource.mutant_transcript"></a>[MutantTranscriptSource](api_phasing.md#varcode.MutantTranscriptSource)
- <a id="varcode.MolecularPhaseResolver"></a><a id="varcode.MolecularPhaseResolver.has_evidence"></a><a id="varcode.MolecularPhaseResolver.mutant_transcript"></a><a id="varcode.MolecularPhaseResolver.in_cis"></a><a id="varcode.MolecularPhaseResolver.phased_partners"></a>[MolecularPhaseResolver](api_phasing.md#varcode.MolecularPhaseResolver)
- <a id="varcode.ReadPhaseResolver"></a>[ReadPhaseResolver](api_phasing.md#varcode.ReadPhaseResolver)
- <a id="varcode.RNAReadPhasingSource"></a><a id="varcode.RNAReadPhasingSource.close"></a><a id="varcode.RNAReadPhasingSource.register_variants"></a><a id="varcode.RNAReadPhasingSource.register_haplotype"></a><a id="varcode.RNAReadPhasingSource.supports_variant"></a><a id="varcode.RNAReadPhasingSource.has_evidence"></a><a id="varcode.RNAReadPhasingSource.in_cis"></a><a id="varcode.RNAReadPhasingSource.partners_in_cis"></a>[RNAReadPhasingSource](api_phasing.md#varcode.RNAReadPhasingSource)
- <a id="varcode.VCFPhaseResolver"></a><a id="varcode.VCFPhaseResolver.in_cis"></a><a id="varcode.VCFPhaseResolver.phased_partners"></a>[VCFPhaseResolver](api_phasing.md#varcode.VCFPhaseResolver)
- <a id="varcodeapply_phase_resolver_to_effects"></a><a id="varcode.apply_phase_resolver_to_effects"></a>[apply_phase_resolver_to_effects](api_phasing.md#varcode.apply_phase_resolver_to_effects)
- <a id="varcode.GermlineContext"></a><a id="varcode.GermlineContext.from_germline_vcf"></a><a id="varcode.GermlineContext.from_multi_sample_vcf"></a><a id="varcode.GermlineContext.from_variants"></a><a id="varcode.GermlineContext.empty"></a><a id="varcode.GermlineContext.__bool__"></a><a id="varcode.GermlineContext.validate_against"></a><a id="varcode.GermlineContext.variants_in_window"></a>[GermlineContext](api_phasing.md#varcode.GermlineContext)
- <a id="varcode.Completeness"></a>[Completeness](api_phasing.md#varcode.Completeness)
- <a id="varcodepredict_germline_aware_effect"></a><a id="varcode.predict_germline_aware_effect"></a>[predict_germline_aware_effect](api_phasing.md#varcode.predict_germline_aware_effect)
- <a id="varcodeapply_germline_to_transcript"></a><a id="varcode.apply_germline_to_transcript"></a>[apply_germline_to_transcript](api_phasing.md#varcode.apply_germline_to_transcript)
- <a id="varcodeenumerate_phase_hypotheses"></a><a id="varcode.enumerate_phase_hypotheses"></a>[enumerate_phase_hypotheses](api_phasing.md#varcode.enumerate_phase_hypotheses)
- <a id="varcodedetect_loh"></a><a id="varcode.detect_loh"></a>[detect_loh](api_phasing.md#varcode.detect_loh)
- <a id="varcodedefault_germline_window"></a><a id="varcode.default_germline_window"></a>[default_germline_window](api_phasing.md#varcode.default_germline_window)

## RNA and transcripts

<a id="mutant-transcripts"></a><a id="observed-rna-import"></a><a id="varcodemutanttranscript"></a><a id="rna-evidence"></a><a id="varcodernaevidenceresolver"></a><a id="varcodenullrnaevidenceresolver"></a>[RNA and transcripts reference](api_rna.md)

- <a id="varcode.load_exacto_fusions"></a>[load_exacto_fusions](api_rna.md#varcode.load_exacto_fusions)
- <a id="varcode.make_fusion_outcome"></a>[make_fusion_outcome](api_rna.md#varcode.make_fusion_outcome)
- <a id="varcode.RNAEvidence"></a>[RNAEvidence](api_rna.md#varcode.RNAEvidence)
- <a id="varcode.MutantTranscript"></a><a id="varcode.MutantTranscript.reference_transcript"></a><a id="varcode.MutantTranscript.edits"></a><a id="varcode.MutantTranscript.reference_segments"></a><a id="varcode.MutantTranscript.cdna_sequence"></a><a id="varcode.MutantTranscript.mutant_protein_sequence"></a><a id="varcode.MutantTranscript.annotator_name"></a><a id="varcode.MutantTranscript.evidence"></a><a id="varcode.MutantTranscript.is_identical_to_reference"></a><a id="varcode.MutantTranscript.is_structural"></a><a id="varcode.MutantTranscript.total_length_delta"></a><a id="varcode.MutantTranscript.from_sequence"></a>[MutantTranscript](api_rna.md#varcode.MutantTranscript)
- <a id="varcodeapply_variant_to_transcript"></a><a id="varcode.apply_variant_to_transcript"></a>[apply_variant_to_transcript](api_rna.md#varcode.apply_variant_to_transcript)
- <a id="varcodeapply_variants_to_transcript"></a><a id="varcode.apply_variants_to_transcript"></a>[apply_variants_to_transcript](api_rna.md#varcode.apply_variants_to_transcript)
- <a id="varcode.RNAEvidenceResolver"></a><a id="varcode.RNAEvidenceResolver.observed_outcomes"></a>[RNAEvidenceResolver](api_rna.md#varcode.RNAEvidenceResolver)
- <a id="varcode.NullRNAEvidenceResolver"></a>[NullRNAEvidenceResolver](api_rna.md#varcode.NullRNAEvidenceResolver)
- <a id="varcodeapply_rna_evidence_to_effects"></a><a id="varcode.apply_rna_evidence_to_effects"></a>[apply_rna_evidence_to_effects](api_rna.md#varcode.apply_rna_evidence_to_effects)
- <a id="varcodemake_rna_outcome"></a><a id="varcode.make_rna_outcome"></a>[make_rna_outcome](api_rna.md#varcode.make_rna_outcome)

## Annotators

<a id="varcodeeffectannotator"></a><a id="varcodefasteffectannotator"></a><a id="varcodeproteindiffeffectannotator"></a><a id="registry"></a><a id="experimental-transcript-model"></a><a id="varcodetranscriptmodeleffectannotator"></a>[Annotators reference](api_annotators.md)

- <a id="varcode.EffectAnnotator"></a>[EffectAnnotator](api_annotators.md#varcode.EffectAnnotator)
- <a id="varcode.FastEffectAnnotator"></a><a id="varcode.FastEffectAnnotator.version"></a><a id="varcode.FastEffectAnnotator.annotate_on_transcript"></a><a id="varcode.FastEffectAnnotator.annotate_with_context"></a>[FastEffectAnnotator](api_annotators.md#varcode.FastEffectAnnotator)
- <a id="varcode.ProteinDiffEffectAnnotator"></a><a id="varcode.ProteinDiffEffectAnnotator.annotate_with_context"></a><a id="varcode.ProteinDiffEffectAnnotator.annotate_on_transcript"></a>[ProteinDiffEffectAnnotator](api_annotators.md#varcode.ProteinDiffEffectAnnotator)
- <a id="varcode.register_annotator"></a>[register_annotator](api_annotators.md#varcode.register_annotator)
- <a id="varcode.get_annotator"></a>[get_annotator](api_annotators.md#varcode.get_annotator)
- <a id="varcode.get_default_annotator"></a>[get_default_annotator](api_annotators.md#varcode.get_default_annotator)
- <a id="varcode.set_default_annotator"></a>[set_default_annotator](api_annotators.md#varcode.set_default_annotator)
- <a id="varcode.use_annotator"></a>[use_annotator](api_annotators.md#varcode.use_annotator)
- <a id="varcode.TranscriptModelEffectAnnotator"></a><a id="varcode.TranscriptModelEffectAnnotator.annotate_with_context"></a>[TranscriptModelEffectAnnotator](api_annotators.md#varcode.TranscriptModelEffectAnnotator)
- <a id="varcodepredict_transcript_model_effect"></a><a id="varcode.predict_transcript_model_effect"></a>[predict_transcript_model_effect](api_annotators.md#varcode.predict_transcript_model_effect)
