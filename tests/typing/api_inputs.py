"""Static consumer checks; checked by mypy, never executed as tests."""

from collections import OrderedDict
from collections.abc import Mapping, Sequence
from pathlib import Path
from types import NotImplementedType
from typing import assert_type

from pyensembl import Transcript

from varcode import (
    AssembledAllele,
    Candidate, CandidatePath, classification_priority, get_effect_candidates,
    iter_candidate_effects, iter_candidate_paths,
    ContigResolution, ContigResolver, CrypticSplicePolicy,
    DNAEffectAnnotator, DNAHypothesis, DNASearchAnchor, EffectAnnotator, EffectCollection, EffectReference,
    GenomicSpliceJunction, GeneLike, Genotype, InvalidContigError, MutationEffect, MutantTranscript,
    ORFAnnotation, SpanningDeletionRecord, TranscriptLike, UnifiedDNAEffectAnnotator,
    Unresolved, Variant, VariantCollection, Zygosity,
    build_dna_search_anchor, combine_cis_snvs,
    effect_collection_from_json, effect_collection_from_serializable_repr,
    effect_collection_to_json, effect_from_json, effect_from_serializable_repr,
    effect_serialization_classes, effect_to_json,
    get_annotator, get_default_annotator, register_annotator, registered_annotators,
    resolve_annotator, set_default_annotator, use_annotator,
    load_maf, load_maf_dataframe, load_vcf, normalize_contig_name, resolve_contig_name,
)
from varcode.nucleotides import normalize_nucleotide_string
from varcode.vcf import read_vcf_into_dataframe


def check_candidate_traversal(effect: MutationEffect) -> None:
    assert_type(get_effect_candidates(effect), tuple[Candidate, ...])
    for terminal, path in iter_candidate_paths(effect):
        assert_type(terminal, MutationEffect | None)
        assert_type(path, CandidatePath)
        if terminal is not None:
            assert_type(classification_priority(terminal), int)
    for terminal in iter_candidate_effects(effect):
        assert_type(terminal, MutationEffect | None)
    iter_candidate_paths("not an effect")  # type: ignore[arg-type]


def check_inputs(variant: Variant, other: object, path: Path) -> None:
    assert_type(load_vcf(path), VariantCollection)
    assert_type(load_maf(path), VariantCollection)
    load_maf_dataframe(path)
    read_vcf_into_dataframe(path)
    assert_type(normalize_nucleotide_string("ACGT"), str)
    assert_type(variant == other, bool)
    assert_type(variant < other, bool)
    assert_type(normalize_contig_name("chr1"), str)
    assert_type(resolve_contig_name("chr1", ["1"]), ContigResolution)
    assert_type(ContigResolver(["chrM"]).resolve("MT").resolved_name, str | None)
    assert_type(
        ContigResolver(["chrM"], contig_lengths={"chrM": 16569}).resolve("MT").resolved_name,
        str | None)

    # These ignores must remain necessary (warn_unused_ignores is enabled),
    # proving the public exports reject incompatible arguments statically.
    load_vcf(b"variants.vcf")  # type: ignore[arg-type]
    load_maf(42)  # type: ignore[arg-type]
    load_maf_dataframe(b"variants.maf")  # type: ignore[arg-type]
    read_vcf_into_dataframe(42)  # type: ignore[arg-type]
    normalize_nucleotide_string(None)  # type: ignore[arg-type]
    normalize_contig_name(None)  # type: ignore[arg-type]
    resolve_contig_name("chr1", [1])  # type: ignore[list-item]


class DecliningAnnotator:
    name = "declines"

    def annotate_on_transcript(
            self, variant: Variant, transcript: Transcript,
    ) -> NotImplementedType:
        return NotImplemented


class InvalidAnnotator:
    name = "invalid"

    def annotate_on_transcript(self, variant: Variant, transcript: Transcript) -> None:
        return None


def check_annotation_workflow(
        variant: Variant, variants: VariantCollection, transcript: Transcript,
) -> None:
    assert_type(variant.effects(), EffectCollection)
    assert_type(variants.effects(), EffectCollection)
    assert_type(next(iter(variants)), Variant)
    assert_type(variants[0], Variant)
    assert_type(variants[:1], Sequence[Variant])
    assert_type(variants.filter(lambda v: v.start > 0), VariantCollection)
    assert_type(variants.groupby(lambda v: v.contig), dict[str, VariantCollection])
    assert_type(variants.spanning_deletion_records, tuple[SpanningDeletionRecord, ...])
    assert_type(combine_cis_snvs(variants, sample="tumor"), VariantCollection)
    effects = variant.effects(annotator="fast", raise_on_error=False)
    assert_type(next(iter(effects)), MutationEffect)
    assert_type(effects[0], MutationEffect)
    assert_type(effects[:1], Sequence[MutationEffect])
    assert_type(effects.filter(lambda e: e.modifies_protein_sequence is True), EffectCollection)
    assert_type(effects.groupby(lambda e: e.transcript_id), dict[str | None, EffectCollection])
    assert_type(effects.clone_with_new_elements([]), EffectCollection)
    effects[0].not_a_real_field  # type: ignore[attr-defined]
    assert_type(effects.drop_silent_and_noncoding(), EffectCollection)
    assert_type(effects.top_priority_effect(), MutationEffect)
    assert_type(variant.effect_on_transcript(transcript), MutationEffect)
    effect = variant.effect_on_transcript(transcript)
    assert_type(effect.mutant_protein_sequence, str | None)
    assert_type(effect.original_protein_sequence, str | None)
    assert_type(effect.modifies_protein_sequence, bool | None)
    assert_type(effect.modifies_coding_sequence, bool | None)
    assert_type(effect.aa_mutation_start_offset, int | None)
    assert_type(effect.aa_mutation_end_offset, int | None)
    assert_type(effect.mutant_transcript, MutantTranscript | None)
    assert_type(effect.baseline_transcript, MutantTranscript | None)
    assert_type(effect.orf_annotation, ORFAnnotation | None)
    assert_type(effect.gene_id, str | None)
    assert_type(effect.transcript_id, str | None)
    assert_type(effect.short_description, str)
    assert_type(effect.format_short_description(deletion_notation="range"), str)
    if effect.mutant_protein_sequence is not None:
        assert_type(effect.mutant_protein_sequence, str)
    assert_type(InvalidContigError("missing"), InvalidContigError)
    unresolved = Unresolved(
        variant, transcript, "germline_composition", reference_effect=effect)
    assert_type(unresolved.reference_effect, MutationEffect | None)
    assert_type(unresolved.priority_class, type[MutationEffect])
    Unresolved(variant, transcript, "missing", reference_effect="missense")  # type: ignore[arg-type]

    plugin = DecliningAnnotator()
    assert_type(register_annotator(plugin), DecliningAnnotator)
    assert_type(registered_annotators(), tuple[str, ...])
    assert_type(registered_annotators(include_deprecated=True), tuple[str, ...])
    assert_type(get_annotator(plugin.name), EffectAnnotator)
    assert_type(get_default_annotator(), EffectAnnotator)
    assert_type(resolve_annotator(None), EffectAnnotator)
    assert_type(resolve_annotator(plugin), EffectAnnotator)
    assert_type(set_default_annotator("fast"), None)
    with use_annotator(plugin) as bound:
        assert_type(bound, None)
        assert_type(plugin.annotate_on_transcript(variant, transcript), NotImplementedType)

    assert_type(effect_serialization_classes(), Mapping[tuple[str, str], type])
    assert_type(EffectReference(0), EffectReference)
    assert_type(AssembledAllele("ATG"), AssembledAllele)
    assert_type(effect_to_json(effect), str)
    assert_type(effect_from_json("{}"), MutationEffect)
    assert_type(effect_from_serializable_repr({}), MutationEffect)
    assert_type(effect_collection_to_json(effects), str)
    assert_type(effect_collection_from_json("{}"), EffectCollection)
    assert_type(effect_collection_from_serializable_repr({}), EffectCollection)
    assert_type(CrypticSplicePolicy(scan_flank=100, max_sites_per_junction=4), CrypticSplicePolicy)
    dna_ann = DNAEffectAnnotator()
    assert_type(register_annotator(dna_ann), DNAEffectAnnotator)
    unified = UnifiedDNAEffectAnnotator()
    assert_type(register_annotator(unified), DNAEffectAnnotator)
    junction = GenomicSpliceJunction("17", 101, 200, "-", 201, 100)
    assert_type(junction.intron_start, int)
    hyp = DNAHypothesis(effect=effect, novel_junctions=(junction,))
    assert_type(hyp.short_description, str)
    anchor = build_dna_search_anchor(variant, transcript)
    assert_type(anchor, DNASearchAnchor)
    assert_type(anchor.exit_frame, int | None)

    variant.effects(raise_on_error="no")  # type: ignore[arg-type]
    variants.effects(validate_reference="no")  # type: ignore[arg-type]
    get_annotator(42)  # type: ignore[arg-type]
    registered_annotators(include_deprecated="yes")  # type: ignore[arg-type]
    use_annotator(None)  # type: ignore[arg-type]
    register_annotator(InvalidAnnotator())  # type: ignore[type-var]
    variant.effects(annotator=InvalidAnnotator())  # type: ignore[arg-type]
    effect.mutant_protein_sequence.upper()  # type: ignore[union-attr]
    effect.format_short_description(deletion_notation="invalid")  # type: ignore[arg-type]
    effect_to_json("not an effect")  # type: ignore[arg-type]
    effect_from_json(42)  # type: ignore[arg-type]
    effect_collection_to_json(effect)  # type: ignore[arg-type]
    effect_collection_from_json(42)  # type: ignore[arg-type]


class _CustomGeneAdapter:
    id: str = "ENSG_CUSTOM"
    name: str = "CUSTOM"
    contig: str = "1"
    start: int = 100
    end: int = 500
    strand: str = "+"


class _CustomTranscriptAdapter:
    id: str = "ENST_CUSTOM"
    name: str = "CUSTOM-001"
    contig: str = "1"
    start: int = 100
    end: int = 500
    strand: str = "+"
    gene_id: str = "ENSG_CUSTOM"
    protein_sequence: str | None = "MKL"
    complete: bool = True
    is_protein_coding: bool = True
    gene: GeneLike | None = _CustomGeneAdapter()


class _CustomVariantCollection(VariantCollection):
    pass


class _CustomEffectCollection(EffectCollection):
    pass


def check_typed_traversal_and_protocols(
        path: Path,
        custom_vc: _CustomVariantCollection,
        custom_ec: _CustomEffectCollection,
) -> None:
    loaded = load_vcf(path)
    assert_type(loaded, VariantCollection)
    assert_type(loaded.union(loaded), VariantCollection)
    assert_type(loaded.intersection(loaded), VariantCollection)
    assert_type(loaded.difference(loaded), VariantCollection)
    assert_type(loaded.for_sample("tumor"), VariantCollection)
    assert_type(loaded.heterozygous_in("tumor"), VariantCollection)
    assert_type(loaded.homozygous_alt_in("tumor"), VariantCollection)
    assert_type(loaded.samples, list[str])
    assert_type(loaded.has_sample_data(), bool)
    assert_type(loaded.reference_names(), set[str])
    assert_type(loaded.original_reference_names(), set[str])
    assert_type(loaded.gene_counts(), dict[str, int])
    assert_type(loaded.groupby_gene_name(), dict[str, VariantCollection])
    assert_type(loaded.groupby_gene_id(), dict[str, VariantCollection])
    assert_type(loaded.groupby_gene(), dict[GeneLike, VariantCollection])
    assert_type(loaded.filter_by_gene_expression({"ENSG1": 2.0}), VariantCollection)
    assert_type(loaded.filter_by_transcript_expression({"ENST1": 2.0}), VariantCollection)
    assert_type(loaded.clone_without_ucsc_data(), VariantCollection)
    assert_type(loaded.genotype(loaded[0], "tumor"), Genotype | None)
    assert_type(loaded.zygosity(loaded[0], "tumor"), Zygosity)

    # Subclass preservation on VariantCollection and EffectCollection
    assert_type(custom_vc.union(loaded), _CustomVariantCollection)
    assert_type(custom_vc.for_sample("tumor"), _CustomVariantCollection)
    assert_type(custom_vc.groupby_gene_name(), dict[str, _CustomVariantCollection])
    assert_type(custom_ec.drop_silent_and_noncoding(), _CustomEffectCollection)
    assert_type(custom_ec.groupby_gene_id(), dict[str | None, _CustomEffectCollection])

    # End-to-end typed load_vcf -> VariantCollection -> effects -> EffectCollection traversal
    effects = loaded.effects()
    assert_type(effects.groupby_variant(), dict[Variant, EffectCollection])
    assert_type(effects.groupby_transcript(), dict[TranscriptLike | None, EffectCollection])
    assert_type(effects.groupby_transcript_id(), dict[str | None, EffectCollection])
    assert_type(effects.groupby_transcript_name(), dict[str | None, EffectCollection])
    assert_type(effects.groupby_gene(), dict[GeneLike | None, EffectCollection])
    assert_type(effects.groupby_gene_id(), dict[str | None, EffectCollection])
    assert_type(effects.groupby_gene_name(), dict[str | None, EffectCollection])
    assert_type(effects.gene_counts(), dict[str | None, int])
    assert_type(effects.top_priority_effect_per_variant(), OrderedDict[Variant, MutationEffect])
    assert_type(effects.top_priority_effect_per_transcript_id(), OrderedDict[str | None, MutationEffect])
    assert_type(effects.top_priority_effect_per_gene_id(), OrderedDict[str | None, MutationEffect])
    assert_type(effects.effect_expression({"ENST1": 1.5}), OrderedDict[MutationEffect, float])
    top_expr = effects.top_expression_effect({"ENST1": 1.5})
    assert_type(top_expr, MutationEffect | None)
    assert_type(effects.filter_by_gene_expression({"ENSG1": 1.0}), EffectCollection)
    assert_type(effects.filter_by_transcript_expression({"ENST1": 1.0}), EffectCollection)
    assert_type(effects.orf_annotations, tuple[ORFAnnotation, ...])

    # Structural TranscriptLike / GeneLike adapters without inheriting pyensembl
    tx_adapter: TranscriptLike = _CustomTranscriptAdapter()
    gene_adapter: GeneLike = _CustomGeneAdapter()
    assert_type(tx_adapter.id, str)
    assert_type(gene_adapter.name, str)
    custom_effect = loaded[0].effect_on_transcript(tx_adapter)
    assert_type(custom_effect, MutationEffect)
    assert_type(custom_effect.transcript, TranscriptLike | None)
    assert_type(custom_effect.reference_transcript, TranscriptLike | None)
    assert_type(custom_effect.gene, GeneLike | None)
    if custom_effect.transcript is not None:
        assert_type(custom_effect.transcript.id, str)
    if custom_effect.gene is not None:
        assert_type(custom_effect.gene.name, str)

    # Negative checks: unchecked None / wrong types must be rejected statically
    custom_effect.transcript.id  # type: ignore[union-attr]
    custom_effect.gene.name  # type: ignore[union-attr]
    top_expr.short_description  # type: ignore[union-attr]
    loaded.for_sample(123)  # type: ignore[arg-type]
    loaded[0].effect_on_transcript("not a transcript")  # type: ignore[arg-type]
    effects.effect_expression("not a mapping")  # type: ignore[arg-type]
