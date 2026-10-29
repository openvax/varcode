"""Static consumer checks; checked by mypy, never executed as tests."""

from pathlib import Path
from types import NotImplementedType
from typing import assert_type

from pyensembl import Transcript

from varcode import (
    ContigResolution, ContigResolver, EffectAnnotator, EffectCollection,
    InvalidContigError, MutationEffect, MutantTranscript, ORFAnnotation,
    Unresolved, Variant, VariantCollection,
    get_annotator, get_default_annotator, register_annotator, resolve_annotator,
    set_default_annotator, use_annotator,
    load_maf, load_maf_dataframe, load_vcf, normalize_contig_name, resolve_contig_name,
)
from varcode.nucleotides import normalize_nucleotide_string
from varcode.vcf import read_vcf_into_dataframe


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
    effects = variant.effects(annotator="fast", raise_on_error=False)
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
    assert_type(get_annotator(plugin.name), EffectAnnotator)
    assert_type(get_default_annotator(), EffectAnnotator)
    assert_type(resolve_annotator(None), EffectAnnotator)
    assert_type(resolve_annotator(plugin), EffectAnnotator)
    assert_type(set_default_annotator("fast"), None)
    with use_annotator(plugin) as bound:
        assert_type(bound, None)
        assert_type(plugin.annotate_on_transcript(variant, transcript), NotImplementedType)

    variant.effects(raise_on_error="no")  # type: ignore[arg-type]
    variants.effects(validate_reference="no")  # type: ignore[arg-type]
    get_annotator(42)  # type: ignore[arg-type]
    use_annotator(None)  # type: ignore[arg-type]
    register_annotator(InvalidAnnotator())  # type: ignore[type-var]
    variant.effects(annotator=InvalidAnnotator())  # type: ignore[arg-type]
    effect.mutant_protein_sequence.upper()  # type: ignore[union-attr]
    effect.format_short_description(deletion_notation="invalid")  # type: ignore[arg-type]
