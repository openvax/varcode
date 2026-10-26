"""Native coding placements work without a catalogue, including incomplete products."""

from dataclasses import replace
from types import SimpleNamespace

import pytest

from varcode import EffectCollection, ReferenceMismatchError, Variant, annotate_orf_variants, to_proteome
from varcode.effects import Unresolved
from .test_extra_orfs import synthetic
from . import transcript_orf_fixture

orf_genome = transcript_orf_fixture.orf_genome

orfdata = pytest.importorskip("orfdata")


def source_orf(reverse=False, split=False, *, start=3, end=None, phase=0,
               start_complete=True, stop_complete=True):
    tx = synthetic(reverse, split)
    context = orfdata.ReferenceContext.from_genome(tx.genome)
    exons = tuple((e.start, e.end) for e in sorted(tx.exons, key=lambda e: e.start, reverse=reverse))
    rna = orfdata.SourceRNA("source-rna", tx.sequence, context, "test-source", contig=tx.contig,
                           strand=tx.strand, exons=exons, gene_id="gene", gene_name="GENE")
    reference = orfdata.TranscriptORF.at(
        rna, start, end=end, phase=phase, start_complete=start_complete, stop_complete=stop_complete,
        orf_id="source-orf", source="test-source", orf_type="source-defined",
        evidence=orfdata.ORFEvidence(ribosome_studies=("reported-study",)))
    return tx, reference


def assert_roundtrip(annotation):
    assert annotation.orf_id in str(annotation.effect)
    encoded = annotation.effect.to_json()
    restored = type(annotation.effect).from_json(encoded)
    assert restored.to_json() == encoded
    assert restored.orf_annotation.to_dict() == annotation.to_dict()
    assert restored.orf.record is None
    assert restored.orf.evidence == annotation.reference.evidence
    return restored.orf_annotation


@pytest.mark.parametrize("reverse,split", [(False, False), (True, False), (False, True), (True, True)])
def test_source_defined_complete_orf_uses_public_effects_and_exports(reverse, split, tmp_path):
    tx, reference = source_orf(reverse, split)
    variant = tx.variant(6, "A", "G")
    effects = variant.effects(extra_orfs=[reference])
    annotation, = effects.orf_annotations
    assert annotation.record is None and annotation.reference is reference
    assert annotation.effect.short_description == "p.K2E"
    assert annotation.outcomes[0].mutant.protein_sequence == "MEPG"
    assert annotation.reference.transcript.id == "source-rna"
    assert_roundtrip(annotation)
    loaded = EffectCollection.from_json(effects.to_json())
    assert loaded.orf_annotations[0].to_dict() == annotation.to_dict()
    proteome = to_proteome(effects, sample_id="sample-1")
    protein, = proteome
    assert protein.producer.origin == "nc_orf"
    assert protein.producer.orf_kind == "source-defined"
    assert protein.producer.orf_ids == ("source-orf",)
    assert protein.producer.tool == "varcode"
    path = tmp_path / "products.json.gz"
    proteome.save(path)
    assert orfdata.Proteome.load(path).identity == proteome.identity


@pytest.mark.parametrize("reverse,split", [(False, False), (True, False), (False, True), (True, True)])
@pytest.mark.parametrize("boundary", (False, None))
def test_partial_known_frame_supports_local_change_and_keeps_missing_ends(reverse, split, boundary):
    # The first complete codon begins at RNA offset 6; one preceding base is retained.
    tx, reference = source_orf(reverse, split, start=5, end=14, phase=1,
                              start_complete=boundary, stop_complete=boundary)
    assert reference.coding_sequence == "GAAACCCGG"
    assert reference.protein_sequence == "KP"
    annotation, = tx.variant(6, "A", "G").effects(extra_orfs=[reference]).orf_annotations
    assert annotation.effect.short_description == "p.K1E"
    outcome, = annotation.outcomes
    assert outcome.baseline.protein_sequence == "KP" and outcome.mutant.protein_sequence == "EP"
    assert outcome.mutant.start_codon_present is boundary
    assert outcome.mutant.termination == ("sequence_boundary" if boundary is False else "unknown_boundary")
    assert outcome.mutant.evidence["coding_interval"] == (5, 14)
    assert annotation.unavailable_operations["complete_mutant_protein"] == "reference_coding_fragment"
    assert "start_loss" in annotation.unavailable_operations and "stop_loss" in annotation.unavailable_operations
    assert_roundtrip(annotation)
    protein, = to_proteome(orf_annotations=[annotation])
    assert protein.sequence == "EP"
    assert protein.producer.start_complete is boundary and protein.producer.stop_complete is boundary


@pytest.mark.parametrize("reverse", (False, True))
def test_partial_initial_methionine_is_not_assumed_to_be_an_initiator(reverse):
    tx, reference = source_orf(reverse, start=3, end=15, start_complete=False, stop_complete=False)
    annotation = annotate_orf_variants([tx.variant(3, "A", "G")], reference)
    assert annotation.effect.short_description == "p.M1V"
    assert annotation.outcomes[0].mutant.protein_sequence == "VKPG"


@pytest.mark.parametrize("offset,ref,alt,expected,protein", [
    (9, "CCC", "", "Deletion", "KG"),
    (9, "C", "", "FrameShift", "KP"),
    (6, "A", "T", "PrematureStop", ""),
    (8, "A", "G", "Silent", "KPG"),
])
def test_partial_internal_edits_keep_the_supplied_fragment(offset, ref, alt, expected, protein):
    tx, reference = source_orf(start=6, end=15, start_complete=False, stop_complete=False)
    annotation = annotate_orf_variants([tx.variant(offset, ref, alt)], reference)
    assert type(annotation.effect).__name__ == expected
    assert annotation.outcomes[0].mutant.protein_sequence == protein
    assert_roundtrip(annotation)


def test_partial_codon_edit_retains_reference_and_unavailable_outcome():
    tx, reference = source_orf(start=5, end=14, phase=1, start_complete=False, stop_complete=False)
    annotation = annotate_orf_variants([tx.variant(5, "G", "A")], reference)
    assert isinstance(annotation.effect, Unresolved)
    assert annotation.reference is reference
    assert annotation.unresolved_reason == "partial_orf_boundary_or_partial_codon_edit"
    assert annotation.outcomes[0].baseline.protein_sequence == "KP"
    assert annotation.outcomes[0].mutant is None
    assert_roundtrip(annotation)
    exported = to_proteome(orf_annotations=[annotation])
    assert len(exported) == 0 and any(r.status == "unresolved" for r in exported.records)
    with pytest.raises(ReferenceMismatchError):
        annotate_orf_variants([tx.variant(6, "C", "G")], reference)


def test_unmapped_source_is_retained_as_explicitly_unavailable():
    tx, reference = source_orf()
    rna = replace(reference.transcript, contig=None, strand=None, exons=())
    unmapped = orfdata.TranscriptORF.at(rna, 3, orf_id="unmapped", source="test")
    annotation, = tx.variant(6, "A", "G").effects(extra_orfs=[unmapped]).orf_annotations
    assert annotation.reference is unmapped
    assert annotation.unresolved_reason == "orf_genomic_mapping_unavailable"
    assert_roundtrip(annotation)


def native_change(reference):
    """Pick an internal single-codon edit; Biopython supplies an independent expected residue."""
    from Bio.Seq import Seq

    origin = reference.coding_origin
    offset = reference.coding_phase + 30
    codon = reference.coding_sequence[offset:offset+3]
    for base in "ACGT":
        alternate = base + codon[1:]
        amino = str(Seq(alternate).translate())
        if amino not in ("*", str(Seq(codon).translate())):
            break
    ref, alt = codon[0], base
    if reference.strand == "-":
        ref, alt = str(Seq(ref).complement()), str(Seq(alt).complement())
    variant = Variant(reference.contig, reference.coding_positions[offset], ref, alt,
                      genome=reference.transcript.genome)
    index = (offset - origin) // 3
    expected = reference.protein_sequence[:index] + amino + reference.protein_sequence[index+1:]
    return variant, expected


@pytest.mark.parametrize("transcript_id,boundaries", [
    ("ENST00000395826", (True, True)),  # complete POMC control
    ("ENST00000084795", (False, True)),  # RPL18: native 5-prime partial codon
    ("ENST00000376028", (True, False)),  # KAZN: native cds_end_NF with trailing partial codon
    ("ENST00000339995", (True, True)),  # EIF4G2: annotated GTG initiation
])
def test_native_annotated_complete_and_partial_cds_consumer(orf_genome, transcript_id, boundaries):
    reference = orfdata.TranscriptORF.from_annotation(orf_genome.transcript_by_id(transcript_id))
    assert (reference.start_complete, reference.stop_complete) == boundaries
    variant, expected = native_change(reference)
    annotation, = variant.effects(extra_orfs=[reference]).orf_annotations
    assert annotation.record is None
    assert annotation.outcomes[0].mutant.protein_sequence == expected
    assert annotation.effect.reference_transcript is reference.transcript
    assert_roundtrip(annotation)
    protein, = to_proteome(orf_annotations=[annotation])
    assert protein.producer.origin == "reference_cds"
    assert protein.producer.start_complete is boundaries[0]
    assert protein.producer.stop_complete is boundaries[1]


@pytest.mark.parametrize("product_id,host_id,source_gene,host_gene", [
    ("ENST00000577949", "ENST00000390687", "SNURF", "SNRPN"),
    ("ENST00000637854", "ENST00000325301", "MIURF", "MIEF1"),
])
def test_native_bicistronic_models_keep_full_host_rnas_and_product_genes(
        orf_genome, product_id, host_id, source_gene, host_gene):
    source, host = (orf_genome.transcript_by_id(t) for t in (product_id, host_id))
    own = orfdata.TranscriptORF.from_annotation(source)
    carried = orfdata.TranscriptORF.from_annotation(source, on=host)
    assert carried.transcript is host and source.sequence != host.sequence
    assert (carried.source_gene_name, carried.gene_name) == (source_gene, host_gene)
    variant, expected = native_change(carried)
    annotations = variant.effects(extra_orfs=[own, carried]).orf_annotations
    assert len(annotations) == 2
    assert {a.reference.transcript_id for a in annotations} == {source.id, host.id}
    for annotation in annotations:
        assert annotation.outcomes[0].mutant.protein_sequence == expected
        assert_roundtrip(annotation)
    protein, = to_proteome(orf_annotations=[annotations[1]])
    assert protein.producer.origin == "reference_cds"
    assert set(protein.producer.gene_names) == {source_gene, host_gene}


def test_recoded_annotation_remains_available_with_explicit_unsupported_calculation(orf_genome):
    reference = orfdata.TranscriptORF.from_annotation(orf_genome.transcript_by_id("ENST00000419783"))
    variant, _ = native_change(reference)
    annotation, = variant.effects(extra_orfs=[reference]).orf_annotations
    assert annotation.unresolved_reason == "orf_reference_translation_requires_recoding"
    assert annotation.outcomes[0].baseline.protein_sequence == reference.protein_sequence
    assert annotation.outcomes[0].mutant is None
    assert_roundtrip(annotation)


def test_native_alternative_starts_remain_distinct_on_one_rna(orf_genome):
    transcript = orf_genome.transcript_by_id("ENST00000644866")
    source = transcript_orf_fixture.DATA["source_isoforms"]
    references = [orfdata.TranscriptORF.from_annotation(transcript)]
    for identifier, start in (("P09038-1", 179), ("P09038-3", 221)):
        references.append(orfdata.TranscriptORF.at(
            transcript, start, orf_id=identifier, source="UniProt 2026_03",
            protein=source[identifier]["sequence"], orf_type="alternative_initiation"))
    assert [ref.start_codon for ref in references] == ["ATG", "CTG", "CTG"]
    variant, _ = native_change(references[0])
    effects = variant.effects(extra_orfs=references)
    assert len(effects.orf_annotations) == 3
    assert len({a.orf_id for a in effects.orf_annotations}) == 3
    assert len({a.outcomes[0].mutant.protein_sequence for a in effects.orf_annotations}) == 3
    for annotation in effects.orf_annotations:
        assert annotation.reference.transcript is transcript
        assert annotation.outcomes[0].mutant.start_codon_present
        assert_roundtrip(annotation)


def test_low_level_native_entry_point_retains_its_generating_tool():
    tx = synthetic()
    reference = tx.record().bind(tx)
    annotation = annotate_orf_variants([tx.variant(6, "A", "G")], reference)
    protein, = to_proteome(orf_annotations=[annotation])
    assert protein.producer.tool == "varcode" and protein.producer.method == "transcript_model"
    restored = type(annotation.effect).from_json(annotation.effect.to_json())
    replayed, = to_proteome([restored])
    assert replayed.producer.tool == "varcode" and replayed.producer.tool_version == protein.producer.tool_version


def test_annotation_without_boundary_or_phase_keeps_an_explicit_unknown_frame():
    tx = synthetic()
    tx.protein_id = "synthetic-partial"
    tx.protein = SimpleNamespace(protein_version=None)
    tx.protein_sequence = "PMK"
    tx.gene_id, tx.gene_name = "gene", "GENE"
    tx.contains_start_codon = tx.contains_stop_codon = tx.complete = False
    tx.coding_sequence_position_ranges = ((100, 108),)
    reference = orfdata.TranscriptORF.from_annotation(tx)
    assert reference.coding_origin is None and reference.coding_phase is None
    annotation, = tx.variant(6, "A", "G").effects(extra_orfs=[reference]).orf_annotations
    assert annotation.unresolved_reason == "orf_reading_frame_unknown"
    assert annotation.outcomes[0].baseline.protein_sequence == "PMK"
    assert_roundtrip(annotation)
