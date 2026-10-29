"""Native products and comparisons survive shared, offline proteome export."""

from dataclasses import replace
import io
import json
import subprocess
import sys

import pytest

from varcode import (
    EffectCandidate, EffectCollection, MutantTranscript, ReferenceSegment,
    SpliceOutcomeSet, StructuralVariant, Variant, to_proteome,
)
from varcode.effect_hypotheses import (
    ClassifiedOutcome, EffectHypothesis, RealizedTranscriptProduct,
    merge_classified_outcomes,
)
from varcode.effects import (
    CrypticAcceptor, CrypticDonor, Deletion, Failure, GeneFusion, HypothesisLimit, Intergenic, Intragenic, IntronRetention,
    SpliceAcceptor, SpliceDonor, Substitution, Unresolved,
)
from varcode.effects.effect_classes import PhaseCandidateSet
from varcode.germline import PhasePartition
from varcode.orf import ORFAnnotation, ORFOutcome, ORFProduct, annotate_orf_variants
from varcode.orf_serialization import _Feature, _ReferenceGenome, _ReferenceTranscript

orfdata = pytest.importorskip("orfdata")


@pytest.fixture
def transcript():
    genome = _ReferenceGenome(dict(reference_name="fixture", species="homo_sapiens", taxon_id=9606,
                                  annotation_name="synthetic", annotation_version="1"))
    return _ReferenceTranscript(
        dict(id="tx-1", name="TX1", transcript_id="tx-1", contig="1", strand="+", start=1, end=12,
             sequence="ATGGCTAAATAA", protein_sequence="MAK", protein_id="protein-1",
             biotype="protein_coding", is_protein_coding=True, complete=True,
             contains_start_codon=True, contains_stop_codon=True,
             start_codon_positions=(1, 2, 3), start_codon_spliced_offsets=(0, 1, 2),
             stop_codon_positions=(10, 11, 12), stop_codon_spliced_offsets=(9, 10, 11),
             coding_sequence="ATGGCTAAA", coding_sequence_position_ranges=((1, 9),),
             five_prime_utr_sequence="", three_prime_utr_sequence="", on_backward_strand=False,
             on_forward_strand=True, topology="linear"),
        genome, (_Feature(1, 12, "exon-1"),), _Feature(1, 12, "gene-1", "GENE1"))


@pytest.fixture
def effect(transcript):
    result = Substitution(Variant("1", 4, "G", "A", genome=transcript.genome), transcript, 1, "A", "T")
    result.annotator = "transcript_model"
    result.annotator_version = "fixture-1"
    return result


@pytest.fixture
def reference_orf(transcript):
    return orfdata.ORFRecord(
        orf_id="nc/orf|1", gene_id="gene-1", gene_name="GENE1", transcript_id="tx-1",
        contig="1", strand="+", coding_blocks=((1, 12),), protein_sequence="MAK",
        biotype="lncRNA", taxon_id=9606, species="homo_sapiens", source="fixture",
        ribosome_studies=(), peptides=(), genome_build="fixture",
        reference_context=orfdata.ReferenceContext.from_genome(transcript.genome)).bind(transcript)


def comparisons(proteome):
    return [r for r in proteome.records if r.kind == "comparison"]


def graph(proteome):
    return dict(entries=[e.to_dict() for e in proteome],
                records=[r.to_dict() for r in proteome.records],
                links=[link.to_dict() for link in proteome.links])


def product(sequence, *, start=True):
    return RealizedTranscriptProduct(None, None, sequence, start_codon_present=start)


def outcome(effect, baseline="MAK", mutant="MTK", phase=(), index=0):
    return ClassifiedOutcome(EffectHypothesis(phase=phase, enumeration_index=index), effect,
                             product(baseline), product(mutant), change_aa_offset=1)


def test_coding_substitution_and_offline_roundtrip(effect, tmp_path, monkeypatch):
    predicted = to_proteome(EffectCollection([effect]), sample_id="tumor 1", run_id="run/1")
    entry, = predicted
    assert entry.sequence == "MTK"
    assert entry.producer.origin == "reference_cds"
    assert entry.producer.tool == "varcode"
    assert entry.producer.tool_version == "fixture-1"
    assert entry.producer.exporters == ("varcode",)
    assert entry.producer.transcript_ids == ("tx-1",)
    assert entry.producer.gene_names == ("GENE1",)
    assert entry.producer.start_complete is None  # reference completeness does not prove mutant extent
    comparison, = comparisons(predicted)
    assert comparison.exporters == ("varcode",)
    data = comparison.to_dict()["payload"]
    assert data["baseline"]["protein_sequence"] == "MAK"
    assert data["mutant"]["protein_sequence"] == "MTK"
    assert data["intervals"] == {"mutant_aa": [1, 2], "baseline_aa": [1, 2]}
    assert entry.occurrences == () and entry.orf is None and entry.derivation is None
    assert {link.relation for link in predicted.links} == {"produces", "has_outcome", "contains"}
    path = tmp_path / "predicted.json.gz"
    predicted.save(path)
    monkeypatch.setattr(_ReferenceGenome, "transcripts_at_locus", lambda *args: pytest.fail("genome accessed"))
    restored = orfdata.Proteome.load(path)
    assert restored.identity == predicted.identity
    assert graph(restored) == graph(predicted)
    # The embedded graph stays inert JSON; it is never passed to serializable.from_json.
    assert next(r for r in restored.records if r.kind == "event").to_dict()["payload"]["native"]["orf_effect_schema"] == 1
    assert json.loads(json.dumps(graph(restored))) == json.loads(json.dumps(graph(predicted)))


def test_full_products_remain_distinct_when_local_changes_match(effect):
    prefix = "M" + "A" * 75
    first = outcome(effect, baseline=prefix + "Q", mutant=prefix + "R")
    second = outcome(effect, baseline=prefix + "Q", mutant=prefix + "S", index=1)
    candidates = merge_classified_outcomes([first, second])
    assert len(candidates) == 2  # differences beyond the former local merge window
    effect.realized_candidates = candidates
    predicted = to_proteome([effect])
    assert {e.sequence for e in predicted} == {prefix + "R", prefix + "S"}
    assert len(comparisons(predicted)) == 2


def test_equal_sequences_preserve_baselines_phase_and_group_associations(effect):
    background = Variant("1", 7, "A", "G", genome=effect.variant.genome)
    a = outcome(effect, baseline="MAK", phase=((background, "cis"),))
    b = outcome(effect, baseline="MAR", phase=((background, "trans"),), index=1)
    predicted = to_proteome(merge_classified_outcomes([a, b]))
    assert len(predicted.proteins) == len(comparisons(predicted)) == 2
    assert len(predicted.by_sequence()["MTK"]) == 2
    assert {r.payload["baseline"]["protein_sequence"] for r in comparisons(predicted)} == {"MAK", "MAR"}
    assert len(predicted.fasta_index()) == 2
    grouped, = predicted.fasta_index(deduplicate=True).values()
    assert {e.accession for e in grouped} == {e.accession for e in predicted}
    for deduplicate, count in [(False, 2), (True, 1)]:
        fasta = io.StringIO()
        predicted.write_fasta(fasta, deduplicate=deduplicate)
        assert sum(line.startswith(">") for line in fasta.getvalue().splitlines()) == count
    for comparison in comparisons(predicted):
        context, = comparison.payload["contexts"]
        assert context["ordinal_key"] == ()
        assert context["phase_mass"] is not None
        assert "probability" in context
        assert comparison.payload["hypothesis"]


def test_empty_deletion_interval_is_retained(effect):
    deletion = Deletion(effect.variant, effect.transcript, 1, "A")
    comparison, = comparisons(to_proteome([deletion]))
    assert comparison.payload["intervals"]["mutant_aa"] == (1, 1)
    assert comparison.payload["intervals"]["baseline_aa"] == (1, 2)
    assert comparison.payload["mutant"]["protein_sequence"] == "MK"


def test_native_orf_and_failed_binding(transcript, tmp_path):
    record = orfdata.ORFRecord(
        orf_id="nc/orf|1", gene_id="gene-1", gene_name="GENE1", transcript_id="tx-1",
        contig="1", strand="+", coding_blocks=((1, 12),), protein_sequence="MAK",
        biotype="lncRNA", taxon_id=9606, species="homo_sapiens", source="fixture",
        ribosome_studies=(), peptides=(), genome_build="fixture",
        reference_context=orfdata.ReferenceContext.from_genome(transcript.genome))
    reference = record.bind(transcript)
    variant = Variant("1", 4, "G", "A", genome=transcript.genome)
    annotation = annotate_orf_variants([variant], reference)
    failed = ORFAnnotation(replace(record, orf_id="missing", transcript_id="missing"), None,
                           variants=(variant,), unresolved_reason="transcript absent")
    predicted = to_proteome(orf_annotations=[annotation, failed])
    entry, = predicted
    assert entry.sequence == "MTK"
    assert entry.producer.origin == "nc_orf"
    assert entry.producer.orf_ids == ("nc/orf|1",)
    assert entry.producer.provider == "fixture"
    assert entry.producer.stop_complete is True
    assert {r.status for r in comparisons(predicted)} == {"product_available", "failed_binding"}
    missing = next(r for r in comparisons(predicted) if r.status == "failed_binding")
    assert missing.payload["omission_reason"] == "transcript absent"
    path = tmp_path / "orf.json.gz"
    predicted.save(path)
    assert graph(orfdata.Proteome.load(path)) == graph(predicted)
    assert to_proteome([annotation.effect]).proteins[0].producer.origin == "nc_orf"
    # Attached and standalone routes share identity and do not duplicate products.
    annotation.effect.orf_annotation = annotation
    assert graph(to_proteome([annotation.effect], orf_annotations=[annotation])) == graph(to_proteome(orf_annotations=[annotation]))
    other = annotate_orf_variants([variant], replace(record, source="another catalogue").bind(transcript))
    combined = to_proteome([annotation.effect, other.effect])
    assert len(combined.proteins) == 2
    assert {e.producer.provider for e in combined} == {"fixture", "another catalogue"}


def test_available_partial_unknown_empty_failure_and_limit_are_distinct(effect):
    baseline = ORFProduct(None, "MAK", True, "observed_stop", None, None)
    values = [
        ORFOutcome(effect, baseline, ORFProduct(None, "MTK", True, "sequence_boundary", None, None)),
        ORFOutcome(effect, baseline, ORFProduct(None, "MAK", None, "unavailable", None, None)),
        ORFOutcome(effect, baseline, ORFProduct(None, "", False, "unavailable", None, None)),
        ORFOutcome(effect, baseline, None, unresolved_reason="no mapped RNA"),
    ]
    # Distinct outcome contexts are distinct native inputs.
    values = [replace(value, hypothesis=EffectHypothesis(enumeration_index=i)) for i, value in enumerate(values)]
    values += [Failure(effect.variant, effect.transcript, error=ValueError("fixture failure")),
               Unresolved(effect.variant, effect.transcript, "missing_dna", "no introns"),
               HypothesisLimit(effect.variant, effect.transcript, 2, PhasePartition((), (), ()))]
    predicted = to_proteome(values)
    assert {r.status for r in comparisons(predicted)} == {
        "partial_product", "unchanged", "no_product", "sequence_unavailable",
        "annotation_failure", "unresolved", "hypothesis_limit"}
    assert {e.sequence for e in predicted} == {"MTK", "MAK"}
    assert next(r for r in comparisons(predicted) if r.status == "hypothesis_limit").limits["max_hypotheses"] == 2


def test_fusion_preserves_two_sources_without_parent_entries(effect):
    first = effect.transcript
    second = _ReferenceTranscript({**first.data, "id": "tx-2", "transcript_id": "tx-2"},
                                  first.genome, first.exons, _Feature(1, 12, "gene-2", "GENE2"))
    rna = MutantTranscript(reference_segments=(ReferenceSegment(first, 0, 6, label="5p"),
                                              ReferenceSegment(second, 6, 12, label="3p")),
                           cdna_sequence="ATGGCTAAATAA", mutant_protein_sequence="MAK",
                           annotator_name="rna", evidence={"read_ids": ["read-1"], "support": "RNA"})
    sv = StructuralVariant("1", 6, "BND", mate_contig="1", mate_start=7, genome=first.genome)
    fusion = GeneFusion(sv, first, second, mutant_transcript=rna)
    predicted = to_proteome([EffectCandidate(fusion, source="RNA assembler", evidence={"reads": 1})])
    entry, = predicted
    assert entry.producer.transcript_ids == ("tx-1", "tx-2")
    assert entry.producer.gene_ids == ("gene-1", "gene-2")
    assert entry.producer.sample_id is None
    assert entry.producer.tool == "unknown"
    assert entry.producer.provider == "RNA assembler"
    archive = json.dumps(graph(predicted))
    assert "5p" in archive and "3p" in archive and "read-1" in archive
    assert entry.derivation is None


def test_stable_order_idempotence_and_conflicting_source_content(effect):
    other = Substitution(Variant("1", 7, "A", "G", genome=effect.variant.genome), effect.transcript, 2, "K", "R")
    first, second = to_proteome([effect, other]), to_proteome([other, effect])
    assert {e.accession for e in first} == {e.accession for e in second}
    assert graph(first.merge(second)) == graph(first)
    assert graph(to_proteome([effect, effect])) == graph(to_proteome([effect]))
    changed = Substitution(effect.variant, effect.transcript, 1, "A", "S")
    changed.annotator, changed.annotator_version = effect.annotator, effect.annotator_version
    with pytest.raises(ValueError, match="Conflicting producer record"):
        to_proteome([effect, changed])


def test_unknown_extension_fields_survive_and_nonportable_values_fail(effect, tmp_path):
    candidate = EffectCandidate(effect, source="external", evidence={"future_field": {"z": [1, None, False]}})
    predicted = to_proteome([candidate])
    comparison, = comparisons(predicted)
    assert comparison.to_dict()["payload"]["contexts"][0]["evidence"]["future_field"] == {"z": [1, None, False]}
    with pytest.raises((TypeError, ValueError)):
        to_proteome([EffectCandidate(effect, evidence={"unsupported": object()})])


def test_splice_alternatives_and_added_rna_context_survive(effect):
    other = Deletion(effect.variant, effect.transcript, 1, "A")
    unresolved = Unresolved(effect.variant, effect.transcript, "retained_intron", "DNA unavailable")
    splice = SpliceOutcomeSet(effect.variant, effect.transcript, candidates=(
        EffectCandidate(effect, source="varcode_splice", evidence={"rank": 0}),
        EffectCandidate(other, source="varcode_splice", evidence={"rank": 1}),
        EffectCandidate(unresolved, source="varcode_splice", evidence={"rank": 2})),
        disrupted_signal_class=SpliceDonor)
    splice._extra_candidates = (EffectCandidate(other, source="RNA", evidence={"reads": ["read1"]}),)
    predicted = to_proteome([splice])
    assert len(comparisons(predicted)) == 4
    assert {entry.sequence for entry in predicted} == {"MTK", "MK"}
    assert len(predicted.by_sequence()["MK"]) == 2
    assert next(r for r in comparisons(predicted) if r.status == "unresolved").payload["omission_reason"] == "DNA unavailable"
    assert any(r.payload["contexts"][-1].get("source") == "RNA" for r in comparisons(predicted))


def test_cryptic_positions_distinguish_otherwise_identical_comparisons(effect, tmp_path):
    signal = SpliceDonor(effect.variant, effect.transcript, effect.transcript.exons[0], 1)
    candidates = tuple(EffectCandidate(CrypticDonor(
        effect.variant, effect.transcript, signal, effect.transcript.exons[0],
        cryptic_genomic_position=position)) for position in (15, 18))
    splice = SpliceOutcomeSet(effect.variant, effect.transcript, candidates=candidates,
                              disrupted_signal_class=SpliceDonor)
    predicted = to_proteome(splice)
    records = comparisons(predicted)
    assert len(records) == 2
    assert {r.status for r in records} == {"sequence_unavailable"}
    assert {r.payload["effect_description"] for r in records} == {"cryptic-donor (predicted)"}
    assert {r.payload["native"]["nodes"][0][1]["cryptic_genomic_position"] for r in records} == {15, 18}
    reversed_splice = SpliceOutcomeSet(effect.variant, effect.transcript, candidates=candidates[::-1],
                                       disrupted_signal_class=SpliceDonor)
    assert {r.key for r in comparisons(to_proteome(reversed_splice))} == {r.key for r in records}
    path = tmp_path / "cryptic.json.gz"
    predicted.save(path)
    assert graph(orfdata.Proteome.load(path)) == graph(predicted)
    standalone = to_proteome(splice.candidates)
    assert len(comparisons(standalone)) == 2
    assert len({r.key for r in standalone.records if r.kind == "event"}) == 2
    assert graph(standalone.merge(to_proteome(splice.candidates[::-1]))) == graph(standalone)


def test_equal_fusion_sequences_retain_transcript_model_associations(effect, tmp_path):
    transcript = effect.transcript
    sv = StructuralVariant("1", 6, "BND", mate_contig="1", mate_start=7, genome=transcript.genome)
    models = [MutantTranscript(reference_transcript=transcript,
                               cdna_sequence="ATGACTAAATAA", mutant_protein_sequence="MTK",
                               evidence={"transcript_model_id": model}) for model in ("model-1", "model-2")]
    first, second = [GeneFusion(sv, transcript, transcript, mutant_transcript=model) for model in models]
    fusion = GeneFusion(sv, transcript, transcript, primary_effects=(first, second))
    predicted = to_proteome(fusion)
    assert len(comparisons(predicted)) == len(predicted.proteins) == 2
    assert len(predicted.by_sequence()["MTK"]) == 2
    assert len(predicted.fasta_index(deduplicate=True)) == 1
    assert len([link for link in predicted.links if link.relation == "produces"]) == 2
    for model in ("model-1", "model-2"):
        assert sum(model in json.dumps(r.to_dict()["payload"]["native"]) for r in comparisons(predicted)) == 1
    # Reordering sibling candidates cannot change outcome IDs.
    fusion._primary_effects = (second, first)
    assert {r.key for r in comparisons(to_proteome(fusion))} == {r.key for r in comparisons(predicted)}
    path = tmp_path / "fusion-models.json.gz"
    predicted.save(path)
    assert graph(orfdata.Proteome.load(path)) == graph(predicted)
    standalone = to_proteome([first, second])
    assert len(standalone.proteins) == len(comparisons(standalone)) == 2
    conflicting = GeneFusion(sv, transcript, transcript, mutant_transcript=replace(models[0], mutant_protein_sequence="MTR"))
    with pytest.raises(ValueError, match="Conflicting producer record"):
        to_proteome([first, conflicting])
    # from_sequence also stores cDNA in an assembled segment's provider.
    # Those sequence bytes must not turn changed output into a new source ID.
    observed = [GeneFusion(sv, transcript, transcript, mutant_transcript=MutantTranscript.from_sequence(
        sequence, reference_transcript=transcript, mutant_protein_sequence="MTK",
        evidence={"transcript_model_id": "same-observation"}))
        for sequence in ("ATGACTAAATAA", "ATGACCAAATAA")]
    with pytest.raises(ValueError, match="Conflicting producer record"):
        to_proteome(observed)


@pytest.mark.parametrize("mechanism,coordinates", [
    (CrypticAcceptor, ({"cryptic_genomic_position": 15}, {"cryptic_genomic_position": 18})),
    (IntronRetention, ({"retained_intron_start": 13, "retained_intron_end": 15},
                      {"retained_intron_start": 13, "retained_intron_end": 18})),
])
def test_standalone_mechanism_coordinates_do_not_hide_output_conflicts(effect, mechanism, coordinates):
    signal_class = SpliceAcceptor if mechanism is CrypticAcceptor else SpliceDonor
    signal = signal_class(effect.variant, effect.transcript, effect.transcript.exons[0], 1)
    extra = {"affected_exon": effect.transcript.exons[0]} if mechanism is CrypticAcceptor else {}

    def candidate(position, protein):
        return mechanism(effect.variant, effect.transcript, signal, **extra, **position,
                         mutant_protein_sequence=protein)

    first, second = [candidate(position, "MTK") for position in coordinates]
    predicted = to_proteome([first, second])
    assert len(predicted.proteins) == len(comparisons(predicted)) == 2
    assert len(predicted.by_sequence()["MTK"]) == 2
    with pytest.raises(ValueError, match="Conflicting producer record"):
        to_proteome([first, candidate(coordinates[0], "MTR")])


@pytest.mark.parametrize("in_cis", [False, True])
def test_single_germline_hypothesis_preserves_provenance(in_cis, tmp_path):
    from varcode import GermlineContext, predict_germline_aware_effect
    from varcode.annotators.registry import get_default_annotator
    from .test_germline_effects import _cftr, _somatic, _same_codon_germline, _StubPhaseResolver

    variant, transcript, background = _somatic(), _cftr(), _same_codon_germline()
    annotator = get_default_annotator()
    plain = annotator.annotate_on_transcript(variant, transcript)
    phased = predict_germline_aware_effect(
        variant, transcript, GermlineContext.from_variants([background]), annotator,
        phase_resolver=_StubPhaseResolver(in_cis_answer=in_cis))
    assert phased.germline_phase_state == "phased"
    assert phased.germline_variants_in_window == (background,)
    if not in_cis:
        assert phased.mutant_protein_sequence == plain.mutant_protein_sequence
    predicted = to_proteome(phased)
    assert graph(predicted) != graph(to_proteome(plain))
    assert len(predicted.merge(to_proteome(plain)).proteins) == 2
    native = type(phased).from_json(phased.to_json())
    assert native.germline_phase_state == phased.germline_phase_state
    assert native.germline_variants_in_window == phased.germline_variants_in_window
    assert graph(to_proteome(native)) == graph(predicted)
    assert len(EffectCollection([phased])) == 1  # graph serialization remains hashable
    archive = json.dumps(graph(predicted))
    assert "germline_phase_state" in archive and "germline_variants_in_window" in archive
    path = tmp_path / "single-phase.json.gz"
    predicted.save(path)
    assert graph(orfdata.Proteome.load(path)) == graph(predicted)


def test_native_fusion_roundtrip_preserves_source_identity(tmp_path):
    from .test_germline_effects import _cftr

    first = _cftr()
    second = first.genome.transcript_by_id("ENST00000357654")
    variant = StructuralVariant(first.contig, first.start, "BND", genome=first.genome,
                                mate_contig=second.contig, mate_start=second.start)
    model = MutantTranscript(reference_segments=(ReferenceSegment(first, 0, 6), ReferenceSegment(second, 6, 12)),
                             cdna_sequence="ATGACTAAATAA", mutant_protein_sequence="MTK",
                             evidence={"transcript_model_id": "fusion-model"})
    fusion = GeneFusion(variant, first, second, mutant_transcript=model)
    restored = GeneFusion.from_json(fusion.to_json())
    original, reloaded = to_proteome(fusion), to_proteome(restored)
    assert {r.key for r in original.records} == {r.key for r in reloaded.records}
    assert {p.accession for p in original} == {p.accession for p in reloaded}
    assert graph(original.merge(reloaded)) == graph(original)
    path = tmp_path / "native-fusion.json.gz"
    original.merge(reloaded).save(path)
    assert graph(orfdata.Proteome.load(path)) == graph(original)


def test_fusion_aliases_preserve_reference_scope_and_roles(effect):
    first = effect.transcript
    second = _ReferenceTranscript({**first.data, "id": "tx-2", "transcript_id": "tx-2"},
                                  first.genome, first.exons, _Feature(1, 12, "gene-2", "GENE2"))
    variant = StructuralVariant("1", 6, "BND", mate_contig="1", mate_start=7, genome=first.genome)
    forward = GeneFusion(variant, first, second)
    reversed_roles = GeneFusion(variant, first, second, five_prime_transcript=second, three_prime_transcript=first)
    assert len(comparisons(to_proteome([forward, reversed_roles]))) == 2
    models = [MutantTranscript(
        reference_segments=(ReferenceSegment(first, offset, offset + 6), ReferenceSegment(second, 6, 12)),
        mutant_protein_sequence="MTK", evidence={"support": "same-source"}) for offset in (0, 3)]
    alternatives = [GeneFusion(variant, first, second, mutant_transcript=model) for model in models]
    assert len(to_proteome(alternatives).proteins) == 2
    # Matching transcript IDs must not collapse genuinely different references.
    other_reference = _ReferenceGenome({**first.genome.metadata, "reference_name": "other-assembly"})
    alias = _ReferenceTranscript(dict(first.data), other_reference, first.exons, first.gene)
    with pytest.raises(ValueError, match="reference_name mismatch"):
        to_proteome(GeneFusion(variant, first, alias))


@pytest.mark.parametrize("coding", [False, True])
def test_dog_reference_uses_orfdata_taxonomy(effect, coding):
    from pyensembl import cached_release

    genome = cached_release(99, species="dog")
    assert genome.taxon_id == 9615 and genome.species_taxon_id == 9612
    effect.transcript.genome = effect.variant.genome = genome
    predicted = to_proteome(effect if coding else Intergenic(effect.variant))
    reference = orfdata.ReferenceContext.from_genome(genome)
    assert predicted.reference == reference
    assert comparisons(predicted)[0].payload["source_taxon_id"] == 9615
    if coding:
        assert predicted.proteins[0].taxon_id == 9615
        assert predicted.proteins[0].producer.reference == reference


def test_strain_reference_uses_orfdata_species_resolution(effect, tmp_path):
    genome = _ReferenceGenome(dict(reference_name="R64-1-1", annotation_name="synthetic",
                                  species="saccharomyces_cerevisiae", taxon_id=559292,
                                  species_taxon_id=4932, strain="S288C", strain_taxon_id=559292))
    effect.transcript.genome = effect.variant.genome = genome
    predicted = to_proteome(effect)
    assert predicted.proteins[0].taxon_id == 4932
    assert predicted.reference == orfdata.ReferenceContext.from_genome(genome)
    assert predicted.reference.strain_taxon_id == 559292
    path = tmp_path / "strain.json.gz"
    predicted.save(path)
    assert graph(orfdata.Proteome.load(path)) == graph(predicted)


def test_partial_source_context_is_validated_against_later_reference(effect):
    effect.transcript.genome = _ReferenceGenome(dict(reference_name="different-assembly"))
    with pytest.raises(ValueError, match="reference_name mismatch"):
        to_proteome(effect)


def test_reference_and_isovar_shaped_composition_and_conflicts(effect, tmp_path):
    predicted = to_proteome([effect], sample_id="sample1")
    reference = predicted.reference
    assert reference is not None
    # A reference baseline can be added later; the prediction never required
    # it to be present. This is the same protein API used by from_genome.
    baseline = orfdata.ProteinSequence("protein-1", "MAK", "annotated", 9606, "synthetic")
    reference_proteome = orfdata.Proteome(reference, [baseline])
    producer = orfdata.ProducerInfo(
        origin="reference_cds", method="rna_assembly", tool="isovar", native_id="rna-product-1",
        sample_id="sample1", transcript_ids=("tx-1",), reference=reference,
        start_complete=None, stop_complete=None)
    rna_entry = orfdata.ProteinSequence.from_producer("MTK", taxon_id=9606, producer=producer)
    rna_record = orfdata.ProducerRecord(
        namespace="isovar", schema_version="isovar.fixture.v1", kind="evidence", native_id="reads-1",
        support="RNA", scope={"sample_id": "sample1", "reads": ["read1"]},
        payload={"translation_support": "not_measured", "future_extension": {"a": 1}})
    rna = orfdata.Proteome.from_proteins([rna_entry], reference=reference, records=[rna_record], links=[
        orfdata.ProducerLink(rna_record.key, rna_entry.accession, "protein", "supports_RNA_template")])
    merged = reference_proteome.merge(predicted, rna)
    assert len(merged) == 3 and len(merged.by_sequence()["MTK"]) == 2
    assert merged.merge(predicted, rna).identity == merged.identity
    path = tmp_path / "composed.json.gz"
    merged.save(path)
    assert graph(orfdata.Proteome.load(path)) == graph(merged)
    with pytest.raises(ValueError, match="different reference"):
        predicted.merge(orfdata.Proteome(replace(reference, annotation_version="2"), []))
    bad = replace(rna_entry, sequence="MKK")
    with pytest.raises(ValueError, match="Conflicting protein identity"):
        merged.merge(orfdata.Proteome.from_proteins([bad], reference=reference))
    bad_record = replace(rna_record, payload={"conflicting": True})
    with pytest.raises(ValueError, match="Conflicting producer record"):
        merged.merge(orfdata.Proteome.from_proteins([], reference=reference, records=[bad_record]))


def test_different_known_source_taxa_and_assemblies_are_rejected(effect):
    # A bad pairing must be caught before assigning the reference's taxon.
    effect.variant.genome = _ReferenceGenome(dict(reference_name="dog-reference", species="canis_lupus_familiaris",
                                                taxon_id=9615, annotation_name="synthetic", annotation_version="1"))
    with pytest.raises(ValueError, match="taxon|species"):
        to_proteome([effect])
    effect.variant.genome = _ReferenceGenome({**effect.transcript.genome.metadata, "reference_name": "other"})
    with pytest.raises(ValueError, match="reference_name mismatch"):
        to_proteome([effect])


def test_missing_gene_names_and_terminal_stop_are_preserved(effect):
    effect.transcript.gene = None
    effect.gene = None
    predicted = to_proteome([outcome(effect, mutant="MTK*")])
    entry, = predicted
    assert entry.sequence == "MTK*"
    assert entry.producer.gene_ids == entry.producer.gene_names == ()
    assert "aa=4" in entry.fasta_token()
    assert entry.producer.stop_complete is None


def test_wrapper_export_does_not_split_native_source_identity(effect):
    predicted = to_proteome([effect])
    entry, = predicted
    wrapper_entry = orfdata.ProteinSequence.from_producer(
        entry.sequence, taxon_id=entry.taxon_id,
        producer=replace(entry.producer, exporter=None, exporters=("wrapper",)))
    wrapped = orfdata.Proteome.from_proteins([wrapper_entry], reference=predicted.reference)
    assert wrapper_entry.accession == entry.accession
    combined = predicted.merge(wrapped)
    assert len(combined) == 1
    assert "wrapper" in json.dumps(graph(combined)) and "varcode" in json.dumps(graph(combined))
    assert combined.merge(predicted, wrapped).identity == combined.identity


def test_optional_dependency_does_not_load_at_import():
    subprocess.run([sys.executable, "-c", "import sys; import varcode; assert 'orfdata' not in sys.modules"], check=True)


def test_collection_version_is_not_assigned_to_a_different_producer(effect):
    effect.annotator = "external_annotator"
    effect.annotator_version = None
    entry, = to_proteome(EffectCollection([effect], annotator="fast", annotator_version="99"))
    assert entry.producer.tool == "unknown"
    assert entry.producer.tool_version is None
    assert entry.producer.method == "external_annotator"


def test_collection_provenance_is_retained_without_splitting_products(effect):
    first = EffectCollection([effect], sources={"a.vcf", "b.vcf"}, annotated_at="2026-01-01T00:00:00Z")
    second = EffectCollection([effect], sources={"c.vcf"}, annotated_at="2026-01-02T00:00:00Z")
    predicted = to_proteome(first).merge(to_proteome(second))
    assert len(predicted) == 1
    collections = [record for record in predicted.records if record.kind == "collection"]
    assert {tuple(record.payload["sources"]) for record in collections} == {("a.vcf", "b.vcf"), ("c.vcf",)}
    assert {record.payload["annotated_at"] for record in collections} == {first.annotated_at, second.annotated_at}
    assert len([link for link in predicted.links if link.relation == "contains"]) == 2


def test_unknown_taxon_is_not_guessed(effect):
    unknown = _ReferenceGenome(dict(reference_name="fixture", annotation_name="synthetic", annotation_version="1"))
    effect.transcript.genome = unknown
    effect.variant.genome = unknown
    predicted = to_proteome([effect])
    entry, = predicted
    assert entry.taxon_id is None and entry.producer.reference is None
    assert predicted.reference is None


def test_fasta_preserves_available_producer_metadata(effect):
    entry, = to_proteome([effect])
    header = entry.fasta_header()
    for value in ("fixture-1", "transcript_model", "GENE1", "tx-1"):
        assert value in header


def test_native_ensembl_effect_collection_with_extra_orf_roundtrips(tmp_path):
    from pyensembl import cached_release

    genome = cached_release(95)
    record = orfdata.load_catalog().get("c7riboseqorf56")
    reference = record.bind(genome.transcript_by_id(record.transcript_id))
    variant = Variant("7", 66205486, "C", "G", genome=genome)
    effects = variant.effects(extra_orfs=[reference])
    predicted = effects.to_proteome(sample_id="fixture")
    assert graph(predicted) == graph(to_proteome(effects, sample_id="fixture"))
    orf_entry, = [entry for entry in predicted if entry.producer.orf_ids]
    assert orf_entry.sequence == effects.orf_annotations[0].outcomes[0].mutant.protein_sequence
    assert orf_entry.sequence.startswith("MASRRR")
    assert orf_entry.producer.origin == "nc_orf"
    assert len([r for r in predicted.records if r.kind == "event"]) == len(effects)
    path = tmp_path / "native.json.gz"
    predicted.save(path)
    assert graph(orfdata.Proteome.load(path)) == graph(predicted)
    # Native ORF archives snapshot each gene through its parent transcript;
    # the absence of a second gene.genome pointer must not change source IDs.
    extra = next(e for e in effects if e.orf_annotation is not None)
    archived = type(extra).from_json(extra.to_json())
    assert {e.accession for e in to_proteome(extra)} == {e.accession for e in to_proteome(archived)}


def test_exporter_version_does_not_change_native_prediction_identity(effect, monkeypatch):
    import varcode.proteome as adapter

    first = to_proteome([effect])
    monkeypatch.setattr(adapter, "_EXPORTER", "varcode/future")
    second = to_proteome([effect])
    # Both source entries and graph records keep identity independently of
    # which adapter version serialized the same native result (ORFdata #124).
    assert first.proteins[0].producer.native_id == second.proteins[0].producer.native_id
    assert first.proteins[0].accession == second.proteins[0].accession
    assert {r.key for r in first.records} == {r.key for r in second.records}
    assert len(first.merge(second).proteins) == 1


def test_default_germline_phase_candidates_export_and_archive(tmp_path):
    from varcode import GermlineContext
    from varcode.effects.structural_serialization import structural_effect_from_dict, structural_effect_to_dict
    from .test_germline_effects import _somatic, _same_codon_germline, CFTR_ID

    variant = _somatic()
    germline = GermlineContext.from_variants([_same_codon_germline()])
    effect, = [e for e in variant.effects(germline=germline) if e.transcript_id == CFTR_ID]
    assert isinstance(effect, PhaseCandidateSet)
    predicted = to_proteome([effect])
    assert {r.payload["effect_description"] for r in comparisons(predicted)} == {"p.L159M", "p.S159T"}
    state = structural_effect_to_dict(effect, schema="orf_effect_schema")
    restored = structural_effect_from_dict(PhaseCandidateSet, state, schema="orf_effect_schema")
    assert [c.short_description for c in restored._candidates_raw] == [c.short_description for c in effect._candidates_raw]
    assert [c.evidence for c in restored.candidates] == [c.evidence for c in effect.candidates]
    path = tmp_path / "phase.json.gz"
    predicted.save(path)
    assert graph(orfdata.Proteome.load(path)) == graph(predicted)


def test_transcript_free_gene_effects_have_distinct_identities(effect, tmp_path):
    first = Intragenic(effect.variant, _Feature(1, 12, "gene-1", "GENE1", "1", "+"))
    second = Intragenic(effect.variant, _Feature(1, 12, "gene-2", "GENE2", "1", "+"))
    predicted = to_proteome([first, second])
    assert len(predicted) == 0 and len(comparisons(predicted)) == 2
    assert len({r.key for r in predicted.records if r.kind == "event"}) == 2
    assert {r.payload["gene_ids"] for r in comparisons(predicted)} == {("gene-1",), ("gene-2",)}
    assert to_proteome([second, first]).merge(predicted).summary()["producer_records"] == 4
    path = tmp_path / "intragenic.json.gz"
    predicted.save(path)
    assert graph(orfdata.Proteome.load(path)) == graph(predicted)


def test_bound_orf_hypothesis_limit_preserves_status_and_cap(reference_orf, effect):
    partner = Variant("1", 5, "C", "G", genome=effect.variant.genome)
    annotation = annotate_orf_variants([effect.variant], reference_orf, germline_variants=[partner], max_hypotheses=1)
    assert isinstance(annotation.effect, HypothesisLimit) and annotation.outcomes == ()
    predicted = to_proteome(orf_annotations=[annotation])
    assert not predicted.proteins
    for record in predicted.records:
        assert record.status == "hypothesis_limit"
        assert record.complete is False
        assert record.limits["max_hypotheses"] == 1
    comparison, = comparisons(predicted)
    assert comparison.payload["omission_reason"] == annotation.effect.reason


@pytest.mark.parametrize("kind", ["failure", "unresolved"])
def test_bound_orf_without_products_retains_underlying_reason(reference_orf, effect, kind):
    native = (Failure(effect.variant, effect.transcript, error="native annotation failed") if kind == "failure" else
              Unresolved(effect.variant, effect.transcript, "missing_sequence", "reference DNA absent"))
    annotation = ORFAnnotation(reference_orf.record, native, reference_orf, variants=(effect.variant,))
    comparison, = comparisons(to_proteome(orf_annotations=[annotation]))
    assert comparison.status == ("annotation_failure" if kind == "failure" else "unresolved")
    assert comparison.payload["omission_reason"] == (native.error if kind == "failure" else native.reason)


@pytest.mark.parametrize("standalone", [False, True])
def test_orf_outcome_sources_remain_distinct(reference_orf, effect, standalone, tmp_path):
    annotation = annotate_orf_variants([effect.variant], reference_orf)
    outcomes = tuple(replace(annotation.outcomes[0], source=source) for source in ("source-A", "source-B"))
    annotated = replace(annotation, outcomes=outcomes)
    predicted = to_proteome(outcomes) if standalone else to_proteome(orf_annotations=[annotated])
    assert len(predicted) == len(comparisons(predicted)) == 2
    assert len(predicted.by_sequence()["MTK"]) == 2
    assert {r.payload["contexts"][-1]["source"] for r in comparisons(predicted)} == {"source-A", "source-B"}
    assert {e.producer.provider for e in predicted} == {"source-A", "source-B"}
    assert len(predicted.merge(predicted)) == 2
    path = tmp_path / "sources.json.gz"
    predicted.save(path)
    assert graph(orfdata.Proteome.load(path)) == graph(predicted)


def test_single_and_multiple_orf_results_use_the_same_export_api(reference_orf, effect):
    annotation = annotate_orf_variants([effect.variant], reference_orf)
    failed = ORFAnnotation(reference_orf.record, None, variants=(effect.variant,), unresolved_reason="unavailable model")
    for result in (annotation, failed):
        assert graph(to_proteome(result)) == graph(to_proteome(orf_annotations=[result]))
    assert graph(to_proteome([annotation, failed])) == graph(to_proteome(orf_annotations=[annotation, failed]))
    assert graph(to_proteome(effect)) == graph(to_proteome([effect]))


def test_ncorf_selector_and_export_delegate_to_orfdata(monkeypatch, tmp_path):
    from pyensembl import cached_release

    genome = cached_release(95)
    catalog = orfdata.load_catalog()
    reference = catalog.get("c7riboseqorf56").bind(genome.transcript_by_id("ENST00000304842"))
    calls = []

    class PlacedORFs:
        def __iter__(self):
            return iter([reference])

    placed = PlacedORFs()
    placed.reference, placed.catalog, placed.excluded = reference.context, catalog, ()

    def load(genome, evidence):
        calls.append((genome, evidence))
        return placed

    monkeypatch.setattr(orfdata, "extra_orfs", load)
    monkeypatch.setattr("varcode.extra_orfs._RESOLVERS", {})
    variant = Variant("7", 66205486, "C", "G", genome=genome)
    first = variant.effects(extra_orfs="all")
    second = variant.effects(extra_orfs="all")
    assert calls == [(genome, None)]
    assert len(first.orf_annotations) == len(second.orf_annotations) == 1
    proteome = first.to_proteome(sample_id="tumor-1")
    assert isinstance(proteome, orfdata.Proteome)
    entry, = [e for e in proteome if e.producer.orf_ids]
    assert entry.producer.orf_ids == (reference.orf_id,)
    assert entry.sequence == first.orf_annotations[0].outcomes[0].mutant.protein_sequence
