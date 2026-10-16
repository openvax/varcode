"""Upstream-frame consequences and native evidence remain independent."""

from dataclasses import replace
import json
from types import SimpleNamespace

import pytest
from pyensembl import cached_release

from varcode import Variant, ReferenceMismatchError
from varcode.uorf import (
    ORFTranscript, annotate_uorf_variants, annotate_uorf_product, load_uorf_catalog,
)
from varcode.effect_hypotheses import RealizedTranscriptProduct
from varcode.nucleotides import reverse_complement

orfdata = pytest.importorskip("orfdata")


@pytest.fixture(scope="module")
def catalog():
    return load_uorf_catalog()


@pytest.fixture(scope="module")
def native(catalog):
    record = catalog.get("c7riboseqorf56")
    transcript = cached_release(95).transcript_by_id(record.transcript_id)
    return record, transcript


def annotation(native, position, ref, alt, **kwargs):
    record, transcript = native
    variant = Variant("7", position, ref, alt, genome=transcript.genome)
    return annotate_uorf_variants([variant], record, transcript, genome_build="GRCh38", **kwargs)


def test_tpst1_plus_four_retains_context_and_coding_consequences(native, catalog):
    result = annotation(native, 66205486, "C", "G")
    assert type(result.effect).__name__ == "Substitution"
    assert result.effect.aa_ref == "P" and result.effect.aa_alt == "A"
    data = json.loads(json.dumps(result.to_dict()))
    context_snv, = data["variant_locations"][0]["initiation_snvs"]
    assert context_snv["context_position"] == 4
    assert context_snv["sequence_changes"] == ["kozak_preference_gained"]
    assert "body" in context_snv["roles"]
    outcome, = data["outcomes"]
    assert outcome["baseline"]["initiation_context"]["downstream"].startswith("C")
    assert outcome["mutant"]["initiation_context"]["downstream"].startswith("G")
    assert outcome["mutant"]["reference_evidence"]["reference_ms_support_on_retained_sequence"]
    assert outcome["mutant"]["reference_evidence"]["mutant_translation_evidence"] == "not_assessed"
    variant = result.effect.variant
    assert any(type(e).__name__ == "FivePrimeUTR" for e in variant.effects())
    queried = variant.orf_effects(catalog=catalog)
    assert any(a.orf_id == result.orf_id and type(a.effect).__name__ == "Substitution" for a in queried)


def test_kozak_snv_does_not_claim_coding_or_measured_translation(native):
    result = annotation(native, 66205480, "A", "C")
    assert len(result.outcomes) == 1
    outcome = result.outcomes[0]
    assert outcome.baseline.protein_sequence == outcome.mutant.protein_sequence
    assert outcome.baseline.initiation_context.upstream != outcome.mutant.initiation_context.upstream


def test_ref_mismatch_and_assembly_fail(native):
    with pytest.raises(ReferenceMismatchError):
        annotation(native, 66205486, "A", "G")
    record, transcript = native
    with pytest.raises(ValueError, match="assembly"):
        ORFTranscript(record, transcript, genome_build="GRCh37")
    with pytest.raises(ValueError, match="protein"):
        ORFTranscript(replace(record, protein_sequence="A" + record.protein_sequence[1:]), transcript, genome_build="GRCh38")
    with pytest.raises(ValueError, match="path"):
        ORFTranscript(replace(record, coding_blocks=((66205482, 66205521), (66240325, 66240404))), transcript, genome_build="GRCh38")


class Transcript:
    """Tiny independently specified transcript, with optional spliced start."""
    id = "test"
    name = "test"
    contig = "1"
    genome = SimpleNamespace(reference_name="GRCh38", release=1,
                             species="homo_sapiens", annotation_name="synthetic",
                             annotation_version="1")
    gene = SimpleNamespace(id="gene", name="GENE")
    is_protein_coding = True
    complete = True

    def __init__(self, reverse=False, split=False):
        self.sequence = "CCCATGAAACCCGGGTAAATGAAATAA"
        self.on_backward_strand = reverse
        self.strand = "-" if reverse else "+"
        positions = list(range(100, 100 + len(self.sequence)))
        if split:
            positions = positions[:5] + [p + 100 for p in positions[5:]]
        if reverse:
            positions = [400 - p for p in positions]
        self.positions = positions
        groups = []
        for p in sorted(positions):
            if groups and groups[-1][1] + 1 == p:
                groups[-1][1] = p
            else:
                groups.append([p, p])
        self.exons = tuple(SimpleNamespace(start=a, end=b) for a, b in groups)
        self.start, self.end = min(positions), max(positions)

    def spliced_offset(self, position):
        return self.positions.index(position)

    def record(self):
        coding = sorted(self.positions[3:18])
        blocks = []
        for p in coding:
            if blocks and blocks[-1][1] + 1 == p:
                blocks[-1][1] = p
            else:
                blocks.append([p, p])
        return orfdata.UORFRecord("test-orf", "gene", "GENE", "test", "1", self.strand,
                                  tuple(map(tuple, blocks)), "MKPG", "uORF", (), ())

    def variant(self, offset, ref, alt):
        position = self.positions[offset]
        if self.on_backward_strand:
            position -= len(ref) - 1
            ref, alt = reverse_complement(ref), reverse_complement(alt)
        return Variant("1", position, ref, alt, genome=cached_release(95))


@pytest.mark.parametrize("taxon,species,assembly", [
    (9606, "homo_sapiens", "GRCh38"),
    (10090, "mus_musculus", "GRCm39"),
    (9615, "canis_lupus_familiaris", "ROS_Cfam_1.0"),
    (4932, "saccharomyces_cerevisiae", "R64-1-1"),
])
def test_modified_products_keep_reference_scope_across_species(taxon, species, assembly):
    from dataclasses import replace
    from pyensembl.species import find_species_by_name
    from varcode.effect_hypotheses import RealizedTranscriptProduct
    from varcode.orf import annotate_orf_product

    t = Transcript()
    t.genome = SimpleNamespace(reference_name=assembly, species=find_species_by_name(species),
                              annotation_name="synthetic", annotation_version="1")
    record = replace(t.record(), taxon_id=taxon, species=species, genome_build=assembly)
    reference = record.bind(t)
    product = RealizedTranscriptProduct(t, t.sequence, None, evidence={
        "base_origins": [(t.contig, p, t.strand) for p in t.positions],
    })
    result = annotate_orf_product(reference, product)
    assert result["reference"]["reference_key"] == reference.context.key
    assert result["reference"]["reference_context"]["taxon_id"] == taxon
    assert result["product"]["protein_sequence"] == "MKPG"
    assert result["product"]["reference_evidence"]["mutant_translation_evidence"] == "not_assessed"
    t.genome.annotation_version = "2"
    with pytest.raises(ValueError, match="annotation_version"):
        annotate_orf_product(reference, product)


def test_legacy_binding_can_declare_custom_reference_context():
    from orfdata import ReferenceContext
    from varcode.uorf import ORFTranscript

    t = Transcript()
    t.genome = SimpleNamespace(reference_name="GRCh38")
    with pytest.raises(ValueError, match="explicit ReferenceContext"):
        ORFTranscript(t.record(), t, genome_build="GRCh38")
    explicit = ReferenceContext(9606, "human", "GRCh38", "synthetic", "1")
    view = ORFTranscript(t.record(), t, genome_build="GRCh38", reference_context=explicit)
    assert view.reference_orf.context == explicit


def test_variant_cannot_cross_taxa_even_with_a_colliding_assembly_label():
    from dataclasses import replace

    t = Transcript()
    t.genome = SimpleNamespace(reference_name="GRCm38", species="human",
                              annotation_name="synthetic", annotation_version="1")
    record = replace(t.record(), genome_build="GRCm38")
    variant = Variant("1", 108, "A", "G", genome=cached_release(95, species="mouse"))
    with pytest.raises(ValueError, match="taxon mismatch"):
        annotate_uorf_variants([variant], record, t, genome_build="GRCm38")


def test_versioned_source_transcript_lookup_and_duplicate_rejection():
    from dataclasses import replace
    from varcode.orf import predict_orf_effects

    t = Transcript()
    t.versioned_id = "test.2"
    record = replace(t.record(), transcript_id="test.2")
    catalog = orfdata.ORFCatalog([record], dict(genome_build="GRCh38", taxon_id=9606,
                                             species="human"))
    variant = Variant("1", 108, "A", "G", genome=cached_release(95))
    annotation, = predict_orf_effects(variant, catalog=catalog, transcripts=[t])
    assert annotation.reference is not None
    assert annotation.reference.model.transcript_id == "test.2"
    with pytest.raises(ValueError, match="Ambiguous supplied transcript"):
        predict_orf_effects(variant, catalog=catalog, transcripts=[t, t])
    t.versioned_id = "test.3"
    annotation, = predict_orf_effects(variant, catalog=catalog, transcripts=[t])
    assert annotation.data["effect"]["type"] == "Unresolved"
    assert "version mismatch" in annotation.data["effect"]["reason"]
    assert "version mismatch" in annotation.unresolved_reason
    assert annotation.outcomes == ()
    assert annotation.data == annotation.to_dict()
    assert annotation.effect is None and annotation.reference is None
    missing, = predict_orf_effects(variant, catalog=catalog, transcripts=[])
    assert missing.effect is None and missing.reference is None
    assert missing.to_dict()["model"]["transcript_id"] == "test.2"
    assert missing.to_dict()["effect"]["type"] == "Unresolved"


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("split", [False, True])
def test_spliced_and_minus_strand_start_context(reverse, split):
    t = Transcript(reverse, split)
    view = ORFTranscript(t.record(), t, genome_build="GRCh38")
    assert view.start_codon_positions == tuple(t.positions[3:6])
    a = annotate_uorf_variants([t.variant(6, "A", "G")], t.record(), t, genome_build="GRCh38")
    outcomes = a.outcomes
    assert any(x.mutant.protein_sequence == "MEPG" for x in outcomes)


@pytest.mark.parametrize("offset,ref,alt,effect,protein", [
    (3, "A", "C", "AlternateStartCodon", "MKPG"),
    (4, "T", "C", "StartLoss", ""),
    (8, "A", "G", "Silent", "MKPG"),
    (6, "A", "T", "PrematureStop", "M"),
    (15, "T", "C", "StopLoss", "MKPGQMK"),
    (9, "CCC", "", "Deletion", "MKG"),
    (9, "C", "", "FrameShift", "MKPGK"),
])
def test_sequence_consequences(offset, ref, alt, effect, protein):
    t = Transcript()
    a = annotate_uorf_variants([t.variant(offset, ref, alt)], t.record(), t, genome_build="GRCh38")
    assert type(a.effect).__name__ == effect
    assert len(a.outcomes) == 1
    outcome = a.outcomes[0]
    assert outcome.mutant.protein_sequence == protein
    if effect == "FrameShift":
        assert outcome.mutant.termination == "observed_stop"
        assert a.effect.aa_mutation_start_offset < 4


def test_known_cis_and_unknown_germline_phase_export_all_outcomes():
    t = Transcript()
    v = t.variant(6, "A", "G")
    partner = t.variant(7, "A", "C")
    a = annotate_uorf_variants([v, partner], t.record(), t, genome_build="GRCh38")
    assert a.outcomes[0].mutant.protein_sequence == "MAPG"
    uncertain = annotate_uorf_variants([v], t.record(), t, genome_build="GRCh38", germline_variants=[partner])
    outcomes = uncertain.outcomes
    assert {o.mutant.protein_sequence for o in outcomes} == {"MAPG", "MEPG"}
    assert {tuple(state for _, state in o.hypothesis.phase) for o in outcomes} == {("cis",), ("trans",)}


def test_observed_product_missing_start_remains_unresolved(native):
    record, transcript = native
    product = RealizedTranscriptProduct(transcript, "ATGAAATAA", None,
        evidence={"base_origins": [None] * 9, "rna_support": 10})
    result = annotate_uorf_product(record, transcript, product, genome_build="GRCh38")
    assert result["status"] == "unresolved"
    assert result["evidence"]["rna_support"] == 10


def test_sid_dna_fusion_keeps_alternative_acceptors_and_native_ms_scope(catalog):
    from varcode.sv_allele_parser import parse_symbolic_alt
    variant = parse_symbolic_alt("7", 66205522, "G", "G[7:66127704[", genome=cached_release(95))
    result, = variant.orf_effects(catalog=catalog)
    assert type(result.effect).__name__ == "GeneFusion"
    outcomes = result.outcomes
    assert any(o.mutant.protein_sequence == "MPSRRRGGSSLIPGRMPILRFSVTTRYFSY"
               for o in outcomes if o.mutant)
    assert len({o.mutant.protein_sequence for o in outcomes if o.mutant}) >= 3
    for outcome in outcomes:
        if outcome.mutant:
            assert not outcome.mutant.reference_evidence["reference_ms_support_on_retained_sequence"]
    assert json.loads(json.dumps(result.to_dict()))["model"]["orf_id"] == "c7riboseqorf56"


def test_stop_loss_with_no_downstream_sequence_marks_boundary():
    t = Transcript()
    t.sequence = t.sequence[:18]
    t.positions = t.positions[:18]
    t.end = t.positions[-1]
    t.exons = (SimpleNamespace(start=t.start, end=t.end),)
    a = annotate_uorf_variants([t.variant(15, "T", "C")], t.record(), t, genome_build="GRCh38")
    assert len(a.outcomes) == 1
    outcome = a.outcomes[0]
    assert outcome.mutant.protein_sequence == "MKPGQ"
    assert outcome.mutant.termination == "sequence_boundary"


def test_padded_insertion_ref_and_missing_intronic_ref_are_not_ignored():
    t = Transcript()
    with pytest.raises(ReferenceMismatchError):
        annotate_uorf_variants([t.variant(8, "C", "CA")], t.record(), t, genome_build="GRCh38")
    spliced = Transcript(split=True)
    v = Variant("1", 150, "A", "G", genome=cached_release(95))
    a = annotate_uorf_variants([v], spliced.record(), spliced, genome_build="GRCh38")
    assert type(a.effect).__name__ == "Unresolved"


def test_generic_orf_annotation_keeps_noncoding_mouse_parent():
    from varcode import annotate_orf_variants
    from orfdata import ORFRecord, ORFEvidence
    t = Transcript()
    t.biotype = "lncRNA"
    t.is_protein_coding = False
    t.genome = cached_release(95, species="mouse")
    model = replace(t.record().model, taxon_id=10090, species="mus_musculus",
                    genome_build="GRCm38", orf_type="lncRNA", source="synthetic-test")
    record = ORFRecord.from_model(model, ORFEvidence())
    reference = record.bind(t)
    variant = Variant("1", 106, "A", "G", genome=t.genome)
    result = annotate_orf_variants([variant], reference)
    assert result.reference is reference
    assert result.reference.transcript is t and not t.is_protein_coding
    assert result.effect.aa_ref == "K" and result.effect.aa_alt == "E"
    data = result.to_dict()
    assert data["reference"]["model"]["taxon_id"] == 10090
    assert data["reference"]["resolved_transcript"]["biotype"] == "lncRNA"
    assert data["outcomes"][0]["mutant"]["protein_sequence"] == "MEPG"
    assert data["outcomes"][0]["mutant"]["reference_evidence"]["mutant_translation_evidence"] == "not_assessed"


def test_reference_peptide_retention_preserves_all_duplicated_copies():
    from varcode.orf import _peptide_overlay, _model_positions
    t = Transcript()
    peptide = orfdata.PeptideEvidence("proteomics", "KP", ((2, 3),), "synthetic", "test")
    record = replace(t.record(), peptides=(peptide,))
    native = _model_positions(record)
    # Duplicate the native KP codons in-frame; both copies retain their ancestry.
    origins = native[:9] + native[3:9] + native[9:]
    overlay = _peptide_overlay(record, "MKPKPG", origins, 0)
    observation, = overlay["peptide_observations"]
    assert observation["retained_mutant_intervals"] == [(2, 3), (4, 5)]
    assert overlay["mutant_translation_evidence"] == "not_assessed"


def test_primary_source_lncRNA_orf_uses_same_variant_api():
    from varcode import annotate_orf_variants
    record = orfdata.load_catalog().get("c10norep110")
    transcript = cached_release(95).transcript_by_id(record.transcript_id)
    reference = record.bind(transcript)
    assert transcript.biotype == "antisense" and not transcript.is_protein_coding
    assert record.gene_name == "NUTM2B-AS1" and record.has_evidence("ms")
    position = reference.model.genomic_positions[3]
    ref = reference.coding_sequence[3]
    if record.strand == "-":
        ref = reverse_complement(ref)
    alt = next(b for b in "ACGT" if b != ref)
    variant = Variant(record.contig, position, ref, alt, genome=transcript.genome)
    annotation = annotate_orf_variants([variant], reference)
    assert type(annotation.effect).__name__ == "PrematureStop"
    assert annotation.reference.transcript is transcript
    assert annotation.outcomes[0].mutant.protein_sequence == "M"


def test_native_orf_access_does_not_serialize(monkeypatch):
    from varcode import ORFAnnotation, ORFOutcome, ORFProduct, annotate_orf_variants
    from varcode.effect_hypotheses import EffectHypothesis
    from varcode.effects import Substitution

    def unexpected_export(*args, **kwargs):
        raise AssertionError("Annotation access must not serialize")

    t = Transcript()
    reference = t.record().bind(t)
    variant = t.variant(6, "A", "G")
    monkeypatch.setattr(ORFAnnotation, "to_dict", unexpected_export)
    monkeypatch.setattr(orfdata.ReferenceORF, "to_dict", unexpected_export)
    monkeypatch.setattr("varcode.orf._plain", unexpected_export)
    result = annotate_orf_variants([variant], reference)
    assert result.reference is reference
    assert result.variants == (variant,)
    assert len(result.outcomes) == 1
    outcome = result.outcomes[0]
    assert isinstance(outcome, ORFOutcome)
    assert isinstance(outcome.effect, Substitution)
    assert isinstance(outcome.hypothesis, EffectHypothesis)
    assert isinstance(outcome.mutant, ORFProduct)
    assert outcome.baseline.protein_sequence == "MKPG"
    assert outcome.mutant.protein_sequence == "MEPG"
    assert outcome.mutant.initiation_context.start_codon == "ATG"


def test_multiple_orfs_are_iterable_and_named_selection_is_explicit():
    from varcode import annotate_orf_variants

    t = Transcript()
    first = t.record()
    second = replace(first, orf_id="second-orf", coding_blocks=((118, 126),),
                     protein_sequence="MK")
    catalog = orfdata.ORFCatalog([first, second], dict(
        genome_build="GRCh38", taxon_id=9606, species="human"))
    annotations = [annotate_orf_variants([t.variant(6, "A", "G")], reference)
                   for reference in catalog.for_transcript(t)]
    assert {a.orf_id for a in annotations} == {"test-orf", "second-orf"}
    selected = catalog.get("second-orf").bind(t)
    assert selected.model.orf_id == "second-orf"
    assert selected.protein_sequence == "MK"
    with pytest.raises(KeyError):
        catalog.get("missing")


@pytest.mark.parametrize("retain_start", [True, False])
def test_native_rna_annotation_and_legacy_json_export(retain_start, monkeypatch):
    from varcode import annotate_orf_product, annotate_orf_rna_product

    t = Transcript()
    reference = t.record().bind(t)
    origins = [(t.contig, p, t.strand) for p in t.positions] if retain_start else [None] * len(t.sequence)
    product = RealizedTranscriptProduct(t, t.sequence, None, evidence={
        "base_origins": origins, "rna_support": 10})
    with monkeypatch.context() as patch:
        def unexpected_export(*args, **kwargs):
            raise AssertionError("RNA annotation must not serialize")
        patch.setattr("varcode.orf._plain", unexpected_export)
        result = annotate_orf_rna_product(reference, product)
    assert result.reference is reference
    assert result.cdna_sequence == t.sequence
    assert result.evidence["rna_support"] == 10
    if retain_start:
        assert result.status == "conditional_translation"
        assert result.product.protein_sequence == "MKPG"
        assert result.product.evidence["rna_support"] == 10
        assert result.effect is not None
    else:
        assert result.status == "unresolved"
        assert result.reason == "reference_start_missing_or_duplicated"
        assert result.product is None and result.effect is None
        assert result.cdna_sequence == t.sequence
        assert result.evidence["rna_support"] == 10
    legacy = annotate_orf_product(reference, product)
    assert isinstance(legacy, dict)
    assert json.loads(json.dumps(legacy)) == json.loads(json.dumps(result.to_dict()))
