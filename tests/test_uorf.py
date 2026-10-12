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

oncoref = pytest.importorskip("oncoref")


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
    outcome, = result.to_dict()["outcomes"]
    assert outcome["baseline"]["protein_sequence"] == outcome["mutant"]["protein_sequence"]
    assert outcome["baseline"]["initiation_context"]["upstream"] != outcome["mutant"]["initiation_context"]["upstream"]


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
    genome = SimpleNamespace(reference_name="GRCh38", release=1)
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
        return oncoref.UORFRecord("test-orf", "gene", "GENE", "test", "1", self.strand,
                                  tuple(map(tuple, blocks)), "MKPG", "uORF", (), ())

    def variant(self, offset, ref, alt):
        position = self.positions[offset]
        if self.on_backward_strand:
            position -= len(ref) - 1
            ref, alt = reverse_complement(ref), reverse_complement(alt)
        return Variant("1", position, ref, alt, genome=cached_release(95))


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("split", [False, True])
def test_spliced_and_minus_strand_start_context(reverse, split):
    t = Transcript(reverse, split)
    view = ORFTranscript(t.record(), t, genome_build="GRCh38")
    assert view.start_codon_positions == tuple(t.positions[3:6])
    a = annotate_uorf_variants([t.variant(6, "A", "G")], t.record(), t, genome_build="GRCh38")
    outcomes = a.to_dict()["outcomes"]
    assert any(x["mutant"]["protein_sequence"] == "MEPG" for x in outcomes)


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
    outcome, = a.to_dict()["outcomes"]
    assert outcome["mutant"]["protein_sequence"] == protein
    if effect == "FrameShift":
        assert outcome["mutant"]["termination"] == "observed_stop"
        assert a.effect.aa_mutation_start_offset < 4


def test_known_cis_and_unknown_germline_phase_export_all_outcomes():
    t = Transcript()
    v = t.variant(6, "A", "G")
    partner = t.variant(7, "A", "C")
    a = annotate_uorf_variants([v, partner], t.record(), t, genome_build="GRCh38")
    assert a.to_dict()["outcomes"][0]["mutant"]["protein_sequence"] == "MAPG"
    uncertain = annotate_uorf_variants([v], t.record(), t, genome_build="GRCh38", germline_variants=[partner])
    outcomes = uncertain.to_dict()["outcomes"]
    assert {o["mutant"]["protein_sequence"] for o in outcomes} == {"MAPG", "MEPG"}
    assert {tuple(state for _, state in o["hypothesis"]["phase"]) for o in outcomes} == {("cis",), ("trans",)}


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
    outcomes = result.to_dict()["outcomes"]
    assert any(o["mutant"]["protein_sequence"] == "MPSRRRGGSSLIPGRMPILRFSVTTRYFSY"
               for o in outcomes if o["mutant"])
    assert len({o["mutant"]["protein_sequence"] for o in outcomes if o["mutant"]}) >= 3
    for outcome in outcomes:
        if outcome["mutant"]:
            assert not outcome["mutant"]["reference_evidence"]["reference_ms_support_on_retained_sequence"]
    assert json.loads(json.dumps(result.to_dict()))["model"]["orf_id"] == "c7riboseqorf56"


def test_stop_loss_with_no_downstream_sequence_marks_boundary():
    t = Transcript()
    t.sequence = t.sequence[:18]
    t.positions = t.positions[:18]
    t.end = t.positions[-1]
    t.exons = (SimpleNamespace(start=t.start, end=t.end),)
    a = annotate_uorf_variants([t.variant(15, "T", "C")], t.record(), t, genome_build="GRCh38")
    outcome, = a.to_dict()["outcomes"]
    assert outcome["mutant"]["protein_sequence"] == "MKPGQ"
    assert outcome["mutant"]["termination"] == "sequence_boundary"


def test_padded_insertion_ref_and_missing_intronic_ref_are_not_ignored():
    t = Transcript()
    with pytest.raises(ReferenceMismatchError):
        annotate_uorf_variants([t.variant(8, "C", "CA")], t.record(), t, genome_build="GRCh38")
    spliced = Transcript(split=True)
    v = Variant("1", 150, "A", "G", genome=cached_release(95))
    a = annotate_uorf_variants([v], spliced.record(), spliced, genome_build="GRCh38")
    assert type(a.effect).__name__ == "Unresolved"
