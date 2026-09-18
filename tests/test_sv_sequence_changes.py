"""SV predictions must survive filtering without calling unknowns unchanged."""

from dataclasses import replace

import pytest
from pyensembl import cached_release

from varcode import (
    EffectCandidate, MutantTranscript, ReferenceSegment, StructuralVariant,
    make_fusion_outcome,
)
from varcode.effects import (
    CrypticExonCandidate, EffectCollection, GeneFusion, LargeDeletion,
    StructuralVariantEffect, TranslocationToIntergenic, Unresolved,
)
from varcode.effects import structural


@pytest.fixture
def cftr():
    return cached_release(81).transcript_by_id("ENST00000003084")


@pytest.fixture
def isolate_span(monkeypatch):
    monkeypatch.setattr(structural, "_fusion_partners", lambda *args: ())
    monkeypatch.setattr(structural, "_enumerate_and_attach_cryptics", lambda *args: None)
    monkeypatch.setattr(structural, "_enumerate_and_attach_splice_outcomes", lambda *args: None)


def _effect(transcript, cdna=None, protein=None, segments=None, evidence=None):
    variant = StructuralVariant(transcript.contig, transcript.start, "BND",
                                genome=transcript.genome)
    model = MutantTranscript(
        reference_transcript=transcript, reference_segments=segments,
        cdna_sequence=cdna, mutant_protein_sequence=protein, evidence=evidence)
    return StructuralVariantEffect(variant, transcript, mutant_transcript=model)


def _identity_segments(transcript):
    return (ReferenceSegment(transcript, 0, len(transcript.sequence)),)


def test_seven_audited_fusion_proteins_survive_filter():
    genome = cached_release(95)
    transcript = genome.transcript_by_id("ENST00000538197")
    variant = StructuralVariant("4", 15012987, "BND", alt="N[4:2664478[",
                                mate_contig="4", mate_start=2664478, genome=genome)
    effect = variant.effect_on_transcript(transcript)
    fusions = [c.effect for c in effect.candidates if isinstance(c.effect, GeneFusion)]
    assert len(fusions) == 7
    assert all(f.mutant_protein_sequence != transcript.protein_sequence for f in fusions)
    assert all(f.modifies_protein_sequence is True for f in fusions)
    assert all(f.modifies_coding_sequence is True for f in fusions)
    assert len(EffectCollection(fusions).drop_silent_and_noncoding()) == 7
    assert list(EffectCollection([effect]).drop_silent_and_noncoding()) == [effect]


@pytest.mark.parametrize("transcript_id", ["ENST00000003084", "ENST00000357654"])
@pytest.mark.parametrize("region", ["coding", "start", "whole", "5utr", "3utr"])
def test_deletion_flags_on_both_strands(transcript_id, region, isolate_span):
    transcript = cached_release(81).transcript_by_id(transcript_id)
    exons = transcript.exons
    if region == "whole":
        start, end = transcript.start, transcript.end
    elif region == "coding":
        start, end = exons[4].start, exons[4].end
    elif region == "start":
        start, end = min(transcript.start_codon_positions), max(transcript.start_codon_positions)
    else:
        exon = exons[0] if region == "5utr" else exons[-1]
        at_low_end = (region == "5utr") == (transcript.strand == "+")
        start, end = ((exon.start, exon.start + 9) if at_low_end
                      else (exon.end - 9, exon.end))
    variant = StructuralVariant(transcript.contig, start, "DEL", end=end,
                                genome=transcript.genome)
    effect = variant.effect_on_transcript(transcript)
    assert isinstance(effect, LargeDeletion)
    expected = region in ("coding", "start", "whole")
    assert effect.modifies_coding_sequence is expected
    assert effect.modifies_protein_sequence is expected
    assert len(EffectCollection([effect]).drop_silent_and_noncoding()) == int(expected)


@pytest.mark.parametrize("change", ["none", "utr", "synonymous", "missense", "stop"])
def test_cds_and_protein_are_distinct_comparisons(cftr, change):
    cdna = cftr.sequence
    start = min(cftr.start_codon_spliced_offsets)
    if change == "utr":
        cdna = ("A" if cdna[0] != "A" else "C") + cdna[1:]
    elif change == "synonymous":
        # A GCT alanine codon changed to GCC.
        pos = next(i for i in range(start + 3, start + 300, 3) if cdna[i:i+3] == "GCT")
        cdna = cdna[:pos] + "GCC" + cdna[pos+3:]
    elif change == "missense":
        cdna = cdna[:start+3] + "TGG" + cdna[start+6:]
    elif change == "stop":
        pos = min(cftr.stop_codon_spliced_offsets)
        cdna = cdna[:pos] + "TAA" + cdna[pos+3:]
    effect = _effect(cftr, cdna=cdna, segments=_identity_segments(cftr))
    assert effect.modifies_coding_sequence is (change in ("synonymous", "missense", "stop"))
    assert effect.modifies_protein_sequence is (change == "missense")


@pytest.mark.parametrize("kind", ["DUP", "INV", "BND"])
def test_unresolved_models_are_not_silent(cftr, kind, isolate_span):
    variant = StructuralVariant(cftr.contig, cftr.start - 100, kind,
                                end=cftr.end - 50, genome=cftr.genome)
    if kind == "BND":
        effect = TranslocationToIntergenic(variant, cftr)
    else:
        effect = variant.effect_on_transcript(cftr)
    assert effect.modifies_coding_sequence is None
    assert effect.modifies_protein_sequence is None
    collection = EffectCollection([effect], annotator="fast", annotator_version="test")
    retained = collection.drop_silent_and_noncoding()
    assert list(retained) == [effect]
    assert retained.annotator == "fast"
    assert retained.annotator_version == "test"
    assert len(collection.drop_silent_and_noncoding(keep_unresolved=False)) == 0


def test_partial_fragment_and_unmapped_assembly_are_unknown(cftr):
    partial = _effect(cftr, segments=(ReferenceSegment(cftr, 0, 500),))
    assembly = _effect(cftr, cdna=cftr.sequence)
    for effect in (partial, assembly):
        assert effect.modifies_protein_sequence is None
        assert effect.modifies_coding_sequence is None


@pytest.mark.parametrize("primary_status", [False, None])
def test_later_changed_candidate_controls_flags_and_filter(cftr, primary_status):
    primary = _effect(cftr, protein=cftr.protein_sequence if primary_status is False else None)
    changed = _effect(cftr, protein=cftr.protein_sequence + "A")
    primary._extra_candidates = (EffectCandidate(changed, source="rna", evidence={"reads": 3}),)
    before = primary.candidates
    assert primary.modifies_coding_sequence is True
    assert primary.modifies_protein_sequence is True
    assert list(EffectCollection([primary]).drop_silent_and_noncoding(False)) == [primary]
    assert primary.candidates == before
    assert primary.candidates[1].evidence == {"reads": 3}
    # Self references and back edges must not recurse indefinitely.
    changed._attach_primary_effects([primary])
    assert primary.modifies_protein_sequence is True


def test_unchanged_plus_unknown_is_unknown(cftr):
    unchanged = _effect(cftr, protein=cftr.protein_sequence)
    unknown = _effect(cftr)
    unchanged._attach_primary_effects([unknown])
    assert unchanged.modifies_protein_sequence is None
    assert len(EffectCollection([unchanged]).drop_silent_and_noncoding()) == 1
    assert len(EffectCollection([unchanged]).drop_silent_and_noncoding(False)) == 0


def test_explicit_protein_without_cds_does_not_prove_cds_unchanged(cftr):
    effect = _effect(cftr, protein=cftr.protein_sequence)
    assert effect.modifies_protein_sequence is False
    assert effect.modifies_coding_sequence is None
    assert len(EffectCollection([effect]).drop_silent_and_noncoding()) == 0
    effect.mutant_transcript = replace(effect.mutant_transcript, mutant_protein_sequence="")
    assert effect.modifies_protein_sequence is True


def test_local_duplication_is_changed(cftr, isolate_span):
    exon = cftr.exons[4]
    variant = StructuralVariant(cftr.contig, exon.start, "DUP", end=exon.end,
                                genome=cftr.genome)
    effect = variant.effect_on_transcript(cftr)
    assert effect.modifies_coding_sequence is True
    assert effect.modifies_protein_sequence is True


def test_replacement_primary_effects_do_not_include_local_unknown(cftr):
    wrapper = _effect(cftr)
    unchanged = _effect(cftr, protein=cftr.protein_sequence)
    wrapper._primary_effects = (unchanged,)
    assert wrapper.modifies_protein_sequence is False
    assert len(EffectCollection([wrapper]).drop_silent_and_noncoding()) == 0


@pytest.mark.parametrize("synonymous", [False, True])
def test_imported_rna_has_explicit_cds_boundaries(cftr, synonymous):
    variant = StructuralVariant(cftr.contig, cftr.start, "BND", genome=cftr.genome)
    start = min(cftr.start_codon_spliced_offsets)
    cdna = cftr.sequence
    if synonymous:
        pos = next(i for i in range(start + 3, start + 300, 3) if cdna[i:i+3] == "GCT")
        cdna = cdna[:pos] + "GCC" + cdna[pos+3:]
    candidate = make_fusion_outcome(
        variant, cftr, sequence=cdna, cds_start=start,
        transcript_model_id="observed-model")
    assert candidate.effect.modifies_coding_sequence is synonymous
    assert candidate.effect.modifies_protein_sequence is False
    assert candidate.evidence["protein_status"] == "predicted_from_observed_rna"


def test_fusion_after_stop_is_unchanged_for_five_prime_partner(cftr):
    partner = cached_release(81).transcript_by_id("ENST00000357654")
    model = structural._build_fusion_mutant_transcript(
        cftr, cftr, cftr.end - 10, partner, partner.end - 10)
    assert model.mutant_protein_sequence == cftr.protein_sequence
    variant = StructuralVariant(cftr.contig, cftr.end - 10, "BND", genome=cftr.genome)
    effect = GeneFusion(variant, cftr, partner, mutant_transcript=model)
    assert effect.modifies_coding_sequence is False
    assert effect.modifies_protein_sequence is False
    assert len(EffectCollection([effect]).drop_silent_and_noncoding()) == 0
    # On the opposite end, compare the same allele with that transcript's
    # reference, but still use the 5' partner's annotated start for the CDS.
    other = GeneFusion(variant, partner, cftr, mutant_transcript=model,
                       five_prime_transcript=cftr, three_prime_transcript=partner)
    assert other.modifies_coding_sequence is True
    assert other.modifies_protein_sequence is True


@pytest.mark.parametrize("region", ["5utr", "3utr"])
def test_local_utr_duplication_does_not_change_cds(cftr, region, isolate_span):
    start = cftr.start + 10 if region == "5utr" else cftr.end - 20
    variant = StructuralVariant(cftr.contig, start, "DUP", end=start + 5,
                                genome=cftr.genome)
    effect = variant.effect_on_transcript(cftr)
    assert effect.modifies_coding_sequence is False
    assert effect.modifies_protein_sequence is False


def test_mitochondrial_translation_uses_correct_table():
    transcript = cached_release(81).transcript_by_id("ENST00000361624")
    effect = _effect(transcript, cdna=transcript.sequence,
                     segments=_identity_segments(transcript))
    assert effect.modifies_coding_sequence is False
    assert effect.modifies_protein_sequence is False


def test_ambiguous_coding_sequence_remains_unknown(cftr):
    start = min(cftr.start_codon_spliced_offsets)
    cdna = cftr.sequence[:start+6] + "N" + cftr.sequence[start+7:]
    effect = _effect(cftr, cdna=cdna, segments=_identity_segments(cftr))
    assert effect.modifies_coding_sequence is None
    assert effect.modifies_protein_sequence is None


def test_split_start_codon_and_deleted_start_base(cftr):
    start = min(cftr.start_codon_spliced_offsets)
    segments = (ReferenceSegment(cftr, 0, start + 1),
                ReferenceSegment(cftr, start + 1, len(cftr.sequence)))
    effect = _effect(cftr, cdna=cftr.sequence, segments=segments)
    assert effect.modifies_protein_sequence is False
    cut_segments = (segments[0], replace(segments[1], start=start + 2))
    cut_cdna = cftr.sequence[:start+1] + cftr.sequence[start+2:]
    cut = _effect(cftr, cdna=cut_cdna, segments=cut_segments)
    # Do not translate from one retained base of a destroyed initiator.
    assert cut.modifies_protein_sequence is None


def test_unresolved_and_cryptic_alternatives_are_not_known_unchanged(cftr):
    primary = _effect(cftr, protein=cftr.protein_sequence)
    unknown = Unresolved(primary.variant, cftr, mechanism="missing_sequence")
    cryptic = CrypticExonCandidate(primary.variant, cftr.contig, cftr.start, cftr.start + 20)
    for alternative in (unknown, cryptic):
        assert alternative.modifies_coding_sequence is None
        assert alternative.modifies_protein_sequence is None
        primary._extra_candidates = (EffectCandidate(alternative, source="test"),)
        assert primary.modifies_protein_sequence is None
        assert len(EffectCollection([primary]).drop_silent_and_noncoding()) == 1
        assert len(EffectCollection([primary]).drop_silent_and_noncoding(False)) == 0


def test_annotation_without_start_codon_remains_unresolved(cftr, monkeypatch):
    exon = cftr.exons[4]
    variant = StructuralVariant(cftr.contig, exon.start, "DEL", end=exon.end,
                                genome=cftr.genome)
    monkeypatch.setattr(type(cftr), "contains_start_codon", property(lambda self: False))
    effect = LargeDeletion(variant, cftr, affected_exons=[exon])
    assert effect.modifies_coding_sequence is True
    assert effect.modifies_protein_sequence is None


def test_annotation_without_cds_features_is_unknown(cftr, monkeypatch):
    def missing_cds(self):
        raise ValueError("No CDS features")

    monkeypatch.setattr(type(cftr), "coding_sequence_position_ranges", property(missing_cds))
    variant = StructuralVariant(cftr.contig, cftr.start, "DEL", end=cftr.end,
                                genome=cftr.genome)
    effect = LargeDeletion(variant, cftr, affected_exons=cftr.exons)
    assert effect.modifies_coding_sequence is None
    assert effect.modifies_protein_sequence is None


def test_empty_and_cyclic_candidate_sets_are_unknown(cftr):
    first, second = _effect(cftr), _effect(cftr)
    first._primary_effects = ()
    assert first.modifies_protein_sequence is None
    first._primary_effects = (second,)
    second._primary_effects = (first,)
    assert first.modifies_protein_sequence is None
