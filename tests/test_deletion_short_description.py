"""One-based deletion labels retain zero-based sequence coordinates (#643)."""

from types import SimpleNamespace

import pytest
from pyensembl import cached_release

from varcode import Variant
from varcode.effects import (
    Deletion, EffectCollection, FrameShift, FrameShiftTruncation, Insertion,
    Intergenic, LargeDeletion, NormalSplicing, PrematureStop,
    StructuralVariantEffect, Substitution,
)
from varcode.effects.classify import classify_from_protein_diff
from varcode.effects.effect_classes import PhaseCandidateSet
from varcode.effect_candidates import EffectCandidate
from varcode.germline import PhaseHypothesis
from varcode.splice_outcomes import SpliceOutcomeSet


def transcript_for(protein):
    return SimpleNamespace(protein_sequence=protein, biotype="protein_coding",
                           gene=SimpleNamespace(biotype="protein_coding"))


@pytest.mark.parametrize("protein,offset,deleted,description,range_description,mutant", [
    ("MALK", 0, "M", "p.M1del", "p.M1del", "ALK"),
    ("MALK", 1, "A", "p.A2del", "p.A2del", "MLK"),
    ("MALK", 3, "K", "p.K4del", "p.K4del", "MAL"),
    ("MALK", 1, "AL", "p.AL2del", "p.A2_L3del", "MK"),
    ("MALKK", 3, "KK", "p.KK4del", "p.K4_K5del", "MAL"),
])
def test_deletion_description_preserves_sequence_coordinates(
        protein, offset, deleted, description, range_description, mutant):
    effect = Deletion(None, transcript_for(protein), offset, deleted)
    assert effect.short_description == description
    assert effect.format_short_description() == description
    assert effect.format_short_description(deletion_notation="range") == range_description
    assert effect.format_short_description(deletion_notation="sequence") == description
    assert effect.aa_mutation_start_offset == effect.aa_mutation_end_offset == offset
    assert effect.aa_ref == deleted
    assert effect.aa_alt == ""
    assert effect.mutant_protein_sequence == mutant


@pytest.mark.parametrize("reference,mutant,offset,deleted,description", [
    ("MLWWE", "MLWE", 3, "W", "p.W4del"),
    ("MACACACT", "MACACT", 5, "AC", "p.AC6del"),
])
def test_repeat_normalization_keeps_the_most_c_terminal_equivalent_deletion(
        reference, mutant, offset, deleted, description):
    effect = classify_from_protein_diff(
        None, transcript_for(reference), reference, mutant,
        length_delta=3 * (len(mutant) - len(reference)))
    assert type(effect) is Deletion
    assert effect.aa_mutation_start_offset == offset
    assert effect.aa_ref == deleted
    assert effect.short_description == description
    assert effect.mutant_protein_sequence == mutant


@pytest.mark.parametrize("annotator", ["fast", "protein_diff", "transcript_model"])
@pytest.mark.parametrize("tid,contig,pos,ref,alt,offset,deleted,description,range_description", [
    # The GOLGA5 KK repeat from #643; Ensembl 81 retains the same K405
    # deletion as the reported Ensembl 87 model and is provisioned in CI.
    ("ENST00000163416", "14", 92811641, "GAGA", "G", 404, "K", "p.K405del", "p.K405del"),
    ("ENST00000003084", "7", 117480101, "AGGTCGCCT", "", 2, "RSP", "p.RSP3del", "p.R3_P5del"),
    # Minus-strand TP53: delete the most C-terminal E of MEEP, p.E3del.
    ("ENST00000269305", "17", 7676589, "CTC", "", 2, "E", "p.E3del", "p.E3del"),
])
def test_native_deletion_labels_and_archives(
        annotator, tid, contig, pos, ref, alt, offset, deleted, description, range_description):
    genome = cached_release(81)
    transcript = genome.transcript_by_id(tid)
    variant = Variant(contig, pos, ref, alt, genome)
    effect = variant.effect_on_transcript(transcript, annotator=annotator)
    effect = getattr(effect, "effect_if_splicing_unchanged", effect)
    assert type(effect) is Deletion
    assert effect.aa_mutation_start_offset == offset
    assert effect.aa_ref == deleted
    assert transcript.protein_sequence[offset:offset + len(deleted)] == deleted
    assert effect.short_description == description
    assert effect.format_short_description(deletion_notation="range") == range_description
    expected = transcript.protein_sequence[:offset] + transcript.protein_sequence[offset + len(deleted):]
    assert effect.mutant_protein_sequence == expected
    restored = EffectCollection.from_json(EffectCollection([effect]).to_json())[0]
    assert restored.short_description == description
    assert restored.format_short_description(deletion_notation="range") == range_description
    assert restored.aa_mutation_start_offset == offset
    assert restored.mutant_protein_sequence == expected


@pytest.mark.parametrize("make_effect,description", [
    (lambda t: Substitution(None, t, 1, "A", "V"), "p.A2V"),
    (lambda t: Insertion(None, t, 1, "RS"), "p.1insRS"),
    (lambda t: PrematureStop(None, t, 1, "A", ""), "p.A2*"),
    (lambda t: FrameShift(None, t, 1, ""), "p.A2fs"),
    (lambda t: FrameShiftTruncation(None, t, 1, "A"), "p.A2fs*"),
    (lambda t: Intergenic(None), "intergenic"),
])
def test_deletion_notation_does_not_change_other_consequences(make_effect, description):
    effect = make_effect(transcript_for("MALK"))
    assert effect.format_short_description() == description
    assert effect.format_short_description(deletion_notation="range") == description


def test_deletion_notation_propagates_through_delegated_descriptions():
    transcript = transcript_for("MALKT")
    deletion = Deletion(None, transcript, 1, "AL")
    normal = NormalSplicing(None, transcript, None, coding_effect=deletion)
    splice = SpliceOutcomeSet(None, transcript, [EffectCandidate(normal)])
    structural = StructuralVariantEffect(None, transcript, primary_effects=[deletion])
    phase = PhaseCandidateSet(None, transcript, [deletion],
                              [PhaseHypothesis()], ())
    for effect, sequence, endpoints in [
        (normal, "normal-splicing (p.AL2del)", "normal-splicing (p.A2_L3del)"),
        (splice, "splice-set:normal-splicing (p.AL2del)", "splice-set:normal-splicing (p.A2_L3del)"),
        (structural, "p.AL2del", "p.A2_L3del"),
        (phase, "?p.AL2del", "?p.A2_L3del"),
    ]:
        assert effect.short_description == sequence
        assert effect.format_short_description() == sequence
        assert effect.format_short_description(deletion_notation="range") == endpoints
        assert effect.short_description == sequence


def test_empty_and_specialized_wrappers_retain_their_labels():
    transcript = transcript_for("MALK")
    deletion = Deletion(None, transcript, 1, "AL")
    for effect, description in [
        (NormalSplicing(None, transcript, None), "normal-splicing"),
        (StructuralVariantEffect(None, transcript), "unresolved-structural-transcript"),
        (LargeDeletion(None, transcript, affected_exons=[], primary_effects=[deletion]), "sv-deletion"),
    ]:
        assert effect.format_short_description() == description
        assert effect.format_short_description(deletion_notation="range") == description


@pytest.mark.parametrize("include_header", [True, False])
def test_native_deletion_export_notation(tmp_path, include_header):
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084")
    deletion = Variant("7", 117480101, "AGGTCGCCT", "", genome).effect_on_transcript(transcript)
    substitution = Variant("7", 117531100, "T", "A", genome).effect_on_transcript(transcript)
    effects = EffectCollection([deletion, substitution], sort_key=False)
    assert list(effects.to_dataframe()["effect"]) == ["p.RSP3del", "p.L159M"]
    assert list(effects.to_dataframe(deletion_notation="range")["effect"]) == ["p.R3_P5del", "p.L159M"]
    for notation, label in [("sequence", "p.RSP3del"), ("range", "p.R3_P5del")]:
        path = tmp_path / (notation + ".csv")
        effects.to_csv(path, include_header=include_header, deletion_notation=notation)
        assert label in path.read_text()
        restored = EffectCollection.from_csv(path, genome=genome)
        restored_deletion = next(e for e in restored if isinstance(e, Deletion))
        assert restored_deletion.short_description == "p.RSP3del"
        assert restored_deletion.mutant_protein_sequence == deletion.mutant_protein_sequence


@pytest.mark.parametrize("notation", ["hgvs", "", None])
def test_invalid_deletion_notation_is_rejected_even_without_deletions(tmp_path, notation):
    effect = Intergenic(None)
    with pytest.raises(ValueError, match="deletion_notation"):
        effect.format_short_description(deletion_notation=notation)
    with pytest.raises(ValueError, match="deletion_notation"):
        EffectCollection([]).to_dataframe(deletion_notation=notation)
    path = tmp_path / "invalid.csv"
    with pytest.raises(ValueError, match="deletion_notation"):
        EffectCollection([]).to_csv(path, deletion_notation=notation)
    assert not path.exists()
