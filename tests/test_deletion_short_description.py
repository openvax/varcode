"""One-based deletion labels retain zero-based sequence coordinates (#643)."""

from types import SimpleNamespace

import pytest
from pyensembl import cached_release

from varcode import Variant
from varcode.effects import Deletion, EffectCollection
from varcode.effects.classify import classify_from_protein_diff


def transcript_for(protein):
    return SimpleNamespace(protein_sequence=protein, biotype="protein_coding",
                           gene=SimpleNamespace(biotype="protein_coding"))


@pytest.mark.parametrize("protein,offset,deleted,description,mutant", [
    ("MALK", 0, "M", "p.M1del", "ALK"),
    ("MALK", 1, "A", "p.A2del", "MLK"),
    ("MALK", 3, "K", "p.K4del", "MAL"),
    ("MALK", 1, "AL", "p.A2_L3del", "MK"),
    ("MALKK", 3, "KK", "p.K4_K5del", "MAL"),
])
def test_deletion_description_preserves_sequence_coordinates(protein, offset, deleted, description, mutant):
    effect = Deletion(None, transcript_for(protein), offset, deleted)
    assert effect.short_description == description
    assert effect.aa_mutation_start_offset == effect.aa_mutation_end_offset == offset
    assert effect.aa_ref == deleted
    assert effect.aa_alt == ""
    assert effect.mutant_protein_sequence == mutant


@pytest.mark.parametrize("reference,mutant,offset,deleted,description", [
    ("MLWWE", "MLWE", 3, "W", "p.W4del"),
    ("MACACACT", "MACACT", 5, "AC", "p.A6_C7del"),
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
@pytest.mark.parametrize("tid,contig,pos,ref,alt,offset,deleted,description", [
    # The GOLGA5 KK repeat from #643; Ensembl 81 retains the same K405
    # deletion as the reported Ensembl 87 model and is provisioned in CI.
    ("ENST00000163416", "14", 92811641, "GAGA", "G", 404, "K", "p.K405del"),
    ("ENST00000003084", "7", 117480101, "AGGTCGCCT", "", 2, "RSP", "p.R3_P5del"),
    # Minus-strand TP53: delete the most C-terminal E of MEEP, p.E3del.
    ("ENST00000269305", "17", 7676589, "CTC", "", 2, "E", "p.E3del"),
])
def test_native_deletion_labels_and_archives(
        annotator, tid, contig, pos, ref, alt, offset, deleted, description):
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
    expected = transcript.protein_sequence[:offset] + transcript.protein_sequence[offset + len(deleted):]
    assert effect.mutant_protein_sequence == expected
    restored = EffectCollection.from_json(EffectCollection([effect]).to_json())[0]
    assert restored.short_description == description
    assert restored.aa_mutation_start_offset == offset
    assert restored.mutant_protein_sequence == expected
