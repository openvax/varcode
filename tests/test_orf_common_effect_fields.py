"""Published ncORF models expose the same product fields through every candidate."""

from varcode import EffectCollection, GermlineContext, Variant, iter_effect_paths, iter_effects
from . import test_extra_orfs

native = test_extra_orfs.native


def test_native_ncorf_phase_candidates_keep_host_identity_and_actual_baselines(native):
    variant, reference, _ = native
    protein = reference.protein_sequence
    assert protein[:2] == "MP"
    # Genomic first/second C of the source ORF's proline codon: CCN -> GCN,
    # CTN or GTN. These controlled edits test a public reference ORF, not a
    # claim of observed patient variants or observed mutant translation.
    germline = Variant("7", 66205487, "C", "T", genome=variant.genome)
    effects = variant.effects(extra_orfs=[reference],
                              germline=GermlineContext.from_variants([germline]))
    expected = {(protein, "MA" + protein[2:]), ("ML" + protein[2:], "MV" + protein[2:])}
    for collection in (effects, EffectCollection.from_json(effects.to_json())):
        annotation, = collection.orf_annotations
        pairs = set()
        for effect in iter_effects(annotation.effect):
            assert effect.orf.orf_id == reference.orf_id
            assert effect.transcript_id == reference.transcript.id
            assert effect.reference_transcript.id == reference.transcript.id
            assert effect.gene_id == reference.transcript.gene_id
            pairs.add((effect.original_protein_sequence, effect.mutant_protein_sequence))
            assert effect.baseline_transcript.mutant_protein_sequence == effect.original_protein_sequence
            assert effect.mutant_transcript.mutant_protein_sequence == effect.mutant_protein_sequence
        assert pairs == expected
        paths = list(iter_effect_paths(annotation.effect))
        assert any(path for _, path in paths)
        assert all(effect.orf is not None for effect, _ in paths)
        assert all(o.mutant.reference_evidence["mutant_translation_evidence"] == "not_assessed"
                   for o in annotation.outcomes)
