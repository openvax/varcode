"""Saved transcript consequences do not need a gene database to be restored."""

import pytest
from pyensembl import Genome, cached_release

from varcode import EffectCollection, Variant
from varcode.effects import Failure, Substitution


@pytest.mark.parametrize("annotator", ["fast", "transcript_model"])
def test_saved_substitution_defers_gene_lookup(monkeypatch, annotator):
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084")
    variant = Variant("7", 117531100, "T", "A", genome)
    effect = variant.effect_on_transcript(transcript, annotator=annotator)
    serialized = EffectCollection([effect]).to_json()
    expected_gene = transcript.gene

    def unavailable(*args, **kwargs):
        raise RuntimeError("annotation database is unavailable")

    with monkeypatch.context() as offline:
        offline.setattr(Genome, "gene_by_id", unavailable)
        restored = EffectCollection.from_json(serialized)[0]
        assert isinstance(restored, Substitution)
        assert restored.short_description == "p.L159M"
        assert restored.transcript_id == transcript.id
        if annotator == "transcript_model":
            assert restored.candidates[0].effect is restored
            assert restored.candidates[0].outcomes[0].mutant == effect.candidates[0].outcomes[0].mutant
        with pytest.raises(RuntimeError, match="annotation database"):
            _ = restored.gene
    assert restored.gene == expected_gene
    assert restored.gene_name == expected_gene.name
    assert restored.gene_id == expected_gene.id
    # Gene has always been writable; a lazy lookup must preserve that API.
    restored.gene = None
    assert restored.gene is None


@pytest.mark.parametrize("with_transcript", [False, True])
def test_failure_restore_does_not_resolve_gene(monkeypatch, with_transcript):
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084") if with_transcript else None
    variant = Variant("7", 117531100, "T", "A", genome)
    serialized = Failure(variant, transcript, "unresolved input").to_json()

    def unavailable(*args, **kwargs):
        raise AssertionError("Failure restoration must not query annotation")

    monkeypatch.setattr(Genome, "gene_by_id", unavailable)
    restored = Failure.from_json(serialized)
    assert restored.error == "unresolved input"
    assert restored.transcript_id == (transcript.id if with_transcript else None)
    if not with_transcript:
        assert restored.gene is None
