"""Saved transcript consequences do not need a gene database to be restored."""

import shutil

import pytest
from pyensembl import Genome, cached_release

from varcode import (
    EffectCollection,
    StructuralVariant,
    Variant,
    iter_candidate_paths,
)
from varcode.effects import (
    AlternateStartCodon,
    Failure,
    StartLoss,
    Substitution,
    Unresolved,
)
from varcode.effects.effect_classes import HaplotypeEffect, PhaseCandidateSet
from varcode.germline import PhaseHypothesis


def _forbid_db_property():
    @property
    def unavailable(_self):
        raise RuntimeError("annotation database is unavailable")

    return unavailable


@pytest.mark.parametrize("annotator", ["fast", "transcript_model"])
def test_saved_substitution_defers_gene_lookup(monkeypatch, annotator):
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084")
    variant = Variant("7", 117531100, "T", "A", genome)
    effect = variant.effect_on_transcript(transcript, annotator=annotator)
    serialized = EffectCollection([effect]).to_json()
    expected_gene = transcript.gene

    with monkeypatch.context() as offline:
        offline.setattr(Genome, "db", _forbid_db_property())
        offline.setattr(genome, "_genes", {})
        restored_collection = EffectCollection.from_json(serialized)
        restored = restored_collection[0]
        assert isinstance(restored, Substitution)
        assert restored.short_description == "p.L159M"
        assert restored.transcript_id == transcript.id
        assert restored.transcript_name == transcript.name
        assert restored.gene_id == expected_gene.id
        assert list(restored_collection.groupby_gene_id()) == [expected_gene.id]
        assert transcript.id in str(restored)
        if annotator == "transcript_model":
            assert restored.candidates[0].effect is restored
            assert restored.candidates[0].outcomes[0].mutant == effect.candidates[0].outcomes[0].mutant
            assert restored.original_protein_sequence == effect.original_protein_sequence
            assert restored.mutant_protein_sequence == effect.mutant_protein_sequence
        with pytest.raises(RuntimeError, match="annotation database"):
            _ = restored.gene
        with pytest.raises(RuntimeError, match="annotation database"):
            _ = restored.gene_name

    assert restored.gene == expected_gene
    assert restored.gene_name == expected_gene.name
    assert restored.gene_id == expected_gene.id
    # Gene has always been writable; lazy lookup and gene_id must preserve overrides.
    restored.gene = None
    assert restored.gene is None
    assert restored.gene_id is None


@pytest.mark.parametrize("with_transcript", [False, True])
def test_failure_restore_does_not_resolve_gene(monkeypatch, with_transcript):
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084") if with_transcript else None
    expected_gene = transcript.gene if transcript is not None else None
    variant = Variant("7", 117531100, "T", "A", genome)
    serialized = Failure(variant, transcript, "unresolved input").to_json()

    with monkeypatch.context() as offline:
        offline.setattr(Genome, "db", _forbid_db_property())
        offline.setattr(genome, "_genes", {})
        restored = Failure.from_json(serialized)
        assert restored.error == "unresolved input"
        assert restored.transcript_id == (transcript.id if with_transcript else None)
        assert restored.gene_id == (expected_gene.id if expected_gene is not None else None)
        assert "unresolved input" in str(restored)
        if not with_transcript:
            assert restored.gene is None

    assert restored.gene == expected_gene


@pytest.mark.parametrize("annotator", ["fast", "transcript_model"])
def test_saved_start_loss_and_alternate_start_restore_without_db(monkeypatch, annotator):
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084")
    start_pos = min(transcript.start_codon_positions)
    start_loss = Variant(transcript.contig, start_pos, "A", "G", genome).effect_on_transcript(
        transcript, annotator=annotator)
    assert isinstance(start_loss, StartLoss)
    alt_start = AlternateStartCodon(
        Variant(transcript.contig, start_pos, "A", "C", genome),
        transcript,
        "ATG",
        "CTG",
    )
    serialized = EffectCollection([start_loss, alt_start]).to_json()

    with monkeypatch.context() as offline:
        offline.setattr(Genome, "db", _forbid_db_property())
        offline.setattr(genome, "_genes", {})
        restored = EffectCollection.from_json(serialized)
        restored_loss, restored_alt = restored[0], restored[1]
        assert isinstance(restored_loss, StartLoss)
        assert restored_loss.short_description == "p.M1? (start-loss)"
        assert restored_loss.aa_ref == "M"
        assert restored_loss.alt_codon_amino_acid == "V"
        assert restored_loss.gene_id == transcript.gene_id
        assert isinstance(restored_alt, AlternateStartCodon)
        assert restored_alt.short_description == "alternate-start-codon (ATG>CTG)"
        assert restored_alt.aa_ref == "M"
        assert restored_alt.gene_id == transcript.gene_id

    # Archives created before aa_ref was serialized still reconstruct when online.
    legacy_state = start_loss.to_dict()
    if "nodes" in legacy_state:
        cls, fields, attached = legacy_state["nodes"][0]
        legacy_fields = {k: v for k, v in fields.items() if k != "aa_ref"}
        legacy_state = {**legacy_state, "nodes": ((cls, legacy_fields, attached), *legacy_state["nodes"][1:])}
    else:
        legacy_state = {k: v for k, v in legacy_state.items() if k != "aa_ref"}
    assert StartLoss.from_dict(legacy_state).aa_ref == "M"


def test_relocated_unindexed_custom_bundle_roundtrip(tmp_path):
    gtf = tmp_path / "custom.gtf"
    gtf.write_text(
        '1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "G1"; gene_name "GEN1"; gene_biotype "protein_coding";\n'
        '1\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; transcript_name "T1-001"; transcript_biotype "protein_coding";\n'
        '1\ttest\texon\t1\t30\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "1"; exon_id "E1";\n'
        '1\ttest\texon\t51\t100\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "2"; exon_id "E2";\n'
        '1\ttest\tCDS\t10\t30\t.\t+\t0\tgene_id "G1"; transcript_id "T1"; protein_id "P1";\n'
        '1\ttest\tCDS\t51\t65\t.\t+\t0\tgene_id "G1"; transcript_id "T1"; protein_id "P1";\n'
        '1\ttest\tstart_codon\t10\t12\t.\t+\t0\tgene_id "G1"; transcript_id "T1";\n'
        '1\ttest\tstop_codon\t66\t68\t.\t+\t0\tgene_id "G1"; transcript_id "T1";\n'
        '1\ttest\tgene\t200\t250\t.\t+\t.\tgene_id "G2"; gene_name "GEN2"; gene_biotype "lincRNA";\n'
        '1\ttest\ttranscript\t200\t250\t.\t+\t.\tgene_id "G2"; transcript_id "T2"; transcript_name "T2-001"; transcript_biotype "lincRNA";\n'
        '1\ttest\texon\t200\t250\t.\t+\t.\tgene_id "G2"; transcript_id "T2"; exon_number "1"; exon_id "E3";\n'
    )
    seq_t1 = "A" * 9 + "ATGGCTAAAGAACAATTAGGA" + "CCATTAGTTAAAGCA" + "TAA" + "C" * 32
    cdna = tmp_path / "cdna.fa"
    cdna.write_text(f">T1\n{seq_t1}\n>T2\n" + "A" * 51 + "\n")
    prot = tmp_path / "prot.fa"
    prot.write_text(">P1\nMAKEQLGPLVKA\n")

    cache_dir = tmp_path / "cache"
    genome = Genome(
        reference_name="custom1",
        annotation_name="custom_ann",
        annotation_version="1",
        gtf_path_or_url=str(gtf),
        transcript_fasta_paths_or_urls=[str(cdna)],
        protein_fasta_paths_or_urls=[str(prot)],
        cache_directory_path=str(cache_dir),
    )
    genome.index()

    t1 = genome.transcript_by_id("T1")
    v_sub = Variant("1", 13, "G", "A", genome)
    v_sub2 = Variant("1", 16, "A", "G", genome)
    sub_eff = v_sub.effect_on_transcript(t1, annotator="transcript_model")
    sub_eff2 = v_sub2.effect_on_transcript(t1, annotator="transcript_model")
    h1 = PhaseHypothesis(cis=(), trans=(v_sub2,), phase_state="unknown", haplotype="h1")
    h2 = PhaseHypothesis(cis=(v_sub2,), trans=(), phase_state="unknown", haplotype="h2")

    effects = [
        *Variant("1", 10, "A", "G", genome).effects(annotator="transcript_model"),  # StartLoss
        *v_sub.effects(annotator="transcript_model"),                               # Substitution
        *Variant("1", 12, "GGCT", "G", genome).effects(annotator="transcript_model"),  # Deletion
        *Variant("1", 16, "A", "T", genome).effects(annotator="transcript_model"),  # PrematureStop
        *Variant("1", 16, "AA", "A", genome).effects(annotator="transcript_model"),  # FrameShift
        *Variant("1", 66, "T", "C", genome).effects(annotator="transcript_model"),  # StopLoss
        *Variant("1", 31, "G", "A", genome).effects(annotator="fast"),              # SpliceOutcomeSet
        *Variant("1", 210, "A", "G", genome).effects(annotator="transcript_model"), # NoncodingTranscript
        *StructuralVariant(
            contig="1", start=25, end=25, sv_type="BND", alt="G[1:210[", ref="G",
            mate_contig="1", mate_start=210, mate_orientation="[[", genome=genome,
        ).effects(),
        AlternateStartCodon(Variant("1", 10, "A", "C", genome), t1, "ATG", "CTG"),
        HaplotypeEffect((v_sub, v_sub2), t1, sub_eff.mutant_transcript, phase_source="vcf_ps"),
        PhaseCandidateSet(v_sub, t1, candidates=(sub_eff, sub_eff2), hypotheses=(h1, h2), germline_variants=(v_sub2,)),
        Unresolved(v_sub, t1, mechanism="unsupported_annotation", reason="fixture", reference_effect=sub_eff),
        Failure(v_sub, t1, "fixture error"),
    ]
    serialized = EffectCollection(effects).to_json()

    # Wipe the PyEnsembl SQLite/FASTA index directory to simulate an unindexed relocated bundle.
    shutil.rmtree(cache_dir)

    restored = EffectCollection.from_json(serialized)
    assert len(restored) == len(effects)
    for item in restored:
        assert item.short_description
        assert item.format_short_description(deletion_notation="range")
        assert str(item)
        assert item == item
        assert isinstance(hash(item), int)
        _ = item.transcript_id
        _ = item.transcript_name
        _ = item.gene_id
        _ = list(iter_candidate_paths(item))

    assert len(EffectCollection(list(restored), distinct=True)) == len(restored)
    assert set(restored.groupby_gene_id()) == {"G1", "G2"}
    assert len(restored.drop_silent_and_noncoding()) > 0
