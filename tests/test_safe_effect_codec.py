"""Finite reviewed class registry and safe typed JSON codec for effect graphs (#674)."""

from __future__ import annotations

import builtins
from collections.abc import Mapping
import json
from types import SimpleNamespace

import pytest
from pyensembl import cached_release

from varcode import (
    EffectCandidate,
    EffectCollection,
    EffectReference,
    GermlineContext,
    SpliceOutcomeSet,
    StructuralVariant,
    Variant,
    VariantCollection,
    effect_collection_from_json,
    effect_collection_from_serializable_repr,
    effect_collection_to_json,
    effect_from_json,
    effect_from_serializable_repr,
    effect_serialization_classes,
    effect_to_json,
    predict_transcript_model_effect,
    register_annotator,
    registered_annotators,
    use_annotator,
)
from varcode.effects import (
    CrypticExonCandidate,
    GeneFusion,
    MutationEffect,
    SpliceDonor,
    StartLoss,
    StructuralVariantEffect,
    Substitution,
    TranslocationToIntergenic,
)
from varcode.effects.effect_classes import HaplotypeEffect, PhaseCandidateSet
from varcode.mutant_transcript import AssembledAllele
from varcode.orf import ORFTranscriptView
from varcode.orf_serialization import (
    ReferenceGenomeSnapshot,
    ReferenceTranscriptSnapshot,
    TranscriptFeatureSnapshot,
)

from .test_extra_orfs import synthetic as synthetic_orf_transcript


def test_effect_serialization_classes_is_finite_readonly_allowlist():
    registry = effect_serialization_classes()
    assert isinstance(registry, Mapping)
    with pytest.raises(TypeError):
        registry[("evil.mod", "Evil")] = object  # type: ignore[index]

    assert len(registry) > 40
    for key, cls in registry.items():
        assert isinstance(key, tuple) and len(key) == 2
        assert isinstance(key[0], str) and isinstance(key[1], str)
        assert isinstance(cls, type)
        assert key[0] not in ("os", "subprocess", "sys", "importlib", "shutil")

    # Public classes
    assert registry[("varcode.effects.structural_serialization", "EffectReference")] is EffectReference
    assert registry[("varcode.mutant_transcript", "AssembledAllele")] is AssembledAllele
    assert registry[("varcode.orf", "ORFTranscriptView")] is ORFTranscriptView
    assert registry[("varcode.orf_serialization", "ReferenceGenomeSnapshot")] is ReferenceGenomeSnapshot
    assert registry[("varcode.orf_serialization", "TranscriptFeatureSnapshot")] is TranscriptFeatureSnapshot
    assert registry[("varcode.orf_serialization", "ReferenceTranscriptSnapshot")] is ReferenceTranscriptSnapshot
    assert registry[("varcode.effects.effect_collection", "EffectCollection")] is EffectCollection
    assert registry[("varcode.effects.effect_classes", "Substitution")] is Substitution
    assert registry[("varcode.effects.effect_classes", "GeneFusion")] is GeneFusion
    assert registry[("varcode.effects.effect_classes", "HaplotypeEffect")] is HaplotypeEffect
    assert registry[("varcode.effects.effect_classes", "PhaseCandidateSet")] is PhaseCandidateSet
    assert registry[("varcode.splice_outcomes", "SpliceOutcomeSet")] is SpliceOutcomeSet

    # Backwards-compatible legacy serialized aliases
    assert registry[("varcode.effects.structural_serialization", "_EffectReference")] is EffectReference
    assert registry[("varcode.mutant_transcript", "_AssembledAllele")] is AssembledAllele
    assert registry[("varcode.effects.structural", "_AssembledAllele")] is AssembledAllele
    assert registry[("varcode.orf", "_ORFTranscript")] is ORFTranscriptView
    assert registry[("varcode.orf_serialization", "_ReferenceGenome")] is ReferenceGenomeSnapshot
    assert registry[("varcode.orf_serialization", "_Feature")] is TranscriptFeatureSnapshot
    assert registry[("varcode.orf_serialization", "_ReferenceTranscript")] is ReferenceTranscriptSnapshot


def test_ordinary_and_splice_effect_round_trip():
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084")

    fast_sub = Variant("7", 117531100, "T", "A", genome).effect_on_transcript(
        transcript, annotator="fast"
    )
    tm_sub = Variant("7", 117531100, "T", "A", genome).effect_on_transcript(
        transcript, annotator="transcript_model"
    )
    splice_effect = Variant("7", 117531115, "G", "A", genome).effect_on_transcript(
        transcript, annotator="fast"
    )
    assert isinstance(splice_effect, SpliceOutcomeSet)
    assert splice_effect.disrupted_signal_class is SpliceDonor

    for original in (fast_sub, tm_sub, splice_effect):
        encoded = effect_to_json(original)
        restored_json = effect_from_json(encoded)
        restored_bytes = effect_from_json(encoded.encode("utf-8"))
        restored_repr = effect_from_serializable_repr(json.loads(encoded))
        for restored in (restored_json, restored_bytes, restored_repr):
            assert type(restored) is type(original)
            assert restored.short_description == original.short_description
            assert restored.to_dict() == original.to_dict()

    restored_tm = effect_from_json(effect_to_json(tm_sub))
    assert restored_tm.candidates[0].effect is restored_tm

    # Pre-13.2 archives stored effect classes under varcode.effects.effect_classes;
    # ensure they still decode cleanly after splitting into domain modules.
    legacy_json = effect_to_json(fast_sub).replace(
        '"varcode.effects.coding_effects"', '"varcode.effects.effect_classes"'
    )
    assert '"varcode.effects.effect_classes"' in legacy_json
    restored_legacy = effect_from_json(legacy_json)
    assert type(restored_legacy) is Substitution
    assert restored_legacy.short_description == fast_sub.short_description
    assert restored_tm.most_likely_candidate is restored_tm.candidates[0]

    restored_splice = effect_from_json(effect_to_json(splice_effect))
    assert isinstance(restored_splice, SpliceOutcomeSet)
    assert restored_splice.disrupted_signal_class is SpliceDonor

    collection = EffectCollection(
        [fast_sub, tm_sub, splice_effect],
        annotator="transcript_model",
        annotator_version="11.17.0",
        annotated_at="2026-11-19T00:00:00Z",
    )
    coll_json = effect_collection_to_json(collection)
    for restored_coll in (
        effect_collection_from_json(coll_json),
        effect_collection_from_json(coll_json.encode("utf-8")),
        effect_collection_from_serializable_repr(json.loads(coll_json)),
    ):
        assert len(restored_coll) == 3
        assert restored_coll.annotator == "transcript_model"
        assert restored_coll.annotator_version == "11.17.0"
        assert [e.short_description for e in restored_coll] == [
            e.short_description for e in collection
        ]


def test_phased_and_germline_effect_graphs_round_trip():
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084")

    # 1. Unphased germline + somatic splice/coding candidate graph
    somatic = Variant("7", 117531115, "G", "A", genome)
    germline = Variant("7", 117531100, "T", "A", genome)
    uncertain = predict_transcript_model_effect(
        (somatic,), transcript, germline_variants=(germline,)
    )
    restored_uncertain = effect_from_json(effect_to_json(uncertain))
    assert len(restored_uncertain.candidates) == len(uncertain.candidates) == 8
    assert restored_uncertain.candidates[0].effect is restored_uncertain
    assert restored_uncertain.to_dict() == uncertain.to_dict()

    # 2. Fast annotator PhaseCandidateSet with unknown germline phase
    start_pos = sorted(transcript.start_codon_positions)
    g_start = Variant("7", start_pos[0], "A", "C", genome)
    s_start = Variant("7", start_pos[1], "T", "C", genome)
    phase_set = s_start.effect_on_transcript(
        transcript,
        annotator="fast",
        germline=GermlineContext.from_variants([g_start]),
    )
    assert isinstance(phase_set, PhaseCandidateSet)
    restored_phase_set = effect_from_json(effect_to_json(phase_set))
    assert isinstance(restored_phase_set, PhaseCandidateSet)
    assert {(c.effect.ref_codon, c.effect.alt_codon) for c in restored_phase_set.candidates} == {
        ("CTG", "CCG"),
        ("ATG", "ACG"),
    }
    assert all(isinstance(c.effect, StartLoss) for c in restored_phase_set.candidates)

    # 3. Cis phased joint haplotype effect with Variant-keyed metadata
    v1 = Variant("7", 117531100, "T", "A", genome)
    v2 = Variant("7", 117531101, "T", "G", genome)
    cis_resolver = SimpleNamespace(in_cis=lambda *args, **kwargs: True)
    phased_effects = VariantCollection([v1, v2]).effects(
        annotator="fast", phase_resolver=cis_resolver
    )
    haplotype_effects = [e for e in phased_effects if isinstance(e, HaplotypeEffect)]
    assert haplotype_effects
    restored_phased = effect_collection_from_json(effect_collection_to_json(phased_effects))
    restored_hap = [e for e in restored_phased if isinstance(e, HaplotypeEffect)]
    assert len(restored_hap) == len(haplotype_effects)
    assert restored_hap[0].short_description == haplotype_effects[0].short_description
    assert tuple(restored_hap[0].variants) == tuple(haplotype_effects[0].variants)
    assert restored_hap[0].mutant_protein_sequence == haplotype_effects[0].mutant_protein_sequence


def test_structural_and_extra_orf_effect_graphs_round_trip():
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084")
    partner = genome.transcript_by_id("ENST00000357654")

    # Assembled DEL carrying AssembledAllele segment
    exon = transcript.exons[4]
    sv_del = StructuralVariant(
        "7",
        exon.start + 10,
        "DEL",
        end=exon.end - 10,
        genome=genome,
        alt_assembly="ATGAAATAG",
    )
    del_effect = sv_del.effect_on_transcript(transcript, annotator="fast")
    restored_del = effect_from_json(effect_to_json(del_effect))
    assert restored_del.mutant_transcript == del_effect.mutant_transcript
    assert isinstance(
        restored_del.mutant_transcript.reference_segments[0].source,
        AssembledAllele,
    )

    # GeneFusion with self-linked candidate cycle + cryptic + splice candidates
    sv_bnd = StructuralVariant(
        "7",
        exon.end,
        "BND",
        genome=genome,
        mate_contig="17",
        mate_start=partner.exons[4].start,
    )
    from varcode.effects.structural_candidates import cryptic_exon_candidate

    fusion = GeneFusion(sv_bnd, transcript, partner)
    intergenic = TranslocationToIntergenic(sv_bnd, transcript)
    fusion.add_primary_effects((intergenic,))
    fusion.add_candidates((
        cryptic_exon_candidate(
            CrypticExonCandidate(sv_bnd, "7", 123, 145, donor_score=0.8, acceptor_score=0.9),
            base_evidence={"sv_type": "BND"},
        ),
        EffectCandidate(intergenic, "varcode_splice", {"reads": 2}),
        EffectCandidate(fusion, "RNA", {"witnesses": ["r1"]}),
    ))
    restored_fusion = effect_from_json(effect_to_json(fusion))
    assert isinstance(restored_fusion, GeneFusion)
    assert restored_fusion.candidates[0].effect is restored_fusion
    assert restored_fusion.candidates[-1].effect is restored_fusion
    assert restored_fusion.candidates[1].effect is restored_fusion.candidates[-2].effect

    # Extra-ORF effect graph (orf_effect_schema)
    orf_tx = synthetic_orf_transcript()
    orf_ref = orf_tx.record().bind(orf_tx)
    orf_var = orf_tx.variant(6, "A", "G")
    orf_effects = orf_var.effects(extra_orfs=[orf_ref])
    assert len(orf_effects.orf_annotations) == 1
    restored_orf_coll = effect_collection_from_json(effect_collection_to_json(orf_effects))
    assert len(restored_orf_coll.orf_annotations) == 1
    assert (
        restored_orf_coll.orf_annotations[0].to_dict()
        == orf_effects.orf_annotations[0].to_dict()
    )


def test_legacy_private_class_aliases_remain_readable():
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084")
    exon = transcript.exons[4]
    sv = StructuralVariant(
        "7",
        exon.start + 10,
        "DEL",
        end=exon.end - 10,
        genome=genome,
        alt_assembly="ATGAAATAG",
    )
    effect = sv.effect_on_transcript(transcript, annotator="fast")
    legacy_sv_json = (
        effect_to_json(effect)
        .replace('"EffectReference"', '"_EffectReference"')
        .replace(
            '{"__module__": "varcode.mutant_transcript", "__name__": "AssembledAllele"}',
            '{"__module__": "varcode.effects.structural", "__name__": "_AssembledAllele"}',
        )
    )
    restored_sv = effect_from_json(legacy_sv_json)
    assert restored_sv.to_dict() == effect.to_dict()
    assert isinstance(
        restored_sv.mutant_transcript.reference_segments[0].source,
        AssembledAllele,
    )

    orf_tx = synthetic_orf_transcript()
    orf_ref = orf_tx.record().bind(orf_tx)
    orf_effect = orf_tx.variant(6, "A", "G").effects(extra_orfs=[orf_ref])[0]
    legacy_orf_json = (
        effect_to_json(orf_effect)
        .replace('"EffectReference"', '"_EffectReference"')
        .replace('"ORFTranscriptView"', '"_ORFTranscript"')
        .replace('"ReferenceGenomeSnapshot"', '"_ReferenceGenome"')
        .replace('"TranscriptFeatureSnapshot"', '"_Feature"')
        .replace('"ReferenceTranscriptSnapshot"', '"_ReferenceTranscript"')
    )
    restored_orf = effect_from_json(legacy_orf_json)
    assert restored_orf.orf_annotation.to_dict() == orf_effect.orf_annotation.to_dict()

    legacy_coll_json = effect_collection_to_json(EffectCollection([effect])).replace(
        '"default_effect_sort_key"', '"_default_effect_sort_key"'
    )
    restored_coll = effect_collection_from_json(legacy_coll_json)
    assert len(restored_coll) == 1
    assert restored_coll[0].to_dict() == effect.to_dict()


def test_rejects_unreviewed_classes_functions_and_prevents_arbitrary_imports(
    monkeypatch, tmp_path
):
    # Warm the finite registry before spying on __import__.
    effect_serialization_classes()

    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084")
    effect = Variant("7", 117531100, "T", "A", genome).effect_on_transcript(
        transcript, annotator="transcript_model"
    )
    valid_repr = json.loads(effect_to_json(effect))

    imported_modules: list[str] = []
    real_import = builtins.__import__

    def spy_import(name, globals=None, locals=None, fromlist=(), level=0):
        imported_modules.append(name)
        if name.startswith("untrusted_never_imported"):
            raise AssertionError(f"Decoder attempted dynamic import of {name!r}")
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", spy_import)

    # 1. Unknown external module/class in __class__
    bad_class = {
        "__class__": {"__module__": "untrusted_never_imported.pkg", "__name__": "Exploit"},
        "payload": "boom",
    }
    with pytest.raises(ValueError, match="Unsupported class or function metadata"):
        effect_from_json(json.dumps(bad_class))
    with pytest.raises(ValueError, match="Unsupported class or function metadata"):
        effect_from_serializable_repr(bad_class)
    with pytest.raises(ValueError, match="Unsupported class or function metadata"):
        effect_collection_from_json(json.dumps(bad_class))

    # 2. Stdlib command execution classes / functions
    for mod_name, obj_name in [
        ("subprocess", "Popen"),
        ("os", "system"),
        ("builtins", "eval"),
    ]:
        payload = {"__class__": {"__module__": mod_name, "__name__": obj_name}}
        with pytest.raises(ValueError, match="Unsupported class or function metadata"):
            effect_from_json(json.dumps(payload))
        with pytest.raises(ValueError, match="Unsupported class or function metadata"):
            MutationEffect.from_json(json.dumps(payload))
        with pytest.raises(ValueError, match="Unsupported class or function metadata"):
            EffectCollection.from_json(json.dumps(payload))

    # 3. __function__ metadata and callable-in-__class__ rejected explicitly
    fn_payload = {"__function__": "system", "__module__": "os"}
    with pytest.raises(ValueError, match="Function metadata is not permitted"):
        effect_from_json(json.dumps(fn_payload))
    with pytest.raises(ValueError, match="Function metadata is not permitted"):
        effect_from_serializable_repr(fn_payload)
    callable_as_class = {
        "__class__": {
            "__module__": "varcode.effects.effect_collection",
            "__name__": "default_effect_sort_key",
        }
    }
    with pytest.raises(ValueError, match="Expected __class__ to be a reviewed class"):
        effect_from_json(json.dumps(callable_as_class))

    # 4. Standalone class metadata inside graph nodes[0][0]
    tampered_node = json.loads(json.dumps(valid_repr))
    tampered_node["nodes"]["__value__"][0]["__value__"][0] = {
        "__module__": "untrusted_never_imported.node",
        "__name__": "EvilEffect",
    }
    with pytest.raises(ValueError, match="Unsupported class or function metadata"):
        effect_from_json(json.dumps(tampered_node))

    # 5. Unreviewed class hidden inside __serialized_keys__
    tampered_keys = {
        "__serialized_keys__": [
            json.dumps(
                {
                    "__class__": {
                        "__module__": "untrusted_never_imported.keys",
                        "__name__": "KeyGadget",
                    }
                }
            )
        ],
        "__serialized_key_0__": "value",
    }
    with pytest.raises(ValueError, match="Unsupported class or function metadata"):
        effect_from_json(json.dumps(tampered_keys))

    assert not any(m.startswith("untrusted_never_imported") for m in imported_modules)

    # 6. EffectCollection.from_csv with tampered effect_json column also rejects
    start_loss = Variant(
        "7", min(transcript.start_codon_positions), "A", "G", genome
    ).effect_on_transcript(transcript, annotator="fast")
    csv_path = tmp_path / "effects.csv"
    EffectCollection([start_loss]).to_csv(str(csv_path))
    csv_text = csv_path.read_text(encoding="utf-8")
    assert "effect_json" in csv_text
    csv_text = csv_text.replace(
        f'""__module__"": ""{type(start_loss).__module__}""',
        '""__module__"": ""untrusted_never_imported.csv""',
    )
    assert "untrusted_never_imported.csv" in csv_text
    csv_path.write_text(csv_text, encoding="utf-8")
    with pytest.raises(ValueError, match="Unsupported class or function metadata"):
        EffectCollection.from_csv(str(csv_path))


def test_root_type_validation_and_invalid_inputs():
    genome = cached_release(81)
    transcript = genome.transcript_by_id("ENST00000003084")
    effect = Variant("7", 117531100, "T", "A", genome).effect_on_transcript(
        transcript, annotator="fast"
    )
    collection = EffectCollection([effect])

    with pytest.raises(ValueError, match="Expected archive root to decode to MutationEffect"):
        effect_from_json(effect_collection_to_json(collection))
    with pytest.raises(ValueError, match="Expected archive root to decode to MutationEffect"):
        effect_from_serializable_repr(json.loads(effect_collection_to_json(collection)))

    with pytest.raises(ValueError, match="Expected archive root to decode to EffectCollection"):
        effect_collection_from_json(effect_to_json(effect))
    with pytest.raises(ValueError, match="Expected archive root to decode to EffectCollection"):
        effect_collection_from_serializable_repr(json.loads(effect_to_json(effect)))

    with pytest.raises(ValueError, match="Expected archive root to decode to GeneFusion"):
        GeneFusion.from_json(effect_to_json(effect))

    with pytest.raises(TypeError, match="Expected a MutationEffect"):
        effect_to_json(collection)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="Expected an EffectCollection"):
        effect_collection_to_json(effect)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="Expected JSON str or bytes"):
        effect_from_json(123)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="Expected JSON str or bytes"):
        effect_collection_from_json(123)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="Cannot convert"):
        effect_from_serializable_repr((1, 2, 3))


def test_registered_annotators_public_api():
    assert registered_annotators() == (
        "dna",
        "fast",
        "transcript_model",
    )
    assert registered_annotators(include_deprecated=True) == (
        "dna",
        "fast",
        "protein_diff",
        "transcript_model",
    )

    class CustomAnnotator:
        name = "custom_test_annotator"

        def annotate_on_transcript(self, variant, transcript):
            return NotImplemented

    with use_annotator(CustomAnnotator()):
        assert "custom_test_annotator" in registered_annotators()
    assert "custom_test_annotator" not in registered_annotators()

    persisted = register_annotator(CustomAnnotator())
    try:
        assert persisted.name in registered_annotators()
    finally:
        from varcode.annotators.registry import _REGISTRY

        _REGISTRY.pop("custom_test_annotator", None)
    assert "custom_test_annotator" not in registered_annotators()
