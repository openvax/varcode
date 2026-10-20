"""ORF archives preserve computed outcomes and can be read offline."""

from dataclasses import replace
import json

import pytest

from varcode import EffectCollection, VariantCollection
from varcode.effects import Failure
from varcode.germline import GermlineContext
from . import test_extra_orfs
from .test_extra_orfs import synthetic

native = test_extra_orfs.native

orfdata = pytest.importorskip("orfdata")


def assert_annotation(before, after):
    assert after.to_dict() == before.to_dict()
    assert type(after.record) is type(before.record)
    if after.reference is not None:
        assert isinstance(after.reference, orfdata.ReferenceORF)
        assert after.effect.orf is after.reference
        assert after.reference.context.key == before.reference.context.key
        assert after.reference.transcript.sequence == before.reference.transcript.sequence
    for original, restored in zip(before.outcomes, after.outcomes):
        assert type(original.effect) is type(restored.effect)
        assert restored.effect.short_description == original.effect.short_description


def test_single_orf_reads_without_annotation_or_catalogue(native, monkeypatch):
    from pyensembl import Genome
    variant, reference, _ = native
    effect = variant.effects(extra_orfs=[reference]).orf_annotations[0].effect
    payload = effect.to_json()

    def unavailable(*args, **kwargs):
        raise AssertionError("Archive tried to consult external reference data")

    monkeypatch.setattr(orfdata, "load_catalog", unavailable)
    monkeypatch.setattr(Genome, "transcript_by_id", unavailable)
    monkeypatch.setattr(Genome, "gene_by_id", unavailable)
    monkeypatch.setattr(Genome, "download", unavailable)
    restored = type(effect).from_json(payload)
    assert_annotation(effect.orf_annotation, restored.orf_annotation)
    assert restored.orf_annotation.effect is restored
    assert restored.candidates[0].effect is restored
    assert restored.to_json() == payload


def test_source_collisions_and_germline_outcomes_survive():
    transcript = synthetic()
    first = replace(transcript.record(), source="first")
    second = replace(first, source="second")
    variant = transcript.variant(6, "A", "G")
    partner = transcript.variant(7, "A", "C")
    effects = variant.effects(
        extra_orfs=[first.bind(transcript), second.bind(transcript)],
        germline=GermlineContext.from_variants([partner]))
    restored = EffectCollection.from_json(effects.to_json())
    assert len(restored.orf_annotations) == 2
    for before, after in zip(effects.orf_annotations, restored.orf_annotations):
        assert_annotation(before, after)
        assert {o.mutant.protein_sequence for o in after.outcomes} == {"MAPG", "MEPG"}
        assert {tuple(state for _, state in o.hypothesis.phase) for o in after.outcomes} == {("cis",), ("trans",)}


def test_multiple_placements_and_scoped_evidence_survive():
    from hashlib import sha256
    from orfdata.evidence_summary import EvidenceSummary
    from orfdata.reference import ScopedTarget
    first, second = synthetic(), synthetic()
    second.id, second.name = "second-rna", "SECOND"
    record = replace(first.record(), source="test-source")
    context = record.bind(first).context
    evidence = EvidenceSummary(
        source="test-source", release="v1", target=ScopedTarget(context, "orf", record.orf_id, "test:v1"),
        protein_sha256=sha256(record.protein_sequence.encode()).hexdigest(),
        method="proteomics", metric="peptides", value=3, definition="source peptide count", scope="orf",
        source_url="https://example.org/synthetic-test", source_sha256="0" * 64)
    record = replace(record, summaries=(evidence,))
    variants = VariantCollection([first.variant(6, "A", "G"), first.variant(7, "A", "C")])

    class Phase:
        def in_cis(self, first, second, transcript=None):
            return True

    effects = variants.effects(
        extra_orfs=[record.bind(first), record.bind(second, any_transcript=True)], phase_resolver=Phase())
    restored = EffectCollection.from_json(effects.to_json())
    assert {a.reference.transcript.id for a in restored.orf_annotations} == {"test", "second-rna"}
    for before, after in zip(effects.orf_annotations, restored.orf_annotations):
        assert_annotation(before, after)
        assert after.reference.evidence.summaries[0].key == evidence.key
    joint = [a for a in restored.orf_annotations if len(a.variants) == 2]
    assert {a.reference.transcript.id for a in joint} == {"test", "second-rna"}
    for after in joint:
        assert after.outcomes[0].mutant.protein_sequence == "MAPG"


def test_failed_binding_keeps_record_and_reason(native):
    variant, reference, catalog = native
    bad = replace(reference.record, protein_sequence="A" + reference.record.protein_sequence[1:])
    effects = variant.effects(extra_orfs=orfdata.ORFCatalog([bad], catalog.metadata), raise_on_error=False)
    restored = EffectCollection.from_json(effects.to_json())
    before, = effects.orf_annotations
    after, = restored.orf_annotations
    assert_annotation(before, after)
    assert after.reference is None
    assert after.unresolved_reason
    assert any(isinstance(e, Failure) and e.orf_annotation is after for e in restored)


def test_csv_archive_rejects_missing_or_stale_payloads(native, tmp_path):
    import pandas as pd
    variant, reference, _ = native
    effects = variant.effects(extra_orfs=[reference])
    frame = effects.to_dataframe(include_archive=True)
    path = tmp_path / "effects.csv"
    frame.loc[0, "effect_json"] = None
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="Incomplete effect_json"):
        EffectCollection.from_csv(path)
    frame = effects.to_dataframe(include_archive=True)
    frame.loc[0, "orf_id"] = "another-orf"
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="disagrees"):
        EffectCollection.from_csv(path)


def test_orf_graph_unknown_schema_is_rejected(native):
    variant, reference, _ = native
    effect = variant.effects(extra_orfs=[reference]).orf_annotations[0].effect
    state = effect.to_dict()
    with pytest.raises(ValueError, match="Unsupported.*schema"):
        type(effect).from_dict({**state, "orf_effect_schema": 999})


def test_fusion_alternatives_and_partial_partners_round_trip_offline(native, monkeypatch):
    from pyensembl import Genome
    from varcode.sv_allele_parser import parse_symbolic_alt
    variant, reference, _ = native
    fusion = parse_symbolic_alt("7", 66205522, "G", "G[7:66127704[", genome=variant.genome)
    annotation = fusion.effects(extra_orfs=[reference]).orf_annotations[0]
    payload = annotation.effect.to_json()

    def unavailable(*args, **kwargs):
        raise AssertionError("Fusion archive tried to look up a reference feature")

    monkeypatch.setattr(Genome, "gene_by_id", unavailable)
    monkeypatch.setattr(Genome, "transcript_by_id", unavailable)
    restored = type(annotation.effect).from_json(payload)
    assert_annotation(annotation, restored.orf_annotation)
    assert len({outcome.mutant.protein_sequence for outcome in restored.orf_annotation.outcomes
                if outcome.mutant is not None}) >= 3
