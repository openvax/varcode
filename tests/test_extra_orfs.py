"""One effect collection retains ordinary and additional coding consequences."""

from dataclasses import replace
from types import SimpleNamespace

import pytest
from pyensembl import cached_release

from varcode import Genome, Variant, VariantCollection
from varcode.effects import Failure, Substitution
from varcode.germline import GermlineContext
from .test_uorf import Transcript

orfdata = pytest.importorskip("orfdata")


@pytest.fixture(scope="module")
def native():
    catalog = orfdata.load_catalog()
    record = catalog.get("c7riboseqorf56")
    genome = cached_release(95)
    transcript = genome.transcript_by_id(record.transcript_id)
    variant = Variant("7", 66205486, "C", "G", genome=genome)
    return variant, record.bind(transcript), orfdata.ORFCatalog([record], catalog.metadata)


@pytest.mark.parametrize("input_kind", ["catalog", "bound", "iterable"])
def test_native_utr_and_coding_consequences_in_one_call(native, input_kind):
    variant, reference, catalog = native
    supplied = {"catalog": catalog, "bound": reference, "iterable": [reference]}[input_kind]
    ordinary = variant.effects()
    effects = variant.effects(extra_orfs=supplied)
    assert [(type(e), e.transcript_id) for e in effects if e.orf_annotation is None] == [
        (type(e), e.transcript_id) for e in ordinary]
    extra = [e for e in effects if e.orf is not None]
    assert len(extra) == 1
    effect = extra[0]
    assert isinstance(effect, Substitution)
    assert effect.short_description == "p.P2A"
    assert effect.reference_transcript is effect.orf.transcript
    assert effect.reference_transcript.id == reference.transcript.id
    assert effect.transcript_id == reference.transcript.id
    assert "orf_id=c7riboseqorf56" in str(effect)
    assert ":orf:" not in str(effect)
    assert effect.reference_transcript.protein_sequence == reference.transcript.protein_sequence
    assert effects.top_priority_effect() is effect
    assert effects.annotator == ordinary.annotator
    assert effects.annotator_version == ordinary.annotator_version
    assert effect in effects.drop_silent_and_noncoding()
    assert effects.orf_annotations == (effect.orf_annotation,)
    outcome = effect.orf_annotation.outcomes[0]
    assert outcome.mutant.protein_sequence.startswith("MASRRR")
    assert outcome.mutant.reference_evidence["mutant_translation_evidence"] == "not_assessed"
    assert effect.orf_annotation.to_dict()["model"]["orf_id"] == "c7riboseqorf56"
    assert effect.annotator == "transcript_model"
    assert effect in effects.groupby_transcript_id()[reference.transcript.id]
    assert effect in effects.filter_by_transcript_expression({reference.transcript.id: 10})
    assert hash(effect) == hash(variant.effects(extra_orfs=supplied).top_priority_effect())


def synthetic(reverse=False, split=False):
    transcript = Transcript(reverse=reverse, split=split)
    transcript.genome = cached_release(95)
    transcript.biotype = "lncRNA"
    transcript.is_protein_coding = False
    return transcript


@pytest.mark.parametrize("orf_type", ["uORF", "uoORF", "intORF", "dORF", "doORF", "lncRNA", "processed_transcript"])
@pytest.mark.parametrize("reverse,split", [(False, False), (True, False), (False, True), (True, True)])
def test_source_classes_and_noncoding_hosts_share_the_same_path(orf_type, reverse, split):
    transcript = synthetic(reverse, split)
    record = replace(transcript.record(), biotype=orf_type)
    reference = record.bind(transcript)
    variant = transcript.variant(6, "A", "G")
    result = variant.effects(extra_orfs=[reference])
    annotation = result.orf_annotations[0]
    assert annotation.reference is reference
    assert annotation.outcomes[0].mutant.protein_sequence == "MEPG"
    assert not transcript.is_protein_coding
    assert annotation.effect.reference_transcript is transcript


def test_two_overlapping_orfs_and_duplicate_inputs():
    transcript = synthetic()
    transcript.sequence = "CCCATGATGCCCGGGTAAATGAAATAA"
    first = replace(transcript.record(), protein_sequence="MMPG")
    second = replace(first, orf_id="second", coding_blocks=((106, 117),), protein_sequence="MPG")
    references = [first.bind(transcript), second.bind(transcript)]
    result = transcript.variant(9, "C", "G").effects(extra_orfs=[*references, references[0]])
    assert {a.orf_id for a in result.orf_annotations} == {first.orf_id, second.orf_id}
    assert len(result.orf_annotations) == 2
    assert {a.outcomes[0].mutant.protein_sequence for a in result.orf_annotations} == {"MMAG", "MAG"}


def test_collection_reuses_generator_for_every_variant(native):
    variant, reference, _ = native
    partner = Variant("7", 66205487, "C", "T", genome=variant.genome)
    result = VariantCollection([variant, partner]).effects(extra_orfs=iter([reference]))
    assert len(result.orf_annotations) == 2
    assert {a.variants[0] for a in result.orf_annotations} == {variant, partner}
    assert len(result.groupby_variant()) == 2


def test_same_identifier_in_different_sources_keeps_distinct_cis_groups():
    transcript = synthetic()
    original = transcript.record()
    references = [replace(original, source=source).bind(transcript) for source in ("one", "two")]
    variants = VariantCollection([transcript.variant(6, "A", "G"), transcript.variant(7, "A", "C")])
    resolver = SimpleNamespace(in_cis=lambda *args, **kwargs: True)
    effects = variants.effects(extra_orfs=references, phase_resolver=resolver)
    joint = [a for a in effects.orf_annotations if len(a.variants) == 2]
    assert len(joint) == 2
    assert {a.reference.model.source for a in joint} == {"one", "two"}


def test_scope_validation_includes_annotation_and_empty_collections(native):
    variant, reference, catalog = native
    # Releases 81 and 95 also differ in GRCh38 assembly patch. Hold the
    # assembly constant so this specifically exercises annotation validation.
    different_annotation = Genome(variant.genome)
    different_annotation.annotation_version = "different"
    earlier = Variant("7", variant.start, variant.ref, variant.alt, genome=different_annotation)
    with pytest.raises(ValueError, match="annotation_version"):
        earlier.effects(extra_orfs=[reference], raise_on_error=False)
    # A collection's context still matters when it contains no bound ORFs.
    class EmptyScoped:
        def __iter__(self):
            return iter(())
    empty = EmptyScoped()
    empty.reference = reference.context
    with pytest.raises(ValueError, match="annotation_version"):
        earlier.effects(extra_orfs=empty)
    different_assembly = Genome(variant.genome)
    different_assembly.assembly_accession = "GCA_000001405.18"
    other = Variant("7", variant.start, variant.ref, variant.alt, genome=different_assembly)
    with pytest.raises(ValueError, match="assembly_accession"):
        other.effects(extra_orfs=[reference], raise_on_error=False)
    wrong = orfdata.ORFCatalog([], {**catalog.metadata, "genome_build": "GRCh37"})
    with pytest.raises(ValueError, match="assembly"):
        variant.effects(extra_orfs=wrong)


def test_binding_failure_retains_source_record(native):
    variant, reference, catalog = native
    missing = replace(reference.record, transcript_id="ENST00000000000", initiation_context=None)
    supplied = orfdata.ORFCatalog([missing], catalog.metadata)
    with pytest.raises(ValueError):
        variant.effects(extra_orfs=supplied)
    effects = variant.effects(extra_orfs=supplied, raise_on_error=False)
    failures = [e for e in effects if e.orf_annotation is not None]
    assert len(failures) == 1 and isinstance(failures[0], Failure)
    annotation = failures[0].orf_annotation
    assert annotation.record is missing
    assert annotation.unresolved_reason
    assert annotation.reference is None and annotation.effect is None


def test_collection_exclusions_remain_visible(native):
    variant, reference, catalog = native
    class Excluded:
        def __iter__(self):
            return iter(())
    supplied = Excluded()
    supplied.reference = reference.context
    supplied.catalog = catalog
    supplied.excluded = (SimpleNamespace(orf_id=reference.model.orf_id, reason="unavailable path"),)
    result = variant.effects(extra_orfs=supplied, raise_on_error=False)
    assert result.orf_annotations[0].unresolved_reason == "unavailable path"


def test_patient_baseline_keeps_unknown_phase_outcomes():
    transcript = synthetic()
    reference = transcript.record().bind(transcript)
    variant = transcript.variant(6, "A", "G")
    partner = transcript.variant(7, "A", "C")
    result = variant.effects(extra_orfs=[reference], germline=GermlineContext.from_variants([partner]))
    outcomes = result.orf_annotations[0].outcomes
    assert {o.mutant.protein_sequence for o in outcomes} == {"MAPG", "MEPG"}


def test_collection_cis_products_keep_original_rna_identity():
    transcript = synthetic()
    reference = transcript.record().bind(transcript)
    first, second = transcript.variant(6, "A", "G"), transcript.variant(7, "A", "C")
    class Phase:
        phase_source = "fixture"
        def in_cis(self, left, right, transcript=None):
            if transcript is not None:
                assert transcript is reference.transcript
            return True
    result = VariantCollection([first, second]).effects(extra_orfs=[reference], phase_resolver=Phase())
    joint = [a for a in result.orf_annotations if len(a.variants) == 2]
    assert len(joint) == 1
    assert joint[0].outcomes[0].mutant.protein_sequence == "MAPG"
    assert joint[0].effect.orf is reference


def test_observed_canonical_protein_is_not_attached_to_extra_orf(native):
    variant, reference, _ = native
    resolver = SimpleNamespace(mutant_transcript=lambda variant, transcript: object())
    with pytest.raises(NotImplementedError, match="base origins"):
        variant.effects(extra_orfs=[reference], phase_resolver=resolver)
    result = variant.effects(extra_orfs=[reference], phase_resolver=resolver, raise_on_error=False)
    annotation = result.orf_annotations[0]
    assert annotation.unresolved_reason and annotation.reference is reference
    assert annotation.outcomes == ()


def test_no_extra_input_preserves_ordinary_api(native):
    variant, _, _ = native
    assert variant.effects().orf_annotations == ()
    assert variant.effects(extra_orfs=[]).orf_annotations == ()
    assert all(e.orf is None and e.reference_transcript is e.transcript for e in variant.effects())


def test_exports_cannot_silently_turn_an_orf_into_the_parent_cds(native, tmp_path):
    from varcode.effects import EffectCollection
    variant, reference, _ = native
    result = variant.effects(extra_orfs=[reference])
    restored = EffectCollection.from_json(result.to_json())
    assert restored.orf_annotations[0].to_dict() == result.orf_annotations[0].to_dict()
    effect = result.top_priority_effect()
    assert type(effect).from_json(effect.to_json()).orf_annotation.to_dict() == effect.orf_annotation.to_dict()
    frame = result.to_dataframe()
    row = frame[frame.orf_id.notna()].iloc[0]
    assert row.orf_id == "c7riboseqorf56"
    assert row.transcript_id == reference.transcript.id
    path = tmp_path / "effects.csv"
    result.to_csv(path)
    assert EffectCollection.from_csv(path).orf_annotations[0].to_dict() == result.orf_annotations[0].to_dict()
    frame.to_csv(path, index=False)  # A legacy summary still lacks its model.
    with pytest.raises(ValueError, match="reference coding models"):
        EffectCollection.from_csv(path, genome=variant.genome)


@pytest.mark.parametrize("value", ["all", "ms", "riboseq", "immunopeptidomics",
                                   ["proteomics", "immunopeptidomics"], "rule"])
def test_orfdata_orfs_for_the_variants_own_genome(native, value):
    variant, _, _ = native
    if value == "rule":
        value = orfdata.RiboSeq() & orfdata.MS()
    effects = variant.effects(extra_orfs=value)
    found = {(e.orf.orf_id, e.short_description) for e in effects if e.orf is not None}
    assert ("c7riboseqorf56", "p.P2A") in found
    # Ordinary effects are unchanged.
    assert [(type(e), e.transcript_id) for e in effects if e.orf_annotation is None] == [
        (type(e), e.transcript_id) for e in variant.effects()]


def test_placements_are_made_once_per_genome_and_reused(native):
    from varcode.extra_orfs import prepare_extra_orfs

    variant, _, _ = native
    resolver = prepare_extra_orfs("ms")
    # "ms" is MS(): one rule, one set of placements.
    assert prepare_extra_orfs("ms") is resolver and prepare_extra_orfs(orfdata.MS()) is resolver
    assert prepare_extra_orfs(orfdata.RiboSeq()) is not resolver
    with pytest.raises(ValueError, match="Unknown evidence method"):
        prepare_extra_orfs("MS")
    assert resolver.for_genome(variant.genome) is resolver.for_genome(variant.genome)
    assert prepare_extra_orfs(False) is None
    assert prepare_extra_orfs("all").evidence is None and prepare_extra_orfs("all") is not resolver
    with pytest.raises(ValueError, match="doesn't say which ORFs: pass 'all'"):
        variant.effects(extra_orfs=True)
    # A list of names is any of them, sharing placements with the same rule; "riboseq" is "ribosome".
    either = prepare_extra_orfs({"immunopeptidomics", "proteomics"})
    assert either is prepare_extra_orfs(orfdata.MS(methods="proteomics") | orfdata.MS(methods="immunopeptidomics"))
    assert prepare_extra_orfs("riboseq") is prepare_extra_orfs("ribosome")
    assert len(variant.effects(extra_orfs=False)) == len(variant.effects())
    collection = VariantCollection([variant]).effects(extra_orfs="ms")
    assert any(e.orf is not None and e.orf.orf_id == "c7riboseqorf56" for e in collection)


def test_a_declared_near_cognate_start_initiates_its_orf():
    transcript = synthetic()
    transcript.sequence = "CCCGTGAAACCCGGGTAAATGAAATAA"  # GTG, read as methionine, as the source declares.
    reference = replace(transcript.record(), start_codon="GTG").bind(transcript)
    effect = transcript.variant(6, "A", "G").effects(extra_orfs=[reference]).orf_annotations[0].effect
    assert isinstance(effect, Substitution) and effect.short_description == "p.K2E"


def test_orfs_orfdata_could_not_place_never_fail_zero_config_effects(native):
    from varcode.extra_orfs import ExtraORFInput

    variant, reference, catalog = native

    class Placed:
        def __iter__(self):
            return iter(())
    placed = Placed()
    placed.reference, placed.catalog = reference.context, catalog
    placed.excluded = (SimpleNamespace(orf_id=reference.model.orf_id, reason="unavailable path"),)
    # An explicit collection reports it, naming the ORF; ORFdata's own placements leave it out.
    with pytest.raises(ValueError, match="Extra ORF c7riboseqorf56: unavailable path"):
        variant.effects(extra_orfs=placed)
    quiet = variant.effects(extra_orfs=ExtraORFInput(placed, report_excluded=False))
    assert quiet.orf_annotations == () and len(quiet) == len(variant.effects())


def test_annotations_carry_the_catalogue_summary_not_its_record_lists(native):
    from varcode.extra_orfs import _summary

    summary = _summary(orfdata.load_catalog())
    assert summary["source"] and summary["genome_build"] == "GRCh38"
    assert not {"unresolved", "catalog_manifest", "records"} & set(summary)


def test_pickled_and_copied_effects_keep_their_orf_annotation(native):
    import copy
    import pickle

    variant, _, _ = native
    effects = variant.effects(extra_orfs="all")
    for restored in (pickle.loads(pickle.dumps(effects)), copy.deepcopy(effects)):
        annotations = restored.orf_annotations
        assert [a.orf_id for a in annotations] == [a.orf_id for a in effects.orf_annotations]
        assert all(a.effect.orf_annotation is a for a in annotations)
        assert {e.short_description for e in restored if e.orf is not None} >= {"p.P2A"}


def test_a_bad_extra_orfs_value_raises_before_any_annotation(native, monkeypatch):
    import varcode.variant

    variant, _, _ = native

    def annotated(*args, **kwargs):
        raise AssertionError("annotated before extra_orfs was checked")
    monkeypatch.setattr(varcode.variant, "predict_variant_effects", annotated)
    for bad in (True, "MS", "riboseq "):
        with pytest.raises(ValueError):
            variant.effects(extra_orfs=bad)
        with pytest.raises(ValueError):
            VariantCollection([variant]).effects(extra_orfs=bad)


def test_placements_follow_the_annotation_however_it_is_wrapped(native):
    import numpy as np

    from varcode import Genome
    from varcode.extra_orfs import clear_extra_orfs_cache, prepare_extra_orfs

    variant, reference, _ = native
    resolver = prepare_extra_orfs("ms")
    placed = resolver.for_genome(variant.genome)
    assert resolver.for_genome(Genome(variant.genome)) is placed
    clear_extra_orfs_cache()
    assert resolver.for_genome(variant.genome) is not placed
    # An array is never compared to "all"; an array of bound ORFs is a collection, as before.
    effects = variant.effects(extra_orfs=np.array([reference, reference], dtype=object))
    assert {a.orf_id for a in effects.orf_annotations} == {"c7riboseqorf56"}
