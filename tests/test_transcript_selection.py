"""Select transcript work without changing overlap or consequence semantics."""

from types import SimpleNamespace

import pytest
from pyensembl import Genome, cached_release

from varcode import GermlineContext, StructuralVariant, Variant, VariantCollection
from varcode.effects import (
    Failure, IncompleteTranscript, Intergenic, Intragenic, NoncodingTranscript,
    predict_variant_effects,
)


@pytest.fixture
def local_genome(tmp_path):
    # A gene's span can extend beyond the transcripts overlapping the variant.
    rows, genes = [], set()
    for gene, transcript, biotype, start, end in [
        ("g1", "coding", "protein_coding", 100, 200),
        ("g1", "lnc", "lncRNA", 100, 200),
        ("g1", "pseudo", "polymorphic_pseudogene", 100, 200),
        ("g1", "nmd", "nonsense_mediated_decay", 100, 200),
        ("g1", "ig", "IG_C_gene", 100, 200),
        ("g2", "remote", "lncRNA", 800, 900),
    ]:
        attributes = (
            'gene_id "%s"; gene_name "%s"; gene_biotype "protein_coding"; '
            'transcript_id "%s"; transcript_name "%s"; '
            'transcript_biotype "%s"; exon_id "%s.exon"; exon_number "1";'
        ) % (gene, gene, transcript, transcript, biotype, transcript)
        features = [("transcript", start, end), ("exon", start, end)]
        if transcript == "coding":
            # A truncated start codon makes the coding model incomplete.
            features.append(("start_codon", start, start + 1))
        if gene not in genes:
            features.insert(0, ("gene", 100, 1000))
            genes.add(gene)
        for feature, low, high in features:
            rows.append("1\ttest\t%s\t%d\t%d\t.\t+\t.\t%s\n" % (
                feature, low, high, attributes))
    gtf = tmp_path / "selection.gtf"
    gtf.write_text("".join(rows))
    genome = Genome(reference_name="test", annotation_name="selection",
                    gtf_path_or_url=str(gtf), cache_directory_path=str(tmp_path))
    genome.index()
    return genome


@pytest.fixture
def variant(local_genome):
    return Variant("1", 150, "A", "G", genome=local_genome)


@pytest.mark.parametrize("annotator", ["fast", "protein_diff", "transcript_model"])
def test_default_preserves_all_biotypes_and_incomplete_transcripts(variant, annotator):
    effects = variant.effects(annotator=annotator)
    assert {e.transcript_id for e in effects if e.transcript} == {
        "coding", "lnc", "pseudo", "nmd", "ig"}
    assert type(next(e for e in effects if e.transcript_id == "coding")) is IncompleteTranscript
    assert all(type(e) is NoncodingTranscript for e in effects
               if e.transcript_id in {"lnc", "pseudo", "nmd", "ig"})
    assert [(type(e), e.transcript_id) for e in effects] == [
        (type(e), e.transcript_id) for e in variant.effects(
            annotator=annotator, transcript_filter=lambda t: True)]


@pytest.mark.parametrize("entry", ["variant", "collection", "predict"])
def test_select_by_biotype_before_annotator_and_preserve_true_intragenic(variant, entry):
    calls = []
    annotator = SimpleNamespace(
        name="selection-test", version="1",
        annotate_on_transcript=lambda v, t: calls.append(t.id) or NoncodingTranscript(v, t))
    kwargs = dict(annotator=annotator, transcript_filter=lambda t: t.biotype == "lncRNA")
    if entry == "variant":
        result = variant.effects(**kwargs)
    elif entry == "collection":
        result = VariantCollection([variant]).effects(**kwargs)
    else:
        result = predict_variant_effects(variant, **kwargs)
    assert calls == ["lnc"]
    assert len(result) == 2
    assert next(e for e in result if e.transcript).transcript_id == "lnc"
    assert next(e for e in result if type(e) is Intragenic).gene.id == "g2"
    assert result.annotator == "selection-test" and result.annotator_version == "1"


def test_rejecting_every_transcript_does_not_invent_gene_only_effects(variant):
    result = variant.effects(transcript_filter=lambda t: False)
    assert len(result) == 1
    assert isinstance(result[0], Intragenic) and result[0].gene.id == "g2"
    # No genuine transcript-free gene overlap: the right answer is empty.
    variant._gene_ids = ["g1"]
    result = variant.effects(transcript_filter=lambda t: False)
    assert len(result) == 0 and result.annotator == "fast" and result.annotated_at


@pytest.mark.parametrize("position,expected", [(500, Intragenic), (1100, Intergenic)])
def test_non_transcript_overlaps_do_not_invoke_filter(local_genome, position, expected):
    def reject_unexpected_call(t):
        raise AssertionError("No transcript overlaps this variant")
    result = Variant("1", position, "A", "G", genome=local_genome).effects(
        transcript_filter=reject_unexpected_call)
    assert result and all(isinstance(e, expected) for e in result)


@pytest.mark.parametrize("entry", ["variant", "collection", "empty", "predict"])
@pytest.mark.parametrize("bad_filter", [False, "protein_coding", {"coding"}])
def test_non_callable_filter_is_configuration_error(variant, entry, bad_filter):
    calls = {
        "variant": variant.effects,
        "collection": VariantCollection([variant]).effects,
        "empty": VariantCollection([]).effects,
        "predict": lambda **kw: predict_variant_effects(variant, **kw),
    }
    with pytest.raises(TypeError, match="transcript_filter.*callable"):
        calls[entry](transcript_filter=bad_filter, raise_on_error=False)


def test_filter_errors_follow_annotation_error_policy(variant):
    def select(t):
        if t.id == "lnc":
            raise ValueError("selection failed")
        return False
    with pytest.raises(ValueError, match="selection failed"):
        variant.effects(transcript_filter=select)
    result = variant.effects(transcript_filter=select, raise_on_error=False)
    failure, = [e for e in result if isinstance(e, Failure)]
    assert failure.transcript_id == "lnc" and failure.error == "selection failed"


@pytest.mark.parametrize("annotator", ["fast", "protein_diff", "transcript_model"])
def test_pinned_coding_result_and_joint_phase_respect_selection(annotator):
    genome = cached_release(81)
    variants = VariantCollection([
        Variant("7", 117530914, "G", "T", genome),
        Variant("7", 117530915, "T", "A", genome),
    ])
    keep = lambda t: t.id == "ENST00000003084"
    resolver = SimpleNamespace(in_cis=lambda *args, **kwargs: True)
    expected = [e for e in variants.effects(annotator=annotator, phase_resolver=resolver)
                if e.transcript and keep(e.transcript)]
    result = variants.effects(annotator=annotator, phase_resolver=resolver,
                              transcript_filter=keep)
    assert len(result) == len(expected) == 3
    assert [(e.transcript_id, e.short_description, e.mutant_protein_sequence) for e in result] == [
        (e.transcript_id, e.short_description, e.mutant_protein_sequence) for e in expected]
    assert not variants.effects(annotator=annotator, phase_resolver=resolver,
                                transcript_filter=lambda t: False)


def test_structural_variants_use_the_same_selection(variant):
    sv = StructuralVariant("1", 150, "DEL", end=160, genome=variant.genome)
    result = sv.effects(transcript_filter=lambda t: t.id == "lnc")
    assert {e.transcript_id for e in result if e.transcript} == {"lnc"}
    assert not any(type(e) is Intragenic and e.gene.id == "g1" for e in result)


def test_only_selected_transcripts_reach_germline_and_phase_hooks(variant):
    context = GermlineContext.from_variants([
        Variant("1", 151, "A", "C", genome=variant.genome)])
    calls = []

    def annotate(v, t, germline_ctx, phase_resolver=None):
        calls.append((t.id, germline_ctx, phase_resolver))
        return NoncodingTranscript(v, t)

    annotator = SimpleNamespace(name="selection-test", version="1",
                               annotate_with_context=annotate)
    phase_calls = []
    phase = SimpleNamespace(mutant_transcript=lambda v, t: phase_calls.append(t.id))
    result = variant.effects(annotator=annotator, germline=context, phase_resolver=phase,
                              transcript_filter=lambda t: t.id == "lnc")
    assert calls == [("lnc", context, phase)]
    assert phase_calls == ["lnc"]
    assert {e.transcript_id for e in result if e.transcript} == {"lnc"}
