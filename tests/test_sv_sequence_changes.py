"""SV predictions must survive filtering without calling unknowns unchanged."""

from dataclasses import replace

import pytest
from pyensembl import cached_release

from varcode import (
    AssembledAllele, EffectCandidate, MutantTranscript, ReferenceSegment,
    StructuralVariant, make_fusion_outcome,
)
from varcode.effects import (
    CrypticExonCandidate, EffectCollection, GeneFusion, LargeDeletion,
    StructuralVariantEffect, TranslocationToIntergenic, Unresolved,
)
from varcode.effects import structural, structural_candidates, structural_fusions
from varcode.effects.codon_tables import codon_table_for_transcript
from varcode.nucleotides import reverse_complement

# Legacy wrapper contract; the default composition path has separate regression tests.
pytestmark = pytest.mark.usefixtures("legacy_annotator")


@pytest.fixture
def cftr():
    return cached_release(81).transcript_by_id("ENST00000003084")


@pytest.fixture
def isolate_span(monkeypatch):
    monkeypatch.setattr(structural_fusions, "fusion_partners", lambda *args: ())
    monkeypatch.setattr(
        structural_candidates, "enumerate_and_attach_cryptics", lambda *args: None)
    monkeypatch.setattr(
        structural_candidates, "enumerate_and_attach_splice_outcomes", lambda *args: None)


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
    assert effect.affected_exons
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


@pytest.mark.parametrize("completeness", [
    "partial_start", "partial_end", "partial_both", "unknown", None,
])
@pytest.mark.parametrize("protein", ["reference_fragment", "different", ""])
def test_unmapped_partial_proteins_do_not_establish_a_change(cftr, completeness, protein):
    if protein == "reference_fragment":
        protein = cftr.protein_sequence[:100]
    effect = _effect(cftr)
    effect.mutant_transcript = MutantTranscript.from_sequence(
        cftr.sequence[:432], reference_transcript=cftr,
        mutant_protein_sequence=protein,
        evidence={"protein_completeness": completeness})
    assert effect.modifies_protein_sequence is None
    assert effect.modifies_coding_sequence is None
    assert list(EffectCollection([effect]).drop_silent_and_noncoding()) == [effect]
    assert len(EffectCollection([effect]).drop_silent_and_noncoding(False)) == 0
    assert effect.mutant_protein_sequence == protein


@pytest.mark.parametrize("completeness", ["partial_start", "partial_end", "partial_both"])
def test_partial_label_prevents_complete_orf_fallback(cftr, completeness):
    # A complete-looking ORF must not override the producer's explicit partial
    # interpretation (e.g. an internal alternative start in an observed fragment).
    effect = _effect(cftr, cdna=cftr.sequence, segments=_identity_segments(cftr),
                     evidence={"protein_completeness": completeness})
    assert effect.modifies_protein_sequence is None
    assert effect.modifies_coding_sequence is None
    evidence = dict(effect.mutant_transcript.evidence,
                    cds_start=min(cftr.start_codon_spliced_offsets),
                    cds_end=max(cftr.stop_codon_spliced_offsets) + 1)
    effect.mutant_transcript = replace(effect.mutant_transcript, evidence=evidence)
    assert effect.modifies_protein_sequence is None
    assert effect.modifies_coding_sequence is None


def _mapped_partial(transcript, completeness, change="none"):
    """An observed CDS fragment mapped onto its reference, one codon changed.

    Flags compare observed codons, so no protein is supplied.
    """
    table = codon_table_for_transcript(transcript)
    start = min(transcript.start_codon_spliced_offsets)
    stop = max(transcript.stop_codon_spliced_offsets) + 1
    left = start if completeness == "partial_end" else start + 300
    right = stop if completeness == "partial_start" else left + 300
    cdna = transcript.sequence[left:right]
    if change == "stop_loss":
        pos, codon = len(cdna) - 3, "TGG"
    else:
        pos = next(i for i in range(3, len(cdna), 3)
                   if table.forward_table.get(cdna[i:i + 3], "M") not in "MW")
        original = cdna[pos:pos + 3]
        synonyms = [c for c, aa in sorted(table.forward_table.items())
                    if aa == table.forward_table[original] and c != original]
        codon = {"none": original, "synonymous": synonyms[0], "missense": "TGG",
                 "nonsense": "TAA", "ambiguous": original[:2] + "N"}[change]
    cdna = cdna[:pos] + codon + cdna[pos + 3:]
    return _effect(
        transcript, cdna=cdna,
        segments=(ReferenceSegment(transcript, left, right),),
        evidence={"protein_completeness": completeness,
                  "cds_start": 0, "cds_end": len(cdna)}), pos


@pytest.mark.parametrize("transcript_id", ["ENST00000003084", "ENST00000357654"])
@pytest.mark.parametrize("completeness", ["partial_start", "partial_end", "partial_both"])
@pytest.mark.parametrize("change", ["none", "synonymous", "missense", "nonsense", "ambiguous"])
def test_partial_observation_only_establishes_mapped_local_changes(
        transcript_id, completeness, change):
    transcript = cached_release(81).transcript_by_id(transcript_id)
    effect, _ = _mapped_partial(transcript, completeness, change)
    before = effect.mutant_transcript
    assert effect.modifies_coding_sequence is (
        True if change in ("synonymous", "missense", "nonsense") else None)
    assert effect.modifies_protein_sequence is (
        True if change in ("missense", "nonsense") else None)
    assert effect.mutant_transcript == before


@pytest.mark.parametrize("transcript_id", ["ENST00000003084", "ENST00000357654"])
def test_mapped_stop_loss_is_a_protein_change(transcript_id):
    transcript = cached_release(81).transcript_by_id(transcript_id)
    effect, _ = _mapped_partial(transcript, "partial_start", "stop_loss")
    assert effect.modifies_coding_sequence is True
    assert effect.modifies_protein_sequence is True


def test_mapped_reference_fragment_is_not_a_truncation(cftr):
    # The #462 report, with reference coordinates: equal observed codons and a
    # shorter supplied protein still leave the unobserved suffix unknown.
    effect, _ = _mapped_partial(cftr, "partial_end")
    effect.mutant_transcript = replace(
        effect.mutant_transcript, mutant_protein_sequence=cftr.protein_sequence[:100])
    assert effect.modifies_coding_sequence is None
    assert effect.modifies_protein_sequence is None
    assert len(EffectCollection([effect]).drop_silent_and_noncoding(False)) == 0


def test_partial_uses_mitochondrial_table():
    # TGA and TGG both encode tryptophan in vertebrate mitochondria.
    transcript = cached_release(81).transcript_by_id("ENST00000361624")
    effect, _ = _mapped_partial(transcript, "partial_both")
    model = effect.mutant_transcript
    cdna = model.cdna_sequence
    pos = next(i for i in range(0, len(cdna), 3) if cdna[i:i + 3] in ("TGA", "TGG"))
    swapped = "TGG" if cdna[pos:pos + 3] == "TGA" else "TGA"
    effect.mutant_transcript = replace(
        model, cdna_sequence=cdna[:pos] + swapped + cdna[pos + 3:])
    assert effect.modifies_coding_sequence is True
    assert effect.modifies_protein_sequence is None


@pytest.mark.parametrize("upstream", [False, True])
def test_alternative_start_codon_reads_as_met_only_where_the_orf_begins(cftr, upstream):
    # A 3' partner's start codon is internal to a fused ORF, so CTG there is Leu.
    brca1 = cached_release(81).transcript_by_id("ENST00000357654")
    start = min(cftr.start_codon_spliced_offsets)
    other = min(brca1.start_codon_spliced_offsets) + 3
    segments = (ReferenceSegment(cftr, start, start + 300),)
    cdna = "CTG" + cftr.sequence[start + 3:start + 300]
    if upstream:
        segments = (ReferenceSegment(brca1, other, other + 30),) + segments
        cdna = brca1.sequence[other:other + 30] + cdna
    effect = _effect(cftr, cdna=cdna, segments=segments, evidence={
        "protein_completeness": "partial_both", "cds_start": 0, "cds_end": len(cdna)})
    assert effect.modifies_coding_sequence is True
    assert effect.modifies_protein_sequence is (True if upstream else None)


def test_partial_bounds_accept_numpy_integers(cftr):
    import numpy as np

    effect, _ = _mapped_partial(cftr, "partial_both", "missense")
    evidence = effect.mutant_transcript.evidence
    effect.mutant_transcript = replace(effect.mutant_transcript, evidence=dict(
        evidence, cds_start=np.int64(0), cds_end=np.int64(evidence["cds_end"])))
    assert effect.modifies_coding_sequence is True
    assert effect.modifies_protein_sequence is True


def test_partial_mapped_change_survives_split_reference_segments(cftr):
    effect, pos = _mapped_partial(cftr, "partial_both", "missense")
    segment, = effect.mutant_transcript.reference_segments
    effect.mutant_transcript = replace(effect.mutant_transcript, reference_segments=(
        replace(segment, end=segment.start + pos + 1),
        replace(segment, start=segment.start + pos + 1)))
    assert effect.modifies_coding_sequence is True
    assert effect.modifies_protein_sequence is True


@pytest.mark.parametrize("missing", [
    "frame", "float_frame", "offset_frame", "mapping", "reverse", "length", "edits"])
def test_partial_local_change_requires_usable_coordinates(cftr, missing):
    from varcode import TranscriptEdit

    effect, _ = _mapped_partial(cftr, "partial_both", "missense")
    model = effect.mutant_transcript
    if missing == "frame":
        model = replace(model, evidence={"protein_completeness": "partial_both"})
    elif missing == "float_frame":
        model = replace(model, evidence=dict(model.evidence, cds_start=0.0))
    elif missing == "offset_frame":
        # Codons out of frame with the reference CDS have no counterpart.
        model = replace(model, evidence=dict(model.evidence, cds_start=1))
    elif missing == "mapping":
        model = replace(model, reference_segments=None)
    elif missing == "reverse":
        model = replace(model, reference_segments=(replace(model.reference_segments[0], strand="-"),))
    elif missing == "length":
        model = replace(model, cdna_sequence=model.cdna_sequence + "A")
    else:
        model = replace(model, edits=(TranscriptEdit(0, 0, "A"),))
    effect.mutant_transcript = model
    assert effect.modifies_coding_sequence is None
    assert effect.modifies_protein_sequence is None


@pytest.mark.parametrize("completeness", ["partial_start", "partial_end", "partial_both"])
def test_imported_partial_peptide_without_reference_mapping_is_unknown(cftr, completeness):
    # A generic supplied sequence has no native Exacto coordinate rows.
    variant = StructuralVariant(cftr.contig, cftr.start, "BND", genome=cftr.genome)
    start = min(cftr.start_codon_spliced_offsets)
    cdna = cftr.sequence[:start + 3] + "TGG" + cftr.sequence[start + 6:start + 300]
    candidate = make_fusion_outcome(
        variant, cftr, sequence=cdna, transcript_model_id="observed-model", source="exacto")
    model = candidate.effect.mutant_transcript
    candidate.effect.mutant_transcript = replace(
        model, mutant_protein_sequence="MW" + cftr.protein_sequence[2:100],
        evidence=dict(model.evidence or {}, protein_completeness=completeness,
                      cds_start=start, cds_end=start + 300))
    assert candidate.effect.modifies_coding_sequence is None
    assert candidate.effect.modifies_protein_sequence is None


@pytest.mark.parametrize("transcript_id", ["ENST00000003084", "ENST00000357654"])
@pytest.mark.parametrize("supplied_protein", [False, True])
@pytest.mark.parametrize("mapped", [False, True])
def test_unknown_transcript_start_with_complete_internal_orf(
        transcript_id, supplied_protein, mapped):
    transcript = cached_release(81).transcript_by_id(transcript_id)
    i = transcript.protein_sequence.index("M", 1)
    start = min(transcript.start_codon_spliced_offsets) + 3 * i
    stop = max(transcript.stop_codon_spliced_offsets) + 1
    cdna = transcript.sequence[start:stop]
    effect = _effect(
        transcript, cdna=cdna,
        protein=transcript.protein_sequence[i:] if supplied_protein else None,
        segments=(ReferenceSegment(transcript, start, stop),) if mapped else None,
        evidence={"protein_completeness": "start_to_stop", "cds_start": 0,
                  "cds_end": len(cdna), "sequence_status": "observed_model_completeness_unknown"})
    assert effect.modifies_coding_sequence is None
    assert effect.modifies_protein_sequence is None
    assert len(EffectCollection([effect]).drop_silent_and_noncoding(False)) == 0


@pytest.mark.parametrize("retained_start", [False, True])
def test_unknown_completeness_preserves_mapped_changes_and_annotated_start(cftr, retained_start):
    start = min(cftr.start_codon_spliced_offsets)
    stop = max(cftr.stop_codon_spliced_offsets) + 1
    internal = start + 3 * cftr.protein_sequence.index("M", 1)
    if retained_start:
        # The known initiator followed by a deletion can coincidentally produce
        # the same amino acids as an internal-Met reference suffix.
        segments = (ReferenceSegment(cftr, start, start + 3),
                    ReferenceSegment(cftr, internal + 3, stop))
    else:
        segments = (ReferenceSegment(cftr, internal, stop),)
    cdna = "".join(cftr.sequence[s.start:s.end] for s in segments)
    if not retained_start:
        # A mapped synonymous change still establishes a CDS change even
        # though the shorter protein remains compatible with missing coverage.
        pos = next(i for i in range(3, len(cdna), 3) if cdna[i:i + 3] == "GCT")
        cdna = cdna[:pos] + "GCC" + cdna[pos + 3:]
    effect = _effect(cftr, cdna=cdna, segments=segments,
                     protein=cftr.protein_sequence[(internal - start) // 3:],
                     evidence={"protein_completeness": "start_to_stop", "cds_start": 0,
                               "cds_end": len(cdna),
                               "sequence_status": "observed_model_completeness_unknown"})
    assert effect.modifies_coding_sequence is True
    assert effect.modifies_protein_sequence is (True if retained_start else None)


def test_reference_suffix_without_cdna_needs_transcript_completeness(cftr):
    protein = cftr.protein_sequence[cftr.protein_sequence.index("M", 1):]
    complete = _effect(cftr, protein=protein)
    assert complete.modifies_protein_sequence is True
    unknown = _effect(cftr, protein=protein, evidence={
        "sequence_status": "observed_model_completeness_unknown"})
    assert unknown.modifies_protein_sequence is None
    assert unknown.modifies_coding_sequence is None


def test_unrendered_mapping_cannot_establish_an_observed_initiator(cftr):
    i = cftr.protein_sequence.index("M", 1)
    sequence = cftr.coding_sequence[3 * i:]
    start = min(cftr.start_codon_spliced_offsets)
    stop = max(cftr.stop_codon_spliced_offsets) + 1
    effect = _effect(cftr, cdna=sequence, protein=cftr.protein_sequence[i:],
                     # This stale mapping still includes the missing prefix.
                     segments=(ReferenceSegment(cftr, start, stop),), evidence={
                         "sequence_status": "observed_model_completeness_unknown",
                         "cds_start": 0, "cds_end": len(sequence)})
    assert effect.modifies_coding_sequence is None
    assert effect.modifies_protein_sequence is None


def test_reference_cds_suffix_without_reference_protein_stays_unknown():
    from types import SimpleNamespace

    transcript = SimpleNamespace(
        id="missing-protein", name="missing-protein", gene=SimpleNamespace(id="gene"),
        protein_sequence=None, coding_sequence="ATGGCTATGAAATAA", complete=False,
        contig="1", start=100, genome=cached_release(81))
    effect = _effect(transcript, cdna="ATGAAATAA", evidence={
        "sequence_status": "observed_model_completeness_unknown",
        "cds_start": 0, "cds_end": 9})
    assert effect.modifies_protein_sequence is None
    assert effect.modifies_coding_sequence is None


def test_complete_protein_and_partial_alternative_remain_unknown(cftr):
    primary = _effect(cftr, protein=cftr.protein_sequence,
                      evidence={"protein_completeness": "start_to_stop"})
    partial, _ = _mapped_partial(cftr, "partial_both", "none")
    primary._extra_candidates = (EffectCandidate(partial, source="rna"),)
    assert primary.modifies_protein_sequence is None
    assert len(EffectCollection([primary]).drop_silent_and_noncoding()) == 1
    assert len(EffectCollection([primary]).drop_silent_and_noncoding(False)) == 0
    changed, _ = _mapped_partial(cftr, "partial_both", "missense")
    primary._extra_candidates += (EffectCandidate(changed, source="rna"),)
    assert primary.modifies_protein_sequence is True
    assert len(EffectCollection([primary]).drop_silent_and_noncoding(False)) == 1


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
    model = structural.build_fusion_mutant_transcript(
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


def _render_segments(segments):
    pieces = []
    for seg in segments:
        seq = seg.source.sequence[seg.start:seg.end]
        pieces.append(reverse_complement(seq) if seg.strand == "-" else seq)
    return "".join(pieces)


@pytest.mark.parametrize("pair", [
    ("ENST00000003084", "ENST00000357654"),
    ("ENST00000357654", "ENST00000003084"),
])
@pytest.mark.parametrize("tail_kind", ["partner", "unmapped", "split_partner"])
def test_partial_junction_into_other_sequence_establishes_change(pair, tail_kind):
    genome = cached_release(81)
    first = genome.transcript_by_id(pair[0])
    second = genome.transcript_by_id(pair[1])
    start1 = min(first.start_codon_spliced_offsets)
    start2 = min(second.start_codon_spliced_offsets)
    prefix = (
        (ReferenceSegment(first, start1, start1 + 4),
         ReferenceSegment(first, start1 + 4, start1 + 90))
        if tail_kind == "split_partner"
        else (ReferenceSegment(first, start1, start1 + 90),)
    )
    if tail_kind == "unmapped":
        novel = AssembledAllele("TGG" * 10)
        segments = prefix + (ReferenceSegment(novel, 0, len(novel.sequence)),)
    else:
        segments = prefix + (ReferenceSegment(second, start2 + 30, start2 + 90),)
    cdna = _render_segments(segments)
    effect = _effect(first, cdna=cdna, segments=segments, evidence={
        "protein_completeness": "partial_end", "cds_start": 0, "cds_end": len(cdna),
    })
    assert effect.modifies_coding_sequence is True
    assert effect.modifies_protein_sequence is True
    assert len(EffectCollection([effect]).drop_silent_and_noncoding(False)) == 1


def test_partial_junction_with_identical_observed_codons_stays_unresolved(cftr):
    start = min(cftr.start_codon_spliced_offsets)
    same_tail = AssembledAllele(cftr.sequence[start + 30:start + 60])
    segments = (
        ReferenceSegment(cftr, start, start + 30),
        ReferenceSegment(same_tail, 0, len(same_tail.sequence)),
    )
    cdna = _render_segments(segments)
    effect = _effect(cftr, cdna=cdna, segments=segments, evidence={
        "protein_completeness": "partial_end", "cds_start": 0, "cds_end": len(cdna),
    })
    assert effect.modifies_coding_sequence is None
    assert effect.modifies_protein_sequence is None


@pytest.mark.parametrize("transcript_id", ["ENST00000003084", "ENST00000357654"])
@pytest.mark.parametrize("kind", [
    "ins1", "del1", "inframe_del", "pre_anchor_frameshift", "inverted_tail", "stop_del",
])
def test_partial_indels_frameshifts_and_inversions_on_both_strands(transcript_id, kind):
    transcript = cached_release(81).transcript_by_id(transcript_id)
    start = min(transcript.start_codon_spliced_offsets)
    stop = max(transcript.stop_codon_spliced_offsets) + 1
    completeness = "partial_both"
    if kind == "ins1":
        ins = AssembledAllele("A")
        segments = (
            ReferenceSegment(transcript, start + 30, start + 60),
            ReferenceSegment(ins, 0, 1),
            ReferenceSegment(transcript, start + 60, start + 120),
        )
    elif kind == "del1":
        segments = (
            ReferenceSegment(transcript, start + 30, start + 60),
            ReferenceSegment(transcript, start + 61, start + 120),
        )
    elif kind == "inframe_del":
        segments = (
            ReferenceSegment(transcript, start + 30, start + 60),
            ReferenceSegment(transcript, start + 69, start + 120),
        )
    elif kind == "pre_anchor_frameshift":
        ins = AssembledAllele("A")
        segments = (
            ReferenceSegment(transcript, start + 31, start + 60),
            ReferenceSegment(ins, 0, 1),
            ReferenceSegment(transcript, start + 60, start + 90),
        )
    elif kind == "inverted_tail":
        segments = (
            ReferenceSegment(transcript, start + 30, start + 60),
            ReferenceSegment(transcript, start + 60, start + 120, strand="-"),
        )
    else:
        completeness = "partial_start"
        # Delete the reference stop codon so translation reads into non-stop 3' sequence.
        tail = AssembledAllele("TGG" * 5)
        segments = (
            ReferenceSegment(transcript, stop - 33, stop - 3),
            ReferenceSegment(tail, 0, len(tail.sequence)),
        )
    cdna = _render_segments(segments)
    effect = _effect(transcript, cdna=cdna, segments=segments, evidence={
        "protein_completeness": completeness, "cds_start": 0, "cds_end": len(cdna),
    })
    assert effect.modifies_coding_sequence is True
    assert effect.modifies_protein_sequence is True


def test_partial_bases_after_observed_stop_cannot_anchor_frame(cftr):
    start = min(cftr.start_codon_spliced_offsets)
    stop_prefix = AssembledAllele("CCCTAA")
    segments = (
        ReferenceSegment(stop_prefix, 0, len(stop_prefix.sequence)),
        ReferenceSegment(cftr, start + 30, start + 90),
    )
    cdna = _render_segments(segments)
    effect = _effect(cftr, cdna=cdna, segments=segments, evidence={
        "protein_completeness": "partial_both", "cds_start": 0, "cds_end": len(cdna),
    })
    assert effect.modifies_coding_sequence is None
    assert effect.modifies_protein_sequence is None


@pytest.mark.parametrize("transcript_id", ["ENST00000003084", "ENST00000357654"])
def test_partial_upstream_5utr_atg_does_not_establish_initiation(transcript_id):
    # An observed upstream ATG in the 5' UTR (k < 0) in a partial_end fragment
    # does not establish initiation on the same transcript (#466).
    transcript = cached_release(81).transcript_by_id(transcript_id)
    start = min(transcript.start_codon_spliced_offsets)
    upstream_atg = AssembledAllele("ATGGGGCCC")
    segments = (
        ReferenceSegment(upstream_atg, 0, len(upstream_atg.sequence)),
        ReferenceSegment(transcript, start, start + 90),
    )
    cdna = _render_segments(segments)
    effect = _effect(transcript, cdna=cdna, segments=segments, evidence={
        "protein_completeness": "partial_end", "cds_start": 0, "cds_end": len(cdna),
    })
    assert effect.modifies_coding_sequence is None
    assert effect.modifies_protein_sequence is None
    # A real coding change inside the downstream CDS (k >= 0) is still detected.
    mutated = cdna[:12] + "TGG" + cdna[15:]
    effect.mutant_transcript = replace(effect.mutant_transcript, cdna_sequence=mutated)
    assert effect.modifies_coding_sequence is True
    assert effect.modifies_protein_sequence is True


@pytest.mark.parametrize("pair", [
    ("ENST00000003084", "ENST00000357654"),
    ("ENST00000357654", "ENST00000003084"),
])
@pytest.mark.parametrize("junction_kind", [
    "inframe_internal", "inframe_at_start", "out_of_frame", "out_of_frame_split",
])
def test_partial_fusion_three_prime_partner_cases(pair, junction_kind):
    genome = cached_release(81)
    five = genome.transcript_by_id(pair[0])
    three = genome.transcript_by_id(pair[1])
    start5 = min(five.start_codon_spliced_offsets)
    start3 = min(three.start_codon_spliced_offsets)
    if junction_kind == "inframe_internal":
        segments = (
            ReferenceSegment(five, start5, start5 + 30),
            ReferenceSegment(three, start3 + 30, start3 + 90),
        )
    elif junction_kind == "inframe_at_start":
        # 5' partner CDS drives translation directly into unmutated 3' start (k < 0 -> k = 0).
        segments = (
            ReferenceSegment(five, start5, start5 + 30),
            ReferenceSegment(three, start3, start3 + 90),
        )
    elif junction_kind == "out_of_frame":
        # 5' partner contributes 31 bases, putting all 3' partner codons out of frame.
        segments = (
            ReferenceSegment(five, start5, start5 + 31),
            ReferenceSegment(three, start3 + 30, start3 + 90),
        )
    else:
        segments = (
            ReferenceSegment(five, start5, start5 + 4),
            ReferenceSegment(five, start5 + 4, start5 + 31),
            ReferenceSegment(three, start3 + 30, start3 + 90),
        )
    cdna = _render_segments(segments)
    variant = StructuralVariant(
        five.contig, five.start, "BND", mate_contig=three.contig, mate_start=three.start,
        genome=genome)
    model = MutantTranscript(
        reference_transcript=five, reference_segments=segments, cdna_sequence=cdna,
        evidence={"protein_completeness": "partial_end", "cds_start": 0, "cds_end": len(cdna)})
    effect_5p = GeneFusion(variant, five, three, mutant_transcript=model,
                           five_prime_transcript=five, three_prime_transcript=three)
    effect_3p = GeneFusion(variant, three, five, mutant_transcript=model,
                           five_prime_transcript=five, three_prime_transcript=three)
    for effect in (effect_5p, effect_3p):
        assert effect.modifies_coding_sequence is True
        assert effect.modifies_protein_sequence is True
        assert len(EffectCollection([effect]).drop_silent_and_noncoding(False)) == 1


def test_partial_fusion_three_prime_partner_stays_conservative_when_unanchored_or_ambiguous(cftr):
    genome = cached_release(81)
    brca1 = genome.transcript_by_id("ENST00000357654")
    start5 = min(brca1.start_codon_spliced_offsets)
    start3 = min(cftr.start_codon_spliced_offsets)
    variant = StructuralVariant(
        brca1.contig, brca1.start, "BND", mate_contig=cftr.contig, mate_start=cftr.start,
        genome=genome)

    # Shifted frame out of frame on BOTH partners (brca1 at 0 mod 3, cftr at 1 mod 3 -> cds_start=2).
    segments = (
        ReferenceSegment(brca1, start5, start5 + 31),
        ReferenceSegment(cftr, start3 + 30, start3 + 90),
    )
    cdna = _render_segments(segments)
    unanchored = MutantTranscript(
        reference_transcript=brca1, reference_segments=segments, cdna_sequence=cdna,
        evidence={"protein_completeness": "partial_both", "cds_start": 2, "cds_end": len(cdna)})
    eff_unanchored = GeneFusion(variant, cftr, brca1, mutant_transcript=unanchored,
                                five_prime_transcript=brca1, three_prime_transcript=cftr)
    assert eff_unanchored.modifies_coding_sequence is None
    assert eff_unanchored.modifies_protein_sequence is None

    # Ambiguous N-containing codon upstream of k=0 cannot establish N-terminal extension.
    amb_seg = (
        ReferenceSegment(brca1, start5, start5 + 27),
        ReferenceSegment(AssembledAllele("NNN"), 0, 3),
        ReferenceSegment(cftr, start3, start3 + 90),
    )
    amb_cdna = _render_segments(amb_seg)
    amb_model = MutantTranscript(
        reference_transcript=brca1, reference_segments=amb_seg, cdna_sequence=amb_cdna,
        evidence={"protein_completeness": "partial_end", "cds_start": 0, "cds_end": len(amb_cdna)})
    eff_amb = GeneFusion(variant, cftr, brca1, mutant_transcript=amb_model,
                         five_prime_transcript=brca1, three_prime_transcript=cftr)
    assert eff_amb.modifies_coding_sequence is None
    assert eff_amb.modifies_protein_sequence is None

    # Out-of-frame fusion starting with an ambiguous codon stays unresolved.
    amb_oof_seg = (
        ReferenceSegment(AssembledAllele("NNN"), 0, 3),
        ReferenceSegment(brca1, start5 + 3, start5 + 31),
        ReferenceSegment(cftr, start3 + 30, start3 + 90),
    )
    amb_oof_cdna = _render_segments(amb_oof_seg)
    amb_oof_model = MutantTranscript(
        reference_transcript=brca1, reference_segments=amb_oof_seg, cdna_sequence=amb_oof_cdna,
        evidence={"protein_completeness": "partial_both", "cds_start": 0,
                  "cds_end": len(amb_oof_cdna)})
    eff_amb_oof = GeneFusion(variant, cftr, brca1, mutant_transcript=amb_oof_model,
                             five_prime_transcript=brca1, three_prime_transcript=cftr)
    assert eff_amb_oof.modifies_coding_sequence is None
    assert eff_amb_oof.modifies_protein_sequence is None

    # Short 5'-anchored fragment matching a codon inside the 3' partner stays unresolved.
    short_seg = (ReferenceSegment(brca1, start5, start5 + 3),)  # ATG, also in CFTR
    short_cdna = _render_segments(short_seg)
    short_model = MutantTranscript(
        reference_transcript=brca1, reference_segments=short_seg, cdna_sequence=short_cdna,
        evidence={"protein_completeness": "partial_end", "cds_start": 0, "cds_end": 3})
    eff_short = GeneFusion(variant, cftr, brca1, mutant_transcript=short_model,
                           five_prime_transcript=brca1, three_prime_transcript=cftr)
    assert eff_short.modifies_coding_sequence is None
    assert eff_short.modifies_protein_sequence is None
