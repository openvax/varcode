"""Completion assumptions never replace an RNA observation or its ORF."""

from copy import deepcopy
from dataclasses import replace
from types import SimpleNamespace

import pytest

from varcode import EffectCandidate, ReferenceSegment, reference_completion_hypotheses
from varcode.effects import GeneFusion
from varcode.mutant_transcript import _AssembledAllele

from . import test_exacto_reference_mapping as mapping
from .test_exacto_reference_mapping import load, rows_for, transcript

context = mapping.context


def fragment(context, start=6, end=18, translate=False):
    _, first, second = context
    return load(context, rows_for([
        (first, start, 9, "match", None), (second, 9, end, "match", None)]), translate=translate)


@pytest.mark.parametrize("ends,expected,span,completeness", [
    (("five_prime", "three_prime"), "ATGAAACCCGGGAAATTTCCCTAA", [6, 18], "start_to_stop"),
    (("five_prime",), "ATGAAACCCGGGAAATTT", [6, 18], "partial_end"),
    (("three_prime",), "CCCGGGAAATTTCCCTAA", [0, 12], "unknown"),
])
def test_completion_preserves_observed_bases_and_marks_assumptions(context, ends, expected, span, completeness):
    candidate = fragment(context, translate=True)
    original = candidate.effect.mutant_transcript
    evidence = deepcopy(candidate.evidence)
    hypothesis, = reference_completion_hypotheses(candidate, ends=ends)
    model = hypothesis.effect.mutant_transcript
    assert model.cdna_sequence == expected
    assert model.cdna_sequence[slice(*span)] == original.cdna_sequence
    assert hypothesis.evidence["reference_completion"]["observed_span"] == span
    assert hypothesis.evidence["protein_completeness"] == completeness
    assert hypothesis.source == model.annotator_name == "varcode_reference_completion"
    assert hypothesis.evidence == model.evidence
    assert candidate.effect.mutant_transcript is original
    assert candidate.evidence == evidence
    assert original.cdna_sequence == "CCCGGGAAATTT"
    assert original.mutant_protein_sequence == "PGKF"
    assert model.mutant_protein_sequence == {
        "start_to_stop": "MKPGKFP", "partial_end": "MKPGKF", "unknown": None}[completeness]
    if completeness != "unknown":
        assert hypothesis.evidence["observed_orf_frame"] == "in_frame"
        assert hypothesis.evidence["initiation_status"] == "assumed_reference_start"
    for assumed in hypothesis.evidence["reference_completion"]["assumed_spans"]:
        assert assumed["end"] <= span[0] or assumed["start"] >= span[1]
    assert "".join(s.source.sequence[s.start:s.end] for s in model.reference_segments) == expected


def test_no_added_bases_or_requested_ends_means_no_duplicate(context):
    _, first, second = context
    candidate = load(context, rows_for([(first, 0, 9, "match", None),
                                      (second, 9, 24, "match", None)]))
    assert reference_completion_hypotheses(candidate) == ()
    assert reference_completion_hypotheses(fragment(context), ends=()) == ()


@pytest.mark.parametrize("end", ["five_prime", "three_prime"])
def test_real_rna_end_evidence_prevents_extension_but_read_ends_do_not(context, end):
    candidate = fragment(context)
    candidate = replace(candidate, evidence=dict(candidate.evidence, **{
        "rna_" + end + "_complete": True, "read_count": 12, "read_end": 200}))
    result, = reference_completion_hypotheses(candidate)
    assumed = result.evidence["reference_completion"]["assumed_spans"]
    assert len(assumed) == 1
    assert assumed[0]["label"] != "assumed_reference_" + end
    assert "read_count" not in result.evidence
    assert result.evidence["observed_evidence"]["read_count"] == 12
    result.evidence["observed_evidence"]["exacto_structure"][0]["sequence"] = "changed"
    assert candidate.evidence["exacto_structure"][0]["sequence"] != "changed"


def test_model_end_evidence_is_respected(context):
    candidate = fragment(context)
    model = candidate.effect.mutant_transcript
    candidate.effect.mutant_transcript = replace(model, evidence=dict(
        model.evidence, rna_five_prime_complete=True, rna_three_prime_complete=True))
    assert reference_completion_hypotheses(candidate) == ()


def alternate(source, identifier, prepend="", append=""):
    """Explicit extra terminal exons; preserve the observed exon coordinates."""
    positions = list(source.positions)
    if prepend:
        positions = (list(range(20, 20 + len(prepend))) if source.strand == "+" else
                     list(range(400 + len(prepend) - 1, 399, -1))) + positions
    if append:
        positions += (list(range(400, 400 + len(append))) if source.strand == "+" else
                      list(range(20 + len(append) - 1, 19, -1)))
    return SimpleNamespace(**dict(vars(source), id=identifier,
        sequence=prepend + source.sequence + append,
        positions=positions, spliced_offset=positions.index,
        start_codon_spliced_offsets=[p + len(prepend) for p in source.start_codon_spliced_offsets],
        stop_codon_spliced_offsets=[p + len(prepend) for p in source.stop_codon_spliced_offsets]))


def test_isoform_pairs_remain_distinct_even_with_identical_proteins(context):
    _, first, second = context
    upstream = alternate(first, "T1_alt", prepend="CCA")
    downstream = alternate(second, "T2_alt", append="GGG")
    results = reference_completion_hypotheses(
        fragment(context), five_prime_transcripts=[first, upstream, upstream],
        three_prime_transcripts=[second, downstream])
    assert len(results) == 4
    assert {r.effect.mutant_protein_sequence for r in results} == {"MKPGKFP"}
    assert {(r.effect.transcript.id, r.effect.partner_transcript.id) for r in results} == {
        ("T1", "T2"), ("T1_alt", "T2"), ("T1", "T2_alt"), ("T1_alt", "T2_alt")}
    for result in results:
        model = result.effect.mutant_transcript
        assert "".join(s.source.sequence[s.start:s.end] for s in model.reference_segments) == model.cdna_sequence
        assert result.evidence["reference_completion"]["observed_span"] == (
            [9, 21] if result.effect.transcript.id == "T1_alt" else [6, 18])


@pytest.mark.parametrize("bad", ["gene", "strand", "contig", "genome", "sequence", "splice", "missing"])
def test_incompatible_alternatives_are_not_completed(context, bad):
    _, first, _ = context
    other = alternate(first, "bad")
    if bad in ("gene", "strand", "contig", "genome"):
        setattr(other, "gene_id" if bad == "gene" else bad, "wrong")
    elif bad == "sequence":
        other.sequence = first.sequence[:7] + "T" + first.sequence[8:]
    elif bad == "splice":
        # Add an exon between two observed runs; terminal matching alone
        # would accept it, but the observed splice path rules it out.
        other.spliced_offset = lambda p: first.spliced_offset(p) + (3 if first.spliced_offset(p) >= 7 else 0)
    else:
        other.sequence = None
    assert reference_completion_hypotheses(fragment(context), five_prime_transcripts=[other]) == ()


@pytest.mark.parametrize("side", ["five_prime", "three_prime"])
def test_unmapped_terminal_bases_prevent_only_that_end_completion(context, side):
    candidate = fragment(context)
    model = candidate.effect.mutant_transcript
    segments = list(model.reference_segments)
    index = 0 if side == "five_prime" else -1
    old = segments[index]
    sequence = old.source.sequence[old.start:old.end]
    segments[index] = ReferenceSegment(_AssembledAllele(sequence), 0, len(sequence))
    candidate.effect.mutant_transcript = replace(model, reference_segments=tuple(segments))
    result, = reference_completion_hypotheses(candidate)
    spans = result.evidence["reference_completion"]["assumed_spans"]
    assert len(spans) == 1 and spans[0]["label"] != "assumed_reference_" + side


def test_split_annotated_start_codon_is_retained(context):
    candidate = fragment(context, start=1)
    result, = reference_completion_hypotheses(candidate)
    assert result.effect.mutant_protein_sequence == "MKPGKFP"
    assert result.evidence["cds_start"] == 0


def test_out_of_frame_internal_atg_is_a_separate_producer_hypothesis(context):
    from varcode import load_exacto_fusions
    from .test_exacto_fusions import link, table
    _, first, second = context
    first.sequence = first.coding_sequence = "ATGAAACATGAATTTCCCAAATAA"
    first.protein_sequence = "MKHEFPK"
    rows = rows_for([(first, 6, 12, "match", None), (second, 12, 18, "match", None)])
    sequence = "CATGAAAAATTT"
    primary = []
    origins = [r["index"] for r in rows for _ in r["sequence"]]
    for position in range(1, 10):
        offset = position - 1
        primary.append(dict(peptide_id="internal", primary_structure_index=str(offset),
            type="base", amino_acid="MKN"[offset // 3], amino_acid_index=str(offset // 3),
            codon_index=str(offset % 3), nucleotide=sequence[position], transcript_model_id="1",
            reference_transcript_ids="T1,T2", transcript_structure_index=origins[position],
            read_start=str(position), read_end=str(position)))
    candidate, = load_exacto_fusions(table(rows), table([link()]),
        variants_by_id={"D1": context[0]}, primary_structures_path=table(primary)).candidates
    result, = reference_completion_hypotheses(candidate)
    assert candidate.effect.mutant_protein_sequence == "MKN"
    assert result.evidence["observed_protein_sequence"] == "MKN"
    assert result.effect.mutant_protein_sequence == "MKHEKFP"
    assert result.evidence["observed_orf_frame"] == "out_of_frame"


def test_unchanged_completed_sequence_does_not_reclassify_observed_fragment(context):
    _, first, second = context
    first.sequence = first.coding_sequence = second.sequence
    first.protein_sequence = second.protein_sequence
    candidate = fragment(context, translate=True)
    assert candidate.effect.modifies_protein_sequence is None
    result, = reference_completion_hypotheses(candidate)
    assert result.effect.modifies_protein_sequence is False
    assert result.effect.modifies_coding_sequence is False
    assert candidate.effect.modifies_protein_sequence is None


def test_inframe_internal_atg_suffix_remains_an_observation_alongside_completion(context):
    _, first, second = context
    for t in (first, second):
        t.sequence = t.coding_sequence = "ATGAAAATGGGGTTTCCCAAATAA"
        t.protein_sequence = "MKMGFPK"
    candidate = fragment(context, end=24, translate=True)
    assert candidate.effect.mutant_protein_sequence == "MGFPK"
    assert candidate.evidence["protein_completeness"] == "start_to_stop"
    assert candidate.effect.modifies_protein_sequence is None
    result, = reference_completion_hypotheses(candidate)
    assert result.effect.mutant_protein_sequence == "MKMGFPK"
    assert result.effect.modifies_protein_sequence is False
    assert result.evidence["observed_orf_frame"] == "in_frame"
    assert candidate.effect.mutant_protein_sequence == "MGFPK"
    assert candidate.effect.modifies_protein_sequence is None


def test_observed_insertion_is_never_filled_or_removed(context):
    _, first, second = context
    candidate = load(context, rows_for([(first, 6, 9, "match", None),
        (first, 9, 9, "insertion", "A"), (first, 9, 15, "match", None),
        (second, 15, 18, "match", None)]))
    result, = reference_completion_hypotheses(candidate)
    assert result.effect.mutant_transcript.cdna_sequence == "ATGAAACCCAGGGTTTTTTCCCTAA"
    # ATG AAA CCC AGG GTT TTT TCC CTA, followed by one incomplete base.
    assert result.effect.mutant_protein_sequence == "MKPRVFSL"
    assert result.evidence["protein_completeness"] == "partial_end"
    assert result.effect.modifies_protein_sequence is True


@pytest.mark.parametrize("mutation", ["layout", "end_flag", "recursive", "ends"])
def test_invalid_or_recursive_completion_is_rejected(context, mutation):
    candidate = fragment(context)
    options = {}
    if mutation == "layout":
        candidate.effect.mutant_transcript = replace(candidate.effect.mutant_transcript, cdna_sequence="ACGT")
    elif mutation == "end_flag":
        candidate = replace(candidate, evidence=dict(candidate.evidence, rna_five_prime_complete="yes"))
    elif mutation == "recursive":
        candidate, = reference_completion_hypotheses(candidate)
    else:
        options["ends"] = "five_prime"
    with pytest.raises(ValueError):
        reference_completion_hypotheses(candidate, **options)


def test_selenocysteine_uses_retained_reference_utr_in_the_hypothesis():
    # UGA is Sec only under the existing mapped-UTR policy. The assumed UTR
    # is visible as reference completion, never claimed as observed coverage.
    first = transcript("T1", "+", "1", "ATGAAACCCGGGTTTCCCAAATAA")
    second = transcript("T2", "+", "2", "ATGTGAAAATAACCC")
    second.protein_sequence = "MUK"
    second.coding_sequence = second.sequence[:12]
    second.stop_codon_spliced_offsets = [9, 10, 11]
    from varcode import MutantTranscript
    segments = (ReferenceSegment(first, 6, 9), ReferenceSegment(second, 3, 9))
    model = MutantTranscript(reference_transcript=first, reference_segments=segments,
                             cdna_sequence="CCCTGAAAA")
    candidate = EffectCandidate(GeneFusion(None, first, second, mutant_transcript=model), source="rna")
    result, = reference_completion_hypotheses(candidate)
    assert result.effect.mutant_protein_sequence == "MKPUK"
    assert result.evidence["selenocysteine_decoded_offsets"] == [9]
    result, = reference_completion_hypotheses(candidate, ends=("five_prime",))
    assert result.effect.mutant_protein_sequence == "MKP"
    partial_utr = replace(model,
        reference_segments=(segments[0], ReferenceSegment(second, 3, 13)),
        cdna_sequence="CCCTGAAAATAAC")
    partial = replace(candidate, effect=GeneFusion(None, first, second, mutant_transcript=partial_utr))
    result, = reference_completion_hypotheses(partial, ends=("five_prime",))
    assert result.effect.mutant_protein_sequence == "MKPUK"
    assert result.evidence["selenocysteine_uncertain_offsets"] == [9]


def test_unannotated_partner_still_allows_five_prime_completion(context):
    from varcode.effects import TranslocationToIntergenic
    candidate = fragment(context)
    model = candidate.effect.mutant_transcript
    candidate = replace(candidate, effect=TranslocationToIntergenic(
        candidate.effect.variant, context[1], mutant_transcript=model))
    result, = reference_completion_hypotheses(candidate)
    assert isinstance(result.effect, TranslocationToIntergenic)
    assert result.effect.mutant_transcript.cdna_sequence == "ATGAAACCCGGGAAATTT"
    assert result.evidence["protein_completeness"] == "partial_end"
    assert result.evidence["reference_completion"]["three_prime_transcript_id"] is None


def test_ambiguous_observed_bases_remain_ambiguous(context):
    _, first, second = context
    candidate = load(context, rows_for([(first, 6, 9, "match", None),
        (first, 9, 12, "mismatch", "NNN"), (second, 12, 18, "match", None)]))
    result, = reference_completion_hypotheses(candidate)
    assert result.effect.mutant_transcript.cdna_sequence == "ATGAAACCCNNNAAATTTCCCTAA"
    assert result.effect.mutant_protein_sequence is None
    assert result.evidence["protein_status"] == "unresolved_ambiguous_coding_sequence"
    assert result.effect.modifies_protein_sequence is None


def test_partial_end_request_accepts_an_iterator(context):
    result, = reference_completion_hypotheses(fragment(context), ends=iter(["five_prime"]))
    assert result.effect.mutant_protein_sequence == "MKPGKF"


def test_original_three_prime_partner_perspective_uses_the_five_prime_start(context):
    candidate = fragment(context)
    _, first, second = context
    candidate = replace(candidate, effect=GeneFusion(context[0], second, first,
        mutant_transcript=candidate.effect.mutant_transcript,
        five_prime_transcript=first, three_prime_transcript=second))
    result, = reference_completion_hypotheses(candidate)
    assert result.effect.five_prime_transcript is first
    assert result.effect.three_prime_transcript is second
    assert result.effect.mutant_protein_sequence == "MKPGKFP"


def test_isoform_with_a_different_upstream_cds_generates_its_own_protein(context):
    _, first, _ = context
    upstream = alternate(first, "T1_long", prepend="ATG")
    upstream.start_codon_spliced_offsets = [0, 1, 2]
    upstream.coding_sequence = upstream.sequence
    upstream.protein_sequence = "M" + first.protein_sequence
    results = reference_completion_hypotheses(fragment(context), five_prime_transcripts=[first, upstream])
    assert [r.effect.mutant_protein_sequence for r in results] == ["MKPGKFP", "MMKPGKFP"]


@pytest.mark.parametrize("starts", [[], [0, 1], [0, 2, 3]])
def test_missing_or_partial_start_annotation_does_not_guess_an_orf(context, starts):
    candidate = fragment(context)
    context[1].start_codon_spliced_offsets = starts
    result, = reference_completion_hypotheses(candidate)
    assert result.effect.mutant_protein_sequence is None
    assert result.evidence["protein_completeness"] == "unknown"


def test_pyensembl_missing_start_annotation_is_an_unresolved_protein(context):
    class Noncoding(SimpleNamespace):
        @property
        def start_codon_spliced_offsets(self):
            raise ValueError("No start_codon features found")

    _, first, _ = context
    candidate = fragment(context)
    noncoding = Noncoding(**{k: v for k, v in vars(first).items() if k != "start_codon_spliced_offsets"})
    result, = reference_completion_hypotheses(candidate, five_prime_transcripts=[noncoding])
    assert result.effect.mutant_protein_sequence is None


def test_observed_stop_is_not_overridden_by_three_prime_completion(context):
    _, first, second = context
    candidate = load(context, rows_for([(first, 6, 9, "match", None),
        (first, 9, 12, "mismatch", "TAA"), (second, 12, 18, "match", None)]))
    result, = reference_completion_hypotheses(candidate)
    assert result.effect.mutant_protein_sequence == "MKP"
    assert result.evidence["cds_end"] == 12
    assert result.evidence["protein_completeness"] == "start_to_stop"


def test_real_transcripts_and_candidate_json_roundtrip():
    from pyensembl import cached_release
    from varcode import MutantTranscript, StructuralVariant
    genome = cached_release(81)
    first = genome.transcript_by_id("ENST00000003084")
    second = genome.transcript_by_id("ENST00000357654")
    left = min(first.start_codon_spliced_offsets) + 6
    right = min(second.start_codon_spliced_offsets) + 9
    observed = first.sequence[left:left + 9] + second.sequence[right:right + 9]
    model = MutantTranscript(reference_transcript=first,
        reference_segments=(ReferenceSegment(first, left, left + 9), ReferenceSegment(second, right, right + 9)),
        cdna_sequence=observed, evidence={"read_count": 7})
    variant = StructuralVariant(first.contig, first.start, "BND", genome=genome,
                                mate_contig=second.contig, mate_start=second.start)
    candidate = EffectCandidate(GeneFusion(variant, first, second, mutant_transcript=model), source="rna")
    result, = reference_completion_hypotheses(candidate)
    restored = EffectCandidate.from_json(result.to_json())
    assert restored.source == result.source
    assert restored.evidence == result.evidence
    assert restored.effect.mutant_protein_sequence == result.effect.mutant_protein_sequence
    completed = restored.effect.mutant_transcript
    assert completed.cdna_sequence == first.sequence[:left + 9] + second.sequence[right:]
    assert completed.cdna_sequence[slice(*restored.evidence["reference_completion"]["observed_span"])] == observed
