"""Sequence access must describe the products used for classification."""

import pickle
from types import SimpleNamespace

import pytest
from Bio.Seq import Seq
from pyensembl import cached_release

from varcode import EffectCollection, GermlineContext, Variant
from varcode.effects import Unresolved


@pytest.fixture
def cftr():
    return cached_release(81).transcript_by_id("ENST00000003084")


def translated(cdna, transcript):
    coding = cdna[min(transcript.start_codon_spliced_offsets):]
    return str(Seq(coding[:len(coding) // 3 * 3]).translate(to_stop=True))


def edited(sequence, transcript, variant):
    offset = transcript.spliced_offset(variant.start)
    assert sequence[offset:offset + len(variant.ref)] == variant.ref
    return sequence[:offset] + variant.alt + sequence[offset + len(variant.ref):]


@pytest.mark.parametrize("annotator", ["fast", "protein_diff", "transcript_model"])
def test_patient_sequences_keep_inherited_residue_outside_somatic_codon(cftr, annotator):
    somatic = Variant("7", 117531100, "T", "A", cftr.genome)
    germline = Variant("7", 117531097, "A", "C", cftr.genome)
    baseline = edited(cftr.sequence, cftr, germline)
    mutant = edited(baseline, cftr, somatic)
    result = somatic.effect_on_transcript(
        cftr, annotator=annotator,
        germline=GermlineContext.from_variants([germline]),
        phase_resolver=SimpleNamespace(in_cis=lambda *args: True))
    assert result.short_description == "p.L159M"
    assert result.original_protein_sequence == translated(baseline, cftr)
    assert result.mutant_protein_sequence == translated(mutant, cftr)
    assert result.baseline_transcript.cdna_sequence == baseline
    assert result.mutant_transcript.cdna_sequence == mutant
    assert result.mutant_transcript.mutant_protein_sequence == result.mutant_protein_sequence
    effects = EffectCollection([result])
    for restored in (EffectCollection.from_json(effects.to_json()),
                     pickle.loads(pickle.dumps(effects))):
        assert restored[0].baseline_transcript.cdna_sequence == baseline
        assert restored[0].mutant_transcript.cdna_sequence == mutant
        assert restored[0].original_protein_sequence == translated(baseline, cftr)
        assert restored[0].mutant_protein_sequence == translated(mutant, cftr)


def test_unavailable_splice_product_does_not_borrow_canonical_sequence(cftr):
    somatic = Variant("7", 117531115, "G", "A", cftr.genome)
    result = somatic.effect_on_transcript(cftr, annotator="transcript_model")
    unavailable = [outcome for candidate in result.candidates
                   for outcome in candidate.outcomes if isinstance(outcome.effect, Unresolved)]
    assert {outcome.effect.mechanism for outcome in unavailable} == {"intron_retention", "cryptic"}
    for outcome in unavailable:
        assert outcome.mutant.cdna_sequence is None
        assert outcome.mutant.protein_sequence is None
        assert outcome.mutant.start_codon_present is None
        assert outcome.baseline.cdna_sequence == cftr.sequence
        assert outcome.effect.mutant_protein_sequence is None
        assert outcome.effect.mutant_transcript.cdna_sequence is None


def test_resolved_splice_and_normal_products_use_standard_sequence_access(cftr):
    somatic = Variant("7", 117531115, "G", "A", cftr.genome)
    result = somatic.effect_on_transcript(cftr, annotator="transcript_model")
    for candidate in result.candidates:
        for outcome in candidate.outcomes:
            if isinstance(outcome.effect, Unresolved):
                continue
            assert outcome.effect.baseline_transcript.cdna_sequence == outcome.baseline.cdna_sequence
            assert outcome.effect.mutant_transcript.cdna_sequence == outcome.mutant.cdna_sequence
            assert outcome.effect.mutant_protein_sequence == translated(outcome.mutant.cdna_sequence, cftr)
            assert outcome.effect.mutant_protein_sequence == outcome.mutant.protein_sequence


@pytest.mark.parametrize("cis", [True, False, None])
def test_full_transcript_uses_distant_supplied_germline(cftr, cis):
    somatic = Variant("7", 117531100, "T", "A", cftr.genome)
    germline = Variant("7", 117587799, "A", "C", cftr.genome)
    result = somatic.effect_on_transcript(
        cftr, annotator="transcript_model",
        germline=GermlineContext.from_variants([germline]),
        phase_resolver=SimpleNamespace(in_cis=lambda *args: cis))
    baseline_sequences = ({cftr.sequence, edited(cftr.sequence, cftr, germline)}
                          if cis is None else
                          {edited(cftr.sequence, cftr, germline) if cis else cftr.sequence})
    assert {c.effect.baseline_transcript.cdna_sequence for c in result.candidates} == baseline_sequences
    assert len(result.candidates) == len(baseline_sequences)
    for candidate in result.candidates:
        baseline = candidate.effect.baseline_transcript.cdna_sequence
        mutant = edited(baseline, cftr, somatic)
        assert candidate.effect.original_protein_sequence == translated(baseline, cftr)
        assert candidate.effect.mutant_protein_sequence == translated(mutant, cftr)
        assert candidate.effect.mutant_transcript.cdna_sequence == mutant


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("with_sequence", [False, True])
def test_splice_products_assemble_supplied_bases_on_both_strands(reverse, with_sequence):
    from varcode import predict_transcript_model_effect
    from varcode.genomic_layout import SequenceUnavailable
    from .test_splice_graph import LongIntronTranscript, LONG_FORWARD, long_provider, point

    transcript = LongIntronTranscript(reverse)
    somatic = point(178, "C", "T") if reverse else point(23, "G", "A")
    germline = point(180, "T", "C") if reverse else point(21, "A", "G")
    full_provider = long_provider(reverse)

    def sequence_provider(contig, start, end):
        if not with_sequence and not any(
                exon.start <= start <= end <= exon.end for exon in transcript.exons):
            raise SequenceUnavailable("only annotated exons available")
        return full_provider(contig, start, end)

    result = predict_transcript_model_effect(
        (somatic,), transcript, germline_variants=(germline,),
        phase_resolver=SimpleNamespace(in_cis=lambda *args: True),
        sequence_provider=sequence_provider)
    genomic = LONG_FORWARD[:20] + "G" + LONG_FORWARD[21:22] + "A" + LONG_FORWARD[23:]
    expected = {
        "normal": "ATGCGAGCTTAG",
        "exon_skip": "ATGGCTTAG",
        "intron_retention": "ATG" + genomic[19:122] + "TAG",
        # The fixture's sole strong replacement donor ends at base 31.
        "cryptic": "ATG" + genomic[19:31] + "GCTTAG",
    }
    for candidate in result.candidates:
        for outcome in candidate.outcomes:
            mechanism, = [option.mechanism for _, option in outcome.hypothesis.splice_plan.choices]
            assert outcome.effect.baseline_transcript.cdna_sequence == "ATGCGAGCTTAG"
            if not with_sequence and mechanism in ("intron_retention", "cryptic"):
                assert isinstance(outcome.effect, Unresolved)
                assert outcome.effect.mutant_transcript.cdna_sequence is None
                continue
            cdna = expected[mechanism]
            assert outcome.effect.mutant_transcript.cdna_sequence == cdna
            assert outcome.effect.original_protein_sequence == "MRA"
            assert outcome.effect.mutant_protein_sequence == str(
                Seq(cdna[:len(cdna) // 3 * 3]).translate(to_stop=True))


@pytest.mark.parametrize("mutant_protein,length_delta,expected_class,offset", [
    ("MAAAAA", 0, "PrematureStop", 6),
    ("MAAAAA", -1, "FrameShiftTruncation", 6),
    ("MAAAAAAGH", 1, "FrameShift", 7),
    ("MAAAAAAAAAAGH", 6, "StopLoss", 11),
])
def test_stop_and_frameshift_offsets_use_actual_patient_baseline(
        mutant_protein, length_delta, expected_class, offset):
    from varcode import MutantTranscript
    from varcode.effects.classify import classify_from_protein_diff

    transcript = SimpleNamespace(
        protein_sequence="MQA", start_codon_spliced_offsets=(0, 1, 2),
        gene=SimpleNamespace(id="test-gene", name="SYN"))
    baseline_protein = "MAAAAAAAAAA"
    baseline = MutantTranscript(reference_transcript=transcript,
                                mutant_protein_sequence=baseline_protein)
    mutant = MutantTranscript(reference_transcript=transcript,
                              mutant_protein_sequence=mutant_protein)
    result = classify_from_protein_diff(
        SimpleNamespace(), transcript, baseline_protein, mutant_protein,
        length_delta=length_delta, baseline_transcript=baseline, mutant_transcript=mutant)
    assert type(result).__name__ == expected_class
    assert result.aa_mutation_start_offset == offset
    assert result.original_protein_sequence == baseline_protein
    assert result.mutant_protein_sequence == mutant_protein
