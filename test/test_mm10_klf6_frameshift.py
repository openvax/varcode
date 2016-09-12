from varcode import Variant
from varcode.effects import FrameShift
from varcode.effects.effect_prediction_coding_frameshift import (
    predict_frameshift_coding_effect,
    cdna_codon_sequence_after_insertion_frameshift,
)

from nose.tools import eq_

def validate_effect_values(effect):
    eq_(effect.__class__, FrameShift)
    transcript = effect.transcript
    eq_(transcript.name, "Klf6-201")
    eq_(transcript.spliced_offset(5864876), 462)
    eq_(effect.shifted_sequence, "GEEGGIRTEDFF")

def test_mm10_Klf6_frameshift():
    variant = Variant("chr13", 5864876, "", "G", "GRCm38")
    effects = variant.effects()
    eq_(len(effects), 1)
    validate_effect_values(effects[0])

def test_mm10_Klf6_frameshift_coding_effect_fn():
    variant = Variant("chr13", 5864876, "", "G", "GRCm38")
    eq_(len(variant.transcripts), 1)
    t = variant.transcripts[0]
    eq_(t.name, "Klf6-201")
    # first start codon offset is 150
    # mutation occurs after offset 462
    effect = predict_frameshift_coding_effect(
        trimmed_cdna_ref="",
        trimmed_cdna_alt="G",
        cds_offset=462 - 150,
        sequence_from_start_codon=t.sequence[150:],
        variant=variant,
        transcript=t)
    validate_effect_values(effect)

def test_mm10_Klf6_frameshift_cdna_codon_sequence():
    variant = Variant("chr13", 5864876, "", "G", "GRCm38")
    eq_(len(variant.transcripts), 1)
    t = variant.transcripts[0]
    eq_(t.name, "Klf6-201")
    mutant_codon_index, seq_after_mutated_codon = \
        cdna_codon_sequence_after_insertion_frameshift(
            sequence_from_start_codon=t.sequence[150:],
            cds_offset_before_insertion=462 - 150,
            inserted_nucleotides="G")
    eq_(mutant_codon_index, 104)
    expected_sequence = t.sequence[462] + "G" + t.sequence[463:]

    print("Reference sequence (first 20 bases): %s" % t.sequence[462:482])
    print("Expected sequence (first 20 bases):  %s" % expected_sequence[:20])
    print("Returned sequence (first 20 bases):  %s" % seq_after_mutated_codon[:20])
    eq_(seq_after_mutated_codon, expected_sequence)
