expected_effect_properties = [
    'gene',
    'gene_name',
    'gene_id',
    'transcript',
    'transcript_name',
    'transcript_id',
    'modifies_coding_sequence',
    'modifies_protein_sequence',
    'aa_mutation_start_offset',
    'aa_mutation_end_offset',
    'mutant_protein_sequence',
    'short_description'
]

def check_effect_properties(effect):
    assert effect is not None
    # try accessing all the properties to make sure none crash
    for attribute_name in expected_effect_properties:
        getattr(effect, attribute_name)
    assert effect.short_description is not None
    assert effect.modifies_coding_sequence in {False, True}
    assert effect.modifies_protein_sequence in {False, True}

def expect_effect(
        variant,
        transcript_id,
        effect_class):
    transcript = variant.ensembl.transcript_by_id(transcript_id)
    effect = variant.effect_on_transcript(transcript)
    assert isinstance(effect, effect_class), \
        "Expected %s on %s to be %s, got %s" % (
            variant, transcript, effect_class.__name__, effect)
    assert effect_class.__name__ in str(effect), \
        "Expected string representation of %s to include effect name %s" % (
            effect, effect_class.__name__)
    assert effect.short_description is not None, \
        "Expected effect %s to have a `short_description` property" % (effect,)
    check_effect_properties(effect)
