# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#         http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from six import integer_types

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
    assert len(str(effect)) > 0
    assert len(repr(effect)) > 0
    assert effect.short_description is not None, \
        "Expected effect %s to have a `short_description` property" % (effect,)
    assert len(effect.short_description) > 0
    assert effect.__class__.__name__ in str(effect), \
        "Expected string representation of %s to include effect name %s" % (
            effect, effect.__class__.__name__)

def expect_effect(
        variant,
        transcript_id=None,
        effect_class=None,
        protein_sequence=None,
        **kwargs):
    if transcript_id is None:
        effects = variant.effects()
        effect = effects.top_priority_effect()
    else:
        transcript = variant.ensembl.transcript_by_id(transcript_id)
        effect = variant.effect_on_transcript(transcript)
    check_effect_properties(effect)
    if effect_class is not None:
        assert effect.__class__ is effect_class, \
            "Expected effect class %s but got %s" % (
                effect_class.__name__,
                effect.__class__.__name__)
    if protein_sequence is not None:
        assert effect.mutant_protein_sequence == protein_sequence, \
            "Expected protein sequence %s but got %s" % (
                protein_sequence,
                effect.mutant_protein_sequence)
    for field, expected_value in kwargs.items():
        actual_value = getattr(effect, field)
        if isinstance(expected_value, integer_types):
            format_string = "Expected %s=%d but got %s"
        elif isinstance(expected_value, float):
            format_string = "Expected %s=%f but got %s"
        else:
            format_string = "Expected %s='%s' but got '%s'"
        assert actual_value == expected_value, format_string % (field, expected_value, actual_value)
