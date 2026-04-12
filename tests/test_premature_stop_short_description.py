# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Regression tests for https://github.com/openvax/varcode/issues/216

PrematureStop.short_description didn't distinguish an insertion that
introduces a stop (empty aa_ref) from a substitution or deletion that
ends with a stop. Previously produced ambiguous strings like "p.7*";
insertion cases should read "p.7ins*".
"""

from varcode import Variant
from varcode.effects import PrematureStop


# Reporter's example: chr1:99772782 A>ATAA on ENSMUST00000086738 (GRCm38).
# Protein context: M D S V P R L T S I L (1-indexed positions 1..11)
# The A at 99772782 is the third nucleotide of the R codon AGA.
# Inserting TAA after that A gives: AGA TAA CTG ACC ...
# Translation: AGA = R, TAA = stop.  The insertion creates a stop
# codon right after R without changing R itself.

def test_insertion_of_stop_codon_uses_ins_notation_in_short_description():
    variant = Variant(
        contig=1,
        start=99772782,
        ref="A",
        alt="ATAA",
        genome="GRCm38",
    )
    transcript = variant.ensembl.transcript_by_id("ENSMUST00000086738")
    effect = variant.effect_on_transcript(transcript)
    assert effect.__class__ is PrematureStop, \
        "Expected PrematureStop, got %s" % effect.__class__.__name__
    # aa_ref is empty (the insertion did not replace any amino acid),
    # so the short description should use "ins" notation instead of
    # the ambiguous "p.7*".
    assert effect.aa_ref == "", \
        "Expected aa_ref='', got %r" % effect.aa_ref
    assert effect.short_description == "p.7ins*", \
        "Expected 'p.7ins*', got %r" % effect.short_description


# Direct unit tests on short_description formatting logic.
# These use a minimal fake transcript that just satisfies the
# attributes accessed during __init__.

class _FakeGene:
    biotype = "protein_coding"


class _FakeTranscript:
    biotype = "protein_coding"
    gene = _FakeGene()

    def __init__(self, protein_length=100):
        self.protein_sequence = "X" * protein_length


def _premature_stop(aa_mutation_start_offset, aa_ref, aa_alt):
    return PrematureStop(
        variant=None,
        transcript=_FakeTranscript(),
        aa_mutation_start_offset=aa_mutation_start_offset,
        aa_ref=aa_ref,
        aa_alt=aa_alt,
    )


def test_short_description_pure_insertion_of_stop():
    # Empty aa_ref, empty aa_alt: insertion of just a stop codon.
    e = _premature_stop(aa_mutation_start_offset=6, aa_ref="", aa_alt="")
    assert e.short_description == "p.7ins*"


def test_short_description_insertion_of_aa_and_stop():
    # Empty aa_ref, non-empty aa_alt: insertion of amino acids followed
    # by a stop.
    e = _premature_stop(aa_mutation_start_offset=6, aa_ref="", aa_alt="K")
    assert e.short_description == "p.7insK*"


def test_short_description_substitution_ending_with_stop():
    # Non-empty aa_ref, empty aa_alt: amino acid replaced by stop.
    e = _premature_stop(aa_mutation_start_offset=6, aa_ref="L", aa_alt="")
    assert e.short_description == "p.L7*"


def test_short_description_complex_substitution_ending_with_stop():
    # Both non-empty: existing behaviour preserved.
    e = _premature_stop(aa_mutation_start_offset=6, aa_ref="LT", aa_alt="K")
    assert e.short_description == "p.LT7K*"
