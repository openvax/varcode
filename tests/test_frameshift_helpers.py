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

from varcode.effects.effect_prediction_coding_frameshift import (
    cdna_codon_sequence_after_insertion_frameshift,
    cdna_codon_sequence_after_deletion_or_substitution_frameshift,
)

from .common import eq_

def test_cdna_codon_sequence_after_insertion_frameshift_before_codon():
    # insertion: T_ATGCCCTAG
    i, s = cdna_codon_sequence_after_insertion_frameshift(
        sequence_from_start_codon="ATGCCCTAG",
        cds_offset_before_insertion=-1,
        inserted_nucleotides="T")
    eq_(i, 0)
    eq_(s, "TATGCCCTAG")

def test_cdna_codon_sequence_after_insertion_frameshift_in_middle_of_codon():
    # insertion: A_T_TGCCCTAG
    i, s = cdna_codon_sequence_after_insertion_frameshift(
        sequence_from_start_codon="ATGCCCTAG",
        cds_offset_before_insertion=0,
        inserted_nucleotides="T")
    eq_(i, 0)
    eq_(s, "ATTGCCCTAG")

def test_cdna_codon_sequence_after_insertion_frameshift_at_end_of_codon():
    # insertion: AT_T_GCCCTAG
    i, s = cdna_codon_sequence_after_insertion_frameshift(
        sequence_from_start_codon="ATGCCCTAG",
        cds_offset_before_insertion=1,
        inserted_nucleotides="T")
    eq_(i, 0)
    eq_(s, "ATTGCCCTAG")

def test_cdna_codon_sequence_after_insertion_frameshift_after_codon():
    # insertion: ATG_T_CCCTAG
    i, s = cdna_codon_sequence_after_insertion_frameshift(
        sequence_from_start_codon="ATGCCCTAG",
        cds_offset_before_insertion=2,
        inserted_nucleotides="T")
    eq_(i, 1)
    eq_(s, "TCCCTAG")

def test_cdna_codon_sequence_after_deletion_or_substitution_frameshift_delA():
    i, s = cdna_codon_sequence_after_deletion_or_substitution_frameshift(
        sequence_from_start_codon="ATGCCCTAG",
        cds_offset=0,
        trimmed_cdna_ref="A",
        trimmed_cdna_alt="")
    eq_(i, 0)
    eq_(s, "TGCCCTAG")


def test_cdna_codon_sequence_after_deletion_or_substitution_frameshift_AT_to_C():
    i, s = cdna_codon_sequence_after_deletion_or_substitution_frameshift(
        sequence_from_start_codon="ATGCCCTAG",
        cds_offset=0,
        trimmed_cdna_ref="AT",
        trimmed_cdna_alt="C")
    eq_(i, 0)
    eq_(s, "CGCCCTAG")
