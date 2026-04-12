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
Regression tests for https://github.com/openvax/varcode/issues/208

Silent.aa_pos was off by one for synonymous SNPs — the code returned
aa_mutation_start_offset which had already been incremented past the
shared prefix from trim_shared_flanking_strings.

Test transcript ENSMUST00000086738 (Cntnap5b-201, GRCm38), with the
coding sequence context provided by the reporter:

    ATG GAT TCT GTA CCA AGA CTG ACC AGC ATT TTG
     M   D   S   V   P   R   L   T   S   I   L
     0   1   2   3   4   5   6   7   8   9  10

The R codon AGA starts at CDS position 15 = transcript offset
(start_codon_spliced_offsets[0] + 15). On genome GRCm38 the reporter
specifies position chr1:99772782 = the last nucleotide of the R codon
(AGA, so position offset 2).
"""

from varcode import Variant
from varcode.effects import Silent


TRANSCRIPT_ID = "ENSMUST00000086738"


def _effect(ref, alt, start=99772782):
    variant = Variant(
        contig=1,
        start=start,
        ref=ref,
        alt=alt,
        genome="GRCm38",
    )
    transcript = variant.ensembl.transcript_by_id(TRANSCRIPT_ID)
    return variant.effect_on_transcript(transcript)


def test_synonymous_snv_at_third_codon_position_aa_pos_is_5():
    # AGA -> AGG (R/R), SNV at the third nucleotide of codon 5 (R).
    effect = _effect("A", "G", start=99772782)
    assert effect.__class__ is Silent, \
        "Expected Silent, got %s" % effect.__class__.__name__
    assert effect.aa_pos == 5, \
        "Expected aa_pos=5, got %d" % effect.aa_pos
    assert effect.aa_ref == "R", \
        "Expected aa_ref='R', got %r" % effect.aa_ref


def test_synonymous_snv_at_first_codon_position_aa_pos_is_6():
    # CTG -> TTG (L/L), SNV at the first nucleotide of codon 6 (L).
    effect = _effect("C", "T", start=99772783)
    assert effect.__class__ is Silent, \
        "Expected Silent, got %s" % effect.__class__.__name__
    assert effect.aa_pos == 6, \
        "Expected aa_pos=6, got %d" % effect.aa_pos
    assert effect.aa_ref == "L", \
        "Expected aa_ref='L', got %r" % effect.aa_ref


def test_synonymous_mnv_spanning_two_codons_aa_pos_is_5():
    # AC -> GT, affects third base of codon 5 (R) and first of codon 6 (L).
    # AGA CTG -> AGG TTG (RL/RL), fully synonymous.
    effect = _effect("AC", "GT", start=99772782)
    assert effect.__class__ is Silent, \
        "Expected Silent, got %s" % effect.__class__.__name__
    assert effect.aa_pos == 5, \
        "Expected aa_pos=5, got %d" % effect.aa_pos
    assert effect.aa_ref == "RL", \
        "Expected aa_ref='RL', got %r" % effect.aa_ref
