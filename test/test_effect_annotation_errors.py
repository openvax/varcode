# Copyright (c) 2015-2016. Mount Sinai School of Medicine
#
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
Any time the effect annotation of a variant is erroneous, add it to this
test suite.
"""

from varcode import Variant
from varcode.effects import Silent, PrematureStop

def expect_effect(variant, effect_class=None, protein_sequence=None):
    effects = variant.effects()
    top_effect = effects.top_priority_effect()
    if effect_class is not None:
        assert top_effect.__class__ is effect_class, \
            "Expected effect class %s but got %s" % (
                effect_class.__name__,
                top_effect.__class__.__name__)
    if protein_sequence is not None:
        assert top_effect.mutant_protein_sequence == protein_sequence, \
            "Expected protein sequence %s but got %s" % (
                protein_sequence,
                top_effect.mutant_protein_sequence)

def test_issue172_insertion_after_stop_codon():
    # Issue: https://github.com/hammerlab/varcode/issues/172
    # Insertion immediately after the stop codon
    # TGA ACT ATT -> TAG CAG CAG ACT ATT...
    #
    # ##fileformat=VCFv4.1
    # #CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT
    # chr1    100484701   .   A   ACAGCAG 5000    .   .   .
    variant = Variant("chr1", 100484701, "A", "ACAGCAG", ensembl="GRCm38")
    expect_effect(
        variant=variant,
        effect_class=Silent)

def test_issue167_insertion_of_stop_codon():
    # Issue: https://github.com/hammerlab/varcode/issues/167
    # The coding sequence at chr1/99772765 has transcript ID ENSMUST00000086738
    # and looks like:
    # ATG GAT TCT GTA CCA AGA CTG ACC AGC ATT TTG
    #                       ^ position 99772782
    # The protein is ENSMUSP00000083944 and looks like:
    #  M   D   S   V   P   R   L   T   S   I   L
    # ##fileformat=VCFv4.1
    # CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT
    # chr1    99772782    .   ACTG    ATGA    5000    .   .   .
    variant = Variant("chr1", 99772782, "A", "ATGA", ensembl="GRCm38")
    expect_effect(
        variant=variant,
        effect_class=PrematureStop,
        protein_sequence="MDSVPR")
