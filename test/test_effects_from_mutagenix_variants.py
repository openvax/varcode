# Copyright (c) 2016. Mount Sinai School of Medicine
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
List of variants copied from:
    https://mutagenetix.utsouthwestern.edu
    /incidental/incidental_rec.cfm?
    mid=&so=rb&ac=1&r0=0&nr=100&rn=29&rl=1&scd=IGL01779&mid=153891
"""

from varcode import Variant
from varcode.effects import Substitution

from .common import expect_effect

def test_substitution_Akt1_chr12_112657169_C_T_G286R():
    expect_effect(
        variant=Variant("chr12", 112657169, "C", "T", "GRCm38"),
        effect_class=Substitution,
        aa_mutation_start_offset=285,
        aa_ref="G",
        aa_alt="R")

def test_substitution_Apof_chr10_128269477_A_G_I167V():
    expect_effect(
        variant=Variant("chr10", 128269477, "A", "G", "GRCm38"),
        effect_class=Substitution,
        aa_mutation_start_offset=166,
        aa_ref="I",
        aa_alt="V")

def test_substitution_Csmd3_chr15_47857894_A_T_V1551D():
    expect_effect(
        variant=Variant("chr15", 47857894, "A", "T", "GRCm38"),
        effect_class=Substitution,
        aa_mutation_start_offset=1550,
        aa_ref="V",
        aa_alt="D")

def test_substitution_Pprc1_chr19_46062202_T_A_I130N():
    expect_effect(
        variant=Variant("chr19", 46062202, "T", "A", "GRCm38"),
        effect_class=Substitution,
        aa_mutation_start_offset=129,
        aa_ref="I",
        aa_alt="N")

def test_substitution_Vipr1_chr9_121664630_T_C_F249S():
    expect_effect(
        variant=Variant("chr9", 121664630, "T", "C", "GRCm38"),
        effect_class=Substitution,
        aa_mutation_start_offset=248,
        aa_ref="F",
        aa_alt="S")
