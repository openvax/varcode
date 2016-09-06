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

Total: 28 list
Gene    Ref Var Chr/Loc Mutation    Predicted Effect    Zygosity
Akt1    C   T   12: 112,657,169 G286R   probably damaging   Het
Apof    A   G   10: 128,269,477 I167V   probably benign Het
Arhgap15    A   G   2: 44,065,045   E220G   probably damaging   Het
Clca3a2 T   A   3: 144,819,378  Y31F    possibly damaging   Het
Clmn    T   C   12: 104,782,140 I383V   probably damaging   Het
Col8a1  A   T   16: 57,628,363  H261Q   unknown Het
Csmd3   A   T   15: 47,857,894  V1551D  probably benign Het
Ddx60   G   A   8: 62,017,823   V1450M  probably damaging   Het
Ethe1   A   T   7: 24,595,009   H79L    probably damaging   Het
Fhdc1   G   A   3: 84,444,735   A⇒V possibly damaging   Het
Gm11110 C   T   17: 57,102,087  V72M    unknown Het
Hs1bp3  C   T   12: 8,341,945   T349I   probably benign Het
Ifna16  A   T   4: 88,676,645   I71N    probably damaging   Het
Il18bp  A   G   7: 102,016,795  Y59H    probably benign Het
Kcnt1   T   A   2: 25,900,967   I511N   probably damaging   Het
Lphn3   A   T   5: 81,387,870   I119F   probably damaging   Het
Mlph    A   G   1: 90,942,950   M528V   probably benign Het
Olfr646 A   G   7: 104,106,633  D118G   probably damaging   Het
Pprc1   T   A   19: 46,062,202  I130N   probably damaging   Het
Rfx1    T   A   8: 84,092,662       probably benign Het
Rnf17   A   T   14: 56,462,063  I553F   possibly damaging   Het
Scaper  A   T   9: 55,892,240   H180Q   probably benign Het
Slc26a4 C   T   12: 31,528,854          Het
Slc30a10    A   T   1: 185,464,179  Q346L   possibly damaging   Het
Stambpl1    A   T   19: 34,240,027  H422L   probably damaging   Het
Trim67  G   A   8: 124,828,121  G701R   probably damaging   Het
Vipr1   T   C   9: 121,664,630  F249S   possibly damaging   Het
Vmn2r117    G   T   17: 23,477,241  D397E   probably benign Het
"""

from varcode import Variant
from varcode.effects import Substitution

from .common import expect_effect

def test_substitution_Akt1_chr12_112657169_C_T_G286R():
    expect_effect(
        variant=Variant("chr12", 112657169, "C", "T", "GRCm38"),
        effect_class=Substitution,
        aa_ref="G",
        aa_alt="R")
