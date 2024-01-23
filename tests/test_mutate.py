# Copyright (c) 2015. Mount Sinai School of Medicine
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

from varcode.effects import mutate
from nose.tools import eq_

def test_snp_mutation():
    seq = "AACCTT"
    mutated = mutate.substitute(seq, 1, "A", "G")
    eq_(mutated, "AGCCTT")

def test_deletion_mutation():
    seq = "AACT"
    mutated = mutate.substitute(seq, 1, "ACT", "T")
    eq_(mutated, "AT")

def test_insert_before():
    mutated = mutate.insert_before("AACT", 1, "GG")
    eq_(mutated, "AGGACT")

def test_insert_after():
    mutated = mutate.insert_after("AACT", 1, "GG")
    eq_(mutated, "AAGGCT")
