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

from nose.tools import eq_

from varcode.string_helpers import trim_shared_flanking_strings

def test_trim_shared_string_endings():
    # empty strings
    eq_(trim_shared_flanking_strings("", "A"), ("", "A", "", ""))
    eq_(trim_shared_flanking_strings("A", ""), ("A", "", "", ""))

    # string pairs with shared prefixes
    eq_(trim_shared_flanking_strings("AA", "AA"), ("", "", "AA", ""))
    eq_(trim_shared_flanking_strings("AB", "AA"), ("B", "A", "A", ""))
    eq_(trim_shared_flanking_strings("AA", "AB"), ("A", "B", "A", ""))
    eq_(trim_shared_flanking_strings("AB", "A"), ("B", "", "A", ""))
    eq_(trim_shared_flanking_strings("AB", "A"), ("B", "", "A", ""))
    eq_(trim_shared_flanking_strings("A", "AB"), ("", "B", "A", ""))

    # string pairs with shared suffixes
    eq_(trim_shared_flanking_strings("CCAT", "GT"),
        ("CCA", "G", "", "T"))
    eq_(trim_shared_flanking_strings("CCAT", "GT"),
        ("CCA", "G", "", "T"))

    # string pairs with shared prefixes+suffixes
    eq_(trim_shared_flanking_strings(
        "AATG", "AACG"), ("T", "C", "AA", "G"))
    eq_(trim_shared_flanking_strings(
        "ABG", "AG"), ("B", "", "A", "G"))
    eq_(trim_shared_flanking_strings(
        "AG", "ABG"), ("", "B", "A", "G"))
