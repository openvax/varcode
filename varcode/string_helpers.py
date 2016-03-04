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

from __future__ import print_function, division, absolute_import

def trim_shared_prefix(ref, alt):
    """
    Sometimes mutations are given with a shared prefix between the reference
    and alternate strings. Examples: C>CT (nucleotides) or GYFP>G (amino acids).

    This function trims the common prefix and returns the disjoint ref
    and alt strings, along with the shared prefix.
    """
    n_ref = len(ref)
    n_alt = len(alt)
    n_min = min(n_ref, n_alt)
    i = 0
    while i < n_min and ref[i] == alt[i]:
        i += 1

    # guaranteed that ref and alt agree on all the characters
    # up to i'th position, so it doesn't matter which one we pull
    # the prefix out of
    prefix = ref[:i]
    ref_suffix = ref[i:]
    alt_suffix = alt[i:]
    return ref_suffix, alt_suffix, prefix

def trim_shared_suffix(ref, alt):
    """
    Reuse the `trim_shared_prefix` function above to implement similar
    functionality for string suffixes.

    Given ref='ABC' and alt='BC', we first revese both strings:
        reverse_ref = 'CBA'
        reverse_alt = 'CB'
    and then the result of calling trim_shared_prefix will be:
        ('A', '', 'CB')
    We then reverse all three of the result strings to get back
    the shared suffix and both prefixes leading up to it:
        ('A', '', 'BC')
    """
    n_ref = len(ref)
    n_alt = len(alt)
    n_min = min(n_ref, n_alt)
    i = 0
    while i < n_min and ref[-i - 1] == alt[-i - 1]:
        i += 1

    # i is length of shared suffix.
    if i == 0:
        return (ref, alt, '')
    return (ref[:-i], alt[:-i], ref[-i:])

def trim_shared_flanking_strings(ref, alt):
    """
    Given two nucleotide or amino acid strings, identify
    if they have a common prefix, a common suffix, and return
    their unique components along with the prefix and suffix.

    For example, if the input ref = "SYFFQGR" and alt = "SYMLLFIFQGR"
    then the result will be:
        ("F", "MLLFI", "SY", "FQGR")
    """
    ref, alt, prefix = trim_shared_prefix(ref, alt)
    ref, alt, suffix = trim_shared_suffix(ref, alt)
    return ref, alt, prefix, suffix
