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

from __future__ import print_function, division, absolute_import

import numpy as np

# include all pseudonucleotides encoding repeats and uncertain bases
VALID_NUCLEOTIDES = {'A', 'C', 'T', 'G'}

EXTENDED_NUCLEOTIDES = {
    'A', 'C', 'T', 'G',
    'Y',  # Pyrimidine (C or T)
    'R',  # Purine (A or G)
    'W',  # weak (A or T)
    'S',  # strong (G or C)
    'K',  # keto (T or G)
    'M',  # amino (C or A)
    'D',  # A, G, T (not C)
    'V',  # A, C, G (not T)
    'H',  # A, C, T (not G)
    'B',  # C, G, T (not A)
    'X',  # any base
    'N',  # any base
}

def normalize_nucleotide_string(nucleotides, allow_extended_nucleotides=False):
    """
    Normalizes a nucleotide string by converting various ways of encoding empty
    strings into "", making all letters upper case, and checking to make sure
    all letters in the string are actually nucleotides.

    Parameters
    ----------
    nucleotides : str
        Sequence of nucleotides, e.g. "ACCTG"

    extended_nucleotides : bool
        Allow non-canonical nucleotide characters like 'X' for unknown base
    """
    # some MAF files represent deletions/insertions with NaN ref/alt values
    if isinstance(nucleotides, float) and np.isnan(nucleotides):
        return ""
    # VCF files sometimes have '.' ref or alt for insertions and deletions
    elif nucleotides == ".":
        return ""
    # MAF files sometimes have '-' ref or alt for insertions and deletions
    elif nucleotides == "-":
        return ""

    if not isinstance(nucleotides, str):
        raise TypeError(
                "Expected nucleotide string, got %s : %s" % (
                    nucleotides, type(nucleotides)))

    nucleotides = nucleotides.upper()

    for letter in set(nucleotides):
        if allow_extended_nucleotides:
            valid_nucleotides = EXTENDED_NUCLEOTIDES
        else:
            valid_nucleotides = VALID_NUCLEOTIDES
        if letter not in valid_nucleotides:
            raise ValueError(
                "Invalid character in nucleotide string: %s" % letter)

    return nucleotides
