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

import numpy as np

import typechecks

# include all pseudonucleotides encoding repeats and uncertain bases
STANDARD_NUCLEOTIDES = {'A', 'C', 'T', 'G'}

PURINE_NUCLEOTIDES = {'A', 'G'}

PYRIMIDINE_NUCLEOTIDES = {'C', 'T'}

AMINO_NUCLEOTIDES = {'A', 'C'}

KETO_NUCLEOTIDES = {'T', 'G'}

STRONG_NUCLEOTIDES = {'G', 'C'}

WEAK_NUCLEOTIDES = {'A', 'T'}

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

def is_purine(nucleotide, allow_extended_nucleotides=False):
    """Is the nucleotide a purine"""
    if not allow_extended_nucleotides and nucleotide not in STANDARD_NUCLEOTIDES:
        raise ValueError("{} is a non-standard nucleotide, neither purine or pyrimidine".format(nucleotide))
    return nucleotide in PURINE_NUCLEOTIDES

def all_standard_nucleotides(nucleotides):
    return all(base in STANDARD_NUCLEOTIDES for base in nucleotides)

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

    # VCF files sometimes have '.' ref or alt for insertions and deletions, and
    # MAF files sometimes have '-' ref or alt for insertions and deletions.
    if nucleotides == "." or nucleotides == "-":
        return ""

    typechecks.require_string(nucleotides, "nucleotide string")

    nucleotides = nucleotides.upper()

    if allow_extended_nucleotides:
        valid_nucleotides = EXTENDED_NUCLEOTIDES
    else:
        valid_nucleotides = STANDARD_NUCLEOTIDES

    if not set(nucleotides) <= valid_nucleotides:
        raise ValueError(
            "Invalid character(s) in nucleotide string: %s" % (
                ",".join(set(nucleotides) - valid_nucleotides),))

    return nucleotides
