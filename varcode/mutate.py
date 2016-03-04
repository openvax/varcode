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


def insert_before(sequence, offset, new_residues):
    """Mutate the given sequence by inserting the string `new_residues` before
    `offset`.

    Parameters
    ----------
    sequence : sequence
        String of amino acids or DNA bases

    offset : int
        Base 0 offset from start of sequence, after which we should insert
        `new_residues`.

    new_residues : sequence
    """
    assert 0 < offset <= len(sequence), \
        "Invalid position %d for sequence of length %d" % (
            offset, len(sequence))
    prefix = sequence[:offset]
    suffix = sequence[offset:]
    return prefix + new_residues + suffix

def insert_after(sequence, offset, new_residues):
    """Mutate the given sequence by inserting the string `new_residues` after
    `offset`.

    Parameters
    ----------
    sequence : sequence
        String of amino acids or DNA bases

    offset : int
        Base 0 offset from start of sequence, after which we should insert
        `new_residues`.

    new_residues : sequence
    """
    assert 0 <= offset < len(sequence), \
        "Invalid position %d for sequence of length %d" % (
            offset, len(sequence))
    prefix = sequence[:offset + 1]
    suffix = sequence[offset + 1:]
    return prefix + new_residues + suffix

def substitute(sequence, offset, ref, alt):
    """Mutate a sequence by substituting given `alt` at instead of `ref` at the
    given `position`.

    Parameters
    ----------
    sequence : sequence
        String of amino acids or DNA bases

    offset : int
        Base 0 offset from start of `sequence`

    ref : sequence or str
        What do we expect to find at the position?

    alt : sequence or str
        Alternate sequence to insert
    """
    n_ref = len(ref)
    sequence_ref = sequence[offset:offset + n_ref]
    assert str(sequence_ref) == str(ref), \
        "Reference %s at offset %d != expected reference %s" % \
        (sequence_ref, offset, ref)
    prefix = sequence[:offset]
    suffix = sequence[offset + n_ref:]
    return prefix + alt + suffix
