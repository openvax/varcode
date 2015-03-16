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


def mutate(sequence, position, variant_ref, variant_alt):
    """
    Mutate a sequence by substituting given `alt` at instead of `ref` at the
    given `position`.

    Parameters
    ----------
    sequence : sequence
        String of amino acids or DNA bases

    position : int
        Position in the sequence, starts from 0

    variant_ref : sequence or str
        What do we expect to find at the position?

    variant_alt : sequence or str
        Alternate sequence to insert
    """
    n_variant_ref = len(variant_ref)
    sequence_ref = sequence[position:position + n_variant_ref]
    assert str(sequence_ref) == str(variant_ref), \
        "Reference %s at position %d != expected reference %s" % \
        (sequence_ref, position, variant_ref)
    prefix = sequence[:position]
    suffix = sequence[position + n_variant_ref:]
    return prefix + variant_alt + suffix
