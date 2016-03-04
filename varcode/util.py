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
import random

from Bio.Seq import reverse_complement
from pyensembl import EnsemblRelease
from pyensembl import MAX_ENSEMBL_RELEASE

from .nucleotides import STANDARD_NUCLEOTIDES
from .variant import Variant
from .variant_collection import VariantCollection

# cache lists of all transcript IDs for difference Ensembl releases
_transcript_ids_cache = {}

def random_variants(
        count,
        ensembl_release=MAX_ENSEMBL_RELEASE,
        deletions=True,
        insertions=True,
        random_seed=None):
    """
    Generate a VariantCollection with random variants that overlap
    at least one complete coding transcript.
    """
    rng = random.Random(random_seed)
    ensembl = EnsemblRelease(ensembl_release)

    if ensembl_release in _transcript_ids_cache:
        transcript_ids = _transcript_ids_cache[ensembl_release]
    else:
        transcript_ids = ensembl.transcript_ids()
        _transcript_ids_cache[ensembl_release] = transcript_ids

    variants = []

    # we should finish way before this loop is over but just in case
    # something is wrong with PyEnsembl we want to avoid an infinite loop
    for _ in range(count * 100):
        if len(variants) < count:
            transcript_id = rng.choice(transcript_ids)
            transcript = ensembl.transcript_by_id(transcript_id)

            if not transcript.complete:
                continue

            exon = rng.choice(transcript.exons)
            base1_genomic_position = rng.randint(exon.start, exon.end)
            transcript_offset = transcript.spliced_offset(base1_genomic_position)
            seq = transcript.sequence

            ref = str(seq[transcript_offset])
            if transcript.on_backward_strand:
                ref = reverse_complement(ref)

            alt_nucleotides = [x for x in STANDARD_NUCLEOTIDES if x != ref]

            if insertions:
                nucleotide_pairs = [
                    x + y
                    for x in STANDARD_NUCLEOTIDES
                    for y in STANDARD_NUCLEOTIDES
                ]
                alt_nucleotides.extend(nucleotide_pairs)
            if deletions:
                alt_nucleotides.append("")
            alt = rng.choice(alt_nucleotides)
            variant = Variant(
                transcript.contig,
                base1_genomic_position,
                ref=ref,
                alt=alt,
                ensembl=ensembl)
            variants.append(variant)
        else:
            return VariantCollection(variants)
    raise ValueError(
        ("Unable to generate %d random variants, "
         "there may be a problem with PyEnsembl") % count)
