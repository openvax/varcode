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
import random
import logging

from Bio.Seq import reverse_complement
from pyensembl import EnsemblRelease
from pyensembl.release_info import MAX_ENSEMBL_RELEASE

from .nucleotides import VALID_NUCLEOTIDES
from .variant import Variant
from .variant_collection import VariantCollection

# cache lists of all transcript IDs for difference Ensembl releases
_transcript_ids_cache = {}


def random_variants(
        count,
        ensembl_release=MAX_ENSEMBL_RELEASE,
        deletions=True,
        insertions=True):
    """
    Generate a VariantCollection with random variants that overlap
    at least one complete coding transcript.
    """
    ensembl = EnsemblRelease(ensembl_release)
    if ensembl_release in _transcript_ids_cache:
        transcript_ids = _transcript_ids_cache[ensembl_release]
    else:
        transcript_ids = ensembl.transcript_ids()
        _transcript_ids_cache[ensembl_release] = transcript_ids
    variants = []
    while len(variants) < count:
        transcript_id = random.choice(transcript_ids)
        transcript = ensembl.transcript_by_id(transcript_id)
        if not transcript.complete:
            continue

        exon = random.choice(transcript.exons)
        base1_genomic_position = random.randint(exon.start, exon.end)
        transcript_offset = transcript.spliced_offset(base1_genomic_position)

        try:
            seq = transcript.sequence
        except ValueError as e:
            logging.warn(e)
            # can't get sequence for non-coding transcripts
            continue
        ref = str(seq[transcript_offset])
        if transcript.on_backward_strand:
            ref = reverse_complement(ref)
        alt_nucleotides = [x for x in VALID_NUCLEOTIDES if x != ref]
        if insertions:
            nucleotide_pairs = [
                x + y
                for x in alt_nucleotides
                for y in VALID_NUCLEOTIDES
            ]
            alt_nucleotides.extend(nucleotide_pairs)
        if deletions:
            alt_nucleotides.append("")
        alt = random.choice(alt_nucleotides)
        variant = Variant(
            transcript.contig,
            base1_genomic_position,
            ref=ref,
            alt=alt,
            ensembl=ensembl)
        variants.append(variant)
    return VariantCollection(variants)
