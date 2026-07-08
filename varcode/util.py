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

import random

from pyensembl import genome_for_reference_name

from .nucleotides import STANDARD_NUCLEOTIDES, reverse_complement
from .variant import Variant
from .variant_collection import VariantCollection

# cache lists of all transcript IDs for difference Ensembl releases
_transcript_ids_cache = {}

def random_variants(
        count,
        genome_name="GRCh38",
        deletions=True,
        insertions=True,
        random_seed=None,
        ensembl=None):
    """
    Generate a VariantCollection with random variants that overlap
    at least one complete coding transcript.

    Parameters
    ----------
    ensembl : pyensembl.EnsemblRelease, optional
        Explicit genome to draw transcripts from. When ``None`` (default)
        the genome is resolved from ``genome_name`` via
        :func:`genome_for_reference_name`, which picks pyensembl's *latest*
        release for that assembly. Pass a pinned release (e.g.
        ``cached_release(81)``) when the caller needs deterministic data
        that is actually installed — ``genome_for_reference_name`` can
        resolve to a release that hasn't been downloaded.
    """
    rng = random.Random(random_seed)
    if ensembl is None:
        ensembl = genome_for_reference_name(genome_name)

    if ensembl in _transcript_ids_cache:
        transcript_ids = _transcript_ids_cache[ensembl]
    else:
        transcript_ids = ensembl.transcript_ids()
        _transcript_ids_cache[ensembl] = transcript_ids

    # Only draw from transcripts on contigs the genome considers valid.
    # ``transcript_ids()`` includes transcripts on alternate/patch scaffolds
    # (e.g. 'CHR_HSCHR11_2_CTG1') whose contig is NOT in ``genome.contigs()``;
    # a Variant built there passes construction but raises
    # ``ValueError: Invalid contig name`` lazily, when effect prediction
    # accesses ``.transcripts`` (see Variant._check_that_genome_has_contig).
    # Mirror that validity set here so we never hand back a variant that
    # can't be annotated.
    valid_contigs = set(ensembl.contigs())

    variants = []

    # we should finish way before this loop is over but just in case
    # something is wrong with PyEnsembl we want to avoid an infinite loop
    for _ in range(count * 100):
        if len(variants) < count:
            transcript_id = rng.choice(transcript_ids)
            transcript = ensembl.transcript_by_id(transcript_id)

            if not transcript.complete:
                continue

            if transcript.contig not in valid_contigs:
                # Alternate/patch scaffold — Variant would reject this contig
                # during annotation. Skip and draw another.
                continue

            try:
                exon = rng.choice(transcript.exons)
                base1_genomic_position = rng.randint(exon.start, exon.end)
                transcript_offset = transcript.spliced_offset(
                    base1_genomic_position)
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
            except ValueError:
                # Some transcripts live on alternate/patch contigs (e.g.
                # 'CHR_HSCHR19LRC_COX1_CTG3_1') that Variant rejects as
                # non-standard for the reference, and a few have sequence /
                # offset edge cases. Skip and draw another transcript rather
                # than failing the whole generator on an unlucky pick — this
                # otherwise made the result depend on the (often unseeded)
                # draw order.
                continue
            variants.append(variant)
        else:
            return VariantCollection(variants)
    raise ValueError(
        ("Unable to generate %d random variants, "
         "there may be a problem with PyEnsembl") % count)
