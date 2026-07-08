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

    variants = []

    # we should finish way before this loop is over but just in case
    # something is wrong with PyEnsembl we want to avoid an infinite loop
    for _ in range(count * 100):
        if len(variants) < count:
            transcript_id = rng.choice(transcript_ids)
            transcript = ensembl.transcript_by_id(transcript_id)

            if not transcript.complete:
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
                # Force the lazy contig validation NOW, through the exact path
                # effect prediction uses: Variant.transcripts calls
                # _check_that_genome_has_contig, which reads a process-wide
                # valid-contig cache keyed by reference name. A variant built on
                # an alternate/patch scaffold (e.g. 'CHR_HSCHR6_MHC_MCF_CTG1')
                # constructs fine but raises here; skipping now guarantees we
                # never return a variant that blows up a later .effects().
                # Going through .transcripts -- not a separately computed contig
                # set -- is what makes this consistent with the caller (earlier
                # attempts compared against the wrong set and let scaffolds
                # through).
                overlapping = variant.transcripts
            except ValueError:
                # Alternate/patch-scaffold contig, or a transcript with a
                # sequence/offset edge case: skip and draw another rather than
                # failing the whole generator on an unlucky (often unseeded)
                # pick.
                continue
            if not overlapping:
                continue
            variants.append(variant)
        else:
            return VariantCollection(variants)
    raise ValueError(
        ("Unable to generate %d random variants, "
         "there may be a problem with PyEnsembl") % count)
