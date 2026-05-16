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

"""RNA BAM-backed phasing source.

This module intentionally contains the pysam/CIGAR interpretation code
separately from :mod:`varcode.phasing`, whose job is only to define the
generic phasing protocols and resolvers.
"""

from collections import defaultdict
from typing import Optional, Sequence


_CIGAR_MATCH = 0
_CIGAR_INSERTION = 1
_CIGAR_DELETION = 2
_CIGAR_REF_SKIP = 3
_CIGAR_SOFT_CLIP = 4
_CIGAR_EQUAL = 7
_CIGAR_DIFF = 8

_CONSUMES_REFERENCE = {
    _CIGAR_MATCH,
    _CIGAR_DELETION,
    _CIGAR_REF_SKIP,
    _CIGAR_EQUAL,
    _CIGAR_DIFF,
}
_CONSUMES_QUERY = {
    _CIGAR_MATCH,
    _CIGAR_INSERTION,
    _CIGAR_SOFT_CLIP,
    _CIGAR_EQUAL,
    _CIGAR_DIFF,
}


class RNAReadPhasingSource:
    """BAM-backed phasing source for RNA read co-occurrence.

    This is the lightweight alternative to an Isovar-style assembly
    source. It reads quality-filtered alignments from an RNA-seq BAM and
    answers whether variants are observed on the same read or paired-end
    fragment. It does **not** assemble contigs and does not provide
    ``mutant_transcript``; callers that need observed mutant proteins
    should use an assembly-backed source.

    Usage::

        source = RNAReadPhasingSource("tumor.rna.bam")
        resolver = MolecularPhaseResolver(source)
        effects = variants.effects(phase_resolver=resolver)

    Parameters
    ----------
    bam_path : str
        Coordinate-sorted, indexed BAM path.
    variants : sequence, optional
        Optional universe used by :meth:`partners_in_cis`. Variants seen
        through :meth:`has_evidence` or :meth:`in_cis` are registered
        automatically, so this is mainly a convenience for callers that
        query ``MolecularPhaseResolver.phased_partners`` directly.
    min_mapping_quality : int
        Minimum MAPQ for reads contributing evidence.
    min_base_quality : int
        Minimum base quality for SNV/MNV and insertion allele calls.
    min_alt_reads : int
        Minimum alt-supporting reads/fragments required for
        ``has_evidence`` and minimum co-occurring fragments required for
        a cis/trans call.
    max_distance_from_read_edge : int, optional
        Discard allele calls whose queried bases are closer than this
        many bases to either read edge. Set to ``None`` to disable.
    require_proper_pair : bool
        For paired reads, discard fragments not marked proper pair.
        Unpaired reads are still accepted.
    skip_duplicates, skip_secondary, skip_supplementary : bool
        Standard SAM flag filters.
    """

    source = "rna_reads"

    def __init__(
            self,
            bam_path: str,
            *,
            variants=None,
            min_mapping_quality: int = 20,
            min_base_quality: int = 20,
            min_alt_reads: int = 2,
            max_distance_from_read_edge: Optional[int] = 5,
            require_proper_pair: bool = True,
            skip_duplicates: bool = True,
            skip_secondary: bool = True,
            skip_supplementary: bool = True):
        try:
            import pysam
        except ImportError as e:
            raise ImportError(
                "RNAReadPhasingSource requires pysam. Install with "
                "`pip install varcode[rna]`.") from e
        self._pysam = pysam
        self._bam = pysam.AlignmentFile(bam_path, "rb")
        self.min_mapping_quality = min_mapping_quality
        self.min_base_quality = min_base_quality
        self.min_alt_reads = min_alt_reads
        self.max_distance_from_read_edge = max_distance_from_read_edge
        self.require_proper_pair = require_proper_pair
        self.skip_duplicates = skip_duplicates
        self.skip_secondary = skip_secondary
        self.skip_supplementary = skip_supplementary
        self._known_variants = []
        self._known_variant_keys = set()
        self._support_cache = {}
        self._phase_cache = {}
        self._contig_cache = {}
        if variants is not None:
            self.register_variants(variants)

    def close(self):
        """Close the underlying BAM handle."""
        self._bam.close()

    def register_variants(self, variants):
        """Register variants used by :meth:`partners_in_cis`.

        This is optional for ordinary ``MolecularPhaseResolver.in_cis`` use,
        where both queried variants are registered automatically.
        """
        for variant in variants:
            self._register_variant(variant)

    def _register_variant(self, variant):
        key = self._variant_key(variant)
        if key not in self._known_variant_keys:
            self._known_variant_keys.add(key)
            self._known_variants.append(variant)

    @staticmethod
    def _variant_key(variant):
        return (
            variant.contig,
            variant.start,
            variant.end,
            variant.ref,
            variant.alt,
            getattr(variant, "reference_name", None),
        )

    def _bam_contig(self, variant):
        key = variant.contig
        if key in self._contig_cache:
            return self._contig_cache[key]
        candidates = [variant.contig]
        if variant.contig.startswith("chr"):
            candidates.append(variant.contig[3:])
        else:
            candidates.append("chr" + variant.contig)
        if variant.contig == "MT":
            candidates.extend(["M", "chrM"])
        elif variant.contig in ("M", "chrM"):
            candidates.extend(["MT", "chrMT"])
        references = set(self._bam.references)
        for contig in candidates:
            if contig in references:
                self._contig_cache[key] = contig
                return contig
        self._contig_cache[key] = None
        return None

    def _passes_read_filters(self, read):
        if read.is_unmapped:
            return False
        if read.mapping_quality < self.min_mapping_quality:
            return False
        if self.skip_duplicates and read.is_duplicate:
            return False
        if self.skip_secondary and read.is_secondary:
            return False
        if self.skip_supplementary and read.is_supplementary:
            return False
        if (self.require_proper_pair and read.is_paired and
                not read.is_proper_pair):
            return False
        return True

    def _base_calls_ok(self, read, query_positions):
        if not query_positions:
            return False
        query_length = read.query_length
        if query_length is None:
            query_length = len(read.query_sequence)
        qualities = read.query_qualities
        for query_pos in query_positions:
            if query_pos is None:
                return False
            if self.max_distance_from_read_edge is not None:
                edge_distance = min(query_pos, query_length - query_pos - 1)
                if edge_distance < self.max_distance_from_read_edge:
                    return False
            if self.min_base_quality:
                if qualities is None or qualities[query_pos] < self.min_base_quality:
                    return False
        return True

    @staticmethod
    def _aligned_pairs(read):
        return read.get_aligned_pairs(matches_only=False)

    @staticmethod
    def _cigar_events(read):
        ref_pos = read.reference_start
        query_pos = 0
        for operation, length in read.cigartuples or ():
            ref_start = ref_pos
            query_start = query_pos
            if operation in _CONSUMES_REFERENCE:
                ref_pos += length
            if operation in _CONSUMES_QUERY:
                query_pos += length
            yield operation, ref_start, ref_pos, query_start, query_pos

    def _substitution_allele(self, read, variant):
        start0 = variant.start - 1
        positions = range(start0, start0 + len(variant.ref))
        ref_to_query = {
            ref_pos: query_pos
            for query_pos, ref_pos in self._aligned_pairs(read)
            if ref_pos is not None
        }
        query_positions = [ref_to_query.get(pos) for pos in positions]
        if not self._base_calls_ok(read, query_positions):
            return None
        observed = "".join(read.query_sequence[pos] for pos in query_positions)
        if observed.upper() == variant.alt.upper():
            return "alt"
        if observed.upper() == variant.ref.upper():
            return "ref"
        return None

    def _has_exact_deletion(self, read, start0, end0):
        for operation, ref_start, ref_end, _, _ in self._cigar_events(read):
            if (operation == _CIGAR_DELETION and
                    ref_start == start0 and ref_end == end0):
                return True
        return False

    def _deletion_allele(self, read, variant):
        if variant.alt:
            return None
        start0 = variant.start - 1
        end0 = variant.end
        if self._has_exact_deletion(read, start0, end0):
            return "alt"
        target = set(range(start0, end0))
        ref_to_query = {
            ref_pos: query_pos
            for query_pos, ref_pos in self._aligned_pairs(read)
            if ref_pos is not None
        }
        if not target.issubset(ref_to_query):
            return None
        query_positions = [ref_to_query[pos] for pos in sorted(target)]
        if any(pos is None for pos in query_positions):
            return None
        if not self._base_calls_ok(read, query_positions):
            return None
        observed = "".join(read.query_sequence[pos] for pos in query_positions)
        if observed.upper() == variant.ref.upper():
            return "ref"
        return None

    def _inserted_positions_after_anchor(self, read, anchor0):
        inserted = []
        for operation, ref_start, _, query_start, query_end in self._cigar_events(read):
            if operation == _CIGAR_INSERTION and ref_start == anchor0 + 1:
                inserted.extend(range(query_start, query_end))
        return inserted

    def _insertion_allele(self, read, variant):
        if variant.ref:
            return None
        anchor0 = variant.start - 1
        ref_to_query = {
            ref_pos: query_pos
            for query_pos, ref_pos in self._aligned_pairs(read)
            if ref_pos is not None
        }
        anchor_query = ref_to_query.get(anchor0)
        if anchor_query is None:
            return None
        inserted_positions = self._inserted_positions_after_anchor(read, anchor0)
        if inserted_positions:
            if not self._base_calls_ok(read, inserted_positions):
                return None
            inserted = "".join(read.query_sequence[pos] for pos in inserted_positions)
            if inserted.upper() == variant.alt.upper():
                return "alt"
            return None
        if self._base_calls_ok(read, [anchor_query]):
            return "ref"
        return None

    def _read_allele(self, read, variant):
        if not self._passes_read_filters(read):
            return None
        if variant.is_insertion:
            return self._insertion_allele(read, variant)
        if variant.is_deletion:
            return self._deletion_allele(read, variant)
        if len(variant.ref) == len(variant.alt) and variant.ref != variant.alt:
            return self._substitution_allele(read, variant)
        return None

    def _fetch_reads_for_variant(self, variant):
        contig = self._bam_contig(variant)
        if contig is None:
            return None
        start0 = max(0, variant.start - 1)
        end0 = max(start0 + 1, variant.end)
        return self._bam.fetch(contig, start0, end0)

    def supports_variant(self, variant) -> Optional[int]:
        """Count quality-filtered reads/fragments supporting ``variant.alt``.

        Returns ``None`` when the variant's contig is absent from the BAM.
        """
        self._register_variant(variant)
        key = self._variant_key(variant)
        if key in self._support_cache:
            return self._support_cache[key]
        reads = self._fetch_reads_for_variant(variant)
        if reads is None:
            self._support_cache[key] = None
            return None
        supporting_fragments = set()
        for read in reads:
            if self._read_allele(read, variant) == "alt":
                supporting_fragments.add(read.query_name)
        count = len(supporting_fragments)
        self._support_cache[key] = count
        return count

    def has_evidence(self, variant) -> bool:
        """True if the BAM has enough alt-supporting reads/fragments."""
        count = self.supports_variant(variant)
        return count is not None and count >= self.min_alt_reads

    def _fetch_reads_for_pair(self, v1, v2):
        if v1.contig != v2.contig:
            return None
        contig = self._bam_contig(v1)
        if contig is None:
            return None
        start0 = max(0, min(v1.start, v2.start) - 1)
        end0 = max(v1.end, v2.end)
        return self._bam.fetch(contig, start0, end0)

    def _fragment_alleles_for_pair(self, v1, v2):
        reads = self._fetch_reads_for_pair(v1, v2)
        if reads is None:
            return ()
        grouped = defaultdict(list)
        for read in reads:
            if self._passes_read_filters(read):
                grouped[read.query_name].append(read)
        fragments = []
        for fragment_reads in grouped.values():
            alleles = []
            for variant in (v1, v2):
                calls = [
                    self._read_allele(read, variant)
                    for read in fragment_reads
                ]
                if "alt" in calls:
                    alleles.append("alt")
                elif "ref" in calls:
                    alleles.append("ref")
                else:
                    alleles.append(None)
            fragments.append(tuple(alleles))
        return fragments

    def _phase_counts(self, v1, v2):
        key = tuple(sorted(
            (self._variant_key(v1), self._variant_key(v2))))
        if key in self._phase_cache:
            return self._phase_cache[key]
        self._register_variant(v1)
        self._register_variant(v2)
        both_alt = 0
        mixed = 0
        for a1, a2 in self._fragment_alleles_for_pair(v1, v2):
            if a1 == "alt" and a2 == "alt":
                both_alt += 1
            elif ((a1 == "alt" and a2 == "ref") or
                    (a1 == "ref" and a2 == "alt")):
                mixed += 1
        counts = both_alt, mixed
        self._phase_cache[key] = counts
        return counts

    def in_cis(self, v1, v2, transcript=None) -> Optional[bool]:
        """Return cis/trans from RNA read or fragment co-occurrence.

        ``True`` means enough fragments support both alts. ``False``
        means enough fragments support one alt with the other's reference
        allele. ``None`` means the BAM does not contain enough
        co-covering evidence to decide.
        """
        both_alt, mixed = self._phase_counts(v1, v2)
        if both_alt >= self.min_alt_reads and both_alt > mixed:
            return True
        if mixed >= self.min_alt_reads and mixed > both_alt:
            return False
        return None

    def partners_in_cis(self, variant) -> Sequence:
        """Known registered variants observed in cis with ``variant``."""
        self._register_variant(variant)
        partners = []
        for other in self._known_variants:
            if other == variant:
                continue
            if self.in_cis(variant, other) is True:
                partners.append(other)
        return tuple(partners)
