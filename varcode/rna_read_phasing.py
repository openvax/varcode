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
from math import exp, lgamma, log, log1p
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


def _binomial_tail(count, trials, rate):
    """P(X >= count) for X ~ Binomial(trials, rate), with 0 < rate < 1."""
    if count <= 0:
        return 1.0
    total, scale = 0.0, lgamma(trials + 1)
    for successes in range(count, trials + 1):
        term = exp(scale - lgamma(successes + 1) - lgamma(trials - successes + 1)
                   + successes * log(rate) + (trials - successes) * log1p(-rate))
        total += term
        # Past the mean the terms only shrink.
        if successes > trials * rate and term <= 1e-17 * total:
            break
    return min(total, 1.0)


def four_gamete_phase(both, first, second, neither, *, min_fragments=2,
                      error_rate=0.01, max_p_value=0.05) -> Optional[bool]:
    """Cis (``True``), trans (``False``) or unknown (``None``) for two
    variants, from the fragments covering both loci: ``both`` carry both
    alt alleles, ``first`` and ``second`` one variant's alt allele alone,
    ``neither`` neither.

    A combination is present when it's on at least ``min_fragments``
    fragments, more than reads showing the wrong allele at either locus
    (``error_rate``) leak into it from its two neighbouring combinations
    together, by a one-sided binomial test at ``max_p_value``. The pair is
    cis when both alt alleles are present together and not each alone, and
    trans when each is present alone and not together. Otherwise, including
    when all three are present, it's unknown.

    This is the four-gamete test (Hudson and Kaplan 1985) read as Nik-Zainal
    et al. 2012 (*Cell* 149:994) read phased mutation pairs. A variant that
    arose later on the other's copy is always with it, while the earlier
    one also appears alone, which isn't trans (#527). Reads showing a
    variant's alt allele in error put it on a small fraction of the other's
    alt fragments, which isn't cis (#547). A germline partner is the earlier
    variant, so the same rule applies. It matches ``IsovarReadPhasing`` in
    Isovar 1.39.7.
    """
    def present(count, neighbour, other_neighbour):
        return (count >= min_fragments and _binomial_tail(
            count, count + neighbour + other_neighbour, error_rate) <= max_p_value)

    together = present(both, first, second)
    first_alone = present(first, neither, both)
    second_alone = present(second, neither, both)
    if together and not (first_alone and second_alone):
        return True
    if first_alone and second_alone and not together:
        return False
    return None


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
        ``has_evidence``, and fragments an allele combination needs to
        count toward a cis/trans call.
    phasing_error_rate : float
        Chance that a read shows the wrong allele at a variant's locus,
        from sequencing error, RNA editing or mismapping. Raise it for
        error-prone reads, such as ONT's.
    max_p_value_for_phasing : float
        One-sided binomial p-value at or below which :meth:`in_cis`
        counts an allele combination as more than such errors.
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
            phasing_error_rate: float = 0.01,
            max_p_value_for_phasing: float = 0.05,
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
        for name, value in (("phasing_error_rate", phasing_error_rate),
                            ("max_p_value_for_phasing", max_p_value_for_phasing)):
            if not 0 < value < 1:
                raise ValueError("%s must be between 0 and 1, not %r" % (name, value))
        self._pysam = pysam
        self._bam = pysam.AlignmentFile(bam_path, "rb")
        self.min_mapping_quality = min_mapping_quality
        self.min_base_quality = min_base_quality
        self.min_alt_reads = min_alt_reads
        self.phasing_error_rate = phasing_error_rate
        self.max_p_value_for_phasing = max_p_value_for_phasing
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
        self._haplotypes = []
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

    def register_haplotype(self, variants, *, flanking_bases=5):
        """Test a known local allele combination by its anchored RNA sequence.

        Parameters
        ----------
        variants : sequence of Variant
            Nonoverlapping substitutions/deletions on one genome and contig.
            All alternate alleles form one hypothesis, not an assumed phase.
            Register competing hypotheses separately when necessary.
        flanking_bases : int
            Unchanged reference bases required on each side (default five).

        Notes
        -----
        A match requires the entire sequence, including both flanks, on ONE
        alignment, with the configured base-quality and read-edge filters.
        Equivalent D/N/split-gap encodings can then support the same known
        sequence. This does not establish a genomic deletion from an RNA skip.
        Reference is fetched from the variants' genome; unavailable sequence,
        overlapping edits and reference mismatches raise ValueError. No BAM
        bases are corrected, no missing bases filled, and no mates assembled.
        Registered variants use these full-context matches instead of individual
        CIGAR calls. Nonmatches remain unknown, not reference/trans evidence.
        """
        from .genome_sequence import reference_range

        variants = tuple(sorted(variants, key=lambda v: v.start))
        if not variants or isinstance(flanking_bases, bool) or not isinstance(flanking_bases, int) or flanking_bases < 1:
            raise ValueError("A nonempty haplotype and positive flanking_bases are required")
        first = variants[0]
        if any(v.contig != first.contig or v.genome is not first.genome for v in variants):
            raise ValueError("Haplotype variants must share a genome dataset and contig")
        if any(not v.ref or len(v.alt) > len(v.ref) for v in variants):
            raise ValueError("Anchored haplotypes currently support substitutions and deletions")
        if any(a.end >= b.start for a, b in zip(variants, variants[1:])):
            raise ValueError("Haplotype edits must not overlap")
        left, right = first.start - flanking_bases, variants[-1].end + flanking_bases
        if left < 1:
            raise ValueError("Haplotype lacks a left genomic anchor")
        sequence = reference_range(first.genome, first.contig, left, right).upper()
        if len(sequence) != right - left + 1 or set(sequence) - set("ACGT"):
            raise ValueError("Unambiguous reference sequence is required across the haplotype")
        for variant in reversed(variants):
            start = variant.start - left
            end = start + len(variant.ref)
            if sequence[start:end] != variant.ref.upper():
                raise ValueError("Haplotype reference allele disagrees with its genome")
            sequence = sequence[:start] + variant.alt.upper() + sequence[end:]
        if set(sequence) - set("ACGT"):
            raise ValueError("Unambiguous alternate sequence is required")
        keys = frozenset(self._variant_key(v) for v in variants)
        record = (keys, left - 1, right - 1, sequence)
        if record not in self._haplotypes:
            self._haplotypes.append(record)
            self.register_variants(variants)
            self._support_cache.clear()
            self._phase_cache.clear()

    def _matches_haplotype(self, read, left0, right0, sequence):
        anchors = {}
        for operation, start, end, query, _ in self._cigar_events(read):
            if operation in (_CIGAR_MATCH, _CIGAR_EQUAL, _CIGAR_DIFF):
                for position in (left0, right0):
                    if start <= position < end:
                        anchors[position] = query + position - start
        if left0 not in anchors or right0 not in anchors:
            return False
        left, right = anchors[left0], anchors[right0]
        if right < left or right - left + 1 != len(sequence):
            return False
        return (self._base_calls_ok(read, range(left, right + 1))
                and read.query_sequence[left:right + 1].upper() == sequence)

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
        if not self._passes_read_filters(read) or read.query_sequence is None:
            return None
        key = self._variant_key(variant)
        registered = False
        for keys, left, right, sequence in self._haplotypes:
            if key in keys:
                registered = True
                if self._matches_haplotype(read, left, right, sequence):
                    return "alt"
        if registered:
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
                supporting_fragments.add((
                    read.get_tag("RG") if read.has_tag("RG") else "", read.query_name))
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
                key = (read.get_tag("RG") if read.has_tag("RG") else "", read.query_name)
                grouped[key].append(read)
        fragments = []
        anchored = any(self._variant_key(v) in keys
                       for v in (v1, v2) for keys, _, _, _ in self._haplotypes)
        for fragment_reads in grouped.values():
            if anchored:
                # Never assemble a new anchored haplotype from mates or from
                # different registered combinations through a shared allele.
                pairs = {tuple(self._read_allele(read, v) for v in (v1, v2))
                         for read in fragment_reads}
                fragments.append(("alt", "alt") if ("alt", "alt") in pairs
                                 else (None, None))
                continue
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
        fragments = self._fragment_alleles_for_pair(v1, v2)
        # Both alt, v1 alt alone, v2 alt alone, neither; the phase call is
        # symmetric in the two variants, so a cached pair can be reversed.
        counts = tuple(fragments.count(combination) for combination in (
            ("alt", "alt"), ("alt", "ref"), ("ref", "alt"), ("ref", "ref")))
        self._phase_cache[key] = counts
        return counts

    def in_cis(self, v1, v2, transcript=None) -> Optional[bool]:
        """Return cis/trans from RNA read or fragment co-occurrence.

        Decided by :func:`four_gamete_phase` from the fragments covering
        both loci: ``True`` when both alt alleles are seen together beyond
        read errors and not each alone, ``False`` when each is seen alone
        and not together, ``None`` otherwise.
        """
        return four_gamete_phase(
            *self._phase_counts(v1, v2), min_fragments=self.min_alt_reads,
            error_rate=self.phasing_error_rate,
            max_p_value=self.max_p_value_for_phasing)

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
