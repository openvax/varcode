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


"""In-cis, same-codon SNV combination transform (:func:`combine_cis_snvs`)."""

from __future__ import annotations

import warnings

from ..errors import InvalidContigError, SampleNotFoundError
from ..genome_sequence import reference_range
from ..phasing import VCFPhaseResolver, query_in_cis
from ..variant import Variant
from ..variant_collection import VariantCollection

__all__ = ["combine_cis_snvs"]


_STANDARD_DNA_BASES = frozenset("ACGT")


def _is_combinable_snv(variant):
    """Return True if ``variant`` is an underived single-nucleotide variant."""
    if getattr(variant, "is_structural", False):
        return False
    if not getattr(variant, "is_snv", False):
        return False
    if getattr(variant, "source_variants", ()):
        return False
    return True


def _transcript_codon_placement(transcript, pos):
    """Return ``(region, exon_idx, codon_index)`` for genomic ``pos`` on ``transcript``.

    Codon membership is computed relative to the annotated CDS start
    (``transcript.first_start_codon_spliced_offset``), not raw
    ``spliced_offset // 3``, so nonzero 5' UTR lengths never shift codon
    boundaries. ``exon_idx`` identifies the containing exon so split codons
    crossing an exon-exon junction are never collapsed across an intron.
    """
    exon_idx = None
    for idx, exon in enumerate(transcript.exons):
        if exon.start <= pos <= exon.end:
            exon_idx = idx
            break
    if exon_idx is None:
        return ("intron", None, None)
    spliced_offset = transcript.spliced_offset(pos)
    cds_start = transcript.first_start_codon_spliced_offset
    cds_end = transcript.last_stop_codon_spliced_offset
    if cds_start <= spliced_offset <= cds_end:
        codon_index = (spliced_offset - cds_start) // 3
        return ("cds", exon_idx, codon_index)
    return ("utr", exon_idx, None)


def _classify_sample_member_genotype(vc, variant, sample_name):
    """Classify ``sample_name``'s genotype state at ``variant`` for merge checks."""
    try:
        gt = vc.genotype(variant, sample_name)
    except SampleNotFoundError:
        return ("missing", None, (), None)
    if gt is None:
        return ("missing", None, (), None)
    called = [a for a in gt.alleles if a is not None]
    if not called:
        return ("missing", None, (), gt)
    alt_idx = vc._alt_index_for(variant)
    alt_slots = tuple(i for i, a in enumerate(gt.alleles) if a == alt_idx)
    if not alt_slots:
        return ("absent", None, (), gt)
    if len(alt_slots) == len(called):
        return ("hom_alt", None, alt_slots, gt)
    if gt.phased and gt.phase_set is not None:
        return ("phased", gt.phase_set, alt_slots, gt)
    return ("unphased_het", None, alt_slots, gt)


def _canonical_combined_gt(states):
    """Build a canonical VCF GT string for a compatible combined MNV."""
    kinds = {s[0] for s in states}
    if kinds == {"missing"}:
        return "./."
    if kinds <= {"missing", "absent"}:
        # Preserve ploidy if known.
        sample_gt = next((s[3] for s in states if s[3] is not None), None)
        ploidy = len(sample_gt.alleles) if sample_gt is not None and sample_gt.alleles else 2
        return "/".join(["0"] * ploidy)
    if kinds == {"hom_alt"}:
        sample_gt = states[0][3]
        ploidy = len(sample_gt.alleles) if sample_gt is not None and sample_gt.alleles else 2
        return "/".join(["1"] * ploidy)
    if kinds == {"phased"}:
        sample_gt = states[0][3]
        slots = set(states[0][2])
        ploidy = len(sample_gt.alleles) if sample_gt is not None and sample_gt.alleles else 2
        alleles = ["1" if i in slots else "0" for i in range(ploidy)]
        return "|".join(alleles)
    if kinds == {"unphased_het"}:
        sample_gt = states[0][3]
        ploidy = len(sample_gt.alleles) if sample_gt is not None and sample_gt.alleles else 2
        n_alt = len(states[0][2])
        alleles = ["0"] * max(0, ploidy - n_alt) + ["1"] * n_alt
        return "/".join(alleles)
    return None


def _check_samples_compatible_for_mnv(vc, members, is_vcf_resolver):
    """Return ``(ok, conflict_reason, merged_gt_by_sample)`` for ``members``.

    Ensures combining ``members`` into one MNV does not drop an uncombined
    haplotype (e.g. ``1/1`` + ``1|0`` where the second haplotype still carries
    only one SNV) and does not corrupt another sample in a multi-sample
    collection (e.g. trans, partial carrier, or unphased in a second sample).
    """
    sample_names = set()
    for by_variant in vc.source_to_metadata_dict.values():
        for v in members:
            meta = by_variant.get(v)
            if meta and meta.get("sample_info"):
                sample_names.update(meta["sample_info"].keys())
    if not sample_names:
        return True, None, {}

    merged_gt_by_sample = {}
    for sample_name in sorted(sample_names):
        states = [_classify_sample_member_genotype(vc, v, sample_name) for v in members]
        kinds = {s[0] for s in states}
        if kinds <= {"missing", "absent"}:
            merged_gt_by_sample[sample_name] = (_canonical_combined_gt(states), None)
            continue
        if kinds == {"hom_alt"}:
            merged_gt_by_sample[sample_name] = (_canonical_combined_gt(states), None)
            continue
        if kinds == {"phased"}:
            phase_sets = {s[1] for s in states}
            slot_sets = {s[2] for s in states}
            if len(phase_sets) == 1 and len(slot_sets) == 1:
                ps = next(iter(phase_sets))
                merged_gt_by_sample[sample_name] = (_canonical_combined_gt(states), ps)
                continue
            return (
                False,
                "sample %r has conflicting phased haplotypes/phase sets across %s"
                % (sample_name, [v.short_description for v in members]),
                {},
            )
        if kinds == {"unphased_het"} and not is_vcf_resolver:
            # An external molecular/read phase resolver established cis for
            # unphased DNA genotypes; allow if no other sample contradicts.
            merged_gt_by_sample[sample_name] = (_canonical_combined_gt(states), None)
            continue
        return (
            False,
            "sample %r has incompatible genotypes across %s (%s)"
            % (
                sample_name,
                [v.short_description for v in members],
                [str(s[3]) if s[3] is not None else None for s in states],
            ),
            {},
        )
    return True, None, merged_gt_by_sample


def _merge_cis_snv_metadata(members, by_variant, merged_gt_by_sample):
    """Build merged source metadata for a combined cis-SNV MNV."""
    member_metas = [by_variant.get(v) for v in members]
    if all(m is None for m in member_metas):
        return None
    ids = [m.get("id") for m in member_metas if m and m.get("id") is not None]
    quals = [m.get("qual") for m in member_metas if m and m.get("qual") is not None]
    raw_filters = [m.get("filter") for m in member_metas if m is not None]
    if not raw_filters or all(f is None for f in raw_filters):
        merged_filter = None
    else:
        non_pass = set()
        for f in raw_filters:
            if f is None:
                continue
            tokens = f if isinstance(f, (list, tuple, set, frozenset)) else str(f).split(";")
            for tok in tokens:
                if tok and tok not in ("PASS", "."):
                    non_pass.add(str(tok))
        merged_filter = sorted(non_pass)

    first_info = next((m.get("info") for m in member_metas if m and m.get("info") is not None), None)
    merged_info = dict(first_info) if first_info is not None else {}
    merged_info["combined_from_cis_snvs"] = True
    merged_info["original_starts"] = tuple(v.start for v in members)
    merged_info["source_ids"] = tuple(
        m.get("id") if m is not None else None for m in member_metas)
    merged_info["source_metadata"] = tuple(
        dict(m) if m is not None else None for m in member_metas)

    has_sample_info = any(m and m.get("sample_info") for m in member_metas)
    merged_sample_info = None
    if has_sample_info:
        merged_sample_info = {}
        all_samples = set()
        for m in member_metas:
            if m and m.get("sample_info"):
                all_samples.update(m["sample_info"].keys())
        for sample_name in sorted(all_samples):
            fields: dict[str, object] = {}
            for m in member_metas:
                s_info = (m or {}).get("sample_info") or {}
                if sample_name in s_info and s_info[sample_name]:
                    for k, val in s_info[sample_name].items():
                        fields.setdefault(k, val)
            if sample_name in merged_gt_by_sample:
                gt_str, ps = merged_gt_by_sample[sample_name]
                if gt_str is not None:
                    fields["GT"] = gt_str
                if ps is not None:
                    fields["PS"] = ps
            merged_sample_info[sample_name] = fields

    return {
        "id": ids[0] if ids else None,
        "qual": min(quals) if quals else None,
        "filter": merged_filter,
        "info": merged_info,
        "sample_info": merged_sample_info,
        "alt_allele_index": 0,
    }


def _build_combined_mnv(members):
    """Build a reference-validated MNV from cis in-codon SNVs, or ``None``."""
    min_pos = members[0].start
    max_pos = members[-1].start
    span_len = max_pos - min_pos + 1
    genome = members[0].genome
    contig = members[0].contig
    if any(v.genome != genome or v.contig != contig for v in members[1:]):
        warnings.warn(
            "combine_cis_snvs: declining %s because members do not share a genome/contig."
            % ([v.short_description for v in members],),
            stacklevel=3,
        )
        return None
    ref_seq = reference_range(genome, contig, min_pos, max_pos)
    if not ref_seq or len(ref_seq) != span_len:
        warnings.warn(
            "combine_cis_snvs: declining %s because reference bases for %s:%d-%d are unavailable."
            % ([v.short_description for v in members], contig, min_pos, max_pos),
            stacklevel=3,
        )
        return None
    ref_upper = ref_seq.upper()
    if any(base not in _STANDARD_DNA_BASES for base in ref_upper):
        warnings.warn(
            "combine_cis_snvs: declining %s because reference span %s:%d-%d contains "
            "ambiguous base(s) (%r)."
            % ([v.short_description for v in members], contig, min_pos, max_pos, ref_upper),
            stacklevel=3,
        )
        return None
    alt_chars = list(ref_upper)
    for v in members:
        offset = v.start - min_pos
        if v.ref.upper() != ref_upper[offset]:
            warnings.warn(
                "combine_cis_snvs: declining %s because %s ref %r disagrees with "
                "reference base %r at %s:%d."
                % (
                    [m.short_description for m in members],
                    v.short_description,
                    v.ref,
                    ref_upper[offset],
                    contig,
                    v.start,
                ),
                stacklevel=3,
            )
            return None
        alt_chars[offset] = v.alt.upper()
    combined = Variant(
        contig=contig,
        start=min_pos,
        ref=ref_upper,
        alt="".join(alt_chars),
        genome=genome,
        allow_extended_nucleotides=any(
            getattr(v, "allow_extended_nucleotides", False) for v in members),
        normalize_contig_names=False,
        chrm_is_rcrs=any(getattr(v, "chrm_is_rcrs", False) for v in members),
    )
    combined.source_variants = tuple(members)
    return combined


def combine_cis_snvs(
        vc: VariantCollection,
        phase_resolver=None,
        *,
        sample: str | None = None,
) -> VariantCollection:
    """Combine cis in-codon SNVs into MNVs without losing transcript or phase interpretations. (reduces)

    Given an explicit phase scope (``phase_resolver`` or ``sample``), replaces
    groups of 2 or 3 pure SNVs that lie in the same CDS codon within a single
    exon and are mutually in cis with a single combined MNV :class:`Variant`.
    Every combined MNV records its source SNVs in ``source_variants`` (ordered
    by ascending genomic coordinate) and retains source metadata provenance.

    Safety and scientific invariants (#368):

    * **Explicit phase/sample scope**: callers must pass either
      ``phase_resolver`` (e.g. :class:`~varcode.VCFPhaseResolver` or
      :class:`~varcode.MolecularPhaseResolver`) or ``sample`` (which constructs
      a :class:`~varcode.VCFPhaseResolver` for that sample), never neither and
      never both.
    * **CDS-relative codon membership**: codon indices are computed from
      ``spliced_offset - transcript.first_start_codon_spliced_offset`` on each
      complete coding transcript, correctly handling both strands and nonzero
      5' UTR lengths.
    * **No cross-intron collapse**: SNVs in a split codon spanning an exon-exon
      junction are never merged into a genomic MNV across the intervening intron.
    * **No lex-first transcript tie-breaking**: if overlapping complete coding
      transcripts disagree on codon grouping (e.g. differing CDS reading frames
      such as ``CDKN2A`` p16INK4a vs p14ARF, or CDS vs UTR/intron), the
      reduction is declined with a warning and all source SNVs pass through
      unchanged.
    * **Reference-validated gap bases**: non-adjacent in-codon SNVs (positions
      1 and 3) fill the intervening base from :func:`~varcode.reference_range`
      and validate all member ``ref`` bases against the reference; missing or
      ambiguous reference bases decline the merge.
    * **Multi-sample and haplotype preservation**: if combining a group would
      drop an uncombined haplotype (e.g. ``1/1`` + ``1|0``) or conflict with
      another sample's genotypes/phase in a multi-sample collection, the merge
      is declined with a warning.
    * **Same-position alternatives**: multi-allelic SNVs sharing a genomic
      position are left untouched.
    * **Idempotence**: already-combined MNVs are not SNVs and carry
      ``source_variants``, so applying :func:`combine_cis_snvs` twice returns
      an unchanged collection on the second call.

    Parameters
    ----------
    vc : VariantCollection
        Input collection. Non-SNVs (indels, MNVs, complex variants, SVs) and
        ``spanning_deletion_records`` pass through unchanged.
    phase_resolver : object, optional
        Phase resolver implementing ``in_cis(v1, v2, transcript=...)``.
        Mutually exclusive with ``sample``.
    sample : str, optional
        Sample name in ``vc`` to scope VCF ``GT``/``PS`` phasing via
        :class:`~varcode.VCFPhaseResolver`. Mutually exclusive with
        ``phase_resolver``.

    Returns
    -------
    VariantCollection
        A new collection with eligible cis in-codon SNVs replaced by combined
        MNVs, preserving ``sources``, merged metadata, and
        ``spanning_deletion_records``.

    Raises
    ------
    ValueError
        If neither ``phase_resolver`` nor ``sample`` is supplied, or if both
        are supplied.
    SampleNotFoundError
        If ``sample`` is supplied and absent from ``vc.samples`` when sample
        data are present.
    """
    if phase_resolver is not None and sample is not None:
        raise ValueError(
            "combine_cis_snvs: pass either phase_resolver or sample, not both.")
    if phase_resolver is None and sample is None:
        raise ValueError(
            "combine_cis_snvs requires an explicit phase_resolver or sample.")

    if phase_resolver is None:
        if vc.has_sample_data() and sample not in vc.samples:
            raise SampleNotFoundError(
                "Sample %r not found. Available samples: %s" % (sample, vc.samples))
        resolver = VCFPhaseResolver(vc, sample)
        is_vcf_resolver = True
    else:
        resolver = phase_resolver
        is_vcf_resolver = isinstance(phase_resolver, VCFPhaseResolver)

    snvs = [v for v in vc if _is_combinable_snv(v)]
    if len(snvs) < 2:
        return vc

    # Skip same-position multi-allelic alternatives so they are never merged
    # into synthetic MNVs or silently collapsed across alt alleles.
    pos_counts: dict[tuple[str, int], int] = {}
    for v in snvs:
        key = (v.contig, v.start)
        pos_counts[key] = pos_counts.get(key, 0) + 1
    multiallelic_positions = {k for k, count in pos_counts.items() if count > 1}
    eligible_snvs = [
        v for v in snvs if (v.contig, v.start) not in multiallelic_positions
    ]
    if len(eligible_snvs) < 2:
        return vc

    placements: dict[Variant, dict[str, tuple[object, tuple]]] = {}
    transcripts_by_id: dict[str, object] = {}
    codon_buckets: dict[tuple[str, int], list[Variant]] = {}

    for v in eligible_snvs:
        try:
            overlapping = v.transcripts
        except InvalidContigError:
            continue
        v_map: dict[str, tuple[object, tuple]] = {}
        for t in overlapping:
            if not getattr(t, "complete", False):
                continue
            placement = _transcript_codon_placement(t, v.start)
            v_map[t.id] = (t, placement)
            transcripts_by_id[t.id] = t
            if placement[0] == "cds":
                codon_buckets.setdefault((t.id, placement[2]), []).append(v)
        placements[v] = v_map

    proposed_groups: set[frozenset[Variant]] = set()
    declined_variants: set[Variant] = set()

    for (t_id, _codon_idx), bucket in codon_buckets.items():
        if len(bucket) < 2:
            continue
        transcript = transcripts_by_id[t_id]
        # Evaluate pairwise phase via the resolver (never raw PS alone).
        n = len(bucket)
        pair_cis: dict[tuple[int, int], bool | None] = {}
        adj: list[list[int]] = [[] for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                ans = query_in_cis(resolver, bucket[i], bucket[j], transcript=transcript)
                pair_cis[(i, j)] = ans
                if ans is True:
                    adj[i].append(j)
                    adj[j].append(i)

        visited = [False] * n
        for start_idx in range(n):
            if visited[start_idx]:
                continue
            stack = [start_idx]
            comp_indices = []
            while stack:
                cur = stack.pop()
                if visited[cur]:
                    continue
                visited[cur] = True
                comp_indices.append(cur)
                for nxt in adj[cur]:
                    if not visited[nxt]:
                        stack.append(nxt)
            if len(comp_indices) < 2:
                continue
            comp_indices.sort()
            comp_variants = [bucket[idx] for idx in comp_indices]

            # Require an all-pairs cis clique (declines contradictory evidence
            # and ambiguous 1/1 + 1|0 + 0|1 triangles).
            is_clique = True
            for a_pos in range(len(comp_indices)):
                for b_pos in range(a_pos + 1, len(comp_indices)):
                    i, j = comp_indices[a_pos], comp_indices[b_pos]
                    if pair_cis[(i, j)] is not True:
                        is_clique = False
                        break
                if not is_clique:
                    break
            if not is_clique:
                warnings.warn(
                    "combine_cis_snvs: declining %s on %s due to contradictory or "
                    "non-transitive cis/trans phase evidence."
                    % ([v.short_description for v in comp_variants], t_id),
                    stacklevel=2,
                )
                declined_variants.update(comp_variants)
                continue

            if len(comp_variants) > 3 or len({v.start for v in comp_variants}) != len(comp_variants):
                warnings.warn(
                    "combine_cis_snvs: declining %s on %s (expected 2-3 distinct positions)."
                    % ([v.short_description for v in comp_variants], t_id),
                    stacklevel=2,
                )
                declined_variants.update(comp_variants)
                continue

            # Never absorb an intervening intron for a split codon crossing an
            # exon-exon junction.
            exon_indices = {placements[v][t_id][1][1] for v in comp_variants}
            min_start = min(v.start for v in comp_variants)
            max_start = max(v.start for v in comp_variants)
            if len(exon_indices) != 1 or (max_start - min_start) > 2:
                declined_variants.update(comp_variants)
                continue

            proposed_groups.add(frozenset(comp_variants))

    if not proposed_groups:
        return vc

    # Drop any proposal containing a variant declined by any transcript, and
    # drop any overlapping proposals (different transcripts grouping an SNV
    # with different partners — never pick lexically first transcript).
    variant_to_groups: dict[Variant, list[frozenset[Variant]]] = {}
    for group in proposed_groups:
        for v in group:
            variant_to_groups.setdefault(v, []).append(group)

    conflicting_groups: set[frozenset[Variant]] = set()
    for group in proposed_groups:
        if any(v in declined_variants for v in group):
            conflicting_groups.add(group)
            continue
        if any(len(variant_to_groups[v]) > 1 for v in group):
            conflicting_groups.add(group)

    if conflicting_groups:
        overlap_groups = [
            g for g in conflicting_groups
            if any(len(variant_to_groups[v]) > 1 for v in g)
        ]
        if overlap_groups:
            warnings.warn(
                "combine_cis_snvs: declining overlapping transcript codon groupings "
                "(incompatible transcript interpretations); source SNVs preserved.",
                stacklevel=2,
            )

    candidate_groups = [g for g in proposed_groups if g not in conflicting_groups]

    # Verify compatibility against ALL overlapping complete coding transcripts
    # (catches differing CDS reading frames even when only 2 SNVs are present,
    # and exon/intron or UTR/CDS splits on alternative isoforms).
    replacement: dict[Variant, tuple[Variant, dict[str, tuple[str | None, object]]]] = {}
    for group in sorted(candidate_groups, key=lambda g: min(v.start for v in g)):
        members = tuple(sorted(group, key=lambda v: v.start))
        overlapping_tids: set[str] = set()
        for v in members:
            overlapping_tids.update(placements.get(v, {}).keys())

        transcript_conflict = None
        for t_id in sorted(overlapping_tids):
            t_places = [
                placements.get(v, {}).get(t_id, (None, ("outside", None, None)))[1]
                for v in members
            ]
            kinds = {p[0] for p in t_places}
            if "cds" in kinds:
                exons = {p[1] for p in t_places}
                codons = {p[2] for p in t_places}
                if kinds != {"cds"} or len(exons) != 1 or len(codons) != 1:
                    transcript_conflict = (
                        "transcript %s assigns incompatible CDS/codon placements %s"
                        % (t_id, t_places)
                    )
                    break
            else:
                if len(kinds) != 1 or (kinds == {"utr"} and len({p[1] for p in t_places}) != 1):
                    transcript_conflict = (
                        "transcript %s splits %s across region/exon boundaries %s"
                        % (t_id, [v.short_description for v in members], t_places)
                    )
                    break
        if transcript_conflict is not None:
            warnings.warn(
                "combine_cis_snvs: declining %s because %s."
                % ([v.short_description for v in members], transcript_conflict),
                stacklevel=2,
            )
            continue

        samples_ok, conflict_reason, merged_gt_by_sample = _check_samples_compatible_for_mnv(
            vc, members, is_vcf_resolver)
        if not samples_ok:
            warnings.warn(
                "combine_cis_snvs: declining %s because %s."
                % ([v.short_description for v in members], conflict_reason),
                stacklevel=2,
            )
            continue

        combined = _build_combined_mnv(members)
        if combined is None:
            continue
        for v in members:
            replacement[v] = (combined, merged_gt_by_sample)

    if not replacement:
        return vc

    out_variants = []
    emitted_combined: set[int] = set()
    for variant in vc:
        entry = replacement.get(variant)
        if entry is None:
            out_variants.append(variant)
            continue
        combined = entry[0]
        if id(combined) in emitted_combined:
            continue
        emitted_combined.add(id(combined))
        out_variants.append(combined)

    out_metadata: dict[str, dict[Variant, object]] = {
        path: {} for path in vc.source_to_metadata_dict
    }
    for path, by_variant in vc.source_to_metadata_dict.items():
        out_by_variant = out_metadata[path]
        for variant, meta in by_variant.items():
            entry = replacement.get(variant)
            if entry is None:
                out_by_variant[variant] = meta
            else:
                combined, merged_gt_by_sample = entry
                if combined not in out_by_variant:
                    merged_meta = _merge_cis_snv_metadata(
                        combined.source_variants, by_variant, merged_gt_by_sample)
                    if merged_meta is not None:
                        out_by_variant[combined] = merged_meta

    return VariantCollection(
        variants=out_variants,
        sources=vc.sources,
        source_to_metadata_dict=out_metadata,
        spanning_deletion_records=vc.spanning_deletion_records,
    )

