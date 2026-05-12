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

"""``VariantCollection`` -> ``VariantCollection`` transforms.

A transform is a pure function: given a ``VariantCollection`` and
optional auxiliary inputs (reference, phase resolver, ...), it returns
a new ``VariantCollection``. Cardinality may be preserved, reduced,
or increased; the docstring of each transform states which.

Every variant in the output carries a ``source_variants`` attribute
recording its provenance:

* ``()`` — pass-through; the variant was not derived from any other.
* ``(v,)`` — derived from a single source variant
  (e.g. left-aligned position).
* ``(v1, v2, ...)`` — combined from multiple source variants
  (e.g. paired breakends merged into one combined SV row).

Transforms compose by application. The intermediate indices a
transform builds (MATEID -> variant, phase set -> variants, ...) live
inside the function and are dropped on return.

Shipped transforms:

* :func:`pair_breakends` — merges MATEID-paired BND rows into a
  single combined ``StructuralVariant``. (reduces)
* :func:`left_align_indels` — shifts indels to their canonical
  leftmost equivalent position. (preserves)

Planned (separate PRs):

* ``combine_cis_snvs`` — combine adjacent in-codon SNVs into MNVs
  given phase. (reduces)
"""

import warnings

from ..genome_sequence import reference_range as _reference_range
from ..structural_variant import StructuralVariant
from ..variant import Variant
from ..variant_collection import VariantCollection


__all__ = ["pair_breakends", "left_align_indels"]


def _is_breakend(variant):
    return getattr(variant, "sv_type", None) == "BND"


def _mate_reference(variant):
    """Return the MATEID/PARID string for a BND variant, or ``None``.

    Looks in the variant's ``info`` dict for ``mateid`` / ``MATEID`` /
    ``parid`` / ``PARID`` — the SV parser writes the lowercase form;
    callers constructing variants directly may use either.
    """
    info = getattr(variant, "info", None) or {}
    for key in ("mateid", "MATEID", "parid", "PARID"):
        value = info.get(key)
        if value is None or value == ".":
            continue
        if isinstance(value, (list, tuple)):
            if not value:
                continue
            value = value[0]
        return str(value)
    return None


def _build_id_indices(vc):
    """Return ``(variant_to_id, id_to_variants)`` from VC source metadata.

    ``id_to_variants`` is a list because a VCF ID is, in principle, only
    unique within a single VCF — when multiple sources merge, two rows
    can technically share an ID. We keep the multiplicity so pairing
    can detect ambiguity instead of silently picking one.
    """
    variant_to_id = {}
    id_to_variants = {}
    for by_variant in vc.source_to_metadata_dict.values():
        for variant, metadata in by_variant.items():
            row_id = (metadata or {}).get("id")
            if row_id is None or row_id == "." or row_id == "":
                continue
            row_id = str(row_id)
            variant_to_id.setdefault(variant, row_id)
            id_to_variants.setdefault(row_id, []).append(variant)
    return variant_to_id, id_to_variants


def _merge_alt_assembly(a, b, a_id, b_id):
    """Pick the combined variant's ``alt_assembly`` from the two halves.

    None+None -> None; one populated -> that one; both equal -> shared;
    both different -> warn and take A's (deterministic via lex source
    ordering).
    """
    aa, ba = a.alt_assembly, b.alt_assembly
    if aa is None and ba is None:
        return None
    if aa is None:
        return ba
    if ba is None:
        return aa
    if aa == ba:
        return aa
    warnings.warn(
        "pair_breakends: paired BND rows %r and %r carry different "
        "alt_assembly sequences; keeping %r's. Inspect "
        "combined.source_variants for both originals."
        % (a_id, b_id, a_id),
        stacklevel=3)
    return aa


def _merge_info(a, b, a_id, b_id):
    """Merge two BND halves' ``info`` dicts into the combined variant's.

    Per-key rules:
      * Both halves agree -> agreed value.
      * Exactly one half has it -> that value.
      * Both halves disagree -> A's value (log at debug-warning level
        if it's a non-trivial field; agree-or-warn would be
        configurable later if needed).

    ``paired_with`` is added pointing at B's row ID for round-trip
    debugging — same info reachable via ``source_variants``, but cheap
    enough to keep as a first-class hint.
    """
    a_info = dict(a.info or {})
    b_info = dict(b.info or {})
    merged = {}
    all_keys = set(a_info) | set(b_info)
    for key in all_keys:
        if key in ("mateid", "MATEID", "parid", "PARID"):
            # The mate pointers no longer apply to the combined row.
            continue
        if key in a_info and key in b_info:
            if a_info[key] == b_info[key]:
                merged[key] = a_info[key]
            else:
                merged[key] = a_info[key]
        elif key in a_info:
            merged[key] = a_info[key]
        else:
            merged[key] = b_info[key]
    if b_id is not None:
        merged["paired_with"] = b_id
    return merged


def _build_combined(a, b, a_id, b_id):
    """Construct the combined ``StructuralVariant`` for a paired BND.

    A is the lex-earlier source row by VCF ID. The combined variant's
    primary endpoint is A's; the mate hint resolves to B's contig/start
    when A doesn't already carry it (which is the common case — A's
    ``mate_contig`` already points at B's contig).
    """
    mate_contig = a.mate_contig if a.mate_contig is not None else b.contig
    mate_start = a.mate_start if a.mate_start is not None else b.start
    combined = StructuralVariant(
        contig=a.contig,
        start=a.start,
        end=a.start,
        sv_type="BND",
        alt=a.symbolic_alt,
        ref=a.ref,
        mate_contig=mate_contig,
        mate_start=mate_start,
        mate_orientation=a.mate_orientation,
        ci_start=a.ci_start,
        ci_end=a.ci_end,
        alt_assembly=_merge_alt_assembly(a, b, a_id, b_id),
        info=_merge_info(a, b, a_id, b_id),
        genome=a.genome,
        normalize_contig_names=False,
    )
    combined.source_variants = (a, b)
    return combined


def _merge_sample_info(a_sample, b_sample, a_id, b_id):
    """Merge per-sample dicts. GT must match across halves; raises on
    disagreement so caller bugs (asymmetric filter application,
    half-only re-genotyping) surface instead of being silently coerced.

    All non-GT FORMAT fields are taken from A.
    """
    if a_sample is None and b_sample is None:
        return None
    if a_sample is None:
        return dict(b_sample)
    if b_sample is None:
        return dict(a_sample)
    merged = {}
    for sample_name in set(a_sample) | set(b_sample):
        a_fields = a_sample.get(sample_name) or {}
        b_fields = b_sample.get(sample_name) or {}
        a_gt = a_fields.get("GT")
        b_gt = b_fields.get("GT")
        if (a_gt is not None and b_gt is not None
                and a_gt != b_gt):
            raise ValueError(
                "pair_breakends: paired BND rows %r and %r disagree on "
                "genotype for sample %r (%r vs %r). Both halves of a "
                "paired breakend describe the same biological event "
                "and should share the same per-sample GT; the mismatch "
                "indicates a caller bug or asymmetric filter "
                "application. Drop one half before pairing if this is "
                "intentional." % (a_id, b_id, sample_name, a_gt, b_gt))
        out_fields = dict(a_fields)
        for key, value in b_fields.items():
            out_fields.setdefault(key, value)
        merged[sample_name] = out_fields
    return merged


def _merge_filter(a_filter, b_filter):
    """Union the two halves' FILTER sets — the stricter filter wins.

    VCF FILTER comes through as a list, a string, or ``None``/``PASS``.
    Normalize to a sorted tuple; ``PASS``-only stays ``("PASS",)``;
    any non-PASS label demotes the row.
    """
    def _norm(f):
        if f is None:
            return frozenset()
        if isinstance(f, (list, tuple, set, frozenset)):
            tokens = [str(x) for x in f if x]
        else:
            tokens = [t for t in str(f).split(";") if t]
        return frozenset(tokens) - {"PASS", ".", ""}
    union = _norm(a_filter) | _norm(b_filter)
    if not union:
        return ["PASS"]
    return sorted(union)


def _merge_metadata(a_meta, b_meta, a_id, b_id):
    """Build the combined variant's per-source metadata entry."""
    if a_meta is None and b_meta is None:
        return None
    a_meta = a_meta or {}
    b_meta = b_meta or {}
    qual_a = a_meta.get("qual")
    qual_b = b_meta.get("qual")
    if qual_a is not None and qual_b is not None:
        qual = min(qual_a, qual_b)
    else:
        qual = qual_a if qual_a is not None else qual_b
    return {
        "id": a_id,
        "qual": qual,
        "filter": _merge_filter(a_meta.get("filter"), b_meta.get("filter")),
        "info": _merge_info_from_metadata(a_meta.get("info"),
                                          b_meta.get("info"),
                                          a_id, b_id),
        "sample_info": _merge_sample_info(a_meta.get("sample_info"),
                                          b_meta.get("sample_info"),
                                          a_id, b_id),
        "alt_allele_index": a_meta.get("alt_allele_index"),
    }


def _merge_info_from_metadata(a_info, b_info, a_id, b_id):
    """Merge metadata-level INFO dicts. Same rules as ``_merge_info``
    but operating on the dict captured at VCF load time (which may be
    richer than the slim ``StructuralVariant.info``)."""
    if a_info is None and b_info is None:
        return None
    a_info = dict(a_info or {})
    b_info = dict(b_info or {})
    merged = {}
    for key in set(a_info) | set(b_info):
        if key in ("MATEID", "PARID"):
            continue
        if key in a_info and key in b_info:
            merged[key] = a_info[key]
        elif key in a_info:
            merged[key] = a_info[key]
        else:
            merged[key] = b_info[key]
    merged["paired_with"] = b_id
    return merged


def pair_breakends(vc):
    """Merge MATEID-paired BND rows into a single combined
    ``StructuralVariant`` per rearrangement event. (reduces)

    For each pair of :class:`~varcode.StructuralVariant` rows where
    row A's ``MATEID`` references row B's VCF ID (and vice versa),
    emit one combined row carrying both endpoints. The combined
    variant's ``source_variants`` attribute holds the two originals.

    Non-BND variants, single-row TRA, and single-ended BNDs (no
    ``MATEID``) pass through unchanged with ``source_variants=()``.

    Pairing rules:

    * Primary key: ``MATEID`` field on each variant's ``info`` dict
      against the VCF row ID stored in the collection's source
      metadata.
    * Alias: ``PARID`` (used by older GRIDSS) is treated as ``MATEID``.
    * Symmetric: row A's ``MATEID`` must equal row B's ID and row B's
      ``MATEID`` must equal row A's ID. Asymmetric references are
      warned and left unpaired.
    * If a ``MATEID`` points to an ID not present in this collection
      (filtered out, chunked load), a warning is emitted and the
      variant passes through unpaired.
    * If three or more rows share a ``MATEID`` group, the whole group
      is left unpaired with a warning (pairing is ambiguous).
    * Already-paired input (``source_variants`` non-empty) passes
      through unchanged — :func:`pair_breakends` is idempotent.

    Metadata merge for paired rows:

    * Genotype: both halves must agree on the per-sample ``GT``. The
      combined variant inherits the shared genotype. Disagreement
      raises :class:`ValueError`.
    * ``alt_assembly``: if exactly one half carries an assembled
      sequence, the combined variant inherits it. If both differ,
      A's wins (deterministic via lex source-ID ordering) with a
      warning.
    * ``filter``: union of both halves' FILTER tokens; ``PASS`` is
      dropped if any non-PASS label is present (stricter wins).
    * ``qual``: minimum of the two halves' quality scores.
    * Other INFO fields: prefer A's; the other half is reachable
      via ``combined.source_variants``.

    Parameters
    ----------
    vc : VariantCollection
        Input collection. May contain a mix of structural and
        non-structural variants.

    Returns
    -------
    VariantCollection
        A new collection. SV rows that were halves of paired BNDs
        are replaced by combined rows; everything else (including
        SNVs/indels/MNVs and unpaired SVs) passes through.
    """
    variant_to_id, id_to_variants = _build_id_indices(vc)

    # For each BND variant, what does it claim as its mate?
    bnd_to_mateid = {}
    for variant in vc:
        if not _is_breakend(variant):
            continue
        if getattr(variant, "source_variants", ()):
            continue
        mate_id = _mate_reference(variant)
        if mate_id is None:
            continue
        bnd_to_mateid[variant] = mate_id

    # How many rows reference each ID as their mate? An ID with
    # in-degree > 1 indicates ambiguity — a clean pair requires the
    # MATEID pointer to be 1:1.
    mate_in_degree = {}
    for mate_id in bnd_to_mateid.values():
        mate_in_degree[mate_id] = mate_in_degree.get(mate_id, 0) + 1

    # Walk BND variants and resolve each to one of: paired (with a
    # specific mate variant), ambiguous (skip with warning naming the
    # connected component), missing-mate (warn), asymmetric (warn).
    replacement = {}     # original variant -> combined variant
    processed = set()    # variants whose pairing decision is final
    ambiguity_warned = set()   # connected-component ids already warned about

    def _component_ids(seed_id):
        """Connected component of mate references reachable from
        ``seed_id``. Used to produce one warning per ambiguous group
        rather than one per row."""
        visited = set()
        stack = [seed_id]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            # Outgoing: variants whose ID is `node` -- find their mateids.
            for v in id_to_variants.get(node, ()):
                mid = bnd_to_mateid.get(v)
                if mid is not None and mid not in visited:
                    stack.append(mid)
            # Incoming: any mateid pointing AT `node`.
            for v, mid in bnd_to_mateid.items():
                if mid == node:
                    own = variant_to_id.get(v)
                    if own is not None and own not in visited:
                        stack.append(own)
        return visited

    for variant in vc:
        if not _is_breakend(variant):
            continue
        if variant in processed:
            continue
        if getattr(variant, "source_variants", ()):
            continue
        own_id = variant_to_id.get(variant)
        if own_id is None:
            continue
        mate_id = bnd_to_mateid.get(variant)
        if mate_id is None:
            continue
        if mate_id not in id_to_variants:
            warnings.warn(
                "pair_breakends: BND %r references MATEID/PARID %r "
                "which is not in this collection; left unpaired. "
                "This typically means the mate row was filter-dropped "
                "or split across chunked loads."
                % (own_id, mate_id),
                stacklevel=2)
            processed.add(variant)
            continue
        # Ambiguity: this variant or its claimed mate has incoming
        # mate-degree > 1, meaning at least one other row also points
        # at the same target. Whole connected component is unsafe.
        if (mate_in_degree.get(own_id, 0) > 1
                or mate_in_degree.get(mate_id, 0) > 1):
            component = _component_ids(own_id)
            component_key = frozenset(component)
            if component_key not in ambiguity_warned:
                ambiguity_warned.add(component_key)
                warnings.warn(
                    "pair_breakends: BND ID group %r has ambiguous "
                    "MATEID/PARID references (at least one ID is "
                    "referenced by more than one row); left unpaired. "
                    "Each rearrangement should produce exactly two BND "
                    "rows with 1:1 mate pointers."
                    % sorted(component),
                    stacklevel=2)
            # Mark every variant in the component as processed so we
            # don't re-warn from a different vantage point.
            for cid in component:
                for v in id_to_variants.get(cid, ()):
                    processed.add(v)
            continue
        # Identify the specific mate variant. With in-degree 1 there's
        # exactly one BND row at `mate_id`.
        mate_candidates = [
            v for v in id_to_variants[mate_id]
            if _is_breakend(v)
        ]
        if not mate_candidates:
            warnings.warn(
                "pair_breakends: BND %r references MATEID %r but that "
                "ID is held by a non-BND variant; left unpaired."
                % (own_id, mate_id),
                stacklevel=2)
            processed.add(variant)
            continue
        mate = mate_candidates[0]
        # Symmetric check: mate must point back at us.
        if bnd_to_mateid.get(mate) != own_id:
            warnings.warn(
                "pair_breakends: BND %r mate reference is asymmetric "
                "(this row -> %r but mate -> %r); left unpaired."
                % (own_id, mate_id, bnd_to_mateid.get(mate)),
                stacklevel=2)
            processed.add(variant)
            processed.add(mate)
            continue
        # Clean pair.
        a, b = sorted([variant, mate], key=lambda v: variant_to_id[v])
        a_id = variant_to_id[a]
        b_id = variant_to_id[b]
        combined = _build_combined(a, b, a_id, b_id)
        replacement[a] = combined
        replacement[b] = combined
        processed.add(a)
        processed.add(b)

    # Build output VC. Pass-through variants land first-occurrence; paired
    # variants are replaced by their combined variant the first time either
    # half is encountered, then the second half is skipped to preserve
    # ordering and cardinality.
    out_variants = []
    emitted_combined = set()
    for variant in vc:
        combined = replacement.get(variant)
        if combined is None:
            out_variants.append(variant)
            continue
        if id(combined) in emitted_combined:
            continue
        emitted_combined.add(id(combined))
        out_variants.append(combined)

    # Build output metadata. Pass-through variants reuse the existing
    # metadata entries by reference; combined variants get a freshly
    # merged entry written to every source path that contained either
    # half.
    out_metadata = {path: {} for path in vc.source_to_metadata_dict}
    for path, by_variant in vc.source_to_metadata_dict.items():
        out_by_variant = out_metadata[path]
        for variant, meta in by_variant.items():
            combined = replacement.get(variant)
            if combined is None:
                out_by_variant[variant] = meta
            elif combined not in out_by_variant:
                a, b = combined.source_variants
                a_meta = by_variant.get(a)
                b_meta = by_variant.get(b)
                a_id = variant_to_id[a]
                b_id = variant_to_id[b]
                out_by_variant[combined] = _merge_metadata(
                    a_meta, b_meta, a_id, b_id)

    return VariantCollection(
        variants=out_variants,
        sources=vc.sources,
        source_to_metadata_dict=out_metadata,
    )


# ====================================================================
# left_align_indels
# ====================================================================


def left_align_indels(vc):
    """Shift indels to their canonical leftmost equivalent position. (preserves)

    Indels in homopolymer or short-tandem-repeat regions can be
    represented at any of several equivalent positions — ``CTT->T``
    inside a ``CT``-repeat means the same biological event as
    ``CT->_`` two positions to the left. Tools that compare variants
    by ``(contig, start, ref, alt)`` see those representations as
    distinct calls. Left-alignment normalizes to a single canonical
    representation per indel: the leftmost equivalent position.

    The algorithm is the standard variant-normalization left-shift
    used by ``bcftools norm`` and GATK ``LeftAlignAndTrimVariants``,
    applied as an opt-in ``VariantCollection -> VariantCollection``
    transform rather than baked into VCF load.

    Reference sequence is read via the genome the variants carry —
    no explicit ``reference`` parameter. Coverage tiers (see
    :mod:`varcode.genome_sequence`):

    * **Chromosome FASTA attached** (via :class:`varcode.Genome`'s
      ``fasta`` slot): indels everywhere shift to canonical
      positions, including in introns and intergenic regions.
    * **No FASTA** (default pyensembl install): indels fully within
      an exon shift via the transcript cDNA fallback. Intronic and
      intergenic indels pass through unchanged. Indels that *start*
      exonic but would shift across an exon boundary stop at the
      boundary and carry ``info["left_align_partial"] = True``.

    Parameters
    ----------
    vc : VariantCollection
        Input collection. May contain a mix of SNVs, MNVs, indels,
        complex variants, and SVs — only pure indels (length-different
        REF/ALT with one side empty after suffix trimming) are
        considered for shifting. Everything else passes through.

    Returns
    -------
    VariantCollection
        A new collection with indels at their canonical leftmost
        positions. Variants that shifted carry
        ``source_variants=(original,)``; everything else passes through
        with ``source_variants=()``.

    See :doc:`/transforms` for the behavior table covering all
    six (location × FASTA-attached) combinations and the metadata
    fields the transform writes.

    Examples
    --------

    >>> from varcode.transforms import left_align_indels
    >>> normalized = left_align_indels(vc)  # doctest: +SKIP
    """
    # original variant -> (shifted variant, bounded_by_coverage flag).
    # The flag rides with the shifted variant so we don't need a
    # second data structure keyed on id() — the lifecycle of the
    # flag and the lifecycle of the shifted variant are identical.
    replacement = {}

    for variant in vc:
        if not variant.is_indel:
            continue
        shifted, partial = _left_align_one(variant, _reference_range)
        if shifted is variant:
            continue
        replacement[variant] = (shifted, partial)

    if not replacement:
        return vc

    out_variants = [
        replacement[v][0] if v in replacement else v
        for v in vc
    ]

    out_metadata = {}
    for path, by_variant in vc.source_to_metadata_dict.items():
        out_by_variant = {}
        for variant, meta in by_variant.items():
            entry = replacement.get(variant)
            if entry is None:
                out_by_variant[variant] = meta
            else:
                shifted, partial = entry
                new_meta = dict(meta) if meta else {}
                info = dict(new_meta.get("info") or {})
                info["original_start"] = variant.start
                if partial:
                    info["left_align_partial"] = True
                new_meta["info"] = info
                out_by_variant[shifted] = new_meta
        out_metadata[path] = out_by_variant

    return VariantCollection(
        variants=out_variants,
        sources=vc.sources,
        source_to_metadata_dict=out_metadata,
    )


def _left_align_one(variant, reference_range_fn):
    """Shift one indel to its canonical leftmost position.

    Returns ``(variant_or_shifted, partial_shift)``. ``partial_shift``
    is True iff the loop exited because reference coverage ran out
    (Tier 2 without a FASTA), not because the reference disagreed —
    consumers can flag the result so downstream tools know a fuller
    shift would have been possible with a chromosome FASTA.

    varcode normalizes indels in :class:`Variant.__init__` so that
    one side of (ref, alt) is empty after shared prefix/suffix
    trimming — a deletion has ``ref != ''`` and ``alt == ''``, an
    insertion has ``ref == ''`` and ``alt != ''``. The check
    position differs between the two:

    * **Deletion** ``ref=payload, alt='', start=N``: removes
      ``payload`` from positions ``[N, N+len-1]``. Shift left valid
      iff ``reference[N-1] == payload[-1]`` — the base entering on
      the left equals the base leaving on the right.
    * **Insertion** ``ref='', alt=payload, start=N``: inserts
      ``payload`` after position ``N``. Shift left valid iff
      ``reference[N] == payload[-1]`` — the anchor base equals the
      last inserted base, allowing the payload to rotate around it.
    """
    if not variant.is_indel:
        return variant, False

    ref, alt = variant.ref, variant.alt
    if ref and alt:
        # Complex variant (both sides non-empty after normalization).
        # Shouldn't happen for is_indel, but defensively pass through.
        return variant, False

    is_insertion = (not ref) and alt
    payload = alt if is_insertion else ref
    if not payload:
        return variant, False

    start = variant.start
    new_start = start
    new_payload = payload
    partial = False

    while new_start > 1:
        # Insertion checks at new_start (the anchor); deletion checks
        # at new_start - 1 (the base just before the deletion).
        check_pos = new_start if is_insertion else (new_start - 1)
        prev = reference_range_fn(
            variant.genome, variant.contig, check_pos, check_pos)
        if not prev:
            # Reference coverage ran out at this position. If we
            # already moved, this is a partial shift — a chromosome
            # FASTA would have allowed continuing.
            partial = (new_start != start)
            break
        if prev.upper() != new_payload[-1].upper():
            break  # natural stop — reference disagrees
        new_payload = prev.upper() + new_payload[:-1]
        new_start -= 1

    if new_start == start:
        return variant, False

    new_ref = new_payload if not is_insertion else ""
    new_alt = new_payload if is_insertion else ""
    shifted = Variant(
        contig=variant.contig,
        start=new_start,
        ref=new_ref,
        alt=new_alt,
        genome=variant.genome,
        allow_extended_nucleotides=getattr(
            variant, "allow_extended_nucleotides", False),
        normalize_contig_names=False,
    )
    shifted.source_variants = (variant,)
    return shifted, partial
