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


"""Canonical left-alignment transform for indels (:func:`left_align_indels`)."""

from __future__ import annotations

from ..genome_sequence import reference_range
from ..variant import Variant
from ..variant_collection import VariantCollection

__all__ = ["left_align_indels"]


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
        shifted, partial = _left_align_one(variant, reference_range)
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
        spanning_deletion_records=vc.spanning_deletion_records,
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
        chrm_is_rcrs=getattr(variant, "chrm_is_rcrs", False),
    )
    shifted.source_variants = (variant,)
    return shifted, partial


