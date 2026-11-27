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

"""Data model for representing the result of applying one or more
variants to a reference transcript.

:class:`TranscriptEdit` and :class:`MutantTranscript` let effect
annotation materialize the mutant sequence, translate it, and compare
it to the reference protein instead of reasoning about offsets against
the reference. :func:`apply_variant_to_transcript` and
:func:`apply_variants_to_transcript` build them; the protein-diff
annotator, splice outcomes, phasing, and germline-aware annotation
consume them.
"""

from dataclasses import dataclass, field
from typing import Optional, Tuple

from serializable import DataclassSerializable

from .effects.transcript_helpers import has_complete_coding_sequence, poly_a_completion
from .nucleotides import reverse_complement


@dataclass
class AssembledAllele(DataclassSerializable):
    """Sequence provider for an externally assembled :class:`ReferenceSegment`.

    Parameters
    ----------
    sequence : str
        Oriented 5'-to-3' nucleotide sequence for the assembled segment.
    """

    sequence: str


# Backwards-compatible alias for existing serialized/pickled graphs.
_AssembledAllele = AssembledAllele


@dataclass(frozen=True)
class TranscriptEdit(DataclassSerializable):
    """A single edit applied to a transcript's spliced mRNA.

    Coordinates are zero-based and refer to positions in the reference
    transcript's cDNA (i.e. :attr:`pyensembl.Transcript.sequence`).
    ``cdna_end`` is exclusive; an insertion has
    ``cdna_start == cdna_end``.
    """
    cdna_start: int
    cdna_end: int
    alt_bases: str
    """Nucleotides inserted at ``cdna_start`` (may be empty for a
    pure deletion)."""

    source_variant: Optional[object] = None
    """Back-pointer to the :class:`varcode.Variant` that motivated
    this edit, when known. Optional so synthetic edits (e.g. from
    an RNA-evidence import that doesn't carry a DNA variant) are
    representable."""

    def __post_init__(self):
        if self.cdna_start < 0:
            raise ValueError(
                "cdna_start must be non-negative, got %d" % self.cdna_start)
        if self.cdna_end < self.cdna_start:
            raise ValueError(
                "cdna_end (%d) must be >= cdna_start (%d)" % (
                    self.cdna_end, self.cdna_start))

    @property
    def is_insertion(self) -> bool:
        return self.cdna_start == self.cdna_end and len(self.alt_bases) > 0

    @property
    def is_deletion(self) -> bool:
        return self.cdna_end > self.cdna_start and len(self.alt_bases) == 0

    @property
    def length_delta(self) -> int:
        """Net change in transcript length after this edit
        (positive for insertion-net, negative for deletion-net)."""
        return len(self.alt_bases) - (self.cdna_end - self.cdna_start)


@dataclass(frozen=True)
class ReferenceSegment(DataclassSerializable):
    """A contiguous run of reference sequence used to assemble a
    :class:`MutantTranscript` for a structural variant (#252).

    A fusion is described by two segments — one from each transcript
    involved, each spanning the contributing exons. A translocation
    to intergenic space is one transcript segment plus a
    ``GenomicInterval`` segment. A large inversion that inverts a
    middle chunk of a single transcript is three segments
    (upstream-forward, reverse-complement middle, downstream-forward).

    The ``source`` is typed as :class:`object` rather than pinned to
    :class:`pyensembl.Transcript` so callers can plug in any object
    that exposes ``sequence`` / ``start`` / ``end`` / ``strand``. This
    keeps the door open for:

    * Personalized / full-genome assemblies — the source can be a
      patient-specific contig object.
    * Long-read-derived local assemblies — the source can be an
      opaque wrapper around a FASTA fragment.
    * In-memory ``GenomicInterval`` stand-ins for intergenic space.

    Coordinates are 0-based, ``start`` inclusive, ``end`` exclusive,
    matching Python slice semantics. The annotator that materializes
    ``cdna_sequence`` concatenates ``source.sequence[start:end]``
    across segments in order (applying reverse-complement when
    ``strand == "-"``).
    """

    source: object
    """Transcript, ``GenomicInterval``, or caller-supplied sequence
    provider. The only contract is: the segment can be rendered into
    a nucleotide string by the producer of the MutantTranscript."""

    start: int
    """0-based inclusive start offset in ``source``'s own coordinate
    system. For a :class:`pyensembl.Transcript` source, this is a
    cDNA offset; for a ``GenomicInterval`` source, this is a genomic
    position."""

    end: int
    """0-based exclusive end offset. ``end - start`` is the segment's
    length in nucleotides."""

    strand: str = "+"
    """``"+"`` = forward, ``"-"`` = reverse-complement. Inversions
    flip strand; the rest of the segments stay ``"+"``."""

    label: str = ""
    """Optional human-readable tag (e.g. ``"5p_partner"``,
    ``"3p_partner"``, ``"intergenic_readthrough"``). Carried through
    in serialization; ignored by translation logic."""

    def __post_init__(self):
        if self.end < self.start:
            raise ValueError(
                "ReferenceSegment end (%d) must be >= start (%d)"
                % (self.end, self.start))
        if self.strand not in ("+", "-"):
            raise ValueError(
                "ReferenceSegment strand must be '+' or '-', got %r"
                % (self.strand,))

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class MutantTranscript(DataclassSerializable):
    """A reference transcript (or assembled set of reference segments) with
    variant-derived edits applied, optional materialized cDNA/protein
    sequences, splice-junction coordinates, and open rearrangement templates.

    **Two shapes:**

    * **Point-variant shape** (``reference_transcript`` is set,
      ``reference_segments`` is ``None``) — derived from a single reference
      transcript by applying zero or more :class:`TranscriptEdit` objects.
    * **Structural-variant shape** (``reference_segments`` is set,
      ``reference_transcript`` may be ``None`` or the primary / 5'-partner
      transcript) — assembled by concatenating :class:`ReferenceSegment` slices
      in order (fusions, intergenic translocations, inversions, deletions,
      duplications, or supplied ``alt_assembly`` contigs).

    **DNA-to-RNA junction & template properties:**

    * :attr:`is_complete` — ``True`` when a full spliced ``cdna_sequence`` is
      resolved; ``False`` for open rearrangements without a supplied assembly.
    * :attr:`novel_junctions` / :attr:`back_junctions` / :attr:`wildtype_junctions` —
      1-based ascending genomic intron tuples ``(contig, start, end, strand)``
      directly matching BAM CIGAR ``N`` intervals (or tandem-duplication
      chimeric back-junction spans).
    * :attr:`mate_region` / :attr:`exit_frame` / :meth:`template_sequence` —
      oriented mate search window, 5' anchor reading frame, and on-demand
      pre-mRNA template materializer for open rearrangements (large insertions
      and gene-to-intergenic translocations).
    """

    reference_transcript: Optional[object] = None
    """The :class:`pyensembl.Transcript` this mutant is derived from
    (point-variant shape), or the primary / 5'-partner transcript
    (SV shape). ``None`` when the SV has no canonical primary
    transcript (e.g. intergenic-to-intergenic BNDs). Not typed
    tightly here so :mod:`pyensembl` isn't a hard import dependency
    for anyone who just wants the dataclass."""

    edits: Tuple[TranscriptEdit, ...] = field(default_factory=tuple)
    """Edits applied to produce this mutant, sorted by
    :attr:`TranscriptEdit.cdna_start`. Empty tuple means the mutant
    is identical to the reference (or, for SV shape, the assembled
    segments carry the rearrangement directly without further
    point-level edits)."""

    reference_segments: Optional[Tuple[ReferenceSegment, ...]] = None
    """Ordered tuple of :class:`ReferenceSegment` objects that, when
    concatenated in order (applying reverse-complement to ``-``
    strand segments), produce the mutant cDNA. ``None`` for the
    point-variant shape; a fusion's segments would be
    ``(5p_partner_segment, 3p_partner_segment)``. Coordinates are
    in each segment's own reference system. A partial structural model may
    describe only a retained reference fragment: its producer identifies this
    in ``evidence`` and leaves the full ``cdna_sequence`` unknown. Concatenating
    those partial segments does not establish a full allele."""

    cdna_sequence: Optional[str] = None
    """The mutated spliced mRNA, when computed. ``None`` if the
    producer hasn't materialized it yet."""

    mutant_protein_sequence: Optional[str] = None
    """The translated mutant protein, stopping at the first stop codon.
    ``None`` when unavailable or when the edit set does not produce a
    coherent ORF (e.g. start-codon loss or an unresolved open rearrangement)."""

    annotator_name: str = "unknown"
    """Name of the :class:`EffectAnnotator` (or other producer) that
    created this ``MutantTranscript``. Used as provenance in
    serialization and for A/B comparisons."""

    evidence: Optional[dict] = None
    """Optional producer-specific evidence (RNA read counts, Isovar
    fragment ids, SpliceAI scores, long-read assembly metadata).
    Shape is annotator-specific and not part of the stable contract;
    consumers that care about a particular evidence shape should
    type-check it at the call site."""

    def __post_init__(self):
        if self.edits:
            starts = [e.cdna_start for e in self.edits]
            if any(starts[i] > starts[i + 1] for i in range(len(starts) - 1)):
                raise ValueError(
                    "MutantTranscript.edits must be sorted by cdna_start")
        if (self.reference_transcript is None
                and self.reference_segments is None):
            raise ValueError(
                "MutantTranscript requires either reference_transcript "
                "(point-variant shape) or reference_segments (SV shape).")

    @classmethod
    def from_sequence(cls, sequence, *, reference_transcript=None,
                      mutant_protein_sequence=None, annotator_name="unknown",
                      evidence=None, label="observed_sequence"):
        """Wrap an oriented external sequence without reconstructing it.

        Sequence is already in transcript 5'-to-3' order. This constructor
        does not infer an ORF, splice structure, or sequence completeness.
        Producers must describe partial observations in ``evidence``.
        """
        return cls(
            reference_transcript=reference_transcript,
            reference_segments=(ReferenceSegment(
                source=AssembledAllele(sequence), start=0, end=len(sequence),
                label=label),),
            cdna_sequence=sequence,
            mutant_protein_sequence=mutant_protein_sequence,
            annotator_name=annotator_name,
            evidence=evidence)

    @property
    def is_identical_to_reference(self) -> bool:
        """True if no edits were applied AND there are no
        reference-rearranging segments. Does NOT check
        ``cdna_sequence`` / ``mutant_protein_sequence`` — a producer
        can legitimately carry an identical sequence with zero edits
        and a single identity segment."""
        return (len(self.edits) == 0
                and self.reference_segments is None)

    @property
    def is_structural(self) -> bool:
        """True when this mutant was assembled from
        :attr:`reference_segments` (SV shape) rather than applying
        :attr:`edits` to a single reference transcript."""
        return self.reference_segments is not None

    @property
    def total_length_delta(self) -> int:
        """Sum of :attr:`TranscriptEdit.length_delta` across all
        edits — how much longer or shorter the mutant cDNA is than
        the reference (point-variant shape). For SV shape, returns
        0; the length of an assembled cDNA is the sum of segment
        lengths, not a delta against a single reference."""
        return sum(e.length_delta for e in self.edits)

    @property
    def is_complete(self) -> bool:
        """True when a full spliced ``cdna_sequence`` has been resolved."""
        return self.cdna_sequence is not None

    @property
    def novel_junctions(self) -> Tuple[Tuple[str, int, int, str], ...]:
        """Novel forward splice junctions as 1-based ascending genomic intron tuples.

        Each tuple is ``(contig, intron_start, intron_end, strand)`` with
        ``intron_start <= intron_end`` on both ``"+"`` and ``"-"`` strands,
        directly matching BAM CIGAR ``N`` intervals.
        """
        raw = (self.evidence or {}).get("novel_junctions") or ()
        return tuple((str(c), int(s), int(e), str(st)) for c, s, e, st in raw)

    @property
    def back_junctions(self) -> Tuple[Tuple[str, int, int, str], ...]:
        """Tandem-duplication back-junctions as 1-based ascending genomic spans.

        Each tuple is ``(contig, span_start, span_end, strand)`` with
        ``span_start <= span_end`` (e.g. exon 25 -> exon 18 in ``EGFR``-KDD,
        which appears in RNA alignments as chimeric/supplementary reads rather
        than a forward CIGAR ``N`` splice).
        """
        raw = (self.evidence or {}).get("back_junctions") or ()
        return tuple((str(c), int(s), int(e), str(st)) for c, s, e, st in raw)

    @property
    def wildtype_junctions(self) -> Tuple[Tuple[str, int, int, str], ...]:
        """Canonical wild-type introns competing near an open rearrangement breakpoint.

        Each tuple uses the same 1-based ascending ``(contig, start, end, strand)``
        convention as :attr:`novel_junctions`.
        """
        raw = (self.evidence or {}).get("wildtype_junctions") or ()
        return tuple((str(c), int(s), int(e), str(st)) for c, s, e, st in raw)

    @property
    def mate_region(self) -> Optional[Tuple[str, int, int, str]]:
        """Oriented 1-based ascending ``(contig, start, end, strand)`` mate search window.

        Populated for intergenic/breakend translocations so downstream RNA
        tools know where to look for cryptic partner exons; ``None`` otherwise.
        """
        raw = (self.evidence or {}).get("mate_search_region")
        if raw is None:
            return None
        c, s, e, st = raw
        return (str(c), int(s), int(e), str(st))

    @property
    def exit_frame(self) -> Optional[int]:
        """CDS reading frame (``0``, ``1``, or ``2``) exiting the 5' intact anchor, or ``None``."""
        val = (self.evidence or {}).get("exit_frame")
        return None if val is None else int(val)

    def template_sequence(self, sequence_provider) -> str:
        """Materialize the oriented 5'-to-3' mutant pre-mRNA template for RNA alignment.

        For open rearrangements (large insertions and gene-to-intergenic
        translocations), concatenates the ordered genomic and literal segments
        in ``evidence["template_segments"]`` in transcript 5'-to-3' orientation
        (applying reverse-complement to ``"-"`` strand segments).

        Parameters
        ----------
        sequence_provider : callable
            Callable ``(contig: str, start_1based: int, end_1based: int) -> str``
            returning forward-strand reference genomic nucleotides (inclusive).

        Returns
        -------
        str
            Oriented 5'-to-3' mutant pre-mRNA template sequence.

        Raises
        ------
        ValueError
            If this ``MutantTranscript`` carries no template segments.
        """
        ev = self.evidence or {}
        segments = ev.get("template_segments")
        if not segments:
            raise ValueError(
                "MutantTranscript has no template_segments for template materialization"
            )
        inserted = ev.get("inserted_sequence")
        pieces = []
        for label, contig, start, end, strand in segments:
            if label == "inserted_sequence":
                if not inserted:
                    continue
                seq = str(inserted)
                pieces.append(reverse_complement(seq) if strand == "-" else seq)
            else:
                genomic = str(sequence_provider(str(contig), int(start), int(end)))
                pieces.append(reverse_complement(genomic) if strand == "-" else genomic)
        return "".join(pieces)


# ---------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------


def resolve_variant_edit(variant, transcript, full_sequence):
    """Resolve a single :class:`~varcode.Variant` to a :class:`TranscriptEdit`.

    Parameters
    ----------
    variant : Variant
        Point variant or indel to map onto ``transcript``.
    transcript : Transcript-like
        Target spliced transcript exposing ``exons`` and ``on_backward_strand``.
    full_sequence : str
        Full spliced cDNA sequence of ``transcript``.

    Returns
    -------
    tuple of (TranscriptEdit, int) or None
        ``(edit, cdna_anchor_offset)`` where ``cdna_anchor_offset`` is the
        pre-edit cDNA offset used by translation logic to decide whether the
        variant lands after the CDS start, or ``None`` when the variant cannot
        be cleanly resolved (multi-exon span, out-of-range, or reference
        mismatch).
    """
    from .effects.transcript_helpers import interval_offset_on_transcript

    variant_start = variant.trimmed_base1_start
    variant_end = variant.trimmed_base1_end
    # Variant must overlap exactly one exon — we don't handle
    # splice-junction-spanning edits at this stage.
    overlapping = [
        ex for ex in transcript.exons
        if variant_start >= ex.start and variant_end <= ex.end
    ]
    if len(overlapping) != 1:
        return None
    try:
        cdna_offset = interval_offset_on_transcript(
            variant_start, variant_end, transcript)
    except (ValueError, KeyError):
        return None
    genome_ref = variant.trimmed_ref
    genome_alt = variant.trimmed_alt
    if transcript.on_backward_strand:
        cdna_ref = reverse_complement(genome_ref)
        cdna_alt = reverse_complement(genome_alt)
    else:
        cdna_ref = genome_ref
        cdna_alt = genome_alt
    n_ref = len(cdna_ref)
    expected_ref = full_sequence[cdna_offset:cdna_offset + n_ref]
    if cdna_ref != expected_ref:
        return None
    # Pure insertions: VCF convention places the alt AFTER the anchor
    # base in genomic coordinates. Forward strand → insertion at
    # cdna_offset + 1. Reverse strand → insertion at cdna_offset (the
    # cDNA runs 3'→5' in genomic space, so "after" in genomic coords
    # is "before" in cDNA coords).
    if n_ref == 0:
        insert_at = (
            cdna_offset if transcript.on_backward_strand else cdna_offset + 1)
        edit = TranscriptEdit(
            cdna_start=insert_at,
            cdna_end=insert_at,
            alt_bases=cdna_alt,
            source_variant=variant,
        )
    else:
        edit = TranscriptEdit(
            cdna_start=cdna_offset,
            cdna_end=cdna_offset + n_ref,
            alt_bases=cdna_alt,
            source_variant=variant,
        )
    return (edit, cdna_offset)


def mutant_cds_start(transcript, edits=()):
    """Map the reference CDS start past edits wholly upstream of it.

    Parameters
    ----------
    transcript : Transcript-like
        Reference transcript exposing ``start_codon_spliced_offsets``.
    edits : iterable of TranscriptEdit, default ()
        Sorted reference-coordinate cDNA edits.

    Returns
    -------
    int
        Zero-based cDNA offset of the start codon in the edited transcript.
    """
    from .effects.selenocysteine import edit_shift

    reference_start = min(transcript.start_codon_spliced_offsets)
    return reference_start + edit_shift(reference_start, edits)


def retained_cds_start(model, transcript):
    """Map all three annotated start-codon bases across ``model.reference_segments``.

    Parameters
    ----------
    model : MutantTranscript
        Segment-assembled mutant transcript model.
    transcript : Transcript-like
        Reference transcript supplying ``start_codon_spliced_offsets``.

    Returns
    -------
    int or None
        Zero-based offset in ``model.cdna_sequence`` where all three bases of
        ``transcript``'s start codon remain contiguous on the ``+`` segment
        strand (including across split segment boundaries), or ``None`` if lost.
    """
    positions = [set(), set(), set()]
    start = min(transcript.start_codon_spliced_offsets)
    offset = 0
    for segment in model.reference_segments or ():
        if segment.source == transcript and segment.strand == "+":
            for i in range(3):
                if segment.start <= start + i < segment.end:
                    positions[i].add(offset + start + i - segment.start)
        offset += segment.length
    return next((p for p in sorted(positions[0])
                 if p + 1 in positions[1] and p + 2 in positions[2]), None)


def translate_from_cds(mutant_cdna, transcript, edits=()):
    """Translate ``mutant_cdna`` from the canonical CDS start to the first stop.

    Parameters
    ----------
    mutant_cdna : str
        Edited spliced cDNA sequence.
    transcript : Transcript-like
        Reference transcript supplying genetic code, start offsets, and Sec sites.
    edits : iterable of TranscriptEdit, default ()
        Reference-coordinate edits applied to obtain ``mutant_cdna``.

    Returns
    -------
    str or None
        Translated protein sequence (with retained recognized initiators decoded
        as ``M``), or ``None`` if the CDS start lies past ``mutant_cdna`` or
        translation fails.
    """
    from .effects.codon_tables import (
        codon_table_for_transcript,
        initiator_as_met,
        translate_sequence,
    )
    from .effects.selenocysteine import edited_selenocysteine
    cds_start = mutant_cds_start(transcript, edits)
    if cds_start >= len(mutant_cdna):
        return None
    codon_table = codon_table_for_transcript(transcript)
    coding = mutant_cdna[cds_start:]
    # Polyadenylation finishes a partial last codon (MT-ND1's stop is TA + A).
    padded = coding + poly_a_completion(transcript, coding)
    truncated = padded[:(len(padded) // 3) * 3]
    try:
        protein = translate_sequence(
            truncated, codon_table=codon_table, to_stop=True,
            selenocysteine={pos - cds_start for pos in
                            edited_selenocysteine(transcript, edits)})
    except ValueError:
        return None
    reference_start = min(transcript.start_codon_spliced_offsets)
    start_removed = any(
        e.cdna_start <= reference_start < e.cdna_end and e.length_delta
        for e in edits)
    if not start_removed:
        protein = initiator_as_met(protein, coding[:3], codon_table)
    return protein


def apply_variant_to_transcript(variant, transcript):
    """Construct a :class:`MutantTranscript` by applying ``variant``
    to ``transcript``'s spliced cDNA.

    Returns a :class:`MutantTranscript` whose ``cdna_sequence`` is
    populated, plus ``mutant_protein_sequence`` when the variant
    lies after the start codon (so translation from the canonical
    start is well-defined). The codon table is selected from the
    transcript's contig — mitochondrial transcripts use NCBI table
    2 automatically (see :func:`varcode.effects.codon_tables.codon_table_for_transcript`).

    Returns ``None`` when the variant can't be cleanly applied:

    * Transcript is not protein-coding or is incomplete.
    * Variant doesn't overlap the transcript at all.
    * Variant spans more than one exon (splice-junction-crossing
      variants need the splice-aware path; not handled here).
    * Reference allele doesn't match the transcript's cDNA at the
      computed offset.

    A caller annotating only this reference-relative edit may use the fast
    :class:`EffectAnnotator` when the model is unavailable. A failed patient
    haplotype comparison must instead remain unresolved: dropping supplied
    cis germline edits would change the comparison being made.
    """
    from pyensembl import Transcript

    if not isinstance(transcript, Transcript):
        return None
    if not transcript.is_protein_coding or not has_complete_coding_sequence(transcript):
        return None

    full_sequence = str(transcript.sequence)
    resolved = resolve_variant_edit(variant, transcript, full_sequence)
    if resolved is None:
        return None
    edit, cdna_offset = resolved
    mutant_cdna = (
        full_sequence[:edit.cdna_start]
        + edit.alt_bases
        + full_sequence[edit.cdna_end:]
    )
    mutant_protein = None
    cds_start = min(transcript.start_codon_spliced_offsets)
    if cdna_offset >= cds_start:
        mutant_protein = translate_from_cds(mutant_cdna, transcript, (edit,))
    return MutantTranscript(
        reference_transcript=transcript,
        edits=(edit,),
        cdna_sequence=mutant_cdna,
        mutant_protein_sequence=mutant_protein,
        annotator_name="protein_diff",
    )


def apply_variants_to_transcript(variants, transcript):
    """Apply a list of variants to a single transcript, yielding one
    :class:`MutantTranscript` that carries all the resulting edits
    (#269). Used for haplotype-aware joint effect prediction — cis
    variants on the same transcript become one combined mutant
    rather than N independent per-variant mutants.

    Edits are applied in cDNA-coordinate order (highest offset
    first, so earlier offsets aren't shifted) to ``transcript``'s
    spliced cDNA. Returns ``None`` when any of the usual
    single-variant preconditions fail (non-coding, incomplete, etc.),
    or when the provided variants conflict — i.e. claim to edit
    overlapping cDNA ranges. The caller is responsible for falling
    back to per-variant effects in that case.

    ``mutant_protein_sequence`` is populated when at least one edit
    lands after the canonical CDS start; the joint cDNA is translated
    from there to the first stop.

    Order of ``variants`` doesn't matter — edits are sorted by cDNA
    offset internally.
    """
    from pyensembl import Transcript

    if not isinstance(transcript, Transcript):
        return None
    if not transcript.is_protein_coding or not has_complete_coding_sequence(transcript):
        return None

    full_sequence = str(transcript.sequence)
    resolved = []
    for variant in variants:
        r = resolve_variant_edit(variant, transcript, full_sequence)
        if r is None:
            return None
        resolved.append(r)
    # Sort by cDNA start so overlap detection is a simple linear scan
    # and we can apply high-to-low without shifting earlier edits.
    resolved.sort(key=lambda pair: pair[0].cdna_start)
    for i in range(len(resolved) - 1):
        e1 = resolved[i][0]
        e2 = resolved[i + 1][0]
        # Range overlap.
        if e1.cdna_end > e2.cdna_start:
            return None
        # An insertion at the same cdna_start as the next edit is
        # order-dependent (whether the inserted bases survive a
        # following deletion, or where they sit relative to a
        # substitution at the same offset, depends on apply order).
        # In contrast, adjacent edits at e1.cdna_end == e2.cdna_start
        # where e1 is NOT an insertion are unambiguous and allowed:
        # high-to-low application applies e2 first, leaving the slice
        # below e2.cdna_start untouched, so e1 then operates on its
        # original range and the two edits compose to a single
        # deterministic result regardless of input order.
        if e1.cdna_end == e2.cdna_start and e1.cdna_start == e1.cdna_end:
            return None
    # Apply edits from the highest cDNA offset down so earlier offsets
    # aren't affected by upstream edits' length changes.
    mutant_cdna = full_sequence
    for edit, _ in sorted(
            resolved, key=lambda pair: pair[0].cdna_start, reverse=True):
        mutant_cdna = (
            mutant_cdna[:edit.cdna_start]
            + edit.alt_bases
            + mutant_cdna[edit.cdna_end:]
        )
    mutant_protein = None
    cds_start = min(transcript.start_codon_spliced_offsets)
    if any(anchor >= cds_start for _, anchor in resolved):
        mutant_protein = translate_from_cds(
            mutant_cdna, transcript, [edit for edit, _ in resolved])
    return MutantTranscript(
        reference_transcript=transcript,
        edits=tuple(edit for edit, _ in resolved),
        cdna_sequence=mutant_cdna,
        mutant_protein_sequence=mutant_protein,
        annotator_name="protein_diff",
    )
