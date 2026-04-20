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
variants to a reference transcript (openvax/varcode#271, stage 1).

:class:`TranscriptEdit` and :class:`MutantTranscript` are the types
used by the forthcoming protein-diff :class:`EffectAnnotator`
(``varcode.annotators``) to reshape effect annotation from
"reason about offsets against the reference" to "materialize the
mutant sequence, translate it, compare to the reference protein."

This module only defines the data shapes — the sequence-construction
logic, fast-path dispatch, and annotator wiring are staged across
later PRs (see the tracking issue). For now, this lets adjacent work
(splice_outcomes rewrite, RNA-evidence ingestion, germline-aware
annotation) reference a stable abstraction.
"""

from dataclasses import dataclass, field
from typing import Optional, Tuple

from serializable import DataclassSerializable


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
    """A reference transcript (or assembled set of reference segments)
    with zero or more variant-derived edits applied, optionally
    carrying the mutated cDNA and protein sequences.

    Producers (the protein-diff annotator, RNA-evidence importers,
    the splice-outcomes rewrite, germline-aware annotation,
    structural-variant annotators) construct this once per
    (transcript-or-segments, variant-set, context) and hand it to
    downstream consumers. Each consumer reads the fields it cares
    about — ``edits`` for provenance, ``cdna_sequence`` /
    ``mutant_protein_sequence`` for protein-level analysis.

    **Two shapes:**

    * **Point-variant shape** (``reference_transcript`` is set,
      ``reference_segments`` is ``None``) — the mutant is derived
      from a single reference transcript by applying zero or more
      :class:`TranscriptEdit` objects. This is the shape used by
      the protein-diff annotator for SNVs, MNVs, and simple indels.

    * **Structural-variant shape** (``reference_segments`` is set,
      ``reference_transcript`` may be ``None`` or the primary /
      5'-partner transcript) — the mutant is assembled by
      concatenating :class:`ReferenceSegment` slices in order. A
      gene fusion is two segments from two transcripts; a
      translocation to intergenic is one transcript segment plus a
      genomic-interval segment; an inversion is three forward/
      reverse/forward segments. ``edits`` may still be populated
      for point-variant edits layered on top of the assembled
      segments, but typically an SV carries no extra edits.

    Sequence fields are ``Optional[str]`` because not every producer
    computes them eagerly. Callers that require the protein check or
    compute it themselves; the protein-diff annotator guarantees it
    for point variants.

    **Forward-looking hooks** (not implemented here; documented so
    new integrations know where to plug in):

    * Personalized / full-genome reference — pass a
      :class:`ReferenceSegment` whose ``source`` is a patient-specific
      contig object. varcode's translation logic reads
      ``source.sequence``; it doesn't care whether that's GRCh38 or
      a custom assembly.
    * Long-read resolution — when an SV has an :attr:`alt_assembly`
      on the :class:`StructuralVariant`, the SV annotator can wrap
      that sequence as a single synthetic segment.
    * SV outcomes ambiguity — a translocation producing many candidate
      ORFs that only RNA can resolve should return
      ``List[MutantTranscript]`` or wrap it in a
      :class:`MultiOutcomeEffect` per #299, each with its own
      ``evidence`` dict capturing the disambiguator.
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
    in each segment's own reference system."""

    cdna_sequence: Optional[str] = None
    """The mutated spliced mRNA, when computed. ``None`` if the
    producer hasn't materialized it yet."""

    mutant_protein_sequence: Optional[str] = None
    """The translated mutant protein, stopping at the first stop
    codon. ``None`` if not yet translated, or if the edit set
    doesn't produce a coherent ORF (e.g. start-codon loss). Callers
    that need a guaranteed-present protein should use the
    protein-diff annotator once it lands."""

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


# ---------------------------------------------------------------------
# Construction (#271 stage 2)
# ---------------------------------------------------------------------


_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


def _reverse_complement(seq: str) -> str:
    return "".join(_COMPLEMENT.get(b, b) for b in reversed(seq.upper()))


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

    Callers that get ``None`` should fall back to the fast
    :class:`EffectAnnotator`. The forthcoming protein-diff annotator
    layers effect classification on top of this builder.
    """
    # Lazy imports to keep module-level deps light.
    from pyensembl import Transcript
    from .effects.codon_tables import (
        codon_table_for_transcript,
        translate_sequence,
    )
    from .effects.transcript_helpers import interval_offset_on_transcript

    if not isinstance(transcript, Transcript):
        return None
    if not transcript.is_protein_coding or not transcript.complete:
        return None

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

    # Strand-flip ref/alt for reverse-strand transcripts.
    genome_ref = variant.trimmed_ref
    genome_alt = variant.trimmed_alt
    if transcript.on_backward_strand:
        cdna_ref = _reverse_complement(genome_ref)
        cdna_alt = _reverse_complement(genome_alt)
    else:
        cdna_ref = genome_ref
        cdna_alt = genome_alt

    # Validate reference matches.
    full_sequence = str(transcript.sequence)
    n_ref = len(cdna_ref)
    expected_ref = full_sequence[cdna_offset:cdna_offset + n_ref]
    if cdna_ref != expected_ref:
        return None

    # For pure insertions (n_ref == 0), VCF convention places the alt
    # AFTER the anchor base in genomic coordinates. On the forward
    # strand that's cdna_offset + 1; on the reverse strand, "after"
    # in genomic coords is "before" in cDNA coords (the cDNA runs
    # 3'→5' in genomic space), so the insertion goes AT cdna_offset.
    if n_ref == 0:
        if transcript.on_backward_strand:
            insert_at = cdna_offset
        else:
            insert_at = cdna_offset + 1
        edit = TranscriptEdit(
            cdna_start=insert_at,
            cdna_end=insert_at,
            alt_bases=cdna_alt,
            source_variant=variant,
        )
        mutant_cdna = (
            full_sequence[:insert_at]
            + cdna_alt
            + full_sequence[insert_at:]
        )
    else:
        edit = TranscriptEdit(
            cdna_start=cdna_offset,
            cdna_end=cdna_offset + n_ref,
            alt_bases=cdna_alt,
            source_variant=variant,
        )
        mutant_cdna = (
            full_sequence[:cdna_offset]
            + cdna_alt
            + full_sequence[cdna_offset + n_ref:]
        )

    mutant_protein = None
    cds_start = min(transcript.start_codon_spliced_offsets)
    if cdna_offset >= cds_start:
        # Variant lands in or after the CDS — translate the new
        # coding sequence from the canonical start to the first
        # stop. The codon table is mt-aware via the transcript's
        # contig.
        codon_table = codon_table_for_transcript(transcript)
        coding = mutant_cdna[cds_start:]
        truncated = coding[:(len(coding) // 3) * 3]
        try:
            mutant_protein = translate_sequence(
                truncated, codon_table=codon_table, to_stop=True)
        except ValueError:
            mutant_protein = None
    # If cdna_offset < cds_start the variant is in the 5' UTR; we
    # don't attempt translation here since it might (depending on the
    # exact variant) preserve, lose, or shift the start codon — that
    # interpretation is the fast annotator's responsibility for now.

    return MutantTranscript(
        reference_transcript=transcript,
        edits=(edit,),
        cdna_sequence=mutant_cdna,
        mutant_protein_sequence=mutant_protein,
        annotator_name="protein_diff",
    )
