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
used by the forthcoming sequence-diff :class:`EffectAnnotator`
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


@dataclass(frozen=True)
class TranscriptEdit:
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
class MutantTranscript:
    """A reference transcript with zero or more variant-derived edits
    applied, optionally carrying the mutated cDNA and protein
    sequences.

    Producers (the sequence-diff annotator, RNA-evidence importers,
    the splice-outcomes rewrite, germline-aware annotation) construct
    this once per (transcript, variant-set, context) and hand it to
    downstream consumers. Each consumer reads the fields it cares
    about — ``edits`` for provenance, ``cdna_sequence`` /
    ``mutant_protein_sequence`` for protein-level analysis.

    Sequence fields are ``Optional[str]`` because not every producer
    computes them eagerly (some stages just need the edit list; full
    translation is lazy). Callers that require the protein must
    check or compute it themselves for now — the sequence-diff
    annotator in stage 2 will guarantee it's populated.
    """

    reference_transcript: object
    """The :class:`pyensembl.Transcript` this mutant is derived from.
    Not typed tightly here so :mod:`pyensembl` isn't a hard import
    dependency for anyone who just wants the dataclass."""

    edits: Tuple[TranscriptEdit, ...] = field(default_factory=tuple)
    """Edits applied to produce this mutant, sorted by
    :attr:`TranscriptEdit.cdna_start`. Empty tuple means the mutant
    is identical to the reference."""

    cdna_sequence: Optional[str] = None
    """The mutated spliced mRNA, when computed. ``None`` if the
    producer hasn't materialized it yet."""

    mutant_protein_sequence: Optional[str] = None
    """The translated mutant protein, stopping at the first stop
    codon. ``None`` if not yet translated, or if the edit set
    doesn't produce a coherent ORF (e.g. start-codon loss). Callers
    that need a guaranteed-present protein should use the
    sequence-diff annotator once it lands."""

    annotator_name: str = "unknown"
    """Name of the :class:`EffectAnnotator` (or other producer) that
    created this ``MutantTranscript``. Used as provenance in
    serialization and for A/B comparisons."""

    evidence: Optional[dict] = None
    """Optional producer-specific evidence (RNA read counts, Isovar
    fragment ids, SpliceAI scores). Shape is annotator-specific and
    not part of the stable contract; consumers that care about a
    particular evidence shape should type-check it at the call site."""

    def __post_init__(self):
        if self.edits:
            starts = [e.cdna_start for e in self.edits]
            if any(starts[i] > starts[i + 1] for i in range(len(starts) - 1)):
                raise ValueError(
                    "MutantTranscript.edits must be sorted by cdna_start")

    @property
    def is_identical_to_reference(self) -> bool:
        """True if no edits were applied. Does NOT check
        ``cdna_sequence`` / ``mutant_protein_sequence`` — a producer
        can legitimately carry an identical sequence with zero edits."""
        return len(self.edits) == 0

    @property
    def total_length_delta(self) -> int:
        """Sum of :attr:`TranscriptEdit.length_delta` across all
        edits — how much longer or shorter the mutant cDNA is than
        the reference."""
        return sum(e.length_delta for e in self.edits)
