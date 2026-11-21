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

"""Structural protocols for gene and transcript views accepted by Varcode (#653).

Varcode annotates ordinary PyEnsembl transcripts, ORF coding adapters
(:class:`~varcode.orf.ORFTranscriptView`), offline archive snapshots
(:class:`~varcode.orf_serialization.ReferenceTranscriptSnapshot`), and
caller-supplied RNA/transcript adapters without requiring a single concrete
base class. These protocols capture the structural surface that public effect
and collection APIs rely on.
"""

from __future__ import annotations

from typing import Protocol


class GeneLike(Protocol):
    """Structural gene view shared by PyEnsembl genes and offline snapshots.

    Attributes
    ----------
    id : str
        Stable gene identifier (e.g. ``"ENSG00000141510"``).
    name : str
        Human-readable gene symbol (e.g. ``"TP53"``).
    contig : str
        Chromosome or scaffold name in the reference annotation's vocabulary.
    start : int
        1-based inclusive genomic start coordinate.
    end : int
        1-based inclusive genomic end coordinate.
    strand : str
        Genomic strand (``"+"`` or ``"-"``).
    """

    @property
    def id(self) -> str: ...

    @property
    def name(self) -> str: ...

    @property
    def contig(self) -> str: ...

    @property
    def start(self) -> int: ...

    @property
    def end(self) -> int: ...

    @property
    def strand(self) -> str: ...


class TranscriptLike(Protocol):
    """Structural transcript view shared by PyEnsembl transcripts, ORF adapters,
    offline snapshots, and external RNA transcript adapters.

    Attributes
    ----------
    id : str
        Stable transcript identifier (e.g. ``"ENST00000269305"``).
    name : str
        Human-readable transcript name (e.g. ``"TP53-001"``).
    contig : str
        Chromosome or scaffold name in the reference annotation's vocabulary.
    start : int
        1-based inclusive genomic start coordinate.
    end : int
        1-based inclusive genomic end coordinate.
    strand : str
        Genomic strand (``"+"`` or ``"-"``).
    gene_id : str
        Identifier of the parent gene.
    protein_sequence : str or None
        Reference translated protein sequence, or ``None`` for non-coding or
        incomplete transcripts without a protein sequence.
    complete : bool
        True when the transcript has an annotated start codon, stop codon, and
        valid coding sequence.
    is_protein_coding : bool
        True when the transcript represents a protein-coding locus.
    gene : GeneLike or None
        Parent gene object when available, or ``None``.
    """

    @property
    def id(self) -> str: ...

    @property
    def name(self) -> str: ...

    @property
    def contig(self) -> str: ...

    @property
    def start(self) -> int: ...

    @property
    def end(self) -> int: ...

    @property
    def strand(self) -> str: ...

    @property
    def gene_id(self) -> str: ...

    @property
    def protein_sequence(self) -> str | None: ...

    @property
    def complete(self) -> bool: ...

    @property
    def is_protein_coding(self) -> bool: ...

    @property
    def gene(self) -> GeneLike | None: ...


__all__ = ["GeneLike", "TranscriptLike"]
