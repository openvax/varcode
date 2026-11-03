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

"""Transcripts annotated across the origin of a circular contig.

A gene that crosses the origin of a circular genome, such as hepatitis B
polymerase, is written as a GenBank join (``join(2307..3215,1..1623)``). It
reaches PyEnsembl as exons in transcript order (``exon_number``) that step
back once, onto base 1. The step is the origin: contiguous sequence, not a
splice junction.

Varcode asks geometric questions about a transcript in its *frame*:
coordinates in which the transcript runs monotonically along its strand.
On a linear contig the frame is the genome's own coordinates. For a
transcript joined across an origin, the frame reads the stretch that starts
at base 1 one contig length further on, which makes the transcript an
ordinary linear interval again. A genomic interval then has up to two
images in the frame, as an arc on a circle can meet a transcript at both of
its ends. The stretch between the transcript's 3' and 5' ends is not on the
transcript, although its linear span ``transcript.start..transcript.end``
includes it.

Varcode does not know contig lengths. It reads the length from the join
itself: the exon on the far side of base 1 ends at the contig's last base.
"""

import weakref
from dataclasses import dataclass


@dataclass(frozen=True)
class OriginJoin:
    """Where a transcript crosses the origin of its circular contig.

    Attributes
    ----------
    exon_index : int
        Zero-based index, in transcript order, of the exon that reaches the
        origin. The next exon continues from it.
    contig_length : int
        Length of the circular contig: the last base of the exon that meets
        the contig's end.
    low_end : int
        Last base of the transcript's stretch that starts at base 1. The
        frame reads positions on that stretch one ``contig_length`` on.
    """

    exon_index: int
    contig_length: int
    low_end: int


def _origin_join(blocks, backward):
    """The :class:`OriginJoin` of ``(start, end)`` blocks in given order."""
    steps = [
        index for index, ((start, _), (following, _)) in enumerate(
            zip(blocks, blocks[1:]))
        if (following > start if backward else following < start)]
    if len(steps) != 1:
        return None
    index = steps[0]
    before, after = blocks[:index + 1], blocks[index + 1:]
    # The stretch from base 1 is after the origin on the forward strand and
    # before it on the reverse strand, where the transcript reads downward.
    low = before if backward else after
    if min(start for start, _ in low) != 1:
        return None
    contig_length = after[0][1] if backward else before[-1][1]
    return OriginJoin(index, contig_length, max(end for _, end in low))


def in_transcript_order(blocks, backward):
    """Return ``(start, end)`` blocks in transcript order.

    Parameters
    ----------
    blocks : sequence of (int, int)
        Genomic blocks as an annotation lists them, normally by
        ``exon_number``.
    backward : bool
        Whether the transcript is on the reverse strand.

    Returns
    -------
    tuple of (int, int)
        The blocks as given when they cross a circular origin. Otherwise
        the blocks in coordinate order along the strand, which also corrects
        annotations that number reverse-strand exons left to right.
    """
    blocks = tuple(blocks)
    if _origin_join(blocks, backward) is not None:
        return blocks
    return tuple(sorted(blocks, reverse=backward))


class TranscriptFrame:
    """A transcript's exons and coordinates in its own frame.

    Parameters
    ----------
    transcript : pyensembl.Transcript or transcript-like
        Anything with ``exons`` (each with ``start`` and ``end``),
        ``on_backward_strand`` and ``contig``.

    Attributes
    ----------
    join : OriginJoin or None
        Where the transcript crosses a circular origin, or ``None``.
    exons : tuple
        The transcript's exon objects in transcript order.
    exon_intervals : tuple of (int, int)
        The same exons as ``(start, end)`` in the frame.
    start, end : int
        The transcript's span in the frame.
    """

    def __init__(self, transcript):
        exons = tuple(transcript.exons)
        backward = transcript.on_backward_strand
        # An origin join always reaches base 1, so other transcripts need
        # no scan.
        join = None
        if exons and min(exon.start for exon in exons) == 1:
            join = _origin_join(
                [(exon.start, exon.end) for exon in exons], backward)
        if join is None:
            exons = tuple(sorted(
                exons, key=lambda exon: exon.start, reverse=backward))
            intervals = tuple((exon.start, exon.end) for exon in exons)
        else:
            low = (range(join.exon_index + 1) if backward
                   else range(join.exon_index + 1, len(exons)))
            shift = join.contig_length
            intervals = tuple(
                (exon.start + shift, exon.end + shift) if index in low
                else (exon.start, exon.end)
                for index, exon in enumerate(exons))
        self.transcript = transcript
        self.contig = transcript.contig
        self.backward = backward
        self.join = join
        self.exons = exons
        self.exon_intervals = intervals
        if intervals:
            self.start = min(start for start, _ in intervals)
            self.end = max(end for _, end in intervals)
        else:
            self.start, self.end = transcript.start, transcript.end

    @property
    def contig_length(self):
        """The circular contig's length, or ``None`` on a linear contig."""
        return None if self.join is None else self.join.contig_length

    @property
    def exon_number_before_origin(self):
        """One-based number of the exon that reaches the origin, or ``None``."""
        return None if self.join is None else self.join.exon_index + 1

    def images(self, start, end):
        """Images of genomic ``start..end`` in the frame.

        Parameters
        ----------
        start, end : int
            One-based inclusive genomic coordinates. On a circular contig
            ``end`` may run past the contig's length, as a record crossing
            the origin writes it.

        Returns
        -------
        tuple of (int, int)
            On a linear contig, ``((start, end),)``. On a circular one, each
            image that meets the transcript's span, in transcript order; an
            arc that meets both of the transcript's ends has two. An
            interval that misses the transcript keeps its own coordinates.
        """
        if self.join is None:
            return ((start, end),)
        length = self.join.contig_length
        found = [
            (start + lap * length, end + lap * length) for lap in (-1, 0, 1)
            if start + lap * length <= self.end
            and end + lap * length >= self.start]
        return tuple(sorted(found, reverse=self.backward)) or ((start, end),)

    def position(self, position):
        """The frame position of a genomic base, first in transcript order."""
        return self.images(position, position)[0][0]

    def genomic(self, position):
        """The genomic position of a frame position."""
        if self.join is None:
            return position
        return (position - 1) % self.join.contig_length + 1

    def overlaps(self, start, end):
        """Whether genomic ``start..end`` meets the transcript's span."""
        return any(
            low <= self.end and high >= self.start
            for low, high in self.images(start, end))

    def contains(self, start, end):
        """Whether genomic ``start..end`` lies within the transcript's span."""
        return any(
            self.start <= low and high <= self.end
            for low, high in self.images(start, end))

    def genomic_intervals(self, start, end):
        """Genomic ``(start, end)`` pieces of a frame interval."""
        if self.join is None:
            return ((start, end),)
        length = self.join.contig_length
        pieces = []
        position = start
        while position <= end:
            genomic = self.genomic(position)
            stop = min(end, position + length - genomic)
            pieces.append((genomic, genomic + stop - position))
            position = stop + 1
        return tuple(pieces)

    def origin_partners(self, exon_number):
        """One-based numbers of the exons forming one stretch on the circle.

        The exons on either side of an origin are one continuous exon. For
        either of them this returns both; otherwise ``(exon_number,)``.
        """
        before = self.exon_number_before_origin
        if before is not None and exon_number in (before, before + 1):
            return (before, before + 1)
        return (exon_number,)


def feature_overlaps(feature, contig, start, end):
    """Whether genomic ``contig:start..end`` meets a transcript, gene or exon.

    Parameters
    ----------
    feature : transcript, gene or exon
        Anything with ``contig``, ``start`` and ``end``. A transcript (with
        ``exons``) or a gene (with ``transcripts``) that reaches base 1 is
        read on the circle when one of its transcripts crosses the origin.
    contig : str
    start, end : int
        One-based inclusive; ``end`` may run past a circular contig's end.

    Returns
    -------
    bool
    """
    if str(feature.contig) != str(contig):
        return False
    linear = feature.start <= end and start <= feature.end
    # Only a feature reaching base 1 can cross an origin; this avoids
    # querying the transcripts of every other gene.
    if feature.start != 1:
        return linear
    transcripts = getattr(feature, "transcripts", None)
    if transcripts is not None:
        frames = [TranscriptFrame(transcript) for transcript in transcripts]
        if any(frame.join is not None for frame in frames):
            return any(frame.overlaps(start, end) for frame in frames)
        return linear
    if getattr(feature, "exons", None) is not None:
        return TranscriptFrame(feature).overlaps(start, end)
    return linear


def transcript_exons(transcript):
    """Exon objects of ``transcript`` in transcript order (see :func:`in_transcript_order`)."""
    return TranscriptFrame(transcript).exons


def exon_number_before_origin(transcript):
    """One-based number of the exon that reaches a circular origin, or ``None``."""
    if transcript.start != 1:
        return None
    return TranscriptFrame(transcript).exon_number_before_origin


def origin_partners(transcript, exon_number):
    """See :meth:`TranscriptFrame.origin_partners`."""
    if transcript.start != 1:
        return (exon_number,)
    return TranscriptFrame(transcript).origin_partners(exon_number)


# ---------------------------------------------------------------------
# Locus queries
# ---------------------------------------------------------------------


_CIRCULAR_GENES = weakref.WeakKeyDictionary()


def _circular_genes(genome, contig):
    """``{gene_id: (gene, transcripts)}`` for genes on ``contig`` with a
    transcript joined across its origin, cached per genome and contig."""
    try:
        cache = _CIRCULAR_GENES.setdefault(genome, {})
    except TypeError:
        cache = {}
    if contig not in cache:
        genes = {}
        try:
            # Every origin-joined transcript reaches base 1.
            at_origin = genome.transcripts_at_locus(contig, 1, 1)
        except ValueError:
            at_origin = []
        for transcript in at_origin:
            if (transcript.gene_id not in genes
                    and TranscriptFrame(transcript).join is not None):
                gene = genome.gene_by_id(transcript.gene_id)
                genes[transcript.gene_id] = (gene, tuple(gene.transcripts))
        cache[contig] = genes
    return cache[contig]


def _gene_overlaps(transcripts, start, end):
    return any(
        TranscriptFrame(transcript).overlaps(start, end)
        for transcript in transcripts)


def transcripts_at_locus(genome, contig, start, end):
    """Transcripts overlapping ``contig:start..end``, read on the circle.

    ``genome.transcripts_at_locus`` without the origin-joined transcripts
    whose linear span alone reaches the locus.
    """
    transcripts = genome.transcripts_at_locus(contig, start, end)
    if not _circular_genes(genome, contig):
        return transcripts
    return [
        transcript for transcript in transcripts
        if TranscriptFrame(transcript).overlaps(start, end)]


def _genes_off_circle(genome, contig, start, end):
    """IDs and names of genes that only an origin join's linear span reaches."""
    off = {
        gene_id: gene.name
        for gene_id, (gene, transcripts) in _circular_genes(genome, contig).items()
        if not _gene_overlaps(transcripts, start, end)}
    return off


def genes_at_locus(genome, contig, start, end):
    """Genes overlapping ``contig:start..end``, read on the circle."""
    genes = genome.genes_at_locus(contig, start, end)
    off = _genes_off_circle(genome, contig, start, end)
    return [gene for gene in genes if gene.id not in off]


def gene_ids_at_locus(genome, contig, start, end):
    """Gene IDs overlapping ``contig:start..end``, read on the circle."""
    gene_ids = genome.gene_ids_at_locus(contig, start, end)
    off = _genes_off_circle(genome, contig, start, end)
    return [gene_id for gene_id in gene_ids if gene_id not in off]


def gene_names_at_locus(genome, contig, start, end):
    """Gene names overlapping ``contig:start..end``, read on the circle."""
    names = genome.gene_names_at_locus(contig, start, end)
    off = _genes_off_circle(genome, contig, start, end)
    if not off:
        return names
    kept = {gene.name for gene in genes_at_locus(genome, contig, start, end)}
    return [name for name in names if name in kept]
