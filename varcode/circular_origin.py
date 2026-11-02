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
polymerase, is written as a GenBank join (``join(2307..3215,1..1623)``) and
reaches PyEnsembl as exons in transcript order that step back to base 1.
The step is the origin, not an intron, and the stretch between the
transcript's 3' and 5' ends lies outside it, although the linear span
``transcript.start..transcript.end`` includes it.
"""


def origin_join_index(transcript):
    """Return where a transcript crosses its circular contig's origin.

    ``transcript.exons`` is in transcript order (PyEnsembl orders it by
    ``exon_number``). On a linear contig that order only moves along the
    strand.

    Returns
    -------
    int or None
        Index ``i`` such that ``exons[i]`` ends at the origin and
        ``exons[i + 1]`` continues from base 1 (forward strand), or
        ``exons[i]`` ends at base 1 and ``exons[i + 1]`` continues from the
        contig's end (reverse strand). ``None`` when the exons never step
        back.

    Raises
    ------
    ValueError
        When the exons step back more than once, or anywhere but at base 1.
        Without the contig's length, varcode cannot place that sequence.
    """
    exons = transcript.exons
    backward = transcript.on_backward_strand
    steps = [
        index for index, (exon, following) in enumerate(zip(exons, exons[1:]))
        if (following.start > exon.start if backward
            else following.start < exon.start)]
    if not steps:
        return None
    if len(steps) == 1:
        index = steps[0]
        at_origin = exons[index] if backward else exons[index + 1]
        if at_origin.start == 1:
            return index
    raise ValueError(
        "Exons of %s step back along the %s strand after exon(s) %s; "
        "varcode follows exons out of coordinate order only across the "
        "origin of a circular contig, where they meet at base 1" % (
            transcript.id, transcript.strand,
            ", ".join(str(step + 1) for step in steps)))


def outside_origin_joined_transcript(transcript, start, end):
    """Whether ``start..end`` lies between an origin-joined transcript's ends.

    That stretch is within ``transcript.start..transcript.end`` but not on
    the transcript, which runs from its 5' end through the origin to its
    3' end. Always ``False`` for a transcript without an origin join, and
    for one whose exons step back other than at the origin.
    """
    # An origin join reaches base 1, so other transcripts need no exons.
    if transcript.start != 1:
        return False
    try:
        if origin_join_index(transcript) is None:
            return False
    except ValueError:
        return False
    exons = transcript.exons
    if transcript.on_backward_strand:
        return exons[0].end < start and end < exons[-1].start
    return exons[-1].end < start and end < exons[0].start


def origin_join_flanks(transcript, flank):
    """Genomic blocks of an origin-joined transcript in transcript order.

    Returns two ``(start, end)`` blocks, the one before the origin and the
    one after it, each extended by up to ``flank`` bases away from the
    origin. The flanks share the stretch between the transcript's ends
    without overlapping. Returns ``None`` without an origin join.
    """
    index = origin_join_index(transcript)
    if index is None:
        return None
    exons = transcript.exons
    if transcript.on_backward_strand:
        # Sense runs down from exons[0].end to base 1, then down from the
        # contig's end (exons[index + 1].end) to exons[-1].start.
        first_end = min(exons[0].end + flank, exons[-1].start - 1)
        second_start = max(exons[-1].start - flank, first_end + 1)
        return (1, first_end), (second_start, exons[index + 1].end)
    # Sense runs up from exons[0].start to the contig's end
    # (exons[index].end), then up from base 1 to exons[-1].end.
    second_end = min(exons[-1].end + flank, exons[0].start - 1)
    first_start = max(exons[0].start - flank, second_end + 1)
    return (first_start, exons[index].end), (1, second_end)
