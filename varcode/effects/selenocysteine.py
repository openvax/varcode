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

"""Where annotated selenocysteine UGA codons survive in a mutant cDNA.

Ensembl marks selenocysteine (Sec) as U in the reference protein, on an
in-frame UGA. UGA encodes Sec only with a SECIS element in the same mRNA's
3' UTR, and SECIS positions aren't annotated. Translators read these codons
as Sec unless the mutant keeps no selenoprotein 3' UTR at all. The SV change
flags also treat a partly kept 3' UTR as uncertain; one stored protein
follows the reference reading.
"""


def reference_selenocysteine(transcript):
    """cDNA offsets of the first base of each annotated Sec codon."""
    protein = getattr(transcript, "protein_sequence", None)
    if not protein or "U" not in protein or not getattr(transcript, "complete", False):
        return ()
    start = min(transcript.start_codon_spliced_offsets)
    return tuple(start + 3 * i for i, aa in enumerate(protein) if aa == "U")


def _three_prime_utr_start(transcript):
    return max(transcript.stop_codon_spliced_offsets) + 1


def edit_shift(offset, edits):
    """Net length change from edits wholly upstream of reference ``offset``."""
    return sum(len(e.alt_bases) - (e.cdna_end - e.cdna_start)
               for e in edits if e.cdna_end <= offset)


def edited_selenocysteine(transcript, edits):
    """Mutant cDNA offsets of annotated Sec codons that no edit touches.

    For a reference cDNA with point-level edits. Empty when the edits delete
    the transcript's whole 3' UTR, since no SECIS can remain.
    """
    offsets = reference_selenocysteine(transcript)
    if not offsets:
        return set()
    utr_start, utr_end = _three_prime_utr_start(transcript), len(transcript.sequence)
    deleted = sum(max(0, min(e.cdna_end, utr_end) - max(e.cdna_start, utr_start))
                  for e in edits)
    if utr_end > utr_start and deleted == utr_end - utr_start:
        return set()
    return {ref + edit_shift(ref, edits) for ref in offsets
            if not any((e.cdna_start < ref + 3 and ref < e.cdna_end)
                       or ref < e.cdna_start == e.cdna_end < ref + 3 for e in edits)}


def segment_selenocysteine(model):
    """Observed offsets of TGA codons mapped by segments onto annotated Sec codons.

    Return ``(decoded, uncertain, terminating)``: a codon is decoded where the
    model keeps its transcript contiguously from that codon through the 3'
    end. With no selenoprotein 3' UTR sequence at all there is no SECIS, so
    UGA terminates. Otherwise decoding is unknown. None when segments don't
    render the cDNA or point edits are overlaid.
    """
    sequence = model.cdna_sequence
    segments = model.reference_segments or ()
    if (sequence is None or model.edits
            or sum(s.length for s in segments) != len(sequence)):
        return None
    runs, offset = [], 0
    for segment in segments:
        if segment.strand == "+":
            last = runs[-1] if runs else None
            if (last and last[0] == segment.source and last[2] == segment.start
                    and last[3] + last[2] - last[1] == offset):
                last[2] = segment.end  # A split that continues the same reference.
            else:
                runs.append([segment.source, segment.start, segment.end, offset])
        offset += segment.length
    decoded, uncertain = set(), set()
    secis_possible = False
    for source, ref_start, ref_end, observed in runs:
        selenocysteine = reference_selenocysteine(source)
        if selenocysteine and ref_end > _three_prime_utr_start(source):
            secis_possible = True
        for ref in selenocysteine:
            pos = observed + ref - ref_start
            if (ref_start <= ref and ref + 3 <= ref_end
                    and sequence[pos:pos + 3].upper() == "TGA"):
                (decoded if ref_end == len(source.sequence) else uncertain).add(pos)
    if not secis_possible:
        return decoded, set(), uncertain
    return decoded, uncertain, set()


def layout_selenocysteine(transcript, layouts):
    """cDNA offsets of annotated Sec codons in concatenated genomic layouts.

    Each layout exposes ``origins()`` as ``(contig, position, kind)`` per base,
    in transcript orientation. All three bases must map consecutively to the
    annotated codon, including across layout boundaries. Empty when no base
    of an annotated 3' UTR remains.
    """
    offsets = set(reference_selenocysteine(transcript))
    if not offsets:
        return set()
    utr_start = _three_prime_utr_start(transcript)
    found, index = set(), 0
    previous = (None, None)
    utr_present = utr_start >= len(transcript.sequence)
    for layout in layouts:
        for contig, position, kind in layout.origins():
            reference = None
            if contig == transcript.contig and kind != "inverted":
                try:
                    reference = transcript.spliced_offset(position)
                except ValueError:
                    reference = None
            if reference is not None:
                if (reference - 2 in offsets
                        and previous == (reference - 2, reference - 1)):
                    found.add(index - 2)
                if reference >= utr_start:
                    utr_present = True
            previous = (previous[1], reference)
            index += 1
    return found if utr_present else set()
