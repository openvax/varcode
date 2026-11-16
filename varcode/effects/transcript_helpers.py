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


def interval_offset_on_transcript(start, end, transcript):
    """
    Given an interval [start:end] and a particular transcript,
    return the start offset of the interval relative to the
    chromosomal positions of the transcript.
    """
    # ensure that start_pos:end_pos overlap with transcript positions
    if start > end:
        raise ValueError(
            "start_pos %d shouldn't be greater than end_pos %d" % (
                start, end))
    if start > transcript.end:
        raise ValueError(
            "Range %d:%d starts after transcript %s (%d:%d)" % (
                start,
                end,
                transcript,
                transcript.start,
                transcript.end))
    if end < transcript.start:
        raise ValueError(
            "Range %d:%d ends before transcript %s (%d:%d)" % (
                start,
                end,
                transcript,
                transcript.start,
                transcript.end))
    # trim the start position to the beginning of the transcript
    if start < transcript.start:
        start = transcript.start
    # trim the end position to the end of the transcript
    if end > transcript.end:
        end = transcript.end
    # return earliest offset into the spliced transcript
    return min(
        transcript.spliced_offset(start),
        transcript.spliced_offset(end))


def has_complete_coding_sequence(transcript):
    """
    Does the transcript have a start codon and a stop codon, so that variants
    in its coding sequence have a defined protein effect?

    This is PyEnsembl's ``complete`` plus a stop codon that polyadenylation
    completes. Six human mitochondrial transcripts (MT-ND1, MT-ND2, MT-CO3,
    MT-ND3, MT-ND4 and MT-CYB) end their reading frame in ``T`` or ``TA``, and
    the poly(A) tail supplies the rest of a ``TAA`` stop codon (NCBI
    NC_012920.1). Ensembl annotates no stop codon for them, so ``complete`` is
    False. See ``Transcript.stop_codon_completed_by_polyadenylation`` in
    PyEnsembl 2.34.1 and later.
    """
    return bool(
        transcript.complete
        or getattr(transcript, "stop_codon_completed_by_polyadenylation", False))


def polyadenylation_completes_stop(transcript):
    """Is the transcript's stop codon completed by the poly(A) tail?"""
    return bool(getattr(transcript, "stop_codon_completed_by_polyadenylation", False))


def poly_a_completion(transcript, coding_sequence):
    """
    The ``A`` bases that polyadenylation adds to finish the last, partial
    codon of ``coding_sequence`` when the transcript's stop codon is completed
    by the poly(A) tail, else ``""``.

    ``coding_sequence`` runs from the start codon to the end of the (possibly
    mutated) transcript. In MT-ND1 it ends in ``TA``, and one ``A`` makes the
    stop codon ``TAA``; in MT-ND4 it ends in ``T`` and two ``A`` do. After a
    mutation the same rule applies to whatever the last partial codon is, so
    only the bases that finish that codon are added and no longer tail is
    invented.
    """
    if not polyadenylation_completes_stop(transcript):
        return ""
    return "A" * (-len(coding_sequence) % 3)


def three_prime_utr_offset(transcript):
    """
    Spliced offset of the first base after the coding sequence.

    For a transcript whose stop codon is completed by polyadenylation the
    coding sequence runs to the end of the transcript, so there is no 3' UTR.
    """
    if polyadenylation_completes_stop(transcript):
        return len(transcript.sequence)
    return max(transcript.stop_codon_spliced_offsets) + 1
