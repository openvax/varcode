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

"""Shared fast path for trivial single-codon SNVs (openvax/varcode#271,
stage 3c).

Both the legacy and forthcoming protein-diff annotators dispatch
coding-SNV variants through the same short-circuit so they agree on
the common case. The full :func:`predict_in_frame_coding_effect` logic
only runs for variants that actually need it — indels, MNVs,
start-/stop-adjacent substitutions, and anything that might fall out
of the straightforward Silent / Substitution / PrematureStop
classification.

Extracting this helper doesn't change behaviour on its own: legacy
still produces the same Effect classes and same ``short_description``
byte-for-byte. The value comes in stage 3d when the protein-diff
annotator shares this code path, removing it as a source of A/B
divergence between the two annotators.
"""

from .codon_tables import codon_table_for_transcript
from .effect_classes import PrematureStop, Silent, Substitution


def try_fast_path_snv(
        variant,
        transcript,
        trimmed_cdna_ref,
        trimmed_cdna_alt,
        sequence_from_start_codon,
        cds_offset):
    """Classify a single-codon SNV without materializing a mutant
    transcript or running the full in-frame pipeline.

    Returns an :class:`Effect` when the variant qualifies:

    * one base of cDNA reference replaced by one base of alt
      (rules out indels and MNVs),
    * the affected codon sits strictly inside the CDS — not the
      start codon (which can produce StartLoss or
      AlternateStartCodon) and not the stop codon (which can
      produce StopLoss with 3' UTR readthrough).

    Returns ``None`` for anything else; callers fall through to
    :func:`predict_in_frame_coding_effect`.
    """
    if len(trimmed_cdna_ref) != 1 or len(trimmed_cdna_alt) != 1:
        return None

    ref_codon_start_offset = cds_offset // 3
    ref_codon_end_offset = ref_codon_start_offset + 1

    # Start codon — defer to the slow path so StartLoss and
    # AlternateStartCodon classification run.
    if ref_codon_start_offset == 0:
        return None

    # Stop codon or past it — slow path can handle StopLoss with
    # 3'UTR readthrough; we can't.
    if ref_codon_end_offset > len(transcript.protein_sequence):
        return None

    codon_start_in_cds = ref_codon_start_offset * 3
    ref_codon = str(sequence_from_start_codon[
        codon_start_in_cds:codon_start_in_cds + 3])
    if len(ref_codon) != 3:
        # Ran off the end of the CDS — shouldn't normally happen once
        # we've passed the stop-codon guard, but be defensive about
        # it and defer to the slow path.
        return None
    offset_in_codon = cds_offset % 3
    mutant_codon = (
        ref_codon[:offset_in_codon]
        + trimmed_cdna_alt
        + ref_codon[offset_in_codon + 1:])

    codon_table = codon_table_for_transcript(transcript)
    aa_ref = transcript.protein_sequence[ref_codon_start_offset]

    if mutant_codon in codon_table.stop_codons:
        return PrematureStop(
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=ref_codon_start_offset,
            aa_ref=aa_ref,
            aa_alt="",
        )

    aa_alt = codon_table.forward_table[mutant_codon]
    if aa_ref == aa_alt:
        return Silent(
            variant=variant,
            transcript=transcript,
            aa_pos=ref_codon_start_offset,
            aa_ref=aa_ref,
        )
    return Substitution(
        variant=variant,
        transcript=transcript,
        aa_mutation_start_offset=ref_codon_start_offset,
        aa_ref=aa_ref,
        aa_alt=aa_alt,
    )
