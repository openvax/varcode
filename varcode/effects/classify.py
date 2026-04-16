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

"""Classify a MutationEffect from a reference/mutant protein pair.

Shared classifier used by both the splice-outcome builder
(:mod:`varcode.splice_outcomes`, #305) and the forthcoming
:class:`ProteinDiffEffectAnnotator` (#309 / stage 3d). Reduces
the protein pair via :func:`trim_shared_flanking_strings` and
dispatches to the standard Effect classes.

This module deliberately does NOT handle:

* Non-coding / incomplete transcripts — caller gates those before
  reaching the classifier.
* Splice-class effects (SpliceDonor, etc.) — those come from the
  position-based splice classifier, not from a protein diff.
* ``AlternateStartCodon`` — requires inspecting the mutant cDNA's
  first codon, not just the protein strings. The caller (the
  annotator) special-cases that before calling this function.
"""

from ..string_helpers import trim_shared_flanking_strings
from .effect_classes import (
    ComplexSubstitution,
    Deletion,
    FrameShift,
    FrameShiftTruncation,
    Insertion,
    PrematureStop,
    Silent,
    StartLoss,
    StopLoss,
    Substitution,
)


def classify_from_protein_diff(
        variant,
        transcript,
        ref_protein,
        mut_protein,
        length_delta=0):
    """Classify a :class:`MutationEffect` by diffing translated
    proteins.

    Parameters
    ----------
    variant : Variant
    transcript : pyensembl.Transcript
    ref_protein : str
        Reference protein sequence (from ``transcript.protein_sequence``).
    mut_protein : str
        Mutant protein sequence (translated from the mutated cDNA,
        stopping at the first stop codon).
    length_delta : int
        Net nucleotide-level length change of the cDNA edit
        (positive for net insertion, negative for net deletion, zero
        for substitution). Used to detect frameshifts: when
        ``length_delta % 3 != 0`` AND the protein diff shows a
        changed tail, the effect is a :class:`FrameShift`.

    Returns
    -------
    MutationEffect subclass instance
    """
    if ref_protein == mut_protein:
        return Silent(
            variant=variant,
            transcript=transcript,
            aa_pos=0,
            aa_ref=ref_protein[:1] if ref_protein else "")

    # Start loss: mutant protein doesn't begin with M.
    if not mut_protein or (mut_protein[0] != "M" and ref_protein and ref_protein[0] == "M"):
        return StartLoss(variant=variant, transcript=transcript)

    ref_delta, alt_delta, prefix, suffix = trim_shared_flanking_strings(
        ref_protein, mut_protein)

    aa_offset = len(prefix)
    n_ref = len(ref_delta)
    n_alt = len(alt_delta)

    # Frameshift: cDNA length change not divisible by 3.
    if length_delta % 3 != 0:
        if n_alt == 0:
            return FrameShiftTruncation(
                variant=variant,
                transcript=transcript,
                stop_codon_offset=aa_offset)
        return FrameShift(
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_offset,
            shifted_sequence=alt_delta)

    # Premature stop: mutant protein shorter than reference and the
    # change is at the tail (the trimmed alt runs to the end of the
    # mutant protein). Use the single reference residue at the stop-
    # creation point as aa_ref (matching legacy's convention, which
    # shows the codon that became a stop rather than the entire
    # truncated tail).
    if len(mut_protein) < len(ref_protein) and aa_offset + n_alt == len(mut_protein):
        return PrematureStop(
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_offset,
            aa_ref=ref_protein[aa_offset] if aa_offset < len(ref_protein) else ref_delta,
            aa_alt=alt_delta)

    # Stop loss: mutant protein longer than reference and the change
    # extends past the original stop.
    if len(mut_protein) > len(ref_protein) and aa_offset + n_ref >= len(ref_protein):
        return StopLoss(
            variant=variant,
            transcript=transcript,
            aa_ref=ref_delta,
            aa_alt=alt_delta)

    # Silent (after trimming — all changes cancelled out).
    if n_ref == 0 and n_alt == 0:
        return Silent(
            variant=variant,
            transcript=transcript,
            aa_pos=aa_offset,
            aa_ref=prefix + suffix)

    # Deletion.
    if n_alt == 0:
        return Deletion(
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_offset,
            aa_ref=ref_delta)

    # Insertion.
    if n_ref == 0:
        return Insertion(
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_offset,
            aa_alt=alt_delta)

    # Substitution (simple 1:1).
    if n_ref == 1 and n_alt == 1:
        return Substitution(
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_offset,
            aa_ref=ref_delta,
            aa_alt=alt_delta)

    # Complex substitution (multi-residue edit).
    return ComplexSubstitution(
        variant=variant,
        transcript=transcript,
        aa_mutation_start_offset=aa_offset,
        aa_ref=ref_delta,
        aa_alt=alt_delta)
