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

Shared classifier used by the splice-outcome builder
(:mod:`varcode.splice_outcomes`), germline-aware prediction
(:mod:`varcode.germline`), and :class:`ProteinDiffEffectAnnotator`.
Reduces the protein pair to its differing region (shared prefix only
for frameshifts, shared prefix and suffix otherwise) and dispatches to
the standard Effect classes.

This module deliberately does NOT handle:

* Non-coding / incomplete transcripts — caller gates those before
  reaching the classifier.
* Splice-class effects (SpliceDonor, etc.) — those come from the
  position-based splice classifier, not from a protein diff.
* ``AlternateStartCodon`` — requires inspecting the mutant cDNA's
  first codon, not just the protein strings. The caller (the
  annotator) special-cases that before calling this function.
"""

from ..string_helpers import (
    trim_shared_flanking_strings,
    trim_shared_prefix,
)
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


def _attach_transcript_models(effect, mutant_transcript, baseline_transcript=None,
                              ref_protein=None):
    """Keep the products used for classification on the ordinary effect API."""
    effect.baseline_transcript = baseline_transcript
    if mutant_transcript is not None and effect.mutant_transcript is None:
        effect.mutant_transcript = mutant_transcript
    if ref_protein is not None:
        effect._original_protein_sequence = ref_protein
    return effect


def classify_from_protein_diff(
        variant,
        transcript,
        ref_protein,
        mut_protein,
        length_delta=0,
        mutant_transcript=None,
        frame_change_aa_offset=None,
        baseline_transcript=None,
        start_complete=True,
        stop_complete=True,
        coding_start=None,
        stop_observed=None):
    """Classify a :class:`MutationEffect` by diffing translated
    proteins.

    Parameters
    ----------
    variant : Variant
    transcript : pyensembl.Transcript
    ref_protein : str
        Baseline protein sequence: the reference or the patient's inherited
        product. Returned sequence accessors retain this actual comparison.
    mut_protein : str
        Mutant protein sequence (translated from the mutated cDNA,
        stopping at the first stop codon).
    length_delta : int
        Net nucleotide-level length change of the cDNA edit
        (positive for net insertion, negative for net deletion, zero
        for substitution). Used to detect frameshifts: when
        ``length_delta % 3 != 0`` AND the protein diff shows a
        changed tail, the effect is a :class:`FrameShift`.
    mutant_transcript : MutantTranscript or None
        Optional edit-level context for the already-materialized
        mutant transcript, retained on every returned effect. When present,
        the classifier can preserve
        insertion semantics for premature-stop calls that would
        otherwise look like a substitution-plus-truncation in a
        whole-protein diff.

    frame_change_aa_offset : int or None
        First affected coding codon when a realized path supplies nucleotide
        coordinates. Distinguishes a coding frameshift whose amino acids
        coincidentally match the reference from a stop-codon extension.
    baseline_transcript : MutantTranscript or None
        Patient baseline, when the comparison includes cis germline edits.
        Keeps start-loss reference codons and offsets in that baseline.

    start_complete, stop_complete : bool
        Whether the reference establishes these boundaries. A coding fragment
        cannot imply start loss from a changed leading methionine or stop loss
        from a longer translated window.
    coding_start : int or None
        Explicit codon origin for a partial coding product; may precede zero
        when the source protein retains a leading unknown partial residue.
    stop_observed : bool or None
        Whether translation found a stop within the supplied coding interval.
        None preserves the established complete-transcript inference.

    Returns
    -------
    MutationEffect subclass instance
    """
    def finish(effect):
        return _attach_transcript_models(
            effect, mutant_transcript, baseline_transcript, ref_protein)

    if ref_protein == mut_protein:
        return finish(Silent(
            variant=variant,
            transcript=transcript,
            aa_pos=0,
            aa_ref=ref_protein[:1] if ref_protein else ""))

    # Start loss: mutant protein doesn't begin with M.
    if start_complete and (not mut_protein or (mut_protein[0] != "M" and ref_protein and ref_protein[0] == "M")):
        reference = {}
        if baseline_transcript is not None:
            from .start_loss import start_loss_details
            baseline = start_loss_details(None, transcript, baseline_transcript)
            reference = dict(ref_codon=baseline["alt_codon"],
                             ref_codon_offset=baseline["alt_codon_offset"])
        return finish(StartLoss(variant=variant, transcript=transcript,
                         mutant_transcript=mutant_transcript, **reference))

    ref_delta, alt_delta, prefix, suffix = trim_shared_flanking_strings(
        ref_protein, mut_protein)

    aa_offset = len(prefix)
    n_ref = len(ref_delta)
    n_alt = len(alt_delta)

    # Frameshift: cDNA length change not divisible by 3.
    if length_delta % 3 != 0:
        # A reading-frame change produces a NOVEL C-terminus that runs to a new stop
        # codon, so a residue it coincidentally shares with the reference protein's own
        # C-terminus must be kept. Derive the novel tail by trimming the shared PREFIX
        # only -- not `alt_delta`, which also had the shared *suffix* trimmed off and so
        # drops real novel residues (e.g. ATM p.F61fs loses its terminal V because the
        # reference ATM protein also ends in V). See openvax/varcode#397.
        _, fs_shifted, fs_prefix = trim_shared_prefix(ref_protein, mut_protein)
        fs_offset = len(fs_prefix)
        # Guard: if the mutation offset is past the end of the reference protein, the
        # edit is in the stop-codon / 3'UTR region. That's a StopLoss, not a FrameShift.
        if stop_complete and (fs_offset if frame_change_aa_offset is None else frame_change_aa_offset) >= len(ref_protein):
            return finish(StopLoss(
                variant=variant,
                transcript=transcript,
                aa_ref="",
                aa_alt=fs_shifted, original_protein_sequence=ref_protein))
        if frame_change_aa_offset is not None:
            fs_offset = min(fs_offset, frame_change_aa_offset)
            fs_shifted = mut_protein[fs_offset:]
        # No novel residues after the shared prefix -> an immediate stop (true truncation).
        if len(fs_shifted) == 0:
            return finish(FrameShiftTruncation(
                variant=variant,
                transcript=transcript,
                stop_codon_offset=fs_offset, original_protein_sequence=ref_protein))
        return finish(FrameShift(
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=fs_offset,
            shifted_sequence=fs_shifted, original_protein_sequence=ref_protein))

    # A mutant protein that is shorter than the reference and truncated
    # at the tail is only a PrematureStop if the mutant CDS actually
    # contains a stop codon. When an in-frame deletion removes the stop
    # codon of a transcript with no 3' UTR to read into, translation
    # simply runs off the end of the available sequence — there is no
    # stop codon at all, so the shortened protein is a C-terminal
    # Deletion (a stop-loss with no predictable readthrough), not a
    # PrematureStop. Mirror the in-frame predictor, which reports a
    # Deletion here. See #394. (When mutant_transcript is unavailable —
    # e.g. the splice-outcome builder — we can't tell, so fall back to
    # the historical PrematureStop classification.)
    mutant_translation_ran_off_end = stop_observed is False
    if (stop_observed is None and mutant_transcript is not None
            and mutant_transcript.cdna_sequence is not None):
        cds_start = min(transcript.start_codon_spliced_offsets) if coding_start is None else coding_start
        n_mutant_codons = (
            len(mutant_transcript.cdna_sequence) - cds_start) // 3
        # translate() stops *before* the stop codon, so a protein that
        # consumed every available codon never encountered one.
        mutant_translation_ran_off_end = len(mut_protein) >= n_mutant_codons

    # Premature stop: mutant protein shorter than reference and the
    # change is at the tail (the trimmed alt runs to the end of the
    # mutant protein). Use the single reference residue at the stop-
    # creation point as aa_ref (matching fast's convention, which
    # shows the codon that became a stop rather than the entire
    # truncated tail).
    if (len(mut_protein) < len(ref_protein)
            and aa_offset + n_alt == len(mut_protein)
            and not mutant_translation_ran_off_end):
        aa_ref = ref_protein[aa_offset] if aa_offset < len(ref_protein) else ref_delta

        # Whole-protein trimming can't distinguish a codon-aligned
        # insertion that introduces an early stop from a substitution
        # of the next reference residue followed by truncation. Use
        # edit provenance when available so pure insertions keep
        # aa_ref="" (matching the historical #174/#175 behaviour
        # preserved by the fast annotator).
        if mutant_transcript is not None and len(mutant_transcript.edits) == 1:
            edit = mutant_transcript.edits[0]
            cds_start = min(transcript.start_codon_spliced_offsets) if coding_start is None else coding_start
            if edit.is_insertion and edit.cdna_start >= cds_start:
                if (edit.cdna_start - cds_start) % 3 == 0:
                    aa_ref = ""
        return finish(PrematureStop(
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_offset,
            aa_ref=aa_ref,
            aa_alt=alt_delta, original_protein_sequence=ref_protein))

    # Stop loss: mutant protein longer than reference, the change
    # extends past the original stop. Two routes into this branch:
    #   1. Reference residues were replaced (n_ref > 0) — the
    #      mutation clearly disrupted residues at/adjacent to the
    #      stop codon.
    #   2. n_ref == 0 after aa-level shared-flank trimming, but the
    #      underlying cDNA edit is NOT a pure insertion — a
    #      substitution or deletion of the stop codon itself
    #      produces zero aa-level ref because the shared prefix
    #      absorbs the entire reference protein. Fall back to the
    #      cDNA edit kind via ``mutant_transcript`` to distinguish
    #      "stop codon mutated → StopLoss" from "in-frame insertion
    #      before stop → Insertion". Closes #319.
    is_cdna_pure_insertion = (
        mutant_transcript is not None
        and len(mutant_transcript.edits) == 1
        and mutant_transcript.edits[0].is_insertion)
    if (stop_complete and len(mut_protein) > len(ref_protein)
            and aa_offset + n_ref >= len(ref_protein)
            and (n_ref > 0 or not is_cdna_pure_insertion)):
        return finish(StopLoss(
            variant=variant,
            transcript=transcript,
            aa_ref=ref_delta,
            aa_alt=alt_delta, original_protein_sequence=ref_protein))

    # Silent (after trimming — all changes cancelled out).
    if n_ref == 0 and n_alt == 0:
        return finish(Silent(
            variant=variant,
            transcript=transcript,
            aa_pos=aa_offset,
            aa_ref=prefix + suffix))

    # Deletion.
    if n_alt == 0:
        return finish(Deletion(
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_offset,
            aa_ref=ref_delta))

    # Insertion.
    if n_ref == 0:
        return finish(Insertion(
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_offset,
            aa_alt=alt_delta))

    # Substitution (simple 1:1).
    if n_ref == 1 and n_alt == 1:
        return finish(Substitution(
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_offset,
            aa_ref=ref_delta,
            aa_alt=alt_delta))

    # Complex substitution (multi-residue edit).
    return finish(ComplexSubstitution(
        variant=variant,
        transcript=transcript,
        aa_mutation_start_offset=aa_offset,
        aa_ref=ref_delta,
        aa_alt=alt_delta))
