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

"""ProteinDiffEffectAnnotator — classify effects from protein diff
instead of offset arithmetic (openvax/varcode#271 stage 3d, #309).

See :doc:`/effect_annotation` for the user-facing guide covering
how fast, protein-diff, splice outcomes, and the SV roadmap
fit together. This module docstring captures implementation
detail relevant to reviewers of the classifier itself.

Algorithm
---------

1. **Gate**: non-coding / incomplete transcripts short-circuit to
   the standard wrappers.

2. **Splice check**: run fast first. If fast classifies the
   variant as a splice effect (SpliceDonor, SpliceAcceptor,
   IntronicSpliceSite), return that directly — protein-diff
   doesn't own splice classification. For ExonicSpliceSite, run
   dual-dispatch: fast provides the splice class, protein-diff
   provides the ``alternate_effect`` via protein diff.

3. **Slow path**: build a :class:`MutantTranscript` via
   :func:`apply_variant_to_transcript`, compare the translated
   mutant protein to the reference, and classify via
   :func:`classify_from_protein_diff` (shared with the
   splice-outcome builder, #305).

4. **AlternateStartCodon**: non-diff special case (resolved in
   #304). When proteins match but the first codon changed to a
   non-ATG alternate start, emit :class:`AlternateStartCodon`
   instead of :class:`Silent`.
"""

from ..effects.classify import classify_from_protein_diff
from ..effects.codon_tables import codon_table_for_transcript
from ..effects.effect_classes import (
    AlternateStartCodon,
    ExonicSpliceSite,
    FivePrimeUTR,
    IncompleteTranscript,
    Intronic,
    IntronicSpliceSite,
    NoncodingTranscript,
    SpliceAcceptor,
    SpliceDonor,
    ThreePrimeUTR,
)
from ..mutant_transcript import apply_variant_to_transcript
from ..version import __version__ as _varcode_version
from .fast import FastEffectAnnotator


class ProteinDiffEffectAnnotator:
    """Classify effects by diffing translated mutant protein against
    the reference protein.

    Produces byte-for-byte identical output to
    :class:`FastEffectAnnotator` on the common case (trivial
    SNVs and simple indels) because both flow through the same
    :func:`classify_from_protein_diff` classifier. Diverges where
    protein-diff's approach is provably more accurate (boundary
    codons, frameshift realignment). Any divergence must appear in
    the parity harness ``EXPECTED_DIFFS`` with an issue link.
    """

    name = "protein_diff"
    version = _varcode_version
    supports = frozenset({"snv", "indel", "mnv"})

    def __repr__(self):
        return "ProteinDiffEffectAnnotator(name=%r, version=%r)" % (
            self.name, self.version)

    def annotate_on_transcript(self, variant, transcript):
        """Classify the effect of ``variant`` on ``transcript``.

        Runs fast first to detect splice-adjacent variants (which
        stay fast-classified); for everything else, builds a
        :class:`MutantTranscript` and diffs the translated protein.
        """
        from pyensembl import Transcript
        if not isinstance(transcript, Transcript):
            raise TypeError(
                "Expected %s : %s to have type Transcript" % (
                    transcript, type(transcript)))

        if not transcript.is_protein_coding:
            return NoncodingTranscript(variant, transcript)

        if not transcript.complete:
            return IncompleteTranscript(variant, transcript)

        # Run fast to get the splice classification.
        fast_effect = FastEffectAnnotator().annotate_on_transcript(
            variant, transcript)

        # Pure-intronic splice effects: fast only — no protein-
        # level diff to compute.
        if isinstance(fast_effect, (SpliceDonor, SpliceAcceptor)):
            return fast_effect
        if (isinstance(fast_effect, IntronicSpliceSite)
                and not isinstance(fast_effect, (SpliceDonor, SpliceAcceptor))):
            return fast_effect

        # Non-splice location-based classes (UTR, deep intronic): fast's
        # position-based classification is authoritative. The protein may
        # be unchanged but "outside the CDS" carries location semantics a
        # whole-protein diff doesn't. Closes #318.
        if isinstance(fast_effect, (ThreePrimeUTR, FivePrimeUTR, Intronic)):
            return fast_effect

        # ExonicSpliceSite: dual-dispatch. Fast provides the
        # splice class; protein-diff provides the alternate_effect
        # via protein diff.
        if isinstance(fast_effect, ExonicSpliceSite):
            mt = apply_variant_to_transcript(variant, transcript)
            if mt is not None and mt.mutant_protein_sequence is not None:
                alt = classify_from_protein_diff(
                    variant=variant,
                    transcript=transcript,
                    ref_protein=str(transcript.protein_sequence),
                    mut_protein=mt.mutant_protein_sequence,
                    length_delta=mt.total_length_delta,
                    mutant_transcript=mt)
                return ExonicSpliceSite(
                    variant=variant,
                    transcript=transcript,
                    exon=fast_effect.exon,
                    alternate_effect=alt)
            return fast_effect

        # Non-splice: protein-diff slow path.
        mt = apply_variant_to_transcript(variant, transcript)
        if mt is None or mt.mutant_protein_sequence is None:
            # UTR, ref-mismatch, splice-junction-spanning, etc.
            return fast_effect

        ref_protein = str(transcript.protein_sequence)
        mut_protein = mt.mutant_protein_sequence

        # Alternate start codon rewrite: if the first codon changed to
        # another recognised start codon in the transcript's codon
        # table (e.g. ATG→CTG/GTG/TTG, or MT ATG→GTG under table 2),
        # the initiator tRNA still loads Met regardless of what the
        # codon would decode to internally. Rewrite the mutant
        # protein's first residue to 'M' so the shared diff classifier
        # sees the biologically correct protein. Closes #320.
        cds_start = min(transcript.start_codon_spliced_offsets)
        ref_first_codon = str(
            transcript.sequence)[cds_start:cds_start + 3]
        mut_first_codon = mt.cdna_sequence[
            cds_start:cds_start + 3].upper()
        if (mut_first_codon != ref_first_codon
                and mut_protein
                and mut_protein[0] != "M"
                and ref_protein
                and ref_protein[0] == "M"):
            codon_table = codon_table_for_transcript(transcript)
            if mut_first_codon in codon_table.start_codons:
                mut_protein = "M" + mut_protein[1:]

        # Proteins match → Silent or AlternateStartCodon. Handle
        # both here because the shared classifier doesn't have
        # access to the cDNA edit offset for the correct aa_pos.
        if ref_protein == mut_protein:
            if ref_first_codon != mut_first_codon:
                codon_table = codon_table_for_transcript(transcript)
                if mut_first_codon in codon_table.start_codons:
                    return AlternateStartCodon(
                        variant=variant,
                        transcript=transcript,
                        ref_codon=ref_first_codon,
                        alt_codon=mut_first_codon)
            from ..effects.effect_classes import Silent
            edit = mt.edits[0] if mt.edits else None
            aa_pos = (edit.cdna_start - cds_start) // 3 if edit else 0
            aa_ref = (
                ref_protein[aa_pos]
                if 0 <= aa_pos < len(ref_protein)
                else "")
            return Silent(
                variant=variant,
                transcript=transcript,
                aa_pos=aa_pos,
                aa_ref=aa_ref)

        return classify_from_protein_diff(
            variant=variant,
            transcript=transcript,
            ref_protein=ref_protein,
            mut_protein=mut_protein,
            length_delta=mt.total_length_delta,
            mutant_transcript=mt)
