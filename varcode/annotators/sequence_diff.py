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

"""SequenceDiffEffectAnnotator — classify effects from protein diff
instead of offset arithmetic (openvax/varcode#271 stage 3d, #309).

See :doc:`/effect_annotation` for the user-facing guide covering
how legacy, sequence-diff, splice outcomes, and the SV roadmap
fit together. This module docstring captures implementation
detail relevant to reviewers of the classifier itself.

Algorithm
---------

1. **Gate**: non-coding / incomplete transcripts short-circuit to
   the standard wrappers.

2. **Splice check**: run legacy first. If legacy classifies the
   variant as a splice effect (SpliceDonor, SpliceAcceptor,
   IntronicSpliceSite), return that directly — sequence-diff
   doesn't own splice classification. For ExonicSpliceSite, run
   dual-dispatch: legacy provides the splice class, sequence-diff
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
    IncompleteTranscript,
    IntronicSpliceSite,
    NoncodingTranscript,
    SpliceAcceptor,
    SpliceDonor,
)
from ..mutant_transcript import apply_variant_to_transcript
from ..version import __version__ as _varcode_version
from .legacy import LegacyEffectAnnotator


class SequenceDiffEffectAnnotator:
    """Classify effects by diffing translated mutant protein against
    the reference protein.

    Produces byte-for-byte identical output to
    :class:`LegacyEffectAnnotator` on the common case (trivial
    SNVs and simple indels) because both flow through the same
    :func:`classify_from_protein_diff` classifier. Diverges where
    sequence-diff's approach is provably more accurate (boundary
    codons, frameshift realignment). Any divergence must appear in
    the parity harness ``EXPECTED_DIFFS`` with an issue link.
    """

    name = "sequence_diff"
    version = _varcode_version
    supports = frozenset({"snv", "indel", "mnv"})

    def __repr__(self):
        return "SequenceDiffEffectAnnotator(name=%r, version=%r)" % (
            self.name, self.version)

    def annotate_on_transcript(self, variant, transcript):
        """Classify the effect of ``variant`` on ``transcript``.

        Runs legacy first to detect splice-adjacent variants (which
        stay legacy-classified); for everything else, builds a
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

        # Run legacy to get the splice classification.
        legacy_effect = LegacyEffectAnnotator().annotate_on_transcript(
            variant, transcript)

        # Pure-intronic splice effects: legacy only — no protein-
        # level diff to compute.
        if isinstance(legacy_effect, (SpliceDonor, SpliceAcceptor)):
            return legacy_effect
        if (isinstance(legacy_effect, IntronicSpliceSite)
                and not isinstance(legacy_effect, (SpliceDonor, SpliceAcceptor))):
            return legacy_effect

        # ExonicSpliceSite: dual-dispatch. Legacy provides the
        # splice class; sequence-diff provides the alternate_effect
        # via protein diff.
        if isinstance(legacy_effect, ExonicSpliceSite):
            mt = apply_variant_to_transcript(variant, transcript)
            if mt is not None and mt.mutant_protein_sequence is not None:
                alt = classify_from_protein_diff(
                    variant=variant,
                    transcript=transcript,
                    ref_protein=str(transcript.protein_sequence),
                    mut_protein=mt.mutant_protein_sequence,
                    length_delta=mt.total_length_delta)
                return ExonicSpliceSite(
                    variant=variant,
                    transcript=transcript,
                    exon=legacy_effect.exon,
                    alternate_effect=alt)
            return legacy_effect

        # Non-splice: sequence-diff slow path.
        mt = apply_variant_to_transcript(variant, transcript)
        if mt is None or mt.mutant_protein_sequence is None:
            # UTR, ref-mismatch, splice-junction-spanning, etc.
            return legacy_effect

        ref_protein = str(transcript.protein_sequence)
        mut_protein = mt.mutant_protein_sequence

        # AlternateStartCodon: non-diff special case (#304).
        # Proteins match but the first codon changed to a non-ATG
        # alternate start — informative for translation-efficiency
        # analysis, so we emit it rather than collapsing into Silent.
        if ref_protein == mut_protein:
            cds_start = min(transcript.start_codon_spliced_offsets)
            ref_first_codon = str(
                transcript.sequence)[cds_start:cds_start + 3]
            mut_first_codon = mt.cdna_sequence[
                cds_start:cds_start + 3].upper()
            if ref_first_codon != mut_first_codon:
                codon_table = codon_table_for_transcript(transcript)
                if mut_first_codon in codon_table.start_codons:
                    return AlternateStartCodon(
                        variant=variant,
                        transcript=transcript,
                        ref_codon=ref_first_codon,
                        alt_codon=mut_first_codon)

        return classify_from_protein_diff(
            variant=variant,
            transcript=transcript,
            ref_protein=ref_protein,
            mut_protein=mut_protein,
            length_delta=mt.total_length_delta)
