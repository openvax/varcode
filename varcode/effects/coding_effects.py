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

"""Coding-sequence and amino-acid consequence classes."""

from __future__ import annotations

from functools import cached_property

from .base import Exonic, TranscriptMutationEffect


class CodingMutation(Exonic):
    """
    Base class for all mutations which result in a modified coding sequence.
    """
    def __str__(self):
        fields = [
            ("variant", self.variant),
            ("transcript_name", self.transcript_name),
            ("transcript_id", self.transcript_id),
            ("effect_description", self.short_description)
        ]
        if self.orf is not None:
            fields.append(("orf_id", self.orf.orf_id))
        fields_str = ", ".join([
            "%s=%s" % (field_name, field_value)
            for (field_name, field_value)
            in fields
        ])
        return "%s(%s)" % (self.__class__.__name__, fields_str)

    _default_modifies_coding_sequence = True


class Silent(CodingMutation):
    """Mutation to an exon of a coding region which doesn't change the
    amino acid sequence.
    """
    def __init__(
            self,
            variant,
            transcript,
            aa_pos,
            aa_ref):
        """
        Parameters
        ----------
        variant : Variant

        transcript : Transcript

        aa_pos : int
            Offset of first synonymous codon in protein sequence

        aa_ref : str
            Reference amino acid(s) at offset
        """
        CodingMutation.__init__(
            self,
            variant=variant,
            transcript=transcript)
        self.aa_pos = aa_pos
        self.aa_ref = str(aa_ref)

    @property
    def short_description(self):
        # HGVS p.= notation for synonymous changes: "p.{aa_ref}{1-indexed pos}="
        # (http://varnomen.hgvs.org/). If aa_ref is empty (the silent
        # change was fully captured by the shared prefix/suffix and we
        # have nothing to name), fall back to a positional-only form.
        if self.aa_ref:
            return "p.%s%d=" % (self.aa_ref, self.aa_pos + 1)
        return "p.%d=" % (self.aa_pos + 1)


class AlternateStartCodon(Silent):
    """Change to the start codon (e.g. ATG>CTG) but without changing the
    starting amino acid from methionine.
    """
    def __init__(
            self,
            variant,
            transcript,
            ref_codon,
            alt_codon,
            aa_ref: str | None = None):
        Silent.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_pos=0,
            aa_ref=transcript.protein_sequence[0] if aa_ref is None else aa_ref)
        self.ref_codon = str(ref_codon)
        self.alt_codon = str(alt_codon)

    @property
    def short_description(self):
        return "alternate-start-codon (%s>%s)" % (
            self.ref_codon, self.alt_codon)


class NonsilentCodingMutation(CodingMutation):
    """
    All coding mutations other than silent codon substitutions
    """

    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_mutation_end_offset,
            aa_ref):
        """
        variant : Variant

        transcript : Transcript

        aa_mutation_start_offset : int
            Offset of first modified amino acid in protein (starting from 0)

        aa_mutation_end_offset : int
            Offset after last mutated amino acid (half-open coordinates)

        aa_ref : str
            Amino acid string of what used to be at aa_mutation_start_offset
            in the wildtype (unmutated) protein.
        """
        CodingMutation.__init__(
            self,
            variant=variant,
            transcript=transcript)
        self.aa_mutation_start_offset = aa_mutation_start_offset
        self.aa_mutation_end_offset = aa_mutation_end_offset
        self.aa_ref = str(aa_ref)

    _default_modifies_protein_sequence = True


class StartGain(TranscriptMutationEffect):
    """A newly created initiation hypothesis on a paired mutant RNA.

    ``discovered_orf`` describes the mutant frame; ``orf`` remains None because
    no reference ORF is asserted. The conditional peptide is available through
    the usual mutant accessor. Baseline translation and actual initiation are
    unknown, so the protein-change flag and original protein remain None.
    """

    _default_modifies_coding_sequence = True
    _default_modifies_protein_sequence = None

    def __init__(self, variant, transcript, discovered_orf, baseline_transcript,
                 mutant_transcript, discovery_policy, hypothesis=None, source_evidence=None):
        super().__init__(variant, transcript)
        self.discovered_orf = discovered_orf
        self.baseline_transcript = baseline_transcript
        self.mutant_transcript = mutant_transcript
        self.discovery_policy = discovery_policy
        self.hypothesis = hypothesis
        self.source_evidence = source_evidence or {}
        self.annotator = "created_orf"
        from ..version import __version__
        self.annotator_version = __version__

    @property
    def short_description(self) -> str:
        orf = self.discovered_orf
        if orf is None:
            return "start-gain (conditional)"
        return "start-gain %s at RNA offset %d (conditional)" % (
            orf.start_codon, orf.start)


class StartLoss(NonsilentCodingMutation):
    """Loss of the annotated start, with available nucleotide context.

    ``ref_codon`` and ``alt_codon`` are in transcript orientation. Their
    offsets are zero-based in the baseline and mutant cDNA, respectively.
    A removed or unmapped start has no alternate codon/offset. Edits and
    sequence remain on ``mutant_transcript`` when available. Neither these
    details nor a codon's internal amino-acid meaning establish initiation.
    The protein outcome remains unknown; alternative-start discovery is
    outside this effect's scope.
    """
    def __init__(
            self,
            variant,
            transcript,
            mutant_transcript=None,
            ref_codon=None,
            alt_codon=None,
            ref_codon_offset=None,
            alt_codon_offset=None,
            codon_change_status=None,
            aa_ref: str | None = None):
        NonsilentCodingMutation.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=0,
            aa_mutation_end_offset=1,
            aa_ref=transcript.protein_sequence[0] if aa_ref is None else aa_ref)
        if codon_change_status is None:
            from .start_loss import start_loss_details
            details = start_loss_details(variant, transcript, mutant_transcript)
            mutant_transcript = details["mutant_transcript"]
            ref_codon = details["ref_codon"] if ref_codon is None else ref_codon
            alt_codon = details["alt_codon"] if alt_codon is None else alt_codon
            ref_codon_offset = details["ref_codon_offset"] if ref_codon_offset is None else ref_codon_offset
            alt_codon_offset = details["alt_codon_offset"] if alt_codon_offset is None else alt_codon_offset
            codon_change_status = details["codon_change_status"]
        self.mutant_transcript = mutant_transcript
        self.ref_codon = ref_codon
        self.alt_codon = alt_codon
        self.ref_codon_offset = ref_codon_offset
        self.alt_codon_offset = alt_codon_offset
        self.codon_change_status = codon_change_status

    def __hash__(self):
        # Transcript-model evidence can contain dictionaries.
        return hash((type(self), self.variant, self.transcript))

    @property
    def alt_codon_amino_acid(self):
        """Internal decoding of the alternate triplet, not an initiated protein residue."""
        if self.alt_codon is None:
            return None
        from .codon_tables import codon_table_for_transcript
        table = codon_table_for_transcript(self.transcript)
        if self.alt_codon in table.stop_codons:
            return "*"
        return table.forward_table.get(self.alt_codon)

    def _get_mutant_protein_sequence(self) -> str | None:
        return None

    @property
    def short_description(self) -> str:
        return "p.%s1? (start-loss)" % (self.aa_ref,)


class KnownAminoAcidChange(NonsilentCodingMutation):
    """
    Coding mutations in which we can predict what the new/mutant protein
    sequence will be.
    """
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset: int,
            aa_ref,
            aa_alt):
        NonsilentCodingMutation.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_mutation_end_offset=aa_mutation_start_offset + len(aa_alt),
            aa_ref=aa_ref)
        self.aa_alt = str(aa_alt)

    def _required_start_offset(self) -> int:
        offset = self.aa_mutation_start_offset
        if offset is None:
            raise ValueError(
                "%s requires a non-None aa_mutation_start_offset" % type(self).__name__)
        return offset

    def _required_original_protein_sequence(self) -> str:
        original = self.original_protein_sequence
        if original is None:
            raise ValueError(
                "%s requires a non-None original_protein_sequence" % type(self).__name__)
        return original

    @property
    def short_description(self) -> str:
        start = self._required_start_offset()
        if len(self.aa_ref) == 0:
            return "p.%dins%s" % (
                start,
                self.aa_alt)
        elif len(self.aa_alt) == 0:
            # Effect offsets stay zero-based; labels name the one-based start.
            return "p.%s%ddel" % (self.aa_ref, start + 1)
        else:
            return "p.%s%d%s" % (
                self.aa_ref,
                start + 1,
                self.aa_alt)

    @cached_property
    def mutant_protein_sequence(self) -> str | None:
        if self.mutant_transcript is not None:
            return self.mutant_transcript.mutant_protein_sequence
        original = self._required_original_protein_sequence()
        start = self._required_start_offset()
        prefix = original[:start]
        suffix = original[start + len(self.aa_ref):]
        return prefix + self.aa_alt + suffix


class Substitution(KnownAminoAcidChange):
    """
    Single amino acid substitution, e.g. BRAF-001 V600E
    """
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_ref,
            aa_alt):
        if len(aa_ref) != 1:
            raise ValueError(
                "Simple substitution can't have aa_ref='%s'" % (aa_ref,))
        if len(aa_alt) != 1:
            raise ValueError(
                "Simple substitution can't have aa_alt='%s'" % (aa_alt,))
        KnownAminoAcidChange.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt=aa_alt)


class ComplexSubstitution(KnownAminoAcidChange):
    """
    In-frame substitution of multiple amino acids, e.g. TP53-002 p.391FY>QQQ
    Can change the length of the protein sequence but since it has
    non-empty ref and alt strings, is more complicated than an insertion or
    deletion alone.
    """
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_ref,
            aa_alt):
        if len(aa_ref) == 1 and len(aa_alt) == 1:
            raise ValueError(
                "ComplexSubstitution can't have aa_ref='%s' and aa_alt='%s'" % (
                    aa_ref, aa_alt))
        KnownAminoAcidChange.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt=aa_alt)


class Insertion(KnownAminoAcidChange):
    """
    In-frame insertion of one or more amino acids.
    """
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_alt):
        KnownAminoAcidChange.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref="",
            aa_alt=aa_alt)


class Deletion(KnownAminoAcidChange):
    """
    In-frame deletion of one or more amino acids.
    """

    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_ref):
        KnownAminoAcidChange.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt="")

    def _format_short_description(self, deletion_notation):
        if len(self.aa_ref) <= 1:
            return self.short_description
        start = self._required_start_offset() + 1
        return "p.%s%d_%s%ddel" % (
            self.aa_ref[0], start,
            self.aa_ref[-1], start + len(self.aa_ref) - 1)


class PrematureStop(KnownAminoAcidChange):
    """In-frame insertion of codons containing a stop codon. May also involve
    insertion/deletion/substitution of other amino acids preceding the stop."""
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_ref="",
            aa_alt="",
            original_protein_sequence: str | None = None):
        """
        Insertion of premature stop codon, possibly preceded by a substitution
        of `aa_ref` amino acids for `aa_alt` alternative residues.

        ``original_protein_sequence`` supplies the patient baseline when it
        differs from the reference transcript; coordinates refer to that protein.
        """
        if original_protein_sequence is not None:
            self._original_protein_sequence = original_protein_sequence
        if "*" in aa_ref:
            raise ValueError(
                ("Unexpected aa_ref = '%s', should only include amino acids "
                 "before the new stop codon.") % aa_ref)
        if "*" in aa_alt:
            raise ValueError(
                ("Unexpected aa_ref = '%s', should only include amino acids "
                 "before the new stop codon.") % aa_alt)
        KnownAminoAcidChange.__init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt=aa_alt)
        self.stop_codon_offset = aa_mutation_start_offset + len(aa_alt)

        original = self._required_original_protein_sequence()
        if self.stop_codon_offset >= len(original):
            raise ValueError(
                ("Premature stop codon cannot be at position %d"
                 " since the original protein of %s has length %d") % (
                    self.stop_codon_offset,
                    transcript,
                    len(original)))

    @property
    def short_description(self) -> str:
        # When aa_ref is empty, the mutation is an insertion (no residues
        # were replaced) — use HGVS "ins" notation so the description is
        # unambiguous. Otherwise fall through to "{ref}{pos}{alt}*".
        start = self._required_start_offset()
        if len(self.aa_ref) == 0:
            return "p.%dins%s*" % (
                start + 1,
                self.aa_alt)
        return "p.%s%d%s*" % (
            self.aa_ref,
            start + 1,
            self.aa_alt)

    @cached_property
    def mutant_protein_sequence(self) -> str | None:
        if self.mutant_transcript is not None:
            return self.mutant_transcript.mutant_protein_sequence
        prefix = self._required_original_protein_sequence()[:self._required_start_offset()]
        return prefix + self.aa_alt


class StopLoss(KnownAminoAcidChange):
    def __init__(
            self,
            variant,
            transcript,
            aa_ref,
            aa_alt,
            original_protein_sequence: str | None = None):
        """Extend the supplied baseline protein, or the reference when omitted."""
        if original_protein_sequence is not None:
            self._original_protein_sequence = original_protein_sequence
        # StopLoss assumes that we deleted some codons ending with a
        # stop codon
        if "*" in aa_ref:
            raise ValueError(
                "StopLoss aa_ref '%s' should not contain '*'" % (
                    aa_ref,))
        if len(aa_alt) == 0:
            raise ValueError(
                "If no amino acids added by StopLoss then it should be Silent")
        # subtract 1 for the stop codon
        n_ref_amino_acids = len(aa_ref)
        protein_length = len(original_protein_sequence if original_protein_sequence is not None
                             else transcript.protein_sequence)
        aa_mutation_start_offset = protein_length - n_ref_amino_acids
        KnownAminoAcidChange.__init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_alt=aa_alt,
            aa_ref=aa_ref)

    @property
    def short_description(self) -> str:
        return "p.%s*%d%s (stop-loss)" % (
            self.aa_ref,
            self._required_start_offset() + 1,
            self.aa_alt)


class FrameShift(KnownAminoAcidChange):
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            shifted_sequence,
            original_protein_sequence: str | None = None):
        """Frameshift mutation preserves all the amino acids up to
        aa_mutation_start_offset and then replaces the rest of the protein with
        new (frameshifted) sequence. Unlike an insertion, where we denote with
        aa_ref as the chracter before the variant sequence, a frameshift starts
        at aa_ref.

        ``original_protein_sequence`` is the patient baseline when supplied.
        """
        if original_protein_sequence is not None:
            self._original_protein_sequence = original_protein_sequence
        protein = (original_protein_sequence if original_protein_sequence is not None
                   else transcript.protein_sequence)
        aa_ref = protein[aa_mutation_start_offset:]
        KnownAminoAcidChange.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt=shifted_sequence)

    @property
    def shifted_sequence(self):
        return self.aa_alt

    @cached_property
    def mutant_protein_sequence(self) -> str | None:
        if self.mutant_transcript is not None:
            return self.mutant_transcript.mutant_protein_sequence
        original_aa_sequence = self._required_original_protein_sequence()
        prefix = original_aa_sequence[:self._required_start_offset()]
        return prefix + self.shifted_sequence

    @property
    def short_description(self) -> str:
        return "p.%s%dfs" % (
            self.aa_ref[0],
            self._required_start_offset() + 1)


class FrameShiftTruncation(PrematureStop, FrameShift):
    """
    A frame-shift mutation which immediately introduces a stop codon.
    """
    def __init__(
            self,
            variant,
            transcript,
            stop_codon_offset,
            aa_ref="",
            original_protein_sequence: str | None = None):
        PrematureStop.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=stop_codon_offset,
            aa_ref=aa_ref,
            aa_alt="",
            original_protein_sequence=original_protein_sequence)

    @property
    def shifted_sequence(self):
        return self.aa_alt

    @property
    def short_description(self) -> str:
        return "p.%s%dfs*" % (
            self.aa_ref,
            self._required_start_offset() + 1)



__all__ = [
    "CodingMutation",
    "Silent",
    "AlternateStartCodon",
    "NonsilentCodingMutation",
    "StartGain",
    "StartLoss",
    "KnownAminoAcidChange",
    "Substitution",
    "ComplexSubstitution",
    "Insertion",
    "Deletion",
    "PrematureStop",
    "StopLoss",
    "FrameShift",
    "FrameShiftTruncation",
]
