# Copyright (c) 2016. Mount Sinai School of Medicine
#
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

from __future__ import print_function, division, absolute_import

from memoized_property import memoized_property
from serializable import Serializable

from .common import bio_seq_to_str

class MutationEffect(Serializable):
    """
    Base class for mutation effects.
    """

    def __init__(self, variant):
        self.variant = variant

    def __str__(self):
        return "%s(%s)" % (self.__class__.__name__, self.variant)

    def __repr__(self):
        return str(self)

    def __lt__(self, other):
        """
        Effects are ordered by their associated variants, which have
        comparison implement in terms of their chromosomal locations.
        """
        return self.variant < other.variant

    @property
    def short_description(self):
        """
        A short but human-readable description of the effect.
        Defaults to class name for most of the non-coding effects,
        but is more informative for coding ones.
        """
        return self.__class__.__name__.lower()

    transcript = None
    gene = None

    @property
    def original_protein_sequence(self):
        """Amino acid sequence of a coding transcript (without the nucleotide
        variant/mutation)
        """
        if self.transcript:
            return self.transcript.protein_sequence
        else:
            return None

    @property
    def gene_name(self):
        if self.gene:
            return self.gene.name
        else:
            return None

    @property
    def gene_id(self):
        if self.gene:
            return self.gene.id
        else:
            return None

    @property
    def transcript_name(self):
        if self.transcript:
            return self.transcript.name
        else:
            return None

    @property
    def transcript_id(self):
        if self.transcript:
            return self.transcript.id
        else:
            return None

    # It's convenient to have a property which tells us:
    # 1) is this a variant effect overlapping a transcript?
    # 2) does that transcript have a coding sequence?
    # 3) does the variant affect the coding sequence?
    modifies_coding_sequence = False

    # Additionally:
    # 4) Does the change to the coding sequence result in a change to the
    # protein sequence?
    modifies_protein_sequence = False

    mutant_protein_sequence = None

    # which residues in the mutant protein sequence are different from
    # the original?
    aa_mutation_start_offset = None
    aa_mutation_end_offset = None


class Intergenic(MutationEffect):
    """Variant has unknown effect if it occurs between genes"""
    short_description = "intergenic"


class Intragenic(MutationEffect):
    """Variant within boundaries of a gene but does not overlap
    introns or exons of any transcript. This seems very peculiar but
    apparently does happen sometimes, maybe some genes have two distinct sets
    of exons which are never simultaneously expressed?
    """
    short_description = "intragenic"

    def __init__(self, variant, gene):
        MutationEffect.__init__(self, variant)
        self.gene = gene

class TranscriptMutationEffect(Intragenic):
    def __init__(self, variant, transcript):
        Intragenic.__init__(self, variant, gene=transcript.gene)
        self.transcript = transcript

    def __str__(self):
        return "%s(variant=%s, transcript_name=%s, transcript_id=%s)" % (
            self.__class__.__name__,
            self.variant.short_description,
            self.transcript.name,
            self.transcript.id)


class Failure(TranscriptMutationEffect):
    """Special placeholder effect for when we want to suppress errors but still
    need to create a non-empty list of effects for each variant.
    """

class NoncodingTranscript(TranscriptMutationEffect):
    """
    Any mutation to a transcript with a non-coding biotype
    """
    short_description = "non-coding-transcript"

class IncompleteTranscript(TranscriptMutationEffect):
    """
    Any mutation to an incompletely annotated transcript with a coding biotype
    """
    short_description = "incomplete"

class FivePrimeUTR(TranscriptMutationEffect):
    """
    Any mutation to the 5' untranslated region (before the start codon) of
    coding transcript.
    """
    short_description = "5' UTR"

class ThreePrimeUTR(TranscriptMutationEffect):
    """
    Any mutation to the 3' untranslated region (after the stop codon) of
    coding transcript.
    """
    short_description = "3' UTR"


class Intronic(TranscriptMutationEffect):
    """
    Mutation in an intronic region of a coding transcript
    """
    def __init__(self, variant, transcript, nearest_exon, distance_to_exon):
        TranscriptMutationEffect.__init__(self, variant, transcript)
        self.nearest_exon = nearest_exon
        self.distance_to_exon = distance_to_exon

    short_description = "intronic"

class SpliceSite(object):
    """
    Parent class for all splice site mutations.
    """
    pass

class IntronicSpliceSite(Intronic, SpliceSite):
    """
    Mutations near exon boundaries, excluding the first two and last two
    nucleotides in an intron, since those are  known to more confidently
    affect splicing and are given their own effect classes below.
    """
    def __init__(self, *args, **kwargs):
        Intronic.__init__(self, *args, **kwargs)

    short_description = "intronic-splice-site"

class SpliceDonor(IntronicSpliceSite):
    """
    Mutation in the first two intron residues.
    """
    def __init__(self, *args, **kwargs):
        Intronic.__init__(self, *args, **kwargs)

    short_description = "splice-donor"

class SpliceAcceptor(IntronicSpliceSite):
    """
    Mutation in the last two intron residues.
    """
    short_description = "splice-acceptor"

class Exonic(TranscriptMutationEffect):
    """
    Any mutation which affects the contents of an exon (coding region or UTRs)
    """
    pass

class ExonLoss(Exonic):
    """
    Deletion of one or more exons in a transcript.
    """
    def __init__(self, variant, transcript, exons):
        Exonic.__init__(self, variant, transcript)
        self.exons = exons

    def __str__(self):
        return "ExonLoss(%s, %s, %s)" % (
            self.variant,
            self.transcript,
            "+".join(str(exon) for exon in self.exons))

    short_description = "exon-loss"

    @property
    def modifies_protein_sequence(self):
        # TODO: distinguish between exon loss in the CDS and UTRs
        return True

    @property
    def modifies_coding_sequence(self):
        # TODO: distinguish between exon loss in the CDS and UTRs
        return True

class ExonicSpliceSite(Exonic, SpliceSite):
    """
    Mutation in the last three nucleotides before an intron
    or in the first nucleotide after an intron.
    """
    def __init__(self, variant, transcript, exon, alternate_effect):
        Exonic.__init__(self, variant, transcript)
        self.exon = exon
        self.alternate_effect = alternate_effect

    def __str__(self):
        return "ExonicSpliceSite(exon=%s, alternate_effect=%s)" % (
            self.exon,
            self.alternate_effect)

    short_description = "exonic-splice-site"

    @property
    def mutant_protein_sequence(self):
        """
        TODO: determine when exonic splice variants cause exon skipping
        vs. translation of the underlying modified coding sequence.

        For now just pretending like there is no effect on splicing.
        """
        return self.alternate_effect.mutant_protein_sequence

    @property
    def modifies_protein_sequence(self):
        return self.alternate_effect.modifies_protein_sequence

    @property
    def modifies_coding_sequence(self):
        return self.alternate_effect.modifies_coding_sequence

class CodingMutation(Exonic):
    """
    Base class for all mutations which result in a modified coding sequence.
    """
    def __str__(self):
        fields = [
            ("variant", self.variant.short_description),
            ("transcript_name", self.transcript.name),
            ("transcript_id", self.transcript.id),
            ("effect_description", self.short_description)
        ]
        fields_str = ", ".join([
            "%s=%s" % (field_name, field_value)
            for (field_name, field_value)
            in fields
        ])
        return "%s(%s)" % (self.__class__.__name__, fields_str)

    @property
    def modifies_coding_sequence(self):
        return True

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

        aa_ref : str or Bio.Seq
            Reference amino acid(s) at offset
        """
        CodingMutation.__init__(
            self,
            variant=variant,
            transcript=transcript)
        self.aa_pos = aa_pos
        self.aa_ref = bio_seq_to_str(aa_ref)

    @property
    def short_description(self):
        return "silent"

class AlternateStartCodon(Silent):
    """Change to the start codon (e.g. ATG>CTG) but without changing the
    starting amino acid from methionine.
    """
    def __init__(
            self,
            variant,
            transcript,
            ref_codon,
            alt_codon):
        Silent.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_pos=0,
            aa_ref=transcript.protein_sequence[0])
        self.ref_codon = bio_seq_to_str(ref_codon)
        self.alt_codon = bio_seq_to_str(alt_codon)

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
        self.aa_ref = bio_seq_to_str(aa_ref)

    @property
    def modifies_protein_sequence(self):
        return True

class StartLoss(NonsilentCodingMutation):
    """
    When a start codon is lost it's difficult to determine if there is
    an alternative Kozak consensus sequence (either before or after the
    original) from which an alternative start codon can be inferred.
    """
    def __init__(
            self,
            variant,
            transcript):
        NonsilentCodingMutation.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=0,
            aa_mutation_end_offset=1,
            aa_ref=transcript.protein_sequence[0])

    @property
    def mutant_protein_sequence(self):
        # TODO: scan downstream in the cDNA sequence to predict an alternate
        # start codon using a Kozak sequence template
        return None

    @property
    def short_description(self):
        return "p.%s1? (start-loss)" % (self.transcript.protein_sequence[0],)

class KnownAminoAcidChange(NonsilentCodingMutation):
    """
    Coding mutations in which we can predict what the new/mutant protein
    sequence will be.
    """
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_ref,
            aa_alt):
        NonsilentCodingMutation.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_mutation_end_offset=aa_mutation_start_offset + len(aa_alt),
            aa_ref=aa_ref)
        self.aa_alt = bio_seq_to_str(aa_alt)

    @property
    def short_description(self):
        if len(self.aa_ref) == 0:
            return "p.%dins%s" % (
                self.aa_mutation_start_offset,
                self.aa_alt)
        elif len(self.aa_alt) == 0:
            return "p.%s%ddel" % (
                self.aa_ref,
                self.aa_mutation_start_offset)
        else:
            return "p.%s%d%s" % (
                self.aa_ref,
                self.aa_mutation_start_offset + 1,
                self.aa_alt)

    @memoized_property
    def mutant_protein_sequence(self):
        original = self.original_protein_sequence
        prefix = original[:self.aa_mutation_start_offset]
        suffix = original[self.aa_mutation_start_offset + len(self.aa_ref):]
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


class PrematureStop(KnownAminoAcidChange):
    """In-frame insertion of codons containing a stop codon. May also involve
    insertion/deletion/substitution of other amino acids preceding the stop."""
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_ref="",
            aa_alt=""):
        """
        Insertion of premature stop codon, possibly preceded by a substitution
        of `aa_ref` amino acids for `aa_alt` alternative residues.
        """
        assert "*" not in aa_ref, \
            ("Unexpected aa_ref = '%s', should only include amino acids "
             "before the new stop codon.") % aa_ref
        assert "*" not in aa_alt, \
            ("Unexpected aa_ref = '%s', should only include amino acids "
             "before the new stop codon.") % aa_alt
        KnownAminoAcidChange.__init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt=aa_alt)
        self.stop_codon_offset = aa_mutation_start_offset + len(aa_alt)

        assert self.stop_codon_offset < len(transcript.protein_sequence), \
            ("Premature stop codon cannot be at position %d"
             " since the original protein of %s has length %d") % (
                self.stop_codon_offset,
                transcript,
                len(transcript.protein_sequence))

    @property
    def short_description(self):
        return "p.%s%d%s*" % (
            self.aa_ref,
            self.aa_mutation_start_offset + 1,
            self.aa_alt)

    @memoized_property
    def mutant_protein_sequence(self):
        prefix = self.original_protein_sequence[:self.aa_mutation_start_offset]
        return prefix + self.aa_alt


class StopLoss(KnownAminoAcidChange):
    def __init__(
            self,
            variant,
            transcript,
            extended_protein_sequence):
        aa_mutation_start_offset = len(transcript.protein_sequence)
        KnownAminoAcidChange.__init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref="*",
            aa_alt=extended_protein_sequence)

    @property
    def extended_protein_sequence(self):
        return self.aa_alt

    @property
    def short_description(self):
        return "p.*%d%s (stop-loss)" % (
            self.aa_mutation_start_offset + 1,
            self.extended_protein_sequence)


class FrameShift(KnownAminoAcidChange):
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            shifted_sequence):
        """Frameshift mutation preserves all the amino acids up to
        aa_mutation_start_offset and then replaces the rest of the protein with
        new (frameshifted) sequence. Unlike an insertion, where we denote with
        aa_ref as the chracter before the variant sequence, a frameshift starts
        at aa_ref.
        """
        aa_ref = transcript.protein_sequence[aa_mutation_start_offset:]
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

    @memoized_property
    def mutant_protein_sequence(self):
        original_aa_sequence = self.original_protein_sequence
        prefix = original_aa_sequence[:self.aa_mutation_start_offset]
        return prefix + self.shifted_sequence

    @property
    def short_description(self):
        return "p.%s%dfs" % (
            self.aa_ref[0],
            self.aa_mutation_start_offset + 1)

class FrameShiftTruncation(PrematureStop, FrameShift):
    """
    A frame-shift mutation which immediately introduces a stop codon.
    """
    def __init__(
            self,
            variant,
            transcript,
            stop_codon_offset,
            aa_ref=""):
        PrematureStop.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=stop_codon_offset,
            aa_ref=aa_ref,
            aa_alt="")

    @property
    def shifted_sequence(self):
        return self.aa_alt

    @property
    def short_description(self):
        return "p.%s%dfs*" % (
            self.aa_ref,
            self.aa_mutation_start_offset + 1)
