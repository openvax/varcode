# Copyright (c) 2015. Mount Sinai School of Medicine
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


class MutationEffect(object):
    """Base class for mutation effects which don't overlap a transcript"""

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

    def short_description(self):
        raise ValueError(
            "Method short_description() not implemented for %s" % (
                self.__class__.__name__,))

    transcript = None
    gene = None

    @memoized_property
    def original_protein_sequence(self):
        """Amino acid sequence of a coding transcript (without the nucleotide
        variant/mutation)
        """
        if self.transcript:
            return self.transcript.protein_sequence
        else:
            return None

    @memoized_property
    def gene_name(self):
        if self.gene:
            return self.gene.name
        else:
            return None

    @memoized_property
    def gene_id(self):
        if self.gene:
            return self.gene.id
        else:
            return None

    @memoized_property
    def transcript_name(self):
        if self.transcript:
            return self.transcript.name
        else:
            return None

    @memoized_property
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
    pass

class Intragenic(MutationEffect):
    """Variant within boundaries of a gene but does not overlap
    introns or exons of any transcript. This seems very peculiar but
    apparently does happen sometimes, maybe some genes have two distinct sets
    of exons which are never simultaneously expressed?
    """

    def __init__(self, variant, gene):
        MutationEffect.__init__(self, variant)
        self.gene = gene

class TranscriptMutationEffect(Intragenic):
    def __init__(self, variant, transcript):
        Intragenic.__init__(self, variant, gene=transcript.gene)
        self.transcript = transcript

    def __str__(self):
        return "%s(variant=%s, transcript=%s)" % (
            self.__class__.__name__,
            self.variant.short_description(),
            self.transcript.name)


class Failure(TranscriptMutationEffect):
    """Special placeholder effect for when we want to suppress errors but still
    need to create a non-empty list of effects for each variant.
    """

class NoncodingTranscript(TranscriptMutationEffect):
    """
    Any mutation to a transcript with a non-coding biotype
    """
    def short_description(self):
        return "non-coding-transcript"

class IncompleteTranscript(TranscriptMutationEffect):
    """
    Any mutation to an incompletely annotated transcript with a coding biotype
    """
    def short_description(self):
        return "incomplete"

class FivePrimeUTR(TranscriptMutationEffect):
    """
    Any mutation to the 5' untranslated region (before the start codon) of
    coding transcript.
    """
    def short_description(self):
        return "5' UTR"

class ThreePrimeUTR(TranscriptMutationEffect):
    """
    Any mutation to the 3' untranslated region (after the stop codon) of
    coding transcript.
    """
    def short_description(self):
        return "3' UTR"


class Intronic(TranscriptMutationEffect):
    """
    Mutation in an intronic region of a coding transcript
    """
    def __init__(self, variant, transcript, nearest_exon, distance_to_exon):
        TranscriptMutationEffect.__init__(self, variant, transcript)
        self.nearest_exon = nearest_exon
        self.distance_to_exon = distance_to_exon

    def short_description(self):
        return "intronic"

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

    def short_description(self):
        return "intronic-splice-site"

class SpliceDonor(IntronicSpliceSite):
    """
    Mutation in the first two intron residues.
    """
    def __init__(self, *args, **kwargs):
        Intronic.__init__(self, *args, **kwargs)

    def short_description(self):
        return "splice-donor"

class SpliceAcceptor(IntronicSpliceSite):
    """
    Mutation in the last two intron residues.
    """
    def short_description(self):
        return "splice-acceptor"

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
            "+".join(self.exons))

    def short_description(self):
        return "exon-loss"

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

    def short_description(self):
        return "exonic-splice-site"

    @memoized_property
    def mutant_protein_sequence(self):
        """
        TODO: determine when exonic splice variants cause exon skipping
        vs. translation of the underlying modified coding sequence.

        For now just pretending like there is no effect on splicing.
        """
        return self.alternate_effect.mutant_protein_sequence

class CodingMutation(Exonic):
    """
    Base class for all mutations which result in a modified coding sequence.
    """
    def __str__(self):
        return "%s(variant=%s, transcript=%s, effect_description=%s)" % (
            self.__class__.__name__,
            self.variant.short_description(),
            self.transcript.name,
            self.short_description())

    modifies_coding_sequence = True


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
        self.aa_ref = aa_ref

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
            aa_ref,
            ref_codon,
            alt_codon):
        Silent.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_pos=0,
            aa_ref=aa_ref)
        self.ref_codon = ref_codon
        self.alt_codon = alt_codon

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
            Amino acid string of what used to be at mutated_protein_start_offset
            in the wildtype (unmutated) protein.
        """
        CodingMutation.__init__(
            self,
            variant=variant,
            transcript=transcript)
        self.aa_mutation_start_offset = aa_mutation_start_offset
        self.aa_mutation_end_offset = aa_mutation_end_offset
        self.aa_ref = aa_ref

    modifies_protein_sequence = True

class BaseSubstitution(NonsilentCodingMutation):
    """
    Coding mutation which replaces some amino acids into others.
    The total number of amino acids changed must be greater than one on
    either the reference or alternate.
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
        self.aa_alt = aa_alt

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

class Substitution(BaseSubstitution):
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
        BaseSubstitution.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt=aa_alt)

class ComplexSubstitution(BaseSubstitution):
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
        BaseSubstitution.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt=aa_alt)

class Insertion(BaseSubstitution):
    """
    In-frame insertion of one or more amino acids.
    """
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_alt):
        BaseSubstitution.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref="",
            aa_alt=aa_alt)

class Deletion(BaseSubstitution):
    """
    In-frame deletion of one or more amino acids.
    """

    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_ref):
        BaseSubstitution.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt="")


class PrematureStop(BaseSubstitution):
    """In-frame insertion of codons ending with a stop codon. May also involve
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
        BaseSubstitution.__init__(
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

    def short_description(self):
        return "p.%s%d%s*" % (
            self.aa_ref,
            self.aa_mutation_start_offset + 1,
            self.aa_alt)

    @memoized_property
    def mutant_protein_sequence(self):
        prefix = self.original_protein_sequence[:self.aa_mutation_start_offset]
        return prefix + self.aa_alt


class StopLoss(BaseSubstitution):
    def __init__(
            self,
            variant,
            transcript,
            extended_protein_sequence):
        aa_mutation_start_offset = len(transcript.protein_sequence)
        self.extended_protein_sequence = extended_protein_sequence
        BaseSubstitution.__init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref="*",
            aa_alt=extended_protein_sequence)

    def short_description(self):
        return "p.*%d%s (stop-loss)" % (
            self.aa_mutation_start_offset + 1,
            self.extended_protein_sequence)

class StartLoss(BaseSubstitution):
    """
    When a start codon is lost it's difficult to determine if there is
    an alternative Kozak consensus sequence (either before or after the
    original) from which an alternative start codon can be inferred.
    """
    def __init__(
            self,
            variant,
            transcript,
            aa_alt="?"):
        BaseSubstitution.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=0,
            aa_ref="M",
            aa_alt=aa_alt)

    @property
    def mutant_protein_sequence(self):
        return None

    def short_description(self):
        return "p.%d? (start-loss)" % (self.aa_mutation_start_offset,)

class FrameShift(NonsilentCodingMutation):
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_ref,
            shifted_sequence):
        """Frameshift mutation preserves all the amino acids up to
        aa_mutation_start_offset and then replaces the rest of the protein with
        new (frameshifted) sequence. Unlike an insertion, where we denote with
        aa_ref as the chracter before the variant sequence, a frameshift starts
        at aa_ref.
        """
        n_new_amino_acids = len(shifted_sequence)
        NonsilentCodingMutation.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_mutation_end_offset=aa_mutation_start_offset + n_new_amino_acids,
            aa_ref=aa_ref)
        self.shifted_sequence = shifted_sequence

    @memoized_property
    def mutant_protein_sequence(self):
        original_aa_sequence = self.original_protein_sequence
        prefix = original_aa_sequence[:self.mutated_protein_start_offset]
        return prefix + self.shifted_sequence

    def short_description(self):
        return "p.%s%dfs" % (
            self.aa_ref,
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
        self.stop_codon_offset = stop_codon_offset
        self.shifted_sequence = ""
        NonsilentCodingMutation.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=stop_codon_offset,
            aa_mutation_end_offset=stop_codon_offset + 1,
            aa_ref=aa_ref)

    def short_description(self):
        return "p.%s%dfs*" % (
            self.aa_ref,
            self.aa_mutation_start_offset + 1)
