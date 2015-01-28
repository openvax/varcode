import Bio.Seq
from memoized_property import memoized_property

class TranscriptMutationEffect(object):

    def __init__(self, variant, transcript):
        self.variant = variant
        self.transcript = transcript

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "%s(variant=%s, transcript=%s)" % (
            self.__class__.__name__, self.variant, self.transcript)

    def short_description():
        raise ValueError(
            "Method short_description() not implemented for %s" % self)

    is_coding = False

    @memoized_property
    def original_nucleotide_sequence(self):
        return self.transcript.coding_sequence

    @memoized_property
    def original_protein_sequence(self):
        return Bio.Seq.translate(
            str(self.original_nucleotide_sequence),
            to_stop=True,
            cds=True)

    @memoized_property
    def mutant_protein_sequence(self):
        raise ValueError(
            "mutant_protein_sequence not implemented for %s" % (
                self.__class__.__name__,))

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
    def short_description(self):
        return "intronic"

class Exonic(TranscriptMutationEffect):
    pass


class CodingSequenceMutation(Exonic):
    """
    Base class for all mutations which result in a modified coding sequence.
    """
    def __init__(self, variant, transcript, aa_pos, aa_ref):
        """
        Parameters
        ----------
        variant : varcode.Variant

        transcript : pyensembl.Transcript

        aa_pos : int
            Position of first modified amino aicd (starting from 0)

        aa_ref : str
            Amino acid string of what used to be at aa_pos in the
            wildtype (unmutated) protein.
        """
        Exonic.__init__(self, variant, transcript)
        self.aa_pos = aa_pos
        self.aa_ref = aa_ref

    def __str__(self):
        return "%s(variant=%s, transcript=%s, effect_description=%s)" % (
            self.__class__.__name__,
            self.variant,
            self.transcript,
            self.short_description())

    is_coding = True


class Silent(CodingSequenceMutation):
    """
    Mutation to an exon of a coding region which doesn't change the
    amino acid sequence.
    """
    def short_description(self):
        return "silent"


class Substitution(CodingSequenceMutation):
    """
    Coding mutation which replaces some amino acids into others.
    The total number of amino acids changed must be greater than one on
    either the reference or alternate.
    """
    def __init__(
            self,
            variant,
            transcript,
            aa_pos,
            aa_ref,
            aa_alt):

        CodingSequenceMutation.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_pos=aa_pos,
            aa_ref=aa_ref)

        self.aa_alt = aa_alt
        self.mutation_start = aa_pos
        self.mutation_end = aa_pos + len(aa_alt)

    def short_description(self):
        if len(self.aa_ref) == 0:
            return "p.%dins%s" % (self.aa_pos, self.aa_alt)
        elif len(self.aa_alt) == 0:
            return "p.%s%ddel" % (self.aa_ref, self.aa_pos)
        else:
            return "p.%s%d%s" % (
                    self.aa_ref,
                    self.aa_pos + 1,
                    self.aa_alt)

    @memoized_property
    def mutant_protein_sequence(self):
        original = self.original_protein_sequence
        prefix = original[:self.aa_pos]
        suffix = original[self.aa_pos + len(self.aa_ref):]
        return prefix + self.aa_alt + suffix

class Insertion(Substitution):
    """
    In-frame insertion of one or more amino acids.
    """
    def __init__(self, variant, transcript, aa_pos, aa_alt):
        Substitution.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_pos=aa_pos,
            aa_ref="",
            aa_alt=aa_alt)

class Deletion(Substitution):
    """
    In-frame deletion of one or more amino acids.
    """

    def __init__(self, variant, transcript, aa_pos, aa_ref):
        Substitution.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_pos=aa_pos,
            aa_ref=aa_ref,
            aa_alt="")

class PrematureStop(Substitution):
    def __init__(
            self,
            variant,
            transcript,
            aa_pos,
            aa_ref):
        Substitution.__init__(
            self,
            variant,
            transcript,
            aa_pos=aa_pos,
            aa_ref=aa_ref,
            aa_alt="*")

    def short_description(self):
        return "p.%s%d*" % (
            self.aa_ref,
            self.aa_pos + 1)

    @memoized_property
    def mutant_protein_sequence(self):
        return self.original_protein_sequence[:self.aa_pos]


class UnpredictableSubstitution(Substitution):
    """
    Variants for which we can't confidently determine a protein sequence.

    Splice site mutations are unpredictable since they require a model of
    alternative splicing that goes beyond this library. Similarly,
    when a start codon is lost it's difficult to determine if there is
    an alternative Kozak consensus sequence (either before or after the
    original) from which an alternative start codon can be inferred.
    """

    @property
    def mutant_protein_sequence(self):
        raise ValueError("Can't determine the protein sequence of %s" % self)

class StopLoss(UnpredictableSubstitution):
    def short_description(self):
        return "*%d%s (stop-loss)" % (self.aa_pos, self.aa_alt)


class StartLoss(UnpredictableSubstitution):
    def __init__(
            self,
            variant,
            transcript,
            aa_pos,
            aa_alt):
        UnpredictableSubstitution.__init__(
            self,
            variant,
            transcript,
            aa_pos=aa_pos,
            aa_ref="M",
            aa_alt=aa_alt)

    def short_description(self):
        return "p.? (start-loss)" % (self.aa_pos, self.aa_)

class FrameShift(CodingSequenceMutation):
    def __init__(
            self,
            variant,
            transcript,
            aa_pos,
            aa_ref,
            shifted_sequence):
        """
        Unlike an insertion, which we denote with aa_ref as the chracter before
        the variant sequence, a frameshift starts at aa_ref
        """
        CodingSequenceMutation.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_pos=aa_pos,
            aa_ref=aa_ref)
        self.shifted_sequence = shifted_sequence
        self.mutation_start = self.aa_pos
        self.mutation_end = self.aa_pos + len(shifted_sequence)

    @memoized_property
    def mutant_protein_sequence(self):
        original_aa_sequence = self.original_protein_sequence[:self.aa_pos]
        return original_aa_sequence + self.shifted_sequence

    def short_description(self):
        return "p.%s%dfs" % (self.aa_ref, self.aa_pos + 1)

class FrameShiftTruncation(PrematureStop, FrameShift):
    """
    A frame-shift mutation which immediately introduces a stop codon.
    """
    def __init__(self, *args, **kwargs):
        super(PrematureStop, self).__init__(*args, **kwargs)

    @memoized_property
    def mutant_protein_sequence(self):
        return self.original_protein_sequence[:self.aa_pos]

    def short_description(self):
        return "p.%s%dfs*" % (self.aa_ref, self.aa_pos + 1)

class _MultipleSubstitution(object):
    """
    We're ordering mutations by their class names,
    need to create this dummy class to capture the
    difference between simple subsitutions (e.g. V600E)
    and compound subsitutions (e.g. QF34YY)
    """
    pass

def get_class(effect):
    if isinstance(effect, Substitution):
        if len(effect.aa_ref) > 1 or len(effect.aa_alt) > 1:
            return _MultipleSubstitution
    return effect.__class__

variant_effect_priority_list = [
    IncompleteTranscript,
    NoncodingTranscript,
    # TODO: Add SpliceDonor and SpliceReceptor mutations and #
    # place them at higher priority positions in this list
    Intronic,
    ThreePrimeUTR,
    # mutations to the upstream 5' UTR may change the ORF (reading frame),
    # so give 5' UTR mutations higher prioriry
    FivePrimeUTR,
    Silent,
    Substitution,
    Insertion,
    Deletion,
    _MultipleSubstitution,
    StopLoss,
    PrematureStop,
    # frame-shift which creates immediate stop codon, same as PrematureStop
    FrameShiftTruncation,
    StartLoss,
    FrameShift
]

variant_effect_priority_dict = {
    variant_effect_class : priority
    for (priority, variant_effect_class)
    in enumerate(variant_effect_priority_list)
}

def top_priority_variant_effect(variant_effects):
    """
    Given a collection of variant effects, return the top priority object.
    """
    best_effect = None
    best_priority = -1
    for variant_effect in variant_effects:
        priority = variant_effect_priority_dict[get_class(variant_effect)]
        if priority > best_priority:
            best_effect = variant_effect
            best_priority = priority
    return best_effect