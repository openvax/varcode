import Bio.Seq

class MutationEffect(object):

    def __init__(self, variant, transcript):
        self.variant = variant
        self.transcript = transcript

    def __repr__(self):
        return str(self)

    def __str__(self):
        raise ValueError(
            "__str__ not implemented for %s" % self.__class__.__name__)

    @property
    def original_nucleotide_sequence(self):
        return self.transcript.sequence

    @property
    def original_protein_sequence(self):
        return Bio.Seq.translate(
            self.original_nucleotide_sequence,
            to_stop=True,
            cds=True)

    @property
    def mutant_protein_sequence(self):
        raise ValueError(
            "mutant_protein_sequence not implemented for %s" % (
                self.__class__.__name__,)

class NoncodingTranscript(MutationEffect):
    """
    Any mutation to a transcript with a non-coding biotype
    """
    def __str__(self):
        return "non-coding-transcript"

class IncompleteTranscript(MutationEffect):
    """
    Any mutation to an incompletely annotated transcript with a coding biotype
    """
    def __str__(self):
        return "incomplete"

class FivePrimeUTR(MutationEffect):
    """
    Any mutation to the 5' untranslated region (before the start codon) of
    coding transcript.
    """
    def __str__(self):
        retunr "5' UTR"


class ThreePrimeUTR(MutationEffect):
    """
    Any mutation to the 3' untranslated region (after the stop codon) of
    coding transcript.
    """
    def __str__(self):
        retunr "3' UTR"

class Intronic(MutationEffect):
    """
    Mutation in an intronic region of a coding transcript
    """
    def __str__(self):
        return "intronic"

class Exonic(MutationEffect):
    def __init__(self, variant, transcript, cdna_pos):
        MutationEffect.__init__(self, variant, transcript)
        self.cdna_pos = cdna_pos

class Silent(Exonic):
    """
    Mutation to an exon of a coding region which doesn't change the
    amino acid sequence.
    """
    def __str__(self):
        return "silent"

class Coding(Exonic):
    def is_coding(self):
        return True

class Substitution(Coding):
    """
    Coding mutation which removes or inserts amino acids at a locus.
    """
    def __init__(
            self, variant, cdna_pos,
            aa_pos, aa_ref, aa_alt):
        Coding.__init__(self, variant, transcript, cdna_pos)
        self.aa_pos = aa_pos
        self.aa_ref = aa_ref
        self.aa_alt = aa_alt

    def __str__(self):
        return "%s p.%s%d%s" % (
            self.transcript.name,
            self.aa_ref,
            self.aa_pos,
            self.aa_alt)

    @property
    def mutant_protein_sequence(self):
        original = self.original_protein_sequence
        prefix = original[:self.aa_pos]
        suffix = original[self.aa_pos + len(self.aa_ref):]
        return prefix + self.aa_alt + suffix


class FrameShift(Exonic):
    def __init__(self, variant, transcript, cdna_pos, aa_pos, shifted_sequence):
        Coding.__init__(self, variant, transcript, cdna_pos)
        self.aa_pos = aa_pos
        self.shifted_sequence = shifted_sequence

    @property
    def mutant_sequence(self):
        original_aa_sequence = self.original_protein_sequence[:self.aa_pos]
        return original_aa_sequence + self.shifted_sequence