from abc import ABCMeta

class MutationEffect(object):
    __meta__ = ABCMeta

    def __init__(self, variant, transcript, cdna_pos):
        self.variant = variant
        self.transcript = transcript
        self.cdna_pos = cdna_pos

    def __repr__(self):
        return str(self)

class CodingMutation(MutationEffect):
    def __init__(
            self, variant, cdna_pos,
            aa_pos, aa_ref, aa_alt):
        MutationEffect.__init__(self, variant, transcript, cdna_pos)
        self.aa_pos = aa_pos
        self.aa_ref = aa_ref
        self.aa_alt = aa_alt

    def __str__(self):
        "CodingMutation(%s, %s%d%s)" % (
            self.transcript,
            self.aa_ref,
            self.aa_pos,
            self.aa_alt)

    @property
    def sequence(self):
        original = self.transcript.protein_sequence
        prefix = original[:self.aa_pos]
        suffix = original[self.aa_pos + len(self.aa_ref):]
        return prefix + self.aa_alt + suffix

class FrameShift(MutationEffect):
    def __init__(self, variant, transcript, cdna_pos, aa_pos, shifted_sequence):
        MutationEffect.__init__(self, variant, transcript, cdna_pos)
        self.aa_pos = aa_pos
        self.shifted_sequence = shifted_sequence

    @property
    def sequence(self):
        original_aa_sequence = self.transcript.protein_sequence[:self.aa_pos]
        return original_aa_sequence + self.shifted_sequence
