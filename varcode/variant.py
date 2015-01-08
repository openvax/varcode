from nucleotides import normalize_nucleotide_string

from pyensembl.locus import normalize_chromosome

class Variant(object):
    def __init__(self, contig, pos, ref, alt, info=None):
        self.contig = normalize_chromosome(contig)
        self.pos = int(pos)
        self.ref = normalize_nucleotide_string(ref)
        self.alt = normalize_nucleotide_string(alt)
        self.info = {} if info is None else info

    @property
    def end_pos(self):
        return self.pos + len(self.ref)

    def __str__(self):
        return "Variant(contig=%s, pos=%d, ref=%s, alt=%s, info=%s)" % (
            self.contig, self.pos, self.ref, self.alt, self.info)

    def __repr__(self):
        return str(self)