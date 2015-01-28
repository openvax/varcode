from nucleotides import normalize_nucleotide_string

from memoized_property import memoized_property
from pyensembl.locus import normalize_chromosome

class Variant(object):
    def __init__(self, contig, pos, ref, alt, info=None):
        self.contig = normalize_chromosome(contig)
        self.pos = int(pos)
        self.ref = normalize_nucleotide_string(ref)
        self.alt = normalize_nucleotide_string(alt)
        self.info = {} if info is None else info

    @memoized_property
    def end_pos(self):
        return self.pos + len(self.ref) - 1

    def __str__(self):
        return "Variant(contig=%s, pos=%d, ref=%s, alt=%s, info=%s)" % (
            self.contig, self.pos, self.ref, self.alt, self.info)

    def __repr__(self):
        return str(self)

    def short_description(self):
        chrom, pos, ref, alt = self.contig, self.pos, self.ref, self.alt
        if ref == alt:
            return "chr%s g.%d %s=%s" % (chrom, pos, ref, alt)
        elif len(ref) == 0 or alt.startswith(ref):
            return "chr%s g.%d ins%s" % (chrom, pos + len(ref),  alt[len(ref):])
        elif len(alt) == 0 or ref.startswith(alt):
            return "chr%s g.%d_%d del%s" % (
                chrom, pos + len(alt), pos + len(ref), ref[len(alt):])
        else:
            return "chr%s g.%d %s>%s" % (chrom, pos, ref, alt)