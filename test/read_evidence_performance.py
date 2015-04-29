"""
Simple script to measure performance of the read_evidence module.

Usage:
    %(prog)s /path/to/reads.bam <region1> <region2> ... --num-loci-per-region X

Example:
    %(prog)s my.bam 8:101723999-101725999 --num-loci-per-region 2000

"""
import argparse
import sys
import random
import time

import varcode
import varcode.read_evidence
from varcode.locus import Locus
    
PARSER = argparse.ArgumentParser(usage=__doc__)
PARSER.add_argument("bam_path")
PARSER.add_argument("regions", nargs="+")
PARSER.add_argument("--num-loci-per-region", type=int, default=100)

def parse_region(s):
    (contig, rest) = s.split(":")
    (start, end) = rest.split("-")
    return varcode.Locus.from_inclusive_coordinates(
        contig, int(start), int(end))

def go(argv):
    args = PARSER.parse_args(argv)

    loci_regions = [parse_region(s) for s in args.regions]

    loci = []
    for region in loci_regions:
        new_loci = random.sample(
            range(region.start, region.end), args.num_loci_per_region)
        loci.extend(
            Locus.from_inclusive_coordinates(region.contig, locus)
            for locus in new_loci)

    print("Loading pileups for %d loci." % len(loci))

    start = time.time()
    varcode.read_evidence.PileupCollection.from_bam(args.bam_path, loci)
    elapsed = time.time() - start

    print("Read pileups for %d loci in %f.2 seconds = %f locus / sec" % (
        len(loci), elapsed, len(loci) / elapsed))

if __name__ == '__main__':
    go(sys.argv[1:])
