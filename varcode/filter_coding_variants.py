import argparse

import vcf

parser = argparse.ArgumentParser()

parser.add_argument("input",
    help="VCF input file")

parser.add_argument("--output-vcf",
    help="VCF output file (keeping only header and coding variants)")

parser.add_argument("--output-csv",
    help="CSV output with chr/pos/ref/alt from VCF and annotation")

if __name__ == "__main__":
    args = parser.parse_args()
    with open(args.input, 'r') as f:
        vcf_reader = vcf.Reader(f)
        for record in vcf_reader:
            print record
            for sample in record.samples:
                print "  ", sample