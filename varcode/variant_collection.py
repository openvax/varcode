from reference_name import infer_reference_name
from variant import Variant

import vcf

class VariantCollection(object):

    def __init__(self, filename):
        assert isinstance(filename, (str,unicode)), \
            "Expected filename to be str, got %s : %s" % (
                filename, type(filename))

        if filename.endswith(".vcf"):
            vcf_reader = vcf.Reader(filename=filename)
            self.reference = infer_reference_name(
                vcf_reader.metadata['reference'])
            self.records = [
                Variant(x.CHROM, x.POS, x.REF, x.ALT[0].sequence)
                for x in vcf_reader
            ]
        elif filename.endswith(".maf"):
            maf_df = maf.load_maf_dataframe(filename)

            if len(maf_df) == 0:
                raise ValueError("Empty MAF file %s" % filename)

            ncbi_builds = maf_df.NCBI_Build.unique()

            if len(ncbi_builds) == 0:
                raise ValueError("No NCBI builds for MAF file %s" % filename)
            elif len(ncbi_builds) > 1:
                raise ValueError(
                    "Multiple NCBI builds (%s) for MAF file %s" % (ncbi_builds, filename))

            ncbi_build = ncbi_builds[0]
            self.reference = infer_reference_name(ncbi_build)
            self.records = []
            for _, x in maf_df.iterrows():
                start_pos = x.Start_Position
                end_pos = x.End_Position
                contig = x.Chromosome
                ref = normalize_nucleotide_string(x.Reference_Allele)

                if x.Tumor_Seq_Allele1 != ref:
                    alt = x.Tumor_Seq_Allele1
                else:
                    assert x.Tumor_Seq_Allele2 != ref, \
                        "Both tumor alleles agree with reference: %s" % (x,)
                    alt = x.Tumor_Seq_Allele2

                alt = normalize_nucleotide_string(alt)

                record = Variant(contig, start_pos, ref, alt)
                self.records.append(record)

        else:
            raise ValueErrr("Unrecognized file type: %s" % filename)

    def __iter__(self):
        return iter(self.records)