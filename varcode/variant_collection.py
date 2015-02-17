from __future__ import print_function, division, absolute_import

from . import maf
from .effect_ordering import effect_priority
from .reference_name import (
    infer_reference_name,
    ensembl_release_number_for_reference_name
)
from .variant import Variant
from .variant_annotator import VariantAnnotator

import vcf


def _load_vcf(filename):
    """
    Load reference name and Variant objects from the given VCF filename.

    Drop any entries whose FILTER field is not one of "." or "PASS".
    """
    vcf_reader = vcf.Reader(filename=filename)
    raw_reference_name = vcf_reader.metadata['reference']

    records = [
        Variant(
            x.CHROM, x.POS,
            x.REF, x.ALT[0].sequence,
            x.INFO)
        for x in vcf_reader
        if x.FILTER is None or x.FILTER == "PASS"
    ]
    return raw_reference_name, records

def _load_maf(filename):
    """
    Load reference name and Variant objects from MAF filename.
    """
    maf_df = maf.load_maf_dataframe(filename)

    if len(maf_df) == 0:
        raise ValueError("Empty MAF file %s" % filename)

    ncbi_builds = maf_df.NCBI_Build.unique()

    if len(ncbi_builds) == 0:
        raise ValueError("No NCBI builds for MAF file %s" % filename)
    elif len(ncbi_builds) > 1:
        raise ValueError(
            "Multiple NCBI builds (%s) for MAF file %s" % (ncbi_builds, filename))

    raw_reference_name = ncbi_builds[0]
    records = []

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

        records.append(Variant(contig, start_pos, ref, alt))

    return raw_reference_name, records

class VariantCollection(object):

    def __init__(self, filename, drop_duplicates=True):
        assert isinstance(filename, str), \
            "Expected filename to be str, got %s : %s" % (
                filename, type(filename))

        if filename.endswith(".vcf"):
            self.raw_reference_name, self.records = _load_vcf(filename)
        elif filename.endswith(".maf"):
            self.raw_reference_name, self.records = _load_maf(filename)
        else:
            raise ValueError("Unrecognized file type: %s" % (filename,))


        self.filename = filename
        self.reference_name = infer_reference_name(self.raw_reference_name)
        self.ensembl_release = ensembl_release_number_for_reference_name(
            self.reference_name)
        self.annot = VariantAnnotator(ensembl_release=self.ensembl_release)

        if drop_duplicates:
            filtered_records = []
            seen = set()
            for record in self.records:
                key = record.short_description()
                if key not in seen:
                    seen.add(key)
                    filtered_records.append(record)
            self.records = filtered_records

    def __len__(self):
        return len(self.records)

    def __iter__(self):
        return iter(self.records)

    def __str__(self):
        s = "VariantCollection(filename=%s, reference=%s)" % (
            self.filename, self.reference_name)
        for record in self.records:
            s += "\n\t%s" % record
        return s

    def __repr__(self):
        return str(self)

    def variant_effects(self, raise_on_error=True):
        """
        Determine the impact of each variant, return a list of
        VariantEffect objects.

        Parameters
        ----------

        raise_on_error : bool, optional.
            Raise exception if error is encountered while annotating
            transcripts, otherwise track errors in VariantEffect.errors
            dictionary (default=True).
        """
        return [
            self.annot.describe_variant(variant, raise_on_error=raise_on_error)
            for variant
            in self.records
        ]

    def effects_to_string(self):
        """
        Create a long string with all transcript effects for each mutation,
        grouped by gene (if a mutation affects multiple genes).
        """
        lines = []
        for variant_effect in self.variant_effects():
            transcript_effect_count = 0
            lines.append("\n%s" % variant_effect.variant)
            transcript_effect_lists = variant_effect.gene_transcript_effects
            for gene, transcript_effects in transcript_effect_lists.iteritems():
                lines.append("  Gene: %s" % gene)
                # print transcript effects with more significant impact
                # on top (e.g. FrameShift should go before NoncodingTranscript)
                for transcript_effect in sorted(
                        transcript_effects,
                        key=effect_priority,
                        reverse=True):
                    transcript_effect_count += 1
                    lines.append("  -- %s" % transcript_effect)
            # if we only printed one effect for this gene then
            # it's redundant to print it again as the highest priority effect
            if transcript_effect_count > 1:
                best = variant_effect.highest_priority_effect
                lines.append("  Highest Priority Effect: %s" % best)
        return "\n".join(lines)

    def print_effects(self):
        """
        Print all variants and their transcript effects (grouped by gene).
        """
        print(self.effects_to_string())
