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


import pandas
from pyensembl import EnsemblRelease
from typechecks import require_string

from .nucleotides import normalize_nucleotide_string
from .reference_name import (
    infer_reference_name,
    ensembl_release_number_for_reference_name
)
from .variant import Variant
from .variant_collection import VariantCollection

TCGA_PATIENT_ID_LENGTH = 12

MAF_COLUMN_NAMES = [
    'Hugo_Symbol',
    'Entrez_Gene_Id',
    'Center',
    'NCBI_Build',
    'Chromosome',
    'Start_Position',
    'End_Position',
    'Strand',
    'Variant_Classification',
    'Variant_Type',
    'Reference_Allele',
    'Tumor_Seq_Allele1',
    'Tumor_Seq_Allele2',
    'dbSNP_RS',
    'dbSNP_Val_Status',
    'Tumor_Sample_Barcode',
    'Matched_Norm_Sample_Barcode',
    'Match_Norm_Seq_Allele1',
    'Match_Norm_Seq_Allele2',
]


def load_maf_dataframe(path, nrows=None, verbose=False):
    """
    Load the guaranteed columns of a TCGA MAF file into a DataFrame
    """
    require_string(path, "Path to MAF")

    n_basic_columns = len(MAF_COLUMN_NAMES)

    df = pandas.read_csv(
        path,
        comment="#",
        sep="\t",
        low_memory=False,
        skip_blank_lines=True,
        header=0)

    if len(df.columns) < n_basic_columns:
        raise ValueError(
            "Too few columns in MAF file %s, expected %d but got  %d : %s" % (
                path, n_basic_columns, len(df.columns), df.columns))

    # check each pair of expected/actual column names to make sure they match
    for expected, actual in zip(MAF_COLUMN_NAMES, df.columns):
        if expected != actual:
            # MAFs in the wild have capitalization differences in their
            # column names, normalize them to always use the names above
            if expected.lower() == actual.lower():
                # using DataFrame.rename in Python 2.7.x doesn't seem to
                # work for some files, possibly because Pandas treats
                # unicode vs. str columns as different?
                df[expected] = df[actual]
                del df[actual]
            else:
                raise ValueError("Expected column %s but got %s" % (
                    expected, actual))
    return df

def load_maf(path):
    """
    Load reference name and Variant objects from MAF filename.
    """
    maf_df = load_maf_dataframe(path)

    if len(maf_df) == 0:
        raise ValueError("Empty MAF file %s" % path)

    ensembl_objects = {}
    variants = []
    for _, x in maf_df.iterrows():
        contig = x.Chromosome
        start_pos = x.Start_Position
        end_pos = x.End_Position
        ref = x.Reference_Allele

        # it's possible in a MAF file to have multiple Ensembl releases
        # mixed in a single MAF file (the genome assembly is
        # specified by the NCBI_Build column)
        ncbi_build = x.NCBI_Build
        if ncbi_build in ensembl_objects:
            ensembl = ensembl_objects[ncbi_build]
        else:
            if isinstance(ncbi_build, int):
                reference_name = "B%d" % ncbi_build
            else:
                reference_name = str(ncbi_build)

            reference_name = infer_reference_name(reference_name)
            ensembl_release = ensembl_release_number_for_reference_name(
                reference_name)
            ensembl = EnsemblRelease(release=ensembl_release)
            ensembl_objects[ncbi_build] = ensembl

        # have to try both Tumor_Seq_Allele1 and Tumor_Seq_Allele2
        # to figure out which is different from the reference allele
        if x.Tumor_Seq_Allele1 != ref:
            alt = x.Tumor_Seq_Allele1
        else:
            if x.Tumor_Seq_Allele2 == ref:
                raise ValueError(
                    "Both tumor alleles agree with reference %s: %s" % (
                        ref, x,))
            alt = x.Tumor_Seq_Allele2

        # nucleotide sequences get normalized in the Variant constructor
        # but also doing it here so we can check the lengths correctly
        ref = normalize_nucleotide_string(ref)
        alt = normalize_nucleotide_string(alt)
        if len(ref) == 0:
            # Since insertions happen *between* reference coordinates
            # it's not clear what their start/end should be. So,
            # by convention, MAF files make the start the reference position
            # before the inserted nucleotides and the end position comes after
            # the inserted nucleotides
            end_offset = 1
        else:
            end_offset = len(ref) - 1
        expected_end_pos = start_pos + end_offset
        if len(ref) > 0 and end_pos != expected_end_pos:
            # only check for correct ending since the meaning of start/end
            # for insertions is different thn for substitutions
            raise ValueError(
                "Expected variant %s:%s '%s' > '%s' to end at %d but got %d" % (
                    contig,
                    start_pos,
                    ref,
                    alt,
                    expected_end_pos,
                    end_pos))

        # keep metadata about the variant and its TCGA annotation
        info = {
            'Hugo_Symbol': x.Hugo_Symbol,
            'Center': x.Center,
            'Strand': x.Strand,
            'Variant_Classification': x.Variant_Classification,
            'Variant_Type': x.Variant_Type,
            'dbSNP_RS': x.dbSNP_RS,
            'dbSNP_Val_Status': x.dbSNP_Val_Status,
            'Tumor_Sample_Barcode': x.Tumor_Sample_Barcode,
            'Matched_Norm_Sample_Barcode': x.Matched_Norm_Sample_Barcode,
        }

        variant = Variant(
            contig,
            start_pos,
            ref,
            alt,
            ensembl=ensembl,
            info=info)

        variants.append(variant)

    return VariantCollection(
        variants=variants,
        path=path)
