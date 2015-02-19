# Copyright (c) 2014. Mount Sinai School of Medicine
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
import logging

from .nucleotides import normalize_nucleotide_string
from .variant import Variant
from .variant_collection import VariantCollection

import pandas

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


def load_maf_dataframe(filename, nrows=None, verbose=False):
    """
    Load the guaranteed columns of a TCGA MAF file into a DataFrame
    """
    logging.info("Opening %s" % filename)

    # skip comments and optional header
    with open(filename) as f:
        lines_to_skip = 0
        for line in f:
            if line.startswith("#") or line.startswith("Hugo_Symbol"):
                lines_to_skip += 1
            else:
                break
    return pandas.read_csv(
        filename,
        skiprows=lines_to_skip,
        sep="\t",
        usecols=range(len(MAF_COLUMN_NAMES)),
        low_memory=False,
        names=MAF_COLUMN_NAMES)

def load_maf(filename):
    """
    Load reference name and Variant objects from MAF filename.
    """
    maf_df = load_maf_dataframe(filename)

    if len(maf_df) == 0:
        raise ValueError("Empty MAF file %s" % filename)

    ncbi_builds = maf_df.NCBI_Build.unique()

    if len(ncbi_builds) == 0:
        raise ValueError("No NCBI builds for MAF file %s" % filename)
    elif len(ncbi_builds) > 1:
        raise ValueError(
            "Multiple NCBI builds (%s) for MAF file %s" % (ncbi_builds, filename))

    reference_name = ncbi_builds[0]
    variants = []

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

        variants.append(Variant(contig, start_pos, ref, alt))

    return VariantCollection(
        variants=variants,
        original_filename=filename,
        reference_name=reference_name)
