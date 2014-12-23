#!/usr/bin/env python

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

import logging

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


