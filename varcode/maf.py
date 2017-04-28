# Copyright (c) 2016-2017. Mount Sinai School of Medicine
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

import pandas
from typechecks import require_string
from pandas import isnull

from .reference import infer_genome
from .variant import Variant, variant_ascending_position_sort_key
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


def load_maf_dataframe(path, nrows=None, raise_on_error=True, encoding=None):
    """
    Load the guaranteed columns of a TCGA MAF file into a DataFrame

    Parameters
    ----------
    path : str
        Path to MAF file

    nrows : int
        Optional limit to number of rows loaded

    raise_on_error : bool
        Raise an exception upon encountering an error or log an error

    encoding : str, optional
        Encoding to use for UTF when reading MAF file.
    """
    require_string(path, "Path to MAF")

    n_basic_columns = len(MAF_COLUMN_NAMES)

    # pylint: disable=no-member
    # pylint gets confused by read_csv
    df = pandas.read_csv(
        path,
        comment="#",
        sep="\t",
        low_memory=False,
        skip_blank_lines=True,
        header=0,
        encoding=encoding)

    if len(df.columns) < n_basic_columns:
        error_message = (
            "Too few columns in MAF file %s, expected %d but got  %d : %s" % (
                path, n_basic_columns, len(df.columns), df.columns))
        if raise_on_error:
            raise ValueError(error_message)
        else:
            logging.warn(error_message)

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
                error_message = (
                    "Expected column %s but got %s" % (expected, actual))
                if raise_on_error:
                    raise ValueError(error_message)
                else:
                    logging.warn(error_message)

    return df

def load_maf(
        path,
        optional_cols=[],
        sort_key=variant_ascending_position_sort_key,
        distinct=True,
        raise_on_error=True,
        encoding=None):
    """
    Load reference name and Variant objects from MAF filename.

    Parameters
    ----------

    path : str
        Path to MAF (*.maf).

    optional_cols : list, optional
        A list of MAF columns to include as metadata if they are present in the MAF.
        Does not result in an error if those columns are not present.

    sort_key : fn
        Function which maps each element to a sorting criterion.
        Set to None to not to sort the variants.

    distinct : bool
        Don't keep repeated variants

    raise_on_error : bool
        Raise an exception upon encountering an error or just log a warning.

    encoding : str, optional
        Encoding to use for UTF when reading MAF file.
    """
    # pylint: disable=no-member
    # pylint gets confused by read_csv inside load_maf_dataframe
    maf_df = load_maf_dataframe(path, raise_on_error=raise_on_error, encoding=encoding)

    if len(maf_df) == 0 and raise_on_error:
        raise ValueError("Empty MAF file %s" % path)

    ensembl_objects = {}
    variants = []
    metadata = {}
    for _, x in maf_df.iterrows():
        contig = x.Chromosome
        if isnull(contig):
            error_message = "Invalid contig name: %s" % (contig,)
            if raise_on_error:
                raise ValueError(error_message)
            else:
                logging.warn(error_message)
                continue

        start_pos = x.Start_Position
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
            ensembl = infer_genome(reference_name)
            ensembl_objects[ncbi_build] = ensembl

        # have to try both Tumor_Seq_Allele1 and Tumor_Seq_Allele2
        # to figure out which is different from the reference allele
        if x.Tumor_Seq_Allele1 != ref:
            alt = x.Tumor_Seq_Allele1
        else:
            if x.Tumor_Seq_Allele2 == ref:
                error_message = (
                    "Both tumor alleles agree with reference %s: %s" % (
                        ref, x,))
                if raise_on_error:
                    raise ValueError(error_message)
                else:
                    logging.warn(error_message)
                    continue
            alt = x.Tumor_Seq_Allele2

        variant = Variant(
            contig,
            start_pos,
            str(ref),
            str(alt),
            ensembl=ensembl)

        # keep metadata about the variant and its TCGA annotation
        metadata[variant] = {
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
        for optional_col in optional_cols:
            if optional_col in x:
                metadata[variant][optional_col] = x[optional_col]

        variants.append(variant)

    return VariantCollection(
        variants=variants,
        source_to_metadata_dict={path: metadata},
        sort_key=sort_key,
        distinct=distinct)
