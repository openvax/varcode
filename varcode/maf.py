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

from __future__ import annotations

import logging
import os
from typing import TYPE_CHECKING

import pandas
from pandas import isnull

from .reference import infer_genome
from .variant import Variant, variant_ascending_position_sort_key
from .variant_collection import VariantCollection

if TYPE_CHECKING:
    from pyensembl import Genome as EnsemblGenome
    from .genome import Genome

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


# MAF columns kept as each variant's metadata.
_METADATA_COLUMNS = (
    'Hugo_Symbol', 'Center', 'Strand', 'Variant_Classification', 'Variant_Type',
    'dbSNP_RS', 'dbSNP_Val_Status', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode')


def _cell(row, column):
    """A MAF cell, with an empty cell (which pandas reads as NaN) as None.

    NaN can't be written as JSON, so a collection holding it in metadata
    couldn't be saved with ``to_json`` (#704).
    """
    value = row[column]
    return None if isnull(value) else value


def load_maf_dataframe(
        path: str | os.PathLike[str], nrows: int | None = None,
        raise_on_error: bool = True, encoding: str | None = None) -> pandas.DataFrame:
    """
    Load the guaranteed columns of a TCGA MAF file into a DataFrame

    Parameters
    ----------
    path : str or os.PathLike[str]
        Path to MAF file

    nrows : int
        Optional limit to number of rows loaded

    raise_on_error : bool
        Raise ``ValueError`` when the 19 required columns
        (``MAF_COLUMN_NAMES``, matched case-insensitively) are out of order;
        otherwise log a warning. A missing required column always raises
        ``ValueError``, since the file can't be read as a MAF without it.

    encoding : str, optional
        Encoding to use for UTF when reading MAF file.
    """
    path = os.fspath(path)
    if not isinstance(path, str):
        raise TypeError("Path to MAF must resolve to str, got %s" % type(path).__name__)

    df = pandas.read_csv(
        path,
        comment="#",
        sep="\t",
        low_memory=False,
        skip_blank_lines=True,
        header=0,
        nrows=nrows,
        encoding=encoding)

    found = {str(column).lower() for column in df.columns}
    missing = [name for name in MAF_COLUMN_NAMES if name.lower() not in found]
    if missing:
        raise ValueError(
            "MAF file %s is missing %d of the %d required columns: %s (found: %s)" % (
                path, len(missing), len(MAF_COLUMN_NAMES), ", ".join(missing),
                ", ".join(map(str, df.columns))))

    # MAFs in the wild have capitalization differences in their column
    # names; give a required one the spec's spelling wherever it is, unless
    # that spelling is already there (TCGA MAFs have both Strand and strand).
    spelled: dict[str, str] = {}
    for column in df.columns:
        spelled.setdefault(str(column).lower(), column)
    df = df.rename(columns={
        spelled[name.lower()]: name for name in MAF_COLUMN_NAMES if name not in df.columns})
    for expected, actual in zip(MAF_COLUMN_NAMES, df.columns):
        if expected != actual:
            error_message = "Expected column %s but got %s" % (expected, actual)
            if raise_on_error:
                raise ValueError(error_message)
            else:
                logging.warning(error_message)

    return df

def _row_variant(x, genomes):
    """The Variant in MAF row ``x``; ValueError if the row can't be read as one.

    ``genomes`` caches the genome of each ``NCBI_Build`` value, since one MAF
    can mix builds.
    """
    empty = [column for column in ("NCBI_Build", "Chromosome", "Start_Position", "Reference_Allele")
             if _cell(x, column) is None]
    if empty:
        raise ValueError("MAF row has no %s: %s" % (", ".join(empty), x))
    ref = x.Reference_Allele

    ncbi_build = x.NCBI_Build
    if ncbi_build not in genomes:
        reference_name = "B%d" % ncbi_build if isinstance(ncbi_build, int) else str(ncbi_build)
        genomes[ncbi_build], _ = infer_genome(reference_name)

    # The alt is the first tumor allele that differs from the reference.
    # An empty cell (NaN) is no allele, not the alt "nan".
    alts = [allele for allele in (_cell(x, "Tumor_Seq_Allele1"), _cell(x, "Tumor_Seq_Allele2"))
            if allele is not None and allele != ref]
    if not alts:
        raise ValueError("No tumor allele differs from reference %s: %s" % (ref, x))

    # GDC GRCh38 MAFs name contigs "chr7" and "chrM"; Ensembl annotation
    # names them "7" and "MT", as load_vcf converts them. Variant leaves
    # hg19's chrM alone, since it isn't Ensembl's MT (#707).
    return Variant(
        x.Chromosome, x.Start_Position, str(ref), str(alts[0]), genomes[ncbi_build],
        convert_ucsc_contig_names=True)


def load_maf(
        path: str | os.PathLike[str],
        optional_cols=None,
        sort_key=variant_ascending_position_sort_key,
        distinct=True,
        raise_on_error=True,
        encoding=None,
        nrows=None) -> VariantCollection:
    """
    Load reference name and Variant objects from MAF filename.

    The ``NCBI_Build`` column selects each variant's reference. UCSC contig
    names, such as the ``chr7`` and ``chrM`` in GDC GRCh38 MAFs, become
    Ensembl's ``7`` and ``MT``. In a GRCh37 MAF, ``chrM`` and ``M`` keep their
    names and fail annotation: hg19's chrM is not the rCRS that Ensembl names
    MT, and its positions are offset from it. Empty cells become None in each
    variant's metadata.

    Parameters
    ----------

    path : str or os.PathLike[str]
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
        Raise ``ValueError`` for a row that can't be read as a variant (an
        empty build, contig, position or reference allele, no tumor allele
        that differs from the reference, or an invalid allele or build), or
        log a warning and skip it. A missing required column always raises.

    encoding : str, optional
        Encoding to use for UTF when reading MAF file.

    nrows : int, optional
        Limit to number of rows loaded
    """
    path = os.fspath(path)
    if optional_cols is None:
        optional_cols = ()
    maf_df = load_maf_dataframe(
        path,
        nrows=nrows,
        raise_on_error=raise_on_error,
        encoding=encoding)

    if len(maf_df) == 0 and raise_on_error:
        raise ValueError("Empty MAF file %s" % path)

    genomes: dict[str | int, EnsemblGenome | Genome] = {}
    variants = []
    metadata = {}
    for _, x in maf_df.iterrows():
        try:
            variant = _row_variant(x, genomes)
        except ValueError as error:
            if raise_on_error:
                raise
            logging.warning("Skipping MAF row: %s", error)
            continue

        # keep metadata about the variant and its TCGA annotation
        metadata[variant] = {
            column: _cell(x, column) for column in _METADATA_COLUMNS}
        for optional_col in optional_cols:
            if optional_col in x:
                metadata[variant][optional_col] = _cell(x, optional_col)

        variants.append(variant)

    return VariantCollection(
        variants=variants,
        source_to_metadata_dict={path: metadata},
        sort_key=sort_key,
        distinct=distinct)
