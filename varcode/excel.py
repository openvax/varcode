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

import pandas as pd

from .variant_collection import VariantCollection

def _find_column(df, queries, original_path):
    """
    Return column name which matches one of the queries
    """
    column_names = df.columns
    normalized_queries = [
        query.strip().lower()
        for query in queries
    ]
    for col in column_names:
        normalized_col = col.lower().strip()
        for query in normalized_queries:
            if query == normalized_col:
                return df[col]
            elif normalized_col.split()[0] == query:
                return df[col]
    raise ValueError(
        "No column found in '%s' for chromosome, valid names would be %s" % (
            original_path,
            queries))

def load_excel(
        excel_file_path,
        chr_column_names=["chr", "Chromsome"],
        pos_column_names=["pos", "Start"],
        ref_column_names=["ref", "Reference"],
        alt_column_names=["alt", "Variant"]):
    """
    Loads a VariantCollection from an Excel file with some leniency
    in the possible names given to columns.

    Parameters
    ----------
    excel_file_path : str
        Path to an Excel file containing columns for chromosome,
        position, reference nucleotides, and variant nucleotides.

    chr_column_names : list of str
        Possible names for columns which indicate the chromosome
        names for variants.

    pos_column_names : list of str
        Possible names for columns which indicate the positions of
        variants.

    ref_column_names : list of str


    Returns a VariantCollection with as many elements as non-header
    rows in the Excel file.
    """
    df = pd.read_csv(excel_file_path)
    chr_column = _find_column(
        df=df, queries=chr_column_names, original_path=excel_file_path)

    pos_col
    for candidate in chr_column_names:
        if any(name == candidate for name in df.columns):
            chr_column_name = candidate
            break

    for col_name in df.columns:
        for keyword in chr_column_names:
            if True:
                pass
