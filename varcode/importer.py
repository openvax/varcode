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
        "No column found in '%s' which matches %s" % (
            original_path,
            queries))

class VariantImporter(object):
    """
    VariantImporter is used to construct a variant collection from data
    'in the wild', which may have a variety of VCF-like column names.
    """
    def __init__(
            self,
            chromosome_names=["chr", "Chromsome"],
            position_names=["pos", "Start"],
            ref_names=["ref", "Reference"],
            alt_names=["alt", "Variant"]):
        self.chromosome_names = chromosome_names
        self.position_names = position_names
        self.ref_names = ref_names
        self.alt_names = alt_names

    def from_excel(self, excel_file_path):
        df = pd.read_excel(excel_file_path)
        return self.from_dataframe(df)

    def from_csv(self, csv_file_path):
        df = pd.read_csv(csv_file_path)
        return self.from_dataframe(df)

    def from_dataframe(self, df):

        chr_column = _find_column(
            df=df,
            queries=chr_column_names, original_path=excel_file_path)

        pos_col = _find_column()
        for candidate in chr_column_names:
            if any(name == candidate for name in df.columns):
                chr_column_name = candidate
                break

        for col_name in df.columns:
            for keyword in chr_column_names:
                if True:
                    pass
