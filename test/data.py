# Copyright (c) 2015. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#         http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Helper functions and shared datasets for tests
"""

import os
from varcode import Variant, VariantCollection, load_maf
import pandas as pd

def data_path(name):
    """
    Return the absolute path to a file in the varcode/test/data directory.
    The name specified should be relative to varcode/test/data.
    """
    return os.path.join(os.path.dirname(__file__), "data", name)

dbnsp_validation_df = pd.read_csv(data_path('dbnsfp_validation_set.csv'))
tcga_ov_variants = load_maf(data_path("tcga_ov.head.maf"))
ov_wustle_variants = load_maf(data_path("ov.wustle.subset5.maf"))

snp_rs4244285 = Variant(
    contig=10,
    start=94781859,
    ref="G",
    alt="A")
snp_rs1537415 = Variant(
    contig=9,
    start=135637876,
    ref="C",
    alt="G")
snp_rs3892097 = Variant(
    contig=22,
    start=42524947,
    ref="G",
    alt="A")

db_snp_variants = VariantCollection([
    snp_rs4244285,
    snp_rs1537415,
    snp_rs3892097,
])
