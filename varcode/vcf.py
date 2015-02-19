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

# required so that 'import vcf' gets the global PyVCF package,
# rather than our local vcf module
from __future__ import absolute_import
from .variant import Variant
from .variant_collection import VariantCollection

# PyVCF
import vcf

def load_vcf(filename, reference_path_field='reference'):
    """
    Load reference name and Variant objects from the given VCF filename.
    Drop any entries whose FILTER field is not one of "." or "PASS".

    Parameters
    ----------

    filename : str

    reference_path_field : str, optional
        Name of metadata field which contains path to reference FASTA
        file (default = 'reference')
    """
    vcf_reader = vcf.Reader(filename=filename)
    reference_name = vcf_reader.metadata[reference_path_field]

    variants = [
        Variant(
            x.CHROM, x.POS,
            x.REF, x.ALT[0].sequence,
            x.INFO)
        for x in vcf_reader
        if not x.FILTER or x.FILTER == "PASS"
    ]
    return VariantCollection(
        variants=variants,
        original_filename=filename,
        reference_path=reference_name)
