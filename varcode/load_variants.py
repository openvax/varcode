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

from .maf import load_maf
from .vcf import load_vcf

def load_variants(filename):
    """
    Creates VariantCollection from name of either VCF or MAF file
    """
    assert isinstance(filename, str), \
        "Expected filename to be str, got %s : %s" % (
            filename, type(filename))

    if filename.endswith(".vcf"):
        return load_vcf(filename)
    elif filename.endswith(".maf"):
        return load_maf(filename)
    else:
        raise ValueError("Unrecognized file type: %s" % (filename,))
