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

# NCBI builds and hg releases aren't identical
# but the differences are all on chrM and unplaced contigs
reference_alias_dict = {
    "NCBI36": ["hg18", "B36", "NCBI36"],
    "GRCh37": ["hg19", "B37", "NCBI37"],
    "GRCh38": ["hg38", "B38", "NCBI38"],
    "GRCm37": ["mm9"],
    "GRCm38": ["mm10"],
}

def infer_reference_name(reference_name_or_path):
    """
    Given a string containing a reference name (such as a path to
    that reference's FASTA file), return its canonical name
    as used by Ensembl.
    """
    for assembly_name in sorted(reference_alias_dict.keys(), reverse=True):
        candiate_list = [assembly_name] + reference_alias_dict[assembly_name]
        for candidate in candiate_list:
            if candidate.lower() in reference_name_or_path.lower():
                return assembly_name
    raise ValueError(
        "Failed to infer genome assembly name for %s" % reference_name_or_path)
