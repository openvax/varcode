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

def infer_reference_name(path):
    # NCBI builds and hg releases aren't identical
    # but the differences are all on chrM and unplaced contigs
    candidates = {
        'NCBI36': ['hg18', 'B36', 'GRCh36', 'NCBI36'],
        'GRCh37': ['hg19', 'B37', 'GRCh37', 'NCBI37'],
        'GRCh38': ['hg38', 'B38', 'GRCh38', 'NCBI38'],
    }

    for name in sorted(candidates.keys(), reverse=True):
        aliases = candidates[name]
        for alias in aliases:
            if alias in path:
                return name
            if alias.lower() in path:
                return name

    raise ValueError(
        "Failed to infer human genome assembly name for %s" % path)

def ensembl_release_number_for_reference_name(name):
    if name == "NCBI36":
        return 54
    elif name == "GRCh37":
        return 75
    else:
        assert name == "GRCh38", "Unrecognized reference name: %s" % name
        return 78
