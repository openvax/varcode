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

from nose.tools import eq_
import warnings
from varcode.reference import infer_reference_name, reference_alias_dict, _most_recent_assembly

## test cases are given as 
## expected response: list of inputs
reference_test_cases = {
    'NCBI36': ['ncbi36p2.fasta', 'b36.fasta', '##reference=file:///var/lib/cwl/ncbi36/homo_sapiens.d1.vd1.fa'],
    'GRCh38': ['grch38p2.fasta', 
                '##reference=file:///var/lib/cwl/job367935311_index_001zdr/GRCh38.d1.vd1.fa',
                '##reference=file:///var/lib/cwl/job367935311_index_001zdr/GRCh38.job36.d1.vd1.fa',
                ],
}

def test_most_recent_assembly():
    eq_(_most_recent_assembly(['ncbi36', 'grch38']), 'grch38')
    eq_(_most_recent_assembly(['ncbi36', 'grch38', '37mm']), 'grch38')
    eq_(_most_recent_assembly(['ncbi36']), 'ncbi36')
    eq_(_most_recent_assembly(['ncbi36', '35']), 'ncbi36')


def test_infer_reference_name_aliases():
    with warnings.catch_warnings(record=True) as w:
        for assembly_name in reference_alias_dict.keys():
            candidate_list = [assembly_name] + reference_alias_dict[assembly_name]
            for candidate in candidate_list:
                eq_(infer_reference_name(candidate), assembly_name)


def test_infer_reference_name_test_cases():
    with warnings.catch_warnings(record=True) as w:
        for assembly_name in reference_test_cases.keys():
            candidate_list = [assembly_name] + reference_test_cases[assembly_name]
            for candidate in candidate_list:
                eq_(infer_reference_name(candidate), assembly_name)

