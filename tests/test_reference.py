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


import warnings

import pytest 

from varcode.reference import infer_reference_name, ensembl_reference_aliases, most_recent_assembly_name
from .common import eq_

## test cases are given as 
## expected response: list of inputs
reference_test_cases = {
    'NCBI36': [
        'ncbi36p2.fasta', 
        'b36.fasta', 
        '##reference=file:///var/lib/cwl/ncbi36/homo_sapiens.d1.vd1.fa'],
    'GRCh38': [
        'grch38p2.fasta', 
        '##reference=file:///var/lib/cwl/job367935311_index_001zdr/GRCh38.d1.vd1.fa',
        '##reference=file:///var/lib/cwl/job367935311_index_001zdr/GRCh38.job36.d1.vd1.fa',
    ],
}

def test_most_recent_assembly():
    eq_(most_recent_assembly_name(['ncbi36', 'grch38']), 'grch38')
    eq_(most_recent_assembly_name(['ncbi36', 'grch38', '37mm']), 'grch38')
    eq_(most_recent_assembly_name(['ncbi36']), 'ncbi36')
    eq_(most_recent_assembly_name(['ncbi36', '35']), 'ncbi36')
def generate_reference_name_aliases():
    with warnings.catch_warnings(record=True) as w:
        for assembly_name, aliases in ensembl_reference_aliases.items():
            candidate_list = [assembly_name] + list(aliases)
            for candidate in candidate_list:
                yield (                
                    candidate,
                    assembly_name
                )

@pytest.mark.parametrize(['candidate', 'assembly_name'], generate_reference_name_aliases())
def test_infer_reference_name_aliases(candidate, assembly_name):
    eq_(infer_reference_name(candidate), assembly_name)
    
def generate_reference_name_fasta_filenames():
    with warnings.catch_warnings(record=True):
        for assembly_name, aliases in reference_test_cases.items():
            candidate_list = [assembly_name] + list(aliases)
            for candidate in candidate_list:
                yield (
                    candidate,
                    assembly_name
                )

@pytest.mark.parametrize(['candidate', 'assembly_name'], generate_reference_name_fasta_filenames())
def test_reference_name_fasta_filenames(candidate, assembly_name):
    eq_(infer_reference_name(candidate), assembly_name)

