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

from varcode import Variant
from varcode.reference import (
    ensembl_reference_aliases, infer_genome, infer_reference_name, most_recent_assembly_name)
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



@pytest.mark.parametrize(["name", "release", "species", "was_ucsc"], [
    ("GRCh38:93", 93, "homo_sapiens", False),
    (" GRCh38 : 93 ", 93, "homo_sapiens", False),
    ("B37:75", 75, "homo_sapiens", False),
    ("hg19:75", 75, "homo_sapiens", True),
    ("GRCm38:95", 95, "mus_musculus", False),
    ("GRCh38.p13:93", 93, "homo_sapiens", False),
    # Ensembl's file names put the release after a dot; other separators work too.
    ("GRCh38.93", 93, "homo_sapiens", False),
    ("GRCh38_93", 93, "homo_sapiens", False),
    ("GRCh38-93", 93, "homo_sapiens", False),
    ("GRCh38 93", 93, "homo_sapiens", False),
    ("hg38.93", 93, "homo_sapiens", True),
    ("mm10.95", 95, "mus_musculus", True),
])
def test_reference_name_with_release_chooses_that_release(name, release, species, was_ucsc):
    genome, converted = infer_genome(name)
    assert (genome.release, genome.species.latin_name, converted) == (release, species, was_ucsc)


@pytest.mark.parametrize(["name", "message"], [
    ("GRCh38:75", "release 75 of homo_sapiens provides GRCh37, not GRCh38"),
    ("GRCh38.75", "release 75 of homo_sapiens provides GRCh37, not GRCh38"),
    ("GRCh37:93", "release 93 of homo_sapiens provides GRCh38, not GRCh37"),
    ("GRCh38:40", "No genome for homo_sapiens in Ensembl release 40"),
])
def test_release_that_does_not_provide_the_assembly_is_rejected(name, message):
    # Previously the release was dropped and the latest release used instead (#512).
    with pytest.raises(ValueError, match=message):
        infer_genome(name)


def test_variant_keeps_the_chosen_release_through_serialization():
    variant = Variant("7", 140753336, "A", "T", genome="GRCh38:93")
    assert variant.genome.release == 93 and variant.reference_name == "GRCh38"
    assert Variant.from_json(variant.to_json()).genome.release == 93


@pytest.mark.parametrize(["name", "reference_name"], [
    ("GRCh38.p13", "GRCh38"),
    ("grch38.d1.vd1", "GRCh38"),
    ("Felis_catus_9.0", "Felis_catus_9.0"),
])
def test_names_ending_in_digits_are_not_read_as_releases(name, reference_name):
    genome, _ = infer_genome(name)
    assert genome.reference_name == reference_name


def test_paths_are_not_read_as_releases():
    genome, _ = infer_genome("##reference=file:///var/lib/cwl/job367935311_index_001zdr/GRCh38.d1.vd1.fa")
    assert genome.reference_name == "GRCh38"
