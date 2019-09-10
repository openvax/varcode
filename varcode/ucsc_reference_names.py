# Copyright (c) 2019. Mount Sinai School of Medicine
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

# NCBI builds and hg releases aren't identical
# but the differences are all on chrM and unplaced contigs
# Mapping between names copied from:
# https://genome.ucsc.edu/FAQ/FAQreleases.html#release1


ucsc_to_ensembl_reference_names = {
    # mouse
    "mm9": "GRCm37",
    "mm10": "GRCm38",

    # human
    "hg18": "NCBI36",
    "hg19": "GRCh37",
    "hg38": "GRCh38",

    # cat
    "felCat9": "Felis_catus_9.0",
    "felCat8": "Felis_catus_8.0",

    # dog
    "canFam3": "CanFam3.1",
    "canFam2": "CanFam2.0",
    "canFam1": "CanFam1.0",

    # rat
    "rn6": "Rnor_6.0",
    "rn5": "Rnor_5.0",

    # chicken
    "galGal6":	"Gallus_gallus_6.0",
    "galGal5":	"Gallus_gallus_5.0",
    "galGal4": "Gallus_gallus_4.0",

    # TODO: add Zebrafish to PyEnsembl
    #   "danRer11": "GRCz11",
    #   "danRer10": "GRCz10",

    # TODO: add cow to PyEnsembl
    #    "bosTau9": "ARS-UCD1.2",
    #    "bosTau8": "UMD3.1",

    # TODO: add pig to PyEnsembl
    #   "susScr11": "Sscrofa11.1",
    #   "susScr3": "Sscrofa10.2",
    #   "susScr2": "Sscrofa9.2",
}

ucsc_reference_names = set(ucsc_to_ensembl_reference_names.keys())

ensembl_to_ucsc_reference_names = {v: k for (k, v) in ucsc_reference_names}