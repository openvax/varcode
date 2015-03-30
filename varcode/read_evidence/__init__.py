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
'''
This subpackage provides functionality for collecting and filtering aligned
sequencing reads from a BAM file, determining the alleles they suggest at
a locus, and assesing the evidence for particular variants.

In this subpackage, the records stored in the BAM file are referred to as
"alignments," whereas the term "read" may be more familiar. We use the term
"alignment" for consistency with the SAM specification, and since an
individual read from the sequencer may generate any number of alignments in
the case of chimeric alignments and secondary alignments.
'''

from .util import alignment_key, read_key
from .pileup import Pileup
from .pileup_element import PileupElement
from .pileup_collection import PileupCollection

__all__ = [
    "PileupCollection",
    "Pileup",
    "PileupElement",
    "alignment_key",
    "read_key",
]
