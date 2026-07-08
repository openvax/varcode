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

"""Regression test for openvax/varcode#397.

``classify_from_protein_diff`` (used by the default ``protein_diff`` annotator)
derived a frameshift's ``shifted_sequence`` from ``trim_shared_flanking_strings``,
which trims a shared *suffix* as well as a shared prefix. A frameshift produces a
NOVEL C-terminus that runs to a new stop codon, so a residue it coincidentally
shares with the reference protein's own C-terminus was wrongly trimmed, dropping
the last 1-2 amino acids of the reported novel tail.

ATM p.F61fs (GRCh38 chr11:g.108227882delT, MANE Select ENST00000675843) is a clean
example: the novel ORF ends ``...R-K-K-Q-N-V`` (74 aa, terminated by TGA), and the
ATM reference protein also ends in ``V`` -- so the shared-suffix trim previously
yielded 73 aa ending ``...F-R-K-K-Q-N``.
"""

import pytest
from pyensembl import cached_release

from varcode import Variant
from varcode.effects import FrameShift

ensembl_grch38 = cached_release(115)

# This regression is pinned against the exact reported variant (ATM p.F61fs
# on MANE Select ENST00000675843), which only exists in newer Ensembl
# releases. Skip cleanly when release 115 isn't installed -- CI's data mirror
# (openvax/ensembl-data) tops out at GRCh38.95, so 115 is unavailable there.
# The same bug CLASS is exercised on an installed release (81) by
# tests/test_annotator_divergence_scenarios.py (CFTR p.L127fs, BRCA1 p.R71fs)
# and tests/test_protein_diff_parity.py, so CI coverage of #396/#397 does not
# depend on release 115.
try:
    ensembl_grch38.transcript_by_id("ENST00000675843")
except Exception as _exc:  # pyensembl raises if the GTF DB isn't downloaded
    pytest.skip(
        "Ensembl release 115 not installed (%s); ATM regression covered on "
        "release 81 elsewhere." % type(_exc).__name__,
        allow_module_level=True)


def _atm_f61fs_effect(annotator):
    # VCF-style anchored representation of chr11:g.108227882delT
    variant = Variant("11", 108227881, "GT", "G", ensembl_grch38)
    effects = variant.effects(annotator=annotator)
    on_mane = [e for e in effects if e.transcript_id == "ENST00000675843"]
    assert len(on_mane) == 1, "expected one effect on ENST00000675843"
    return on_mane[0]


def test_protein_diff_frameshift_keeps_cterminus_matching_reference_suffix():
    eff = _atm_f61fs_effect("protein_diff")
    assert isinstance(eff, FrameShift)
    assert eff.mutant_protein_sequence.endswith("RKKQNV")
    assert len(eff.mutant_protein_sequence) == 74
    assert eff.shifted_sequence.endswith("V")


def test_fast_and_protein_diff_agree_on_frameshift_tail():
    fast = _atm_f61fs_effect("fast")
    pdiff = _atm_f61fs_effect("protein_diff")
    assert fast.mutant_protein_sequence == pdiff.mutant_protein_sequence
