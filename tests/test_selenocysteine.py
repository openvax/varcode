# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#         http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Regression tests for selenocysteine (Sec / U) handling (issue #42).

Selenoproteins recode an in-frame UGA -- normally a stop codon -- as
selenocysteine. Ensembl's reference protein for such transcripts is
full-length with a ``U`` at the recoded position, and varcode annotates
against that reference protein_sequence rather than re-translating the
coding sequence with the standard table. As a result an internal UGA is
NOT mistaken for a stop, so selenoprotein variants are no longer
mis-annotated as StopLoss / UTR effects (the original bug in #42).

GPX1 (glutathione peroxidase 1, ENSG00000233276) is the canonical test
case: a single selenocysteine at residue 49 (transcript ENST00000419783,
minus strand). Coordinates are GRCh38 and the transcript is pinned to the
installed release used across the suite.
"""

from pyensembl import cached_release

from varcode import Variant
from varcode.effects import Substitution, Deletion

from .common import expect_effect

ensembl_grch38 = cached_release(81)

# GPX1 selenoprotein transcript: minus strand, selenocysteine at residue 49.
GPX1_TRANSCRIPT_ID = "ENST00000419783"


def _gpx1_variant(contig_pos, ref, alt):
    return Variant("3", contig_pos, ref, alt, ensembl=ensembl_grch38)


def test_substitution_at_selenocysteine_reports_U_not_stop():
    # The selenocysteine codon itself (UGA) is reported as reference amino
    # acid ``U``, and a change there is an ordinary missense -- NOT a
    # StopLoss and NOT a UTR effect, which is what #42 originally observed.
    variant = _gpx1_variant(49358132, "T", "A")  # UGA -> UGU, Sec -> Cys
    expect_effect(
        variant,
        transcript_id=GPX1_TRANSCRIPT_ID,
        effect_class=Substitution,
        aa_ref="U",
        aa_alt="C",
        short_description="p.U49C")


def test_substitution_downstream_of_selenocysteine_keeps_numbering():
    # A missense well past the selenocysteine keeps correct residue
    # numbering (100, not truncated at the internal UGA), confirming the
    # reference protein reads through the Sec.
    variant = _gpx1_variant(49357702, "G", "A")
    expect_effect(
        variant,
        transcript_id=GPX1_TRANSCRIPT_ID,
        effect_class=Substitution,
        aa_ref="R",
        aa_alt="W",
        short_description="p.R100W")


def test_inframe_deletion_upstream_of_selenocysteine_reads_through():
    # An in-frame deletion upstream of the selenocysteine exercises the
    # re-translation path. The mutant protein must stay full length and
    # retain the ``U`` rather than truncating at the internal UGA.
    variant = _gpx1_variant(49358144, "CAC", "")  # deletes one upstream codon
    transcript = ensembl_grch38.transcript_by_id(GPX1_TRANSCRIPT_ID)
    effect = variant.effect_on_transcript(transcript)
    assert isinstance(effect, Deletion), \
        "Expected an in-frame Deletion, got %s" % type(effect).__name__
    mutant_protein = effect.mutant_protein_sequence
    reference_protein = transcript.protein_sequence
    assert "U" in mutant_protein, \
        "Selenocysteine dropped from mutant protein: %s" % mutant_protein
    assert len(mutant_protein) == len(reference_protein) - 1, \
        "Expected a single-residue deletion, ref=%d mutant=%d" % (
            len(reference_protein), len(mutant_protein))
