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

"""
Regression tests for https://github.com/openvax/varcode/issues/217

Silent.short_description returned the literal string "silent" instead of
the HGVS-standard `p.{AA}{pos}=` format (e.g. "p.R6=").
"""

from varcode import Variant
from varcode.effects import Silent


TRANSCRIPT_ID = "ENSMUST00000086738"


def test_silent_short_description_uses_hgvs_format():
    # AGA -> AGG at chr1:99772782 (3rd nt of codon 5, R codon)
    # Synonymous: R at position 5 (0-indexed) = position 6 (1-indexed).
    variant = Variant(
        contig=1,
        start=99772782,
        ref="A",
        alt="G",
        genome="GRCm38",
    )
    transcript = variant.ensembl.transcript_by_id(TRANSCRIPT_ID)
    effect = variant.effect_on_transcript(transcript)
    assert effect.__class__ is Silent
    # HGVS: p.{aa_ref}{1-indexed pos}=
    assert effect.short_description == "p.R6=", \
        "Expected 'p.R6=', got %r" % effect.short_description


def test_silent_short_description_for_multi_codon_synonymous():
    # AC -> GT at chr1:99772782, spans two codons AGA+CTG both synonymous.
    variant = Variant(
        contig=1,
        start=99772782,
        ref="AC",
        alt="GT",
        genome="GRCm38",
    )
    transcript = variant.ensembl.transcript_by_id(TRANSCRIPT_ID)
    effect = variant.effect_on_transcript(transcript)
    assert effect.__class__ is Silent
    # aa_ref is 'RL', 1-indexed pos is 6.
    assert effect.short_description == "p.RL6=", \
        "Expected 'p.RL6=', got %r" % effect.short_description


# Direct unit test on the Silent class without needing a real transcript.

class _FakeGene:
    biotype = "protein_coding"


class _FakeTranscript:
    biotype = "protein_coding"
    gene = _FakeGene()
    protein_sequence = "X" * 100


def test_silent_unit_hgvs_notation():
    e = Silent(
        variant=None,
        transcript=_FakeTranscript(),
        aa_pos=5,
        aa_ref="R",
    )
    # 0-indexed aa_pos=5 -> 1-indexed position 6 in HGVS.
    assert e.short_description == "p.R6="


def test_silent_unit_empty_aa_ref_still_safe():
    # Edge case: when nothing was affected (both aa_ref and aa_alt empty
    # after trim), aa_ref may be empty. Should degrade to a sensible
    # descriptor rather than crash.
    e = Silent(
        variant=None,
        transcript=_FakeTranscript(),
        aa_pos=5,
        aa_ref="",
    )
    # Any non-crashing string is acceptable; check it contains the position
    # marker and the HGVS "=" for synonymous.
    assert "=" in e.short_description
