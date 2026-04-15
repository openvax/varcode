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

"""Parity harness for SequenceDiffEffectAnnotator vs LegacyEffectAnnotator
(openvax/varcode#271 stage 3d, #309).

**WIP — skipped in CI.** The harness skeleton is in place so the
validation corpus and assertion shape can be discussed before the
classifier lands. The target gate for 3d merge is:

    for variant in CORPUS:
        for transcript in variant.transcripts:
            legacy = LegacyEffectAnnotator().annotate_on_transcript(
                variant, transcript)
            sdiff = SequenceDiffEffectAnnotator().annotate_on_transcript(
                variant, transcript)
            assert type(legacy) is type(sdiff)
            assert legacy.short_description == sdiff.short_description

Any intentional divergence (a legacy bug that sequence_diff corrects)
must be called out explicitly in ``EXPECTED_DIFFS`` with a comment
linking the issue that documents it.
"""

import pytest


pytest.skip(
    "SequenceDiffEffectAnnotator is scaffolding only; parity harness "
    "activates with #309's classifier PR.",
    allow_module_level=True,
)


# Placeholder corpus — 3d will populate from:
#   * tests/data/somatic_hg19_14muts.vcf
#   * the BRCA1/CFTR fixtures used in test_effect_classes.py
#   * the splice-site corpus from test_splice_site_effects.py
#   * the MT cases from tests/test_mitochondrial.py
CORPUS = []

# Variants where sequence_diff and legacy intentionally differ. Each
# entry is a (variant_key, reason_issue_url) pair. Keep empty at merge
# time — any real disagreements must be triaged, not silently accepted.
EXPECTED_DIFFS = {}


def test_sequence_diff_matches_legacy_on_corpus():
    """Byte-for-byte parity on the validation corpus. Gate for 3d merge."""
    pass  # filled in by #309
