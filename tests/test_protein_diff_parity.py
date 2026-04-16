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

"""Parity harness for ProteinDiffEffectAnnotator vs LegacyEffectAnnotator
(openvax/varcode#271 stage 3d, #309).

**WIP — skipped in CI.** The harness skeleton is in place so the
validation corpus and assertion shape can be discussed before the
classifier lands. The target gate for 3d merge is:

    for variant in CORPUS:
        for transcript in variant.transcripts:
            legacy = LegacyEffectAnnotator().annotate_on_transcript(
                variant, transcript)
            sdiff = ProteinDiffEffectAnnotator().annotate_on_transcript(
                variant, transcript)
            assert type(legacy) is type(sdiff)
            assert legacy.short_description == sdiff.short_description

Rollout policy
--------------

During the **initial corpus build** (before 3d merges), parity
deltas are logged as INFO via ``logging.info`` rather than asserted.
The goal at that stage is to **review the parity contract** — each
delta is either:

  (a) a protein-diff bug that needs fixing before merge, or
  (b) a legacy bug that protein-diff corrects (pinned in
      ``EXPECTED_DIFFS`` with an issue link).

Only after every delta has been triaged does the harness flip to
``assert``. At that point the parity contract is committed: any
new delta is a regression.

Each entry in ``EXPECTED_DIFFS`` is a
``(variant_key, issue_url, reason)`` tuple. Empty at merge
time — any real disagreement must be resolved before the harness
becomes a gate.

Corpus composition
------------------

~500 variants pulled from existing test fixtures plus targeted
edge cases from the PR #310 review:

  Existing fixtures:
    * tests/data/somatic_hg19_14muts.vcf (14 somatic SNVs)
    * BRCA1/CFTR coding SNVs / MNVs from test_effect_classes.py
    * Splice-site corpus from test_splice_site_effects.py
    * MT cases from test_mitochondrial.py

  Edge cases to pin (from vaxrank feedback):
    * AlternateStartCodon + downstream missense in the same
      variant — verifies the "only taken when the only change
      is at codon 0" rule.
    * MT PrematureStop and MT StopLoss — verifies the codon
      table routes through apply_variant_to_transcript correctly
      end-to-end.
    * StopLoss with short downstream stop (the extension hits a
      stop a few residues into the 3'UTR) — verifies the bounded
      readthrough case, not just the unbounded one.
    * Pure insertion of 3 residues that match flanking residues —
      verifies trim_shared_flanking_strings doesn't over-trim.
      Legacy regressed on this pre-#201.
    * Phased MNP where two SNVs hit the same codon — verifies
      the classifier when (ref_delta, alt_delta) are both length
      1 but represent a compound substitution.
    * Variant at last exonic base that also changes an amino
      acid — verifies the splice-adjacency gate's dual-dispatch
      behaviour (see protein_diff._variant_is_splice_adjacent
      docstring for the contract).
    * A variant whose effect_type changed across varcode major
      versions (legacy ComplexSubstitution ↔ Substitution) —
      pinned explicitly so parity deltas are intentional.
"""

import pytest


# The parity harness is expensive and WIP. Mark it so it can be
# selected explicitly (`pytest -m parity`) once populated, without
# slowing the default suite run.
pytestmark = pytest.mark.parity

pytest.skip(
    "ProteinDiffEffectAnnotator is scaffolding only; parity harness "
    "activates with #309's classifier PR.",
    allow_module_level=True,
)


# Placeholder corpus — 3d will populate from the sources documented
# in the module docstring.
CORPUS = []

# Variants where protein_diff and legacy intentionally differ. Each
# entry: (variant_key, reason_issue_url, reason_description).
# Empty at merge time — any real disagreements must be triaged, not
# silently accepted.
EXPECTED_DIFFS = {}


def test_protein_diff_matches_legacy_on_corpus():
    """Byte-for-byte parity on the validation corpus. Gate for 3d merge."""
    pass  # filled in by #309
