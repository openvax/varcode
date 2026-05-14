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

"""Joint cis-variant effect prediction (#269).

HaplotypeEffect + apply_variants_to_transcript: when two or more
variants are in cis according to a phase resolver, varcode builds a
single joint :class:`MutantTranscript` rather than computing
per-variant effects independently. Tests exercise both the
:class:`ReadPhaseResolver` and :class:`VCFPhaseResolver` backends.
"""

import os
import tempfile

import pytest

from varcode import (
    MutantTranscript,
    ReadPhaseResolver,
    VCFPhaseResolver,
    Variant,
    VariantCollection,
    apply_variants_to_transcript,
    load_vcf,
)
from varcode.effects.effect_classes import HaplotypeEffect


# -----------------------------------------------------------------
# apply_variants_to_transcript: multi-edit builder
# -----------------------------------------------------------------


def test_apply_variants_to_transcript_joint_cdna_differs_from_per_variant():
    """Apply two distant SNVs on the same transcript. Joint cDNA has
    both base changes; per-variant cDNAs have only one each."""
    from pyensembl import cached_release
    g = cached_release(81)
    transcript = g.transcript_by_id("ENST00000003084")  # CFTR, forward
    v1 = Variant("7", 117531100, "T", "A", g)
    v2 = Variant("7", 117531114, "G", "T", g)
    joint = apply_variants_to_transcript([v1, v2], transcript)
    assert isinstance(joint, MutantTranscript)
    assert len(joint.edits) == 2
    ref = str(transcript.sequence)
    diffs = [i for i in range(len(ref)) if ref[i] != joint.cdna_sequence[i]]
    # Two base differences, both present in the joint cDNA.
    assert len(diffs) == 2


def test_apply_variants_detects_overlapping_edits():
    """Two variants that claim the same base → build returns None so
    the caller can fall back to per-variant effects."""
    from pyensembl import cached_release
    g = cached_release(81)
    transcript = g.transcript_by_id("ENST00000003084")
    # Deliberately overlapping: v2's deletion range covers v1's
    # single-base edit.
    v1 = Variant("7", 117531100, "T", "A", g)
    v2 = Variant("7", 117531100, "TTG", "", g, allow_extended_nucleotides=True)
    # apply_variants_to_transcript should refuse this combination.
    result = apply_variants_to_transcript([v1, v2], transcript)
    assert result is None


def test_apply_variants_order_independent():
    """Order of variants in the input doesn't change the output cDNA."""
    from pyensembl import cached_release
    g = cached_release(81)
    transcript = g.transcript_by_id("ENST00000003084")
    v1 = Variant("7", 117531100, "T", "A", g)
    v2 = Variant("7", 117531114, "G", "T", g)
    a = apply_variants_to_transcript([v1, v2], transcript)
    b = apply_variants_to_transcript([v2, v1], transcript)
    assert a.cdna_sequence == b.cdna_sequence


# -----------------------------------------------------------------
# Read-phasing-resolver-driven joint effect
# -----------------------------------------------------------------


class _StubReadPhasingSource:
    def __init__(self, phasing=None, mutant_transcripts=None):
        self._phasing = dict(phasing) if phasing else {}
        self._mts = dict(mutant_transcripts) if mutant_transcripts else {}

    def has_evidence(self, variant):
        return variant in self._phasing

    def partners_in_cis(self, variant):
        return self._phasing.get(variant, ())

    def mutant_transcript(self, variant, transcript):
        return self._mts.get((variant, transcript.id))


def test_read_resolver_produces_haplotype_effect():
    """Two cis variants via read-phasing stub → collection contains a
    HaplotypeEffect with both variants and the resolver-provided
    MutantTranscript (not DNA inference)."""
    from pyensembl import cached_release
    g = cached_release(81)
    transcript = g.transcript_by_id("ENST00000003084")
    v1 = Variant("7", 117531100, "T", "A", g)
    v2 = Variant("7", 117531114, "G", "T", g)
    stub_mt = MutantTranscript(
        reference_transcript=transcript,
        cdna_sequence="ACGT" * 50,
        mutant_protein_sequence="MOBSERVED" + "A" * 100,
        annotator_name="upstream_rna_tool",
    )
    source = _StubReadPhasingSource(
        phasing={v1: (v1, v2), v2: (v1, v2)},
        mutant_transcripts={
            (v1, transcript.id): stub_mt,
            (v2, transcript.id): stub_mt,
        },
    )
    resolver = ReadPhaseResolver(source)
    collection = VariantCollection([v1, v2])
    effects = collection.effects(phase_resolver=resolver)
    haplotype_effects = [e for e in effects if isinstance(e, HaplotypeEffect)]
    assert len(haplotype_effects) == 1
    he = haplotype_effects[0]
    assert set(he.variants) == {v1, v2}
    assert he.transcript is transcript
    assert he.phase_source == "read_phasing"
    # Prefers the resolver-provided contig over DNA inference.
    assert he.mutant_transcript is stub_mt
    assert he.mutant_protein_sequence.startswith("MOBSERVED")


def test_read_resolver_haplotype_and_per_variant_coexist():
    """Per-variant effects stay on the collection alongside the
    HaplotypeEffect — additive, not replacement."""
    from pyensembl import cached_release
    g = cached_release(81)
    transcript = g.transcript_by_id("ENST00000003084")
    v1 = Variant("7", 117531100, "T", "A", g)
    v2 = Variant("7", 117531114, "G", "T", g)
    stub_mt = MutantTranscript(
        reference_transcript=transcript,
        cdna_sequence="ACGT" * 50,
        mutant_protein_sequence="M" + "A" * 99,
        annotator_name="upstream_rna_tool",
    )
    source = _StubReadPhasingSource(
        phasing={v1: (v1, v2)},
        mutant_transcripts={(v1, transcript.id): stub_mt},
    )
    resolver = ReadPhaseResolver(source)
    collection = VariantCollection([v1, v2])
    effects = collection.effects(phase_resolver=resolver)
    per_variant_for_v1 = [
        e for e in effects
        if e.variant == v1 and not isinstance(e, HaplotypeEffect)]
    per_variant_for_v2 = [
        e for e in effects
        if e.variant == v2 and not isinstance(e, HaplotypeEffect)]
    assert per_variant_for_v1
    assert per_variant_for_v2


def test_read_resolver_no_haplotype_when_only_one_variant():
    from pyensembl import cached_release
    g = cached_release(81)
    transcript = g.transcript_by_id("ENST00000003084")
    v1 = Variant("7", 117531100, "T", "A", g)
    stub_mt = MutantTranscript(
        reference_transcript=transcript,
        cdna_sequence="ACGT" * 50,
        mutant_protein_sequence="M" + "A" * 99,
        annotator_name="upstream_rna_tool",
    )
    source = _StubReadPhasingSource(
        phasing={v1: (v1,)},
        mutant_transcripts={(v1, transcript.id): stub_mt},
    )
    resolver = ReadPhaseResolver(source)
    effects = VariantCollection([v1]).effects(phase_resolver=resolver)
    assert not any(isinstance(e, HaplotypeEffect) for e in effects)


def test_read_resolver_no_haplotype_when_variants_trans():
    from pyensembl import cached_release
    g = cached_release(81)
    transcript = g.transcript_by_id("ENST00000003084")
    v1 = Variant("7", 117531100, "T", "A", g)
    v2 = Variant("7", 117531114, "G", "T", g)
    mt1 = MutantTranscript(
        reference_transcript=transcript,
        cdna_sequence="ACGT" * 50,
        annotator_name="upstream_rna_tool")
    mt2 = MutantTranscript(
        reference_transcript=transcript,
        cdna_sequence="TGCA" * 50,
        annotator_name="upstream_rna_tool")
    # Each variant has evidence only of itself — they're trans.
    source = _StubReadPhasingSource(
        phasing={v1: (v1,), v2: (v2,)},
        mutant_transcripts={
            (v1, transcript.id): mt1,
            (v2, transcript.id): mt2,
        },
    )
    resolver = ReadPhaseResolver(source)
    effects = VariantCollection([v1, v2]).effects(phase_resolver=resolver)
    assert not any(isinstance(e, HaplotypeEffect) for e in effects)


# -----------------------------------------------------------------
# VCF-resolver-driven joint effect
# -----------------------------------------------------------------


VCF_BODY = """##fileformat=VCFv4.1
##reference=GRCh38
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor
7\t117531100\t.\tT\tA\t100\tPASS\t.\tGT:PS\t0|1:100
7\t117531114\t.\tG\tT\t100\tPASS\t.\tGT:PS\t0|1:100
7\t117531120\t.\tG\tA\t100\tPASS\t.\tGT:PS\t1|0:100
"""


@pytest.fixture
def phased_vcf():
    fd, path = tempfile.mkstemp(suffix=".vcf")
    with os.fdopen(fd, "w") as f:
        f.write(VCF_BODY)
    yield path
    os.unlink(path)


def test_vcf_resolver_produces_haplotype_effect_for_cis_pair(phased_vcf):
    """Two variants in the same phase set on the same haplotype
    slot → joint HaplotypeEffect with DNA-inferred
    MutantTranscript."""
    vc = load_vcf(phased_vcf, genome="GRCh38")
    resolver = VCFPhaseResolver(vc, sample="tumor")
    effects = vc.effects(phase_resolver=resolver, raise_on_error=False)
    haplotype_effects = [e for e in effects if isinstance(e, HaplotypeEffect)]
    assert haplotype_effects
    # All HaplotypeEffects on this sample should be from vcf_ps.
    assert all(he.phase_source == "vcf_ps" for he in haplotype_effects)
    # The cis pair is v1 (117531100) + v2 (117531114), both on slot 1.
    # v3 (117531120) is slot 0 → trans → not in the same haplotype
    # group.
    for he in haplotype_effects:
        positions = sorted(v.start for v in he.variants)
        assert 117531120 not in positions


def test_vcf_resolver_haplotype_mutant_transcript_has_both_edits(phased_vcf):
    """DNA-inferred joint MutantTranscript contains edits for every
    cis variant — one combined protein sequence, not two."""
    vc = load_vcf(phased_vcf, genome="GRCh38")
    resolver = VCFPhaseResolver(vc, sample="tumor")
    effects = vc.effects(phase_resolver=resolver, raise_on_error=False)
    he = next(e for e in effects if isinstance(e, HaplotypeEffect))
    assert he.mutant_transcript is not None
    assert len(he.mutant_transcript.edits) == len(he.variants)
    # DNA-inferred, not from Isovar.
    assert he.mutant_transcript.annotator_name == "protein_diff"


# -----------------------------------------------------------------
# HGVS-style short description
# -----------------------------------------------------------------


def test_haplotype_effect_short_description_uses_bracket_notation():
    """HGVS cis notation: ``[v1;v2]`` — distinguishes from the trans
    form ``[v1];[v2]``."""
    from pyensembl import cached_release
    g = cached_release(81)
    transcript = g.transcript_by_id("ENST00000003084")
    v1 = Variant("7", 117531100, "T", "A", g)
    v2 = Variant("7", 117531114, "G", "T", g)
    mt = apply_variants_to_transcript([v1, v2], transcript)
    he = HaplotypeEffect(
        variants=[v1, v2],
        transcript=transcript,
        mutant_transcript=mt,
        phase_source="vcf_ps",
    )
    desc = he.short_description
    assert desc.startswith("[") and desc.endswith("]")
    # Bracket holds both variant descriptions, semicolon-separated
    # (cis convention).
    assert desc.count(";") == 1
