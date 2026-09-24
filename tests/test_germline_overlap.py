"""Allele-only and genotype-bearing regressions for unsupported LOH claims."""

import pytest
from pyensembl import cached_release

from varcode import (
    Completeness, EffectCollection, GermlineAlleleOverlap, GermlineContext,
    Variant, detect_germline_overlap,
)
from varcode.effects import Unresolved
from varcode.transcript_model import predict_transcript_model_effect


@pytest.fixture
def allele():
    return Variant("7", 117531101, "T", "C", genome=81)


@pytest.mark.parametrize("backend", ["fast", "protein_diff", "transcript_model"])
@pytest.mark.parametrize("normal_gt,tumor_gt", [("0/1", "0/1"), ("1/1", "1/1"), ("0/1", "1/1")])
def test_genotype_calls_do_not_manufacture_loh(allele, backend, normal_gt, tumor_gt, tmp_path):
    # Fixture construction is shipped with the test, including both genotypes.
    path = tmp_path / "paired.vcf"
    path.write_text(
        '##fileformat=VCFv4.2\n##reference=GRCh38\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n'
        f'7\t117531101\t.\tT\tC\t.\tPASS\t.\tGT\t{normal_gt}\t{tumor_gt}\n')
    context = GermlineContext.from_multi_sample_vcf(
        str(path), sample="NORMAL", completeness=Completeness.COMPLETE)
    effect = next(e for e in allele.effects(annotator=backend, germline=context)
                  if e.transcript_id == "ENST00000003084")
    assert isinstance(effect, GermlineAlleleOverlap)
    assert effect.is_germline_overlap
    assert effect.is_loh is None and effect.loh_status == "not_assessed"
    assert effect.modifies_coding_sequence is False
    assert effect.modifies_protein_sequence is False
    assert effect.mutant_protein_sequence is None
    restored = EffectCollection.from_json(EffectCollection([effect]).to_json())[0]
    assert restored.is_germline_overlap and restored.is_loh is None


@pytest.mark.parametrize("backend", ["fast", "protein_diff", "transcript_model"])
def test_missing_genotype_still_only_establishes_overlap(allele, backend):
    context = GermlineContext.from_variants([allele])
    effects = allele.effects(annotator=backend, germline=context)
    assert effects
    assert all(isinstance(e, GermlineAlleleOverlap) for e in effects)
    assert all(e.is_loh is None for e in effects)


def test_overlap_requires_matching_build_and_alleles(allele):
    assert detect_germline_overlap(allele, [allele])
    for kwargs in ({"genome": 75}, {"alt": "G"}, {"ref": "A"}, {"start": 117531102}):
        args = dict(contig="7", start=117531101, ref="T", alt="C", genome=81)
        args.update(kwargs)
        assert not detect_germline_overlap(allele, [Variant(**args)])


def test_mixed_group_does_not_apply_inherited_allele_twice(allele):
    transcript = cached_release(81).transcript_by_id("ENST00000003084")
    novel = Variant("7", 117531100, "T", "A", genome=81)
    effect = predict_transcript_model_effect(
        (novel, allele), transcript, germline_variants=(allele,))
    assert isinstance(effect, Unresolved)
    assert effect.mechanism == "germline_overlap_haplotype"
    assert effect.is_germline_overlap and effect.is_loh is None
