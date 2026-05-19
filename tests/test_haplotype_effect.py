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


def test_apply_variants_rejects_insertion_abutting_deletion():
    """An insertion at offset X plus a deletion starting at X resolve
    to edits with the same cdna_start. Applying them in either order
    yields a different mutant cDNA (the deletion either consumes the
    inserted bases or leaves them in place), so the joint result is
    order-dependent. apply_variants_to_transcript must refuse the
    combination rather than silently picking one ordering."""
    from pyensembl import cached_release
    g = cached_release(81)
    transcript = g.transcript_by_id("ENST00000003084")  # CFTR, forward
    # Insertion of "AAA" after genomic position 117531100 → edit at
    # cdna_start = cdna_offset(117531100) + 1, zero-width.
    v_ins = Variant("7", 117531100, "T", "TAAA", g)
    # Deletion of the 2 bases at 117531101-117531102 → edit anchored
    # at the same cdna_start as the insertion above.
    v_del = Variant("7", 117531100, "TTG", "T", g)
    assert apply_variants_to_transcript([v_ins, v_del], transcript) is None
    # Symmetric: input order doesn't matter — both orderings refused.
    assert apply_variants_to_transcript([v_del, v_ins], transcript) is None


def test_apply_variants_rejects_two_insertions_at_same_offset():
    """Two distinct insertions at the same cDNA offset are
    order-dependent (which one ends up upstream in the mutant cDNA
    depends on apply order). Previously caught by a 4-way equality
    check, now caught by the generalised insertion-at-boundary rule —
    re-pin that the case is still refused after the rule restructure."""
    from pyensembl import cached_release
    g = cached_release(81)
    transcript = g.transcript_by_id("ENST00000003084")  # CFTR, forward
    # Both anchored at genomic 117531100 → both edits resolve to a
    # zero-width insertion at the same cdna_start.
    v_ins_a = Variant("7", 117531100, "T", "TAAA", g)
    v_ins_b = Variant("7", 117531100, "T", "TGGG", g)
    assert apply_variants_to_transcript([v_ins_a, v_ins_b], transcript) is None
    assert apply_variants_to_transcript([v_ins_b, v_ins_a], transcript) is None


# -----------------------------------------------------------------
# Joint-effect biological semantics — these tests verify the headline
# value-add of the haplotype pipeline: cis edits compose, they don't
# get classified independently.
# -----------------------------------------------------------------


def _apply_one(variant, transcript):
    """Convenience: per-variant mutant protein via the singular
    builder, for comparison against joint mutant proteins."""
    from varcode import apply_variant_to_transcript
    return apply_variant_to_transcript(variant, transcript)


def test_joint_frameshift_cancellation_at_cdna_level():
    """A +1 insertion followed by a -1 deletion downstream cancel as a
    pair: the joint cDNA has zero net length delta and restores the
    reading frame past the deletion site. The joint mutant cDNA from
    the deletion onward matches the reference cDNA. Each variant
    alone is a frameshift that produces a different cDNA than the
    joint.

    Note: the joint protein-level outcome depends on whether the
    short shifted region between the ins and del happens to contain
    a stop codon (every ~21 bases on average), so this test pins
    cDNA-level composition — the structural guarantee — rather than
    protein-level identity to the reference. That guarantee IS the
    headline of the joint pipeline: cis edits compose at the cDNA
    level, and the protein is read off the composed cDNA, not
    inferred per-variant."""
    from pyensembl import cached_release
    g = cached_release(81)
    transcript = g.transcript_by_id("ENST00000003084")  # CFTR, forward
    # +1 insertion at genomic 117531094 (cdna 600). Insertion lands at
    # cdna_start = 601 per the forward-strand rule.
    v_ins = Variant("7", 117531094, "T", "TA", g)
    # -1 deletion 11 bp downstream (drops the 'T' at genomic 117531105,
    # cdna 611). Net joint length delta = 0.
    v_del = Variant("7", 117531104, "AT", "A", g)
    joint = apply_variants_to_transcript([v_ins, v_del], transcript)
    mt_ins = _apply_one(v_ins, transcript)
    mt_del = _apply_one(v_del, transcript)
    assert joint is not None
    assert joint.total_length_delta == 0, (
        "ins(+1) + del(-1) must yield zero net length delta")
    # Joint cDNA matches reference length exactly.
    ref_cdna = str(transcript.sequence)
    assert len(joint.cdna_sequence) == len(ref_cdna)
    # Past the deletion site (cdna 612 onwards), joint cDNA is back in
    # register with the reference — the frame is restored.
    assert joint.cdna_sequence[612:] == ref_cdna[612:], (
        "Joint cDNA past the deletion site must equal reference cDNA — "
        "the +1/-1 pair cancels the shift.")
    # The joint cDNA differs from EACH per-variant cDNA in the
    # substituted middle stretch — the joint isn't a relabelling of
    # either per-variant outcome.
    assert joint.cdna_sequence != mt_ins.cdna_sequence
    assert joint.cdna_sequence != mt_del.cdna_sequence


def test_joint_premature_stop_truncates_below_per_variant():
    """Two cis SNVs in the same codon, each missense alone, together
    create a TAA stop. Joint protein is dramatically shorter than
    either per-variant protein. Pins that PrematureStop classification
    is driven by the joint protein, not by per-variant proteins."""
    from pyensembl import cached_release
    g = cached_release(81)
    transcript = g.transcript_by_id("ENST00000003084")  # CFTR, forward
    # CFTR cDNA offset 420 = codon "GTA" (Val). Variants:
    #   v1: G→T at genomic 117530914 (cdna 420) → codon "TTA" Leu
    #   v2: T→A at genomic 117530915 (cdna 421) → codon "GAA" Glu
    #   joint: codon "TAA" → stop
    v1 = Variant("7", 117530914, "G", "T", g)
    v2 = Variant("7", 117530915, "T", "A", g)
    mt_v1 = _apply_one(v1, transcript)
    mt_v2 = _apply_one(v2, transcript)
    joint = apply_variants_to_transcript([v1, v2], transcript)
    assert mt_v1.mutant_protein_sequence is not None
    assert mt_v2.mutant_protein_sequence is not None
    assert joint is not None and joint.mutant_protein_sequence is not None
    # Each per-variant protein retains full reference length (missense).
    ref_len = len(str(transcript.protein_sequence))
    assert len(mt_v1.mutant_protein_sequence) == ref_len
    assert len(mt_v2.mutant_protein_sequence) == ref_len
    # Joint protein truncates at the new stop, well below the
    # per-variant lengths.
    assert len(joint.mutant_protein_sequence) < ref_len - 100, (
        "Joint protein length=%d, per-variant length=%d — joint "
        "should be truncated by the new stop." % (
            len(joint.mutant_protein_sequence), ref_len))


def test_joint_upstream_frameshift_rescues_downstream_stop():
    """v1 alone introduces an in-frame premature stop (codon TAT at
    cdna [612, 615) → TAA via T→A at cdna 614). v2 is an upstream +1
    insertion that shifts the reading frame past v1's site. In the
    joint, v1's would-be stop codon is no longer at a codon boundary,
    AND the change v1 introduces (T→A at orig cdna 614) flips the
    next shifted-frame codon away from a stop — so the joint protein
    extends past where v1-only terminated."""
    from pyensembl import cached_release
    g = cached_release(81)
    transcript = g.transcript_by_id("ENST00000003084")  # CFTR, forward
    # v1: T→A at genomic 117531108 (cdna 614, third base of codon TAT
    # at cdna [612,615)) → in-frame TAA stop.
    v1 = Variant("7", 117531108, "T", "A", g)
    # v2: +1 insertion of "A" at cdna 604 — the cdna anchor at genomic
    # 117531097 has base 'A', so Variant ref="A" alt="AA" trims to a
    # pure insertion of "A". Chosen because the +1 shifted reading
    # from cdna 604 has no stop until after v1's site; v1's T→A change
    # then prevents the shifted codon at [615,618) from being TAA.
    v2 = Variant("7", 117531097, "A", "AA", g)
    mt_v1 = _apply_one(v1, transcript)
    joint = apply_variants_to_transcript([v1, v2], transcript)
    assert mt_v1.mutant_protein_sequence is not None
    assert joint is not None and joint.mutant_protein_sequence is not None
    # The v1-alone protein truncates at the new stop. The joint
    # protein extends past v1's truncation because the upstream
    # frameshift shifted v1's codon out of frame.
    assert len(joint.mutant_protein_sequence) > len(mt_v1.mutant_protein_sequence), (
        "Joint protein length=%d should exceed v1-only length=%d "
        "(upstream frameshift removes v1's in-frame stop)." % (
            len(joint.mutant_protein_sequence),
            len(mt_v1.mutant_protein_sequence)))


# -----------------------------------------------------------------
# Coverage gaps in the existing matrix
# -----------------------------------------------------------------


def test_apply_variants_reverse_strand_composes_per_edit_complement():
    """Every existing multi-edit test uses CFTR (forward strand). Two
    cis SNVs on BRCA1 (reverse strand) must each be reverse-
    complemented individually before composing into the joint cDNA.
    Verifies the per-edit complementation isn't lost in the joint
    builder."""
    from pyensembl import cached_release
    g = cached_release(81)
    transcript = g.transcript_by_id("ENST00000357654")  # BRCA1, reverse
    assert transcript.on_backward_strand
    # Two SNVs in BRCA1 exon 11 (cdna offsets 4316 and 4309 — same
    # exon, ~7 bp apart on the reverse strand).
    # Position 43082563: genomic ref 'T' → cdna 'A' at offset 4316.
    # Position 43082570: genomic ref 'C' → cdna 'G' at offset 4309.
    v1 = Variant("17", 43082563, "T", "A", g)  # cdna A→T at 4316
    v2 = Variant("17", 43082570, "C", "A", g)  # cdna G→T at 4309
    joint = apply_variants_to_transcript([v1, v2], transcript)
    assert joint is not None
    assert len(joint.edits) == 2
    # Both edits' alt_bases are reverse-complemented from their
    # genomic alts: "A" → "T" each.
    for edit in joint.edits:
        assert edit.alt_bases == "T", (
            "Both edits' alt_bases must be reverse-complemented to 'T' "
            "on the negative strand; got %r" % edit.alt_bases)
    # Joint cDNA differs from reference at exactly the two offsets.
    ref = str(transcript.sequence)
    diffs = [i for i in range(len(ref)) if ref[i] != joint.cdna_sequence[i]]
    assert sorted(diffs) == [4309, 4316], (
        "Joint cDNA should differ from reference at cdna offsets 4309 "
        "and 4316; got %r" % diffs)


def test_apply_variants_three_cis_variants_compose():
    """The conflict-detection loop only checks adjacent pairs in the
    sorted edit list — adjacency-only is provably sufficient for the
    sorted scan, but no test exercised N≥3. Three distinct SNVs in
    one transcript must compose into a single joint with all three
    edits surviving."""
    from pyensembl import cached_release
    g = cached_release(81)
    transcript = g.transcript_by_id("ENST00000003084")  # CFTR, forward
    v1 = Variant("7", 117531100, "T", "A", g)  # cdna 606 T→A
    v2 = Variant("7", 117531106, "T", "C", g)  # cdna 612 T→C
    v3 = Variant("7", 117531112, "A", "G", g)  # cdna 618 A→G
    joint = apply_variants_to_transcript([v1, v2, v3], transcript)
    assert joint is not None
    assert len(joint.edits) == 3
    ref = str(transcript.sequence)
    diffs = [i for i in range(len(ref)) if ref[i] != joint.cdna_sequence[i]]
    assert sorted(diffs) == [606, 612, 618]
    # Sanity: input order doesn't matter.
    rotated = apply_variants_to_transcript([v3, v1, v2], transcript)
    assert rotated.cdna_sequence == joint.cdna_sequence


def test_apply_variants_cross_exon_cis_variants_compose():
    """Two cis variants in different exons of the same transcript.
    Pins that cDNA-offset resolution composes correctly when the
    two edits are separated by an intron in genomic space but
    adjacent in cDNA space."""
    from pyensembl import cached_release
    g = cached_release(81)
    transcript = g.transcript_by_id("ENST00000003084")  # CFTR, forward
    # Exon 3 (117530899-117531114) and exon 4 (117534276-117534365)
    # — separated by ~3kb of intron, adjacent in cDNA.
    v_exon3 = Variant("7", 117531100, "T", "A", g)   # cdna 606
    v_exon4 = Variant("7", 117534280, "T", "G", g)   # cdna 625
    joint = apply_variants_to_transcript([v_exon3, v_exon4], transcript)
    assert joint is not None
    assert len(joint.edits) == 2
    cdna_starts = sorted(e.cdna_start for e in joint.edits)
    assert cdna_starts == [606, 625]
    ref = str(transcript.sequence)
    diffs = [i for i in range(len(ref)) if ref[i] != joint.cdna_sequence[i]]
    assert sorted(diffs) == [606, 625]


# -----------------------------------------------------------------
# Defensive / contract
# -----------------------------------------------------------------


def test_apply_variants_rejects_duplicate_variant():
    """Passing the same variant twice resolves to two edits with
    identical cdna ranges. The overlap detector catches this (range
    overlap branch). Pinning the behavior so a future "dedupe input
    variants" change is a deliberate choice, not a silent semantic
    shift."""
    from pyensembl import cached_release
    g = cached_release(81)
    transcript = g.transcript_by_id("ENST00000003084")  # CFTR, forward
    v = Variant("7", 117531100, "T", "A", g)
    assert apply_variants_to_transcript([v, v], transcript) is None


def test_apply_variants_rejects_whole_list_on_any_unresolvable():
    """One variant in the list fails to resolve (wrong reference base).
    The whole call must return None — the joint builder doesn't
    silently apply the well-formed subset and drop the bad variant,
    because that would produce a different effect than the caller
    asked for."""
    from pyensembl import cached_release
    g = cached_release(81)
    transcript = g.transcript_by_id("ENST00000003084")  # CFTR, forward
    v_ok = Variant("7", 117531100, "T", "A", g)
    # Wrong reference base — _resolve_variant_edit will return None
    # due to the ref-match gate.
    v_bad = Variant(
        "7", 117531106, "X", "A", g, allow_extended_nucleotides=True)
    assert apply_variants_to_transcript([v_ok, v_bad], transcript) is None
    # Symmetric: order of input doesn't matter.
    assert apply_variants_to_transcript([v_bad, v_ok], transcript) is None


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
