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
Tests for splice-site effect prediction, covering gaps identified in
https://github.com/openvax/varcode/issues/262

Uses CFTR (ENST00000003084) on chr7 forward strand (GRCh38/Ensembl 81)
as the primary test transcript.

Exon 4: genomic 117530899-117531114, transcript offset 405-620
  - Last 3 bases of exon 4 in transcript: AAG (matches MAG splice signal)
  - First base of exon 4 in transcript: G (purine, matches YAG|R acceptor)

Exon 6: genomic 117535248-117535411, transcript offset 711-874
  - Last 3 bases of exon 6: CAG (matches MAG splice signal)
  - First base of exon 6: G (purine)
"""

import pytest
from pyensembl import cached_release

from varcode import Variant
from varcode.effects import (
    ExonicSpliceSite,
    IntronicSpliceSite,
    Intronic,
    SpliceAcceptor,
    SpliceDonor,
    Substitution,
)
from varcode.effects.effect_prediction import _canonical_splice_base

from .common import expect_effect

ensembl_grch38 = cached_release(81)

CFTR_TRANSCRIPT_ID = "ENST00000003084"
BRCA1_TRANSCRIPT_ID = "ENST00000357654"

# -----------------------------------------------------------------------
# Exonic splice site detection: the changes_exonic_splice_site bug
#
# The bug: end_of_variant_exon is set to end_of_reference_exon without
# applying the mutation, so the function checks the reference against
# itself.  This causes splice-preserving mutations to be misclassified
# as ExonicSpliceSite, and splice-disrupting mutations to sometimes be
# missed (if the reference doesn't match the pattern).
# -----------------------------------------------------------------------


def test_exonic_splice_preserving_mutation_not_classified_as_splice_site():
    """SNV at -3 of exon that preserves the MAG splice signal.

    CFTR exon 4 ends with AAG. A->C at the -3 position gives CAG.
    Both AAG and CAG match the MAG pattern (M = A or C).
    This should be a Substitution, NOT an ExonicSpliceSite.

    This test exposes the bug in changes_exonic_splice_site() where
    the reference sequence is checked against itself.
    """
    # chr7:117531112 is the A in AAG (third-to-last base of exon 4)
    variant = Variant("7", 117531112, "A", "C", ensembl_grch38)
    expect_effect(
        variant,
        transcript_id=CFTR_TRANSCRIPT_ID,
        effect_class=Substitution,
    )


def test_exonic_splice_disrupting_mutation_classified_as_splice_site():
    """SNV at last base of exon that breaks the MAG splice signal.

    CFTR exon 4 ends with AAG. G->T at the last position gives AAT.
    AAT does not match MAG (must end in G). This should be ExonicSpliceSite.
    """
    # chr7:117531114 is the G in AAG (last base of exon 4)
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    expect_effect(
        variant,
        transcript_id=CFTR_TRANSCRIPT_ID,
        effect_class=ExonicSpliceSite,
    )


def test_exonic_splice_disrupting_mutation_has_coding_alternate():
    """ExonicSpliceSite should carry the coding effect as alternate_effect.

    When a variant both disrupts a splice signal and changes coding sequence,
    the ExonicSpliceSite.alternate_effect should be the coding effect (e.g.,
    Substitution), so both effects are accessible.
    """
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effect = variant.effect_on_transcript(transcript)
    assert effect.__class__ is ExonicSpliceSite, \
        "Expected ExonicSpliceSite, got %s" % effect.__class__.__name__
    assert effect.alternate_effect is not None, \
        "ExonicSpliceSite should have an alternate_effect"
    assert effect.alternate_effect.__class__ is Substitution, \
        "Expected alternate_effect to be Substitution, got %s" % (
            effect.alternate_effect.__class__.__name__,)


def test_exonic_splice_site_at_exon_start_acceptor_disrupted():
    """SNV at first base of non-first exon that disrupts the 3' acceptor signal.

    The 3' splice site pattern is YAG|R where R is a purine on the exon side.
    CFTR exon 4 starts with G (purine). Changing G->C (pyrimidine) breaks it.
    """
    # chr7:117530899 is the first base of exon 4
    variant = Variant("7", 117530899, "G", "C", ensembl_grch38)
    expect_effect(
        variant,
        transcript_id=CFTR_TRANSCRIPT_ID,
        effect_class=ExonicSpliceSite,
    )


def test_exonic_splice_site_at_exon_start_acceptor_preserved():
    """SNV at first base of exon that preserves the 3' acceptor signal.

    CFTR exon 4 starts with G (purine). Changing G->A (also purine) preserves
    the YAG|R pattern. This should NOT be ExonicSpliceSite.
    """
    # chr7:117530899, G->A (purine to purine, acceptor signal preserved)
    variant = Variant("7", 117530899, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effect = variant.effect_on_transcript(transcript)
    # Should be a coding effect, not ExonicSpliceSite
    assert effect.__class__ is not ExonicSpliceSite, \
        "Purine-to-purine change at exon start should not be ExonicSpliceSite"


# -----------------------------------------------------------------------
# Intronic splice classification: SpliceDonor and SpliceAcceptor tests
#
# These effect classes had ZERO tests (noted by the TODO comment in
# test_effect_classes.py line 40: "# TODO: SpliceDonor, SpliceReceptor").
# The current classification is purely distance-based and ignores
# sequence content.
# -----------------------------------------------------------------------


def test_splice_donor_at_plus_1():
    """Intronic variant at +1 after exon (first intronic base after donor).

    For a forward-strand gene, this is exon.end + 1.
    Should be classified as SpliceDonor.
    """
    # chr7:117531115 is +1 after CFTR exon 4
    variant = Variant("7", 117531115, "G", "T", ensembl_grch38)
    expect_effect(
        variant,
        transcript_id=CFTR_TRANSCRIPT_ID,
        effect_class=SpliceDonor,
    )


def test_splice_donor_at_plus_2():
    """Intronic variant at +2 after exon (second intronic base after donor).

    Canonical + strand base at donor +2 is T (from GT signal).
    Should be classified as SpliceDonor.
    """
    # chr7:117531116 is +2 after CFTR exon 4, canonical ref=T
    variant = Variant("7", 117531116, "T", "A", ensembl_grch38)
    expect_effect(
        variant,
        transcript_id=CFTR_TRANSCRIPT_ID,
        effect_class=SpliceDonor,
    )


def test_splice_acceptor_at_minus_1():
    """Intronic variant at -1 before exon (last intronic base before acceptor).

    For a forward-strand gene, this is exon.start - 1.
    Should be classified as SpliceAcceptor.
    """
    # chr7:117530898 is -1 before CFTR exon 4
    variant = Variant("7", 117530898, "G", "T", ensembl_grch38)
    expect_effect(
        variant,
        transcript_id=CFTR_TRANSCRIPT_ID,
        effect_class=SpliceAcceptor,
    )


def test_splice_acceptor_at_minus_2():
    """Intronic variant at -2 before exon (second-to-last intronic base).

    Canonical + strand base at acceptor -2 is A (from AG signal).
    Should be classified as SpliceAcceptor.
    """
    # chr7:117530897 is -2 before CFTR exon 4, canonical ref=A
    variant = Variant("7", 117530897, "A", "T", ensembl_grch38)
    expect_effect(
        variant,
        transcript_id=CFTR_TRANSCRIPT_ID,
        effect_class=SpliceAcceptor,
    )


def test_intronic_splice_site_at_plus_3_to_6():
    """Intronic variant at +3 to +6 after exon should be IntronicSpliceSite.

    These positions are implicated in splicing but not as confidently as +1/+2.
    """
    for offset in [3, 4, 5, 6]:
        pos = 117531114 + offset  # exon 4 end + offset
        variant = Variant("7", pos, "G", "T", ensembl_grch38)
        transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
        effect = variant.effect_on_transcript(transcript)
        assert effect.__class__ is IntronicSpliceSite, \
            "Expected IntronicSpliceSite at +%d, got %s" % (
                offset, effect.__class__.__name__)


def test_intronic_splice_site_at_minus_3():
    """Intronic variant at -3 before exon should be IntronicSpliceSite.

    Position -3 is part of the 3' splice motif but allows more degeneracy.
    """
    # chr7:117530896 is -3 before CFTR exon 4
    variant = Variant("7", 117530896, "G", "T", ensembl_grch38)
    expect_effect(
        variant,
        transcript_id=CFTR_TRANSCRIPT_ID,
        effect_class=IntronicSpliceSite,
    )


def test_deep_intronic_variant_is_intronic():
    """Intronic variant far from any exon should be plain Intronic."""
    # 20bp into the intron before CFTR exon 4
    variant = Variant("7", 117530899 - 20, "G", "T", ensembl_grch38)
    expect_effect(
        variant,
        transcript_id=CFTR_TRANSCRIPT_ID,
        effect_class=Intronic,
    )


# -----------------------------------------------------------------------
# ExonicSpliceSite interaction with drop_silent_and_noncoding (issue #136)
# -----------------------------------------------------------------------


def test_exonic_splice_site_not_dropped_by_drop_silent_and_noncoding():
    """ExonicSpliceSite with a coding alternate should survive
    drop_silent_and_noncoding().

    Regression test for https://github.com/openvax/varcode/issues/136
    """
    # Use a variant that produces ExonicSpliceSite with Substitution alternate
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    effects = variant.effects()

    # Verify there's at least one ExonicSpliceSite
    splice_effects = [e for e in effects if e.__class__ is ExonicSpliceSite]
    assert len(splice_effects) > 0, "Expected at least one ExonicSpliceSite"

    # After dropping silent and noncoding, the ExonicSpliceSite with a coding
    # alternate should remain
    remaining = effects.drop_silent_and_noncoding()
    assert len(remaining) > 0, \
        "drop_silent_and_noncoding() should not drop ExonicSpliceSite with coding alternate"


# -----------------------------------------------------------------------
# Reverse strand splice site tests (using BRCA1)
#
# BRCA1-001 (ENST00000357654) on chr17 reverse strand.
# Exon 12: ENSE00003527960, genomic 43082404-43082575
# On reverse strand, the "end" of an exon in transcript order is at the
# genomic start coordinate.
# -----------------------------------------------------------------------


def test_splice_donor_reverse_strand():
    """SpliceDonor on reverse strand gene (BRCA1 exon 12).

    On reverse strand, the donor site (3' end of exon in transcript order)
    is at the genomic start of the exon. So the first intronic bases are
    at genomic positions exon.start - 1 and exon.start - 2.

    The canonical GT donor on the - strand appears as AC on the + strand.
    So at exon.start - 1 the + strand reference is C (complement of G),
    and at exon.start - 2 it is A (complement of T).
    """
    # Exon 12 of BRCA1 starts at 43082404 on the reverse strand.
    # Donor +1 intronic = 43082403, canonical + strand ref = C
    variant = Variant("17", 43082403, "C", "T", ensembl_grch38)
    expect_effect(
        variant,
        transcript_id=BRCA1_TRANSCRIPT_ID,
        effect_class=SpliceDonor,
    )


def test_splice_donor_reverse_strand_plus_2():
    """SpliceDonor at +2 on reverse strand (BRCA1 exon 12).

    At exon.start - 2, the canonical + strand ref is A (complement of T).
    """
    variant = Variant("17", 43082402, "A", "G", ensembl_grch38)
    expect_effect(
        variant,
        transcript_id=BRCA1_TRANSCRIPT_ID,
        effect_class=SpliceDonor,
    )


def test_splice_acceptor_reverse_strand():
    """SpliceAcceptor on reverse strand gene (BRCA1 exon 12).

    On reverse strand, the acceptor site (5' end of exon in transcript order)
    is at the genomic end of the exon. The canonical AG acceptor on the - strand
    appears as CT on the + strand. At exon.end + 1 the + strand ref is C
    (complement of G), and at exon.end + 2 it is T (complement of A).
    """
    # Exon 12 of BRCA1 ends at 43082575.
    # Acceptor -1 intronic = 43082576, canonical + strand ref = C
    variant = Variant("17", 43082576, "C", "A", ensembl_grch38)
    expect_effect(
        variant,
        transcript_id=BRCA1_TRANSCRIPT_ID,
        effect_class=SpliceAcceptor,
    )


def test_splice_acceptor_reverse_strand_minus_2():
    """SpliceAcceptor at -2 on reverse strand (BRCA1 exon 12).

    At exon.end + 2, the canonical + strand ref is T (complement of A).
    """
    variant = Variant("17", 43082577, "T", "G", ensembl_grch38)
    expect_effect(
        variant,
        transcript_id=BRCA1_TRANSCRIPT_ID,
        effect_class=SpliceAcceptor,
    )


# -----------------------------------------------------------------------
# _canonical_splice_base unit tests
# -----------------------------------------------------------------------


class TestCanonicalSpliceBase:
    """Direct tests for the _canonical_splice_base helper function."""

    # Forward strand donor: GT at +1,+2
    def test_forward_donor_plus_1(self):
        assert _canonical_splice_base("+", before_exon=False, distance_to_exon=1) == "G"

    def test_forward_donor_plus_2(self):
        assert _canonical_splice_base("+", before_exon=False, distance_to_exon=2) == "T"

    # Forward strand acceptor: AG at -2,-1
    def test_forward_acceptor_minus_1(self):
        assert _canonical_splice_base("+", before_exon=True, distance_to_exon=1) == "G"

    def test_forward_acceptor_minus_2(self):
        assert _canonical_splice_base("+", before_exon=True, distance_to_exon=2) == "A"

    # Reverse strand donor: GT on - strand = AC on + strand
    def test_reverse_donor_minus_1(self):
        assert _canonical_splice_base("-", before_exon=False, distance_to_exon=1) == "C"

    def test_reverse_donor_minus_2(self):
        assert _canonical_splice_base("-", before_exon=False, distance_to_exon=2) == "A"

    # Reverse strand acceptor: AG on - strand = CT on + strand
    def test_reverse_acceptor_plus_1(self):
        assert _canonical_splice_base("-", before_exon=True, distance_to_exon=1) == "C"

    def test_reverse_acceptor_plus_2(self):
        assert _canonical_splice_base("-", before_exon=True, distance_to_exon=2) == "T"

    # Out of range distances return None
    def test_distance_3_returns_none(self):
        assert _canonical_splice_base("+", before_exon=True, distance_to_exon=3) is None

    def test_distance_0_returns_none(self):
        assert _canonical_splice_base("+", before_exon=True, distance_to_exon=0) is None


# -----------------------------------------------------------------------
# Sequence-aware intronic splice classification
#
# Variants at canonical splice dinucleotide positions (+1/+2 or -1/-2)
# should be classified as SpliceDonor/SpliceAcceptor only if the
# reference base matches the canonical signal. If the reference is
# non-canonical, the site is already unusual and should be downgraded
# to IntronicSpliceSite.
# -----------------------------------------------------------------------


class TestSequenceAwareSpliceDonor:
    """Test that SpliceDonor requires canonical reference base."""

    def test_forward_donor_plus_1_canonical_ref_is_splice_donor(self):
        """Ref=G at donor +1 (canonical GT) → SpliceDonor."""
        variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
        expect_effect(
            variant,
            transcript_id=CFTR_TRANSCRIPT_ID,
            effect_class=SpliceDonor,
        )

    def test_forward_donor_plus_1_noncanonical_ref_is_intronic_splice_site(self):
        """Ref=A at donor +1 (not canonical G) → IntronicSpliceSite."""
        variant = Variant("7", 117531115, "A", "G", ensembl_grch38)
        expect_effect(
            variant,
            transcript_id=CFTR_TRANSCRIPT_ID,
            effect_class=IntronicSpliceSite,
        )

    def test_forward_donor_plus_2_canonical_ref_is_splice_donor(self):
        """Ref=T at donor +2 (canonical GT) → SpliceDonor."""
        variant = Variant("7", 117531116, "T", "A", ensembl_grch38)
        expect_effect(
            variant,
            transcript_id=CFTR_TRANSCRIPT_ID,
            effect_class=SpliceDonor,
        )

    def test_forward_donor_plus_2_noncanonical_ref_is_intronic_splice_site(self):
        """Ref=C at donor +2 (not canonical T) → IntronicSpliceSite."""
        variant = Variant("7", 117531116, "C", "A", ensembl_grch38)
        expect_effect(
            variant,
            transcript_id=CFTR_TRANSCRIPT_ID,
            effect_class=IntronicSpliceSite,
        )

    def test_reverse_donor_minus_1_canonical_ref_is_splice_donor(self):
        """Ref=C at reverse-strand donor -1 (canonical) → SpliceDonor."""
        variant = Variant("17", 43082403, "C", "A", ensembl_grch38)
        expect_effect(
            variant,
            transcript_id=BRCA1_TRANSCRIPT_ID,
            effect_class=SpliceDonor,
        )

    def test_reverse_donor_minus_1_noncanonical_ref_is_intronic_splice_site(self):
        """Ref=G at reverse-strand donor -1 (not canonical C) → IntronicSpliceSite."""
        variant = Variant("17", 43082403, "G", "A", ensembl_grch38)
        expect_effect(
            variant,
            transcript_id=BRCA1_TRANSCRIPT_ID,
            effect_class=IntronicSpliceSite,
        )

    def test_reverse_donor_minus_2_canonical_ref_is_splice_donor(self):
        """Ref=A at reverse-strand donor -2 (canonical) → SpliceDonor."""
        variant = Variant("17", 43082402, "A", "G", ensembl_grch38)
        expect_effect(
            variant,
            transcript_id=BRCA1_TRANSCRIPT_ID,
            effect_class=SpliceDonor,
        )

    def test_reverse_donor_minus_2_noncanonical_ref_is_intronic_splice_site(self):
        """Ref=T at reverse-strand donor -2 (not canonical A) → IntronicSpliceSite."""
        variant = Variant("17", 43082402, "T", "G", ensembl_grch38)
        expect_effect(
            variant,
            transcript_id=BRCA1_TRANSCRIPT_ID,
            effect_class=IntronicSpliceSite,
        )


class TestSequenceAwareSpliceAcceptor:
    """Test that SpliceAcceptor requires canonical reference base."""

    def test_forward_acceptor_minus_1_canonical_ref_is_splice_acceptor(self):
        """Ref=G at acceptor -1 (canonical AG) → SpliceAcceptor."""
        variant = Variant("7", 117530898, "G", "A", ensembl_grch38)
        expect_effect(
            variant,
            transcript_id=CFTR_TRANSCRIPT_ID,
            effect_class=SpliceAcceptor,
        )

    def test_forward_acceptor_minus_1_noncanonical_ref_is_intronic_splice_site(self):
        """Ref=T at acceptor -1 (not canonical G) → IntronicSpliceSite."""
        variant = Variant("7", 117530898, "T", "A", ensembl_grch38)
        expect_effect(
            variant,
            transcript_id=CFTR_TRANSCRIPT_ID,
            effect_class=IntronicSpliceSite,
        )

    def test_forward_acceptor_minus_2_canonical_ref_is_splice_acceptor(self):
        """Ref=A at acceptor -2 (canonical AG) → SpliceAcceptor."""
        variant = Variant("7", 117530897, "A", "G", ensembl_grch38)
        expect_effect(
            variant,
            transcript_id=CFTR_TRANSCRIPT_ID,
            effect_class=SpliceAcceptor,
        )

    def test_forward_acceptor_minus_2_noncanonical_ref_is_intronic_splice_site(self):
        """Ref=C at acceptor -2 (not canonical A) → IntronicSpliceSite."""
        variant = Variant("7", 117530897, "C", "G", ensembl_grch38)
        expect_effect(
            variant,
            transcript_id=CFTR_TRANSCRIPT_ID,
            effect_class=IntronicSpliceSite,
        )

    def test_reverse_acceptor_plus_1_canonical_ref_is_splice_acceptor(self):
        """Ref=C at reverse-strand acceptor +1 (canonical) → SpliceAcceptor."""
        variant = Variant("17", 43082576, "C", "A", ensembl_grch38)
        expect_effect(
            variant,
            transcript_id=BRCA1_TRANSCRIPT_ID,
            effect_class=SpliceAcceptor,
        )

    def test_reverse_acceptor_plus_1_noncanonical_ref_is_intronic_splice_site(self):
        """Ref=G at reverse-strand acceptor +1 (not canonical C) → IntronicSpliceSite."""
        variant = Variant("17", 43082576, "G", "A", ensembl_grch38)
        expect_effect(
            variant,
            transcript_id=BRCA1_TRANSCRIPT_ID,
            effect_class=IntronicSpliceSite,
        )

    def test_reverse_acceptor_plus_2_canonical_ref_is_splice_acceptor(self):
        """Ref=T at reverse-strand acceptor +2 (canonical) → SpliceAcceptor."""
        variant = Variant("17", 43082577, "T", "G", ensembl_grch38)
        expect_effect(
            variant,
            transcript_id=BRCA1_TRANSCRIPT_ID,
            effect_class=SpliceAcceptor,
        )

    def test_reverse_acceptor_plus_2_noncanonical_ref_is_intronic_splice_site(self):
        """Ref=A at reverse-strand acceptor +2 (not canonical T) → IntronicSpliceSite."""
        variant = Variant("17", 43082577, "A", "G", ensembl_grch38)
        expect_effect(
            variant,
            transcript_id=BRCA1_TRANSCRIPT_ID,
            effect_class=IntronicSpliceSite,
        )


class TestMultiBaseSpliceVariants:
    """Multi-base variants at splice sites should still be classified by distance
    (the sequence check only applies to SNVs where we can unambiguously
    determine the reference base at a single position)."""

    def test_deletion_at_donor_is_splice_donor(self):
        """2bp deletion spanning donor +1/+2 → SpliceDonor regardless of bases."""
        variant = Variant("7", 117531115, "GT", "", ensembl_grch38)
        expect_effect(
            variant,
            transcript_id=CFTR_TRANSCRIPT_ID,
            effect_class=SpliceDonor,
        )

    def test_insertion_at_donor_is_splice_donor(self):
        """Insertion at donor +1 → SpliceDonor regardless of bases."""
        variant = Variant("7", 117531115, "", "AAA", ensembl_grch38)
        expect_effect(
            variant,
            transcript_id=CFTR_TRANSCRIPT_ID,
            effect_class=SpliceDonor,
        )
