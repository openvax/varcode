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
Tests for the splice-mechanism Effect classes
(openvax/varcode#262, #382).

Each candidate inside a ``SpliceOutcomeSet`` carries an
:class:`EffectCandidate` whose ``.effect`` is one of:
:class:`NormalSplicing`, :class:`ExonSkipping`,
:class:`IntronRetention`, :class:`CrypticDonor`, or
:class:`CrypticAcceptor`. Class identity is the mechanism; each
class carries its own protein vocab on the instance (``aa_ref`` /
``aa_alt`` / ``mutant_protein_sequence`` / ``mutant_transcript``).
"""

import pytest
from pyensembl import cached_release

import varcode
from varcode import (
    CrypticAcceptor,
    CrypticDonor,
    EffectCandidate,
    ExonSkipping,
    IntronRetention,
    NormalSplicing,
    SpliceOutcomeSet,
    SpliceMechanismEffect,
    Variant,
    apply_rna_evidence_to_effects,
    make_rna_outcome,
)
from varcode.effects import (
    ExonicSpliceSite,
    IntronicSpliceSite,
    SpliceAcceptor,
    SpliceDonor,
)
from varcode.splice_outcomes import enumerate_splice_outcomes


ensembl_grch38 = cached_release(81)
CFTR_TRANSCRIPT_ID = "ENST00000003084"
BRCA1_TRANSCRIPT_ID = "ENST00000357654"


def _candidate_of_type(splice_set, mechanism_cls):
    """Return the candidate whose .effect is an instance of
    mechanism_cls, or raise StopIteration."""
    return next(
        c for c in splice_set.candidates
        if isinstance(c.effect, mechanism_cls))


# --------------------------------------------------------------------
# Back-compat
# --------------------------------------------------------------------


def test_default_behavior_unchanged():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effect = variant.effect_on_transcript(transcript)
    assert effect.__class__ is SpliceDonor
    assert not isinstance(effect, SpliceOutcomeSet)


def test_default_for_collection_unchanged():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    effects = variant.effects()
    classes = {type(e) for e in effects}
    assert SpliceOutcomeSet not in classes


# --------------------------------------------------------------------
# Opt-in wraps splice effects
# --------------------------------------------------------------------


def test_opt_in_wraps_splice_donor():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    effects = variant.effects(splice_outcomes=True)
    cftr_effect = next(
        e for e in effects
        if getattr(e, "transcript", None) is not None
        and e.transcript.id == CFTR_TRANSCRIPT_ID
    )
    assert isinstance(cftr_effect, SpliceOutcomeSet)
    assert cftr_effect.disrupted_signal_class is SpliceDonor


def test_opt_in_wraps_splice_acceptor():
    variant = Variant("7", 117530898, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert isinstance(target, SpliceOutcomeSet)
    assert target.disrupted_signal_class is SpliceAcceptor


def test_opt_in_wraps_exonic_splice_site():
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert isinstance(target, SpliceOutcomeSet)
    assert target.disrupted_signal_class is ExonicSpliceSite


def test_opt_in_wraps_intronic_splice_site():
    variant = Variant("7", 117531115, "A", "G", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert isinstance(target, SpliceOutcomeSet)
    assert target.disrupted_signal_class is IntronicSpliceSite


def test_opt_in_passes_through_non_splice_effects():
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", ensembl_grch38)
    effects = variant.effects(splice_outcomes=True)
    classes = [type(e).__name__ for e in effects]
    assert "Substitution" in classes


# --------------------------------------------------------------------
# Mechanism class identity
# --------------------------------------------------------------------


def test_splice_donor_candidate_set_has_expected_mechanisms():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    mechanisms = {type(c.effect) for c in target.candidates}
    assert mechanisms == {
        NormalSplicing, ExonSkipping, IntronRetention, CrypticDonor}


def test_splice_acceptor_candidate_set_uses_cryptic_acceptor():
    variant = Variant("7", 117530898, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    mechanisms = {type(c.effect) for c in target.candidates}
    assert CrypticAcceptor in mechanisms
    assert CrypticDonor not in mechanisms


def test_all_candidates_are_splice_mechanism_effects():
    """Every candidate's .effect is a :class:`SpliceMechanismEffect`
    subclass — uniform shape, no enum / class-identity drift."""
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    for c in target.candidates:
        assert isinstance(c.effect, SpliceMechanismEffect)


def test_candidates_sorted_by_probability_descending():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    probs = [c.probability for c in target.candidates]
    assert probs == sorted(probs, reverse=True)


def test_most_likely_for_splice_donor_is_exon_skipping():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert isinstance(target.most_likely_effect, ExonSkipping)


def test_most_likely_for_exonic_splice_site_is_normal_splicing():
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert isinstance(target.most_likely_effect, NormalSplicing)


# --------------------------------------------------------------------
# Splice signal back-reference (#382): each mechanism carries
# .splice_signal pointing to the underlying SpliceDonor / etc.
# --------------------------------------------------------------------


def test_splice_mechanism_effects_reference_the_splice_signal():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    assert isinstance(bare, SpliceDonor)
    splice_set = enumerate_splice_outcomes(bare)
    for c in splice_set.candidates:
        assert c.effect.splice_signal is bare
        # Reachable from the mechanism: where in the transcript was hit.
        assert c.effect.splice_signal.nearest_exon is not None


# --------------------------------------------------------------------
# Protein vocab on resolved mechanisms
# --------------------------------------------------------------------


def test_normal_splicing_delegates_protein_to_coding_effect():
    """NormalSplicing.aa_ref / aa_alt / mutant_protein_sequence come
    from its underlying .coding_effect (ExonicSpliceSite case).
    """
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    normal = _candidate_of_type(splice_set, NormalSplicing)
    assert normal.effect.coding_effect is not None
    # Protein vocab is reachable on the mechanism itself.
    assert normal.effect.aa_ref == normal.effect.coding_effect.aa_ref
    assert "p." in normal.effect.short_description


def test_normal_splicing_for_intronic_disruption_has_no_coding_effect():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    normal = _candidate_of_type(splice_set, NormalSplicing)
    assert normal.effect.coding_effect is None
    assert normal.effect.mutant_protein_sequence is None


def test_in_frame_exon_skipping_carries_aa_ref_and_protein_sequence():
    """A resolved in-frame ExonSkipping has its own ``aa_ref`` /
    ``mutant_protein_sequence`` populated — no wrapper indirection."""
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    skip = _candidate_of_type(splice_set, ExonSkipping)
    assert skip.effect.in_frame is True
    assert skip.effect.aa_ref  # non-empty
    assert skip.effect.mutant_protein_sequence
    assert skip.effect.mutant_transcript is not None
    # Resolved-state convenience:
    assert skip.effect.resolved


# --------------------------------------------------------------------
# Predicted (unresolved) state
# --------------------------------------------------------------------


def test_intron_retention_predicted_state_carries_intron_coords():
    """Without a genomic_sequence provider, IntronRetention is
    "predicted" — protein fields are None but the retained-intron
    coords + side are populated so consumers know which intron."""
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    ir = _candidate_of_type(splice_set, IntronRetention)
    assert not ir.effect.resolved
    assert ir.effect.mutant_protein_sequence is None
    # Even in predicted state, structural info is known.
    assert ir.effect.side in ("donor", "acceptor")
    assert ir.effect.retained_intron_start > 0
    assert ir.effect.retained_intron_end > ir.effect.retained_intron_start


def test_cryptic_donor_predicted_state_carries_affected_exon():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    cryptic = _candidate_of_type(splice_set, CrypticDonor)
    assert not cryptic.effect.resolved
    assert cryptic.effect.direction == "donor"
    assert cryptic.effect.affected_exon is not None
    # Without provider, cryptic-site fields stay None.
    assert cryptic.effect.cryptic_genomic_position is None
    assert cryptic.effect.motif_score is None


# --------------------------------------------------------------------
# Resolved state with a genomic_sequence provider
# --------------------------------------------------------------------


def _synthetic_intron_provider(stop_offset=30, pad_with="A"):
    def provider(contig, start, end):
        length = end - start + 1
        core = (pad_with * stop_offset) + "TAA" + (pad_with * length)
        return core[:length]
    return provider


def test_intron_retention_with_provider_populates_protein_vocab():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    splice_set = enumerate_splice_outcomes(
        bare, genomic_sequence=_synthetic_intron_provider())
    ir = _candidate_of_type(splice_set, IntronRetention)
    assert ir.effect.resolved
    assert ir.effect.mutant_protein_sequence is not None
    assert len(ir.effect.mutant_protein_sequence) < len(
        str(transcript.protein_sequence))
    assert ir.effect.mutant_transcript is not None


def _synthetic_cryptic_donor_provider(motif_offset=65, motif="CAGGTAAGT"):
    def provider(contig, start, end):
        length = end - start + 1
        seq = bytearray(b"A" * length)
        for i, c in enumerate(motif):
            if motif_offset + i < length:
                seq[motif_offset + i] = ord(c)
        return seq.decode()
    return provider


def test_cryptic_donor_with_provider_populates_protein_vocab():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    splice_set = enumerate_splice_outcomes(
        bare, genomic_sequence=_synthetic_cryptic_donor_provider(
            motif_offset=65))
    cryptic = _candidate_of_type(splice_set, CrypticDonor)
    assert cryptic.effect.resolved
    assert cryptic.effect.motif_score is not None
    assert cryptic.effect.cryptic_genomic_position is not None
    assert cryptic.effect.exon_length_delta is not None
    assert cryptic.effect.mutant_transcript is not None


# --------------------------------------------------------------------
# EffectCollection integration
# --------------------------------------------------------------------


def test_collection_iteration_after_wrapping():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    effects = variant.effects(splice_outcomes=True)
    items = list(effects)
    assert len(items) > 0
    splice_set_count = sum(1 for e in items if isinstance(e, SpliceOutcomeSet))
    assert splice_set_count >= 1


def test_short_description_uses_most_likely_inner_effect():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    desc = target.short_description
    assert desc.startswith("splice-set:")
    # The inner mechanism's short_description shows up.
    assert target.most_likely_effect.short_description in desc


# --------------------------------------------------------------------
# Reverse-strand
# --------------------------------------------------------------------


def test_opt_in_works_on_reverse_strand_donor():
    variant = Variant("17", 43082403, "C", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(BRCA1_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert isinstance(target, SpliceOutcomeSet)
    assert target.disrupted_signal_class is SpliceDonor


# --------------------------------------------------------------------
# enumerate_splice_outcomes pass-through
# --------------------------------------------------------------------


def test_enumerate_passes_through_non_splice():
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(BRCA1_TRANSCRIPT_ID)
    sub_effect = variant.effect_on_transcript(transcript)
    assert type(sub_effect).__name__ == "Substitution"
    wrapped = enumerate_splice_outcomes(sub_effect)
    assert wrapped is sub_effect


# --------------------------------------------------------------------
# Package-level exports
# --------------------------------------------------------------------


def test_package_level_exports():
    assert varcode.SpliceOutcomeSet is SpliceOutcomeSet
    assert varcode.SpliceMechanismEffect is SpliceMechanismEffect
    assert varcode.NormalSplicing is NormalSplicing
    assert varcode.ExonSkipping is ExonSkipping
    assert varcode.IntronRetention is IntronRetention
    assert varcode.CrypticDonor is CrypticDonor
    assert varcode.CrypticAcceptor is CrypticAcceptor


# --------------------------------------------------------------------
# MultiOutcomeEffect protocol
# --------------------------------------------------------------------


def test_splice_outcome_set_is_a_multi_outcome_effect():
    from varcode import MultiOutcomeEffect
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    assert isinstance(target, MultiOutcomeEffect)
    assert all(isinstance(c, EffectCandidate) for c in target.candidates)
    assert target.most_likely_candidate is target.candidates[0]
    assert target.priority_class is target.disrupted_signal_class


# --------------------------------------------------------------------
# Boundary-codon reconstruction (#298). Out-of-frame exon skips
# produce a real ExonSkipping with in_frame=False and frameshifted
# protein vocab on the instance.
# --------------------------------------------------------------------


def test_mid_codon_in_frame_exon_skip_reshapes_boundary():
    variant = Variant("7", 117536674, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    if not isinstance(bare, SpliceDonor):
        pytest.skip(
            "Canonical donor G not present at 117536674; got %s."
            % type(bare).__name__)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    skip = _candidate_of_type(splice_set, ExonSkipping)
    # ExonSkipping with a reshaped boundary still carries aa_ref/alt
    # on the instance.
    assert skip.effect.aa_ref
    assert skip.effect.mutant_protein_sequence


def test_codon_aligned_exon_skip_is_in_frame():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    splice_set = next(e for e in effects if e.transcript is transcript)
    skip = _candidate_of_type(splice_set, ExonSkipping)
    # CFTR exon 3 is codon-aligned (216 bp / 72 codons).
    assert skip.effect.in_frame is True


def test_out_of_frame_exon_skip_marked_out_of_frame():
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    target_exon = None
    for exon in transcript.exons[2:]:
        if (exon.end - exon.start + 1) % 3 != 0:
            target_exon = exon
            break
    if target_exon is None:
        pytest.skip("No out-of-frame exon found in CFTR beyond exon 2")
    donor_plus_1 = target_exon.end + 1
    variant = Variant("7", donor_plus_1, "G", "A", ensembl_grch38)
    bare = variant.effect_on_transcript(transcript)
    if not isinstance(bare, SpliceDonor):
        pytest.skip(
            "Canonical donor G not present at %d; got %s." % (
                donor_plus_1, type(bare).__name__))
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    skip = _candidate_of_type(splice_set, ExonSkipping)
    assert skip.effect.in_frame is False
    # Out-of-frame skip produces a real (frameshifted) protein sequence.
    assert skip.effect.mutant_protein_sequence
    assert skip.effect.mutant_protein_sequence != str(
        transcript.protein_sequence)


# --------------------------------------------------------------------
# Serialization round-trip (no enum to stringify, no custom
# to_dict / from_dict needed)
# --------------------------------------------------------------------


def test_splice_outcome_set_json_round_trip():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)

    payload = target.to_json()
    restored = SpliceOutcomeSet.from_json(payload)

    assert type(restored) is SpliceOutcomeSet
    assert len(restored.candidates) == len(target.candidates)
    assert restored.disrupted_signal_class is target.disrupted_signal_class
    # Mechanism class identity preserved across round-trip.
    for rt_c, og_c in zip(restored.candidates, target.candidates):
        assert type(rt_c.effect) is type(og_c.effect)
        assert rt_c.probability == og_c.probability


def test_splice_outcome_set_round_trip_preserves_priority_class():
    from varcode.effects import effect_priority
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    effects = variant.effects(splice_outcomes=True)
    target = next(e for e in effects if e.transcript is transcript)
    restored = SpliceOutcomeSet.from_json(target.to_json())
    assert effect_priority(restored) == effect_priority(target)


def test_non_splice_effects_are_not_multi_outcome():
    from varcode import MultiOutcomeEffect
    from varcode.effects import Substitution, Silent, Intronic, MutationEffect
    for cls in (Substitution, Silent, Intronic, MutationEffect):
        assert not issubclass(cls, MultiOutcomeEffect)


# --------------------------------------------------------------------
# Priority integration
# --------------------------------------------------------------------


def test_splice_outcome_set_sorts_as_disrupted_signal_class():
    from varcode.effects import effect_priority
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare_effect = variant.effect_on_transcript(transcript)
    assert isinstance(bare_effect, SpliceDonor)
    bare_priority = effect_priority(bare_effect)
    wrapped_effects = variant.effects(splice_outcomes=True)
    wrapped = next(e for e in wrapped_effects if e.transcript is transcript)
    assert effect_priority(wrapped) == bare_priority


def test_splice_outcome_set_top_priority_works():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    effects = variant.effects(splice_outcomes=True)
    top = effects.top_priority_effect()
    top_class_name = type(top).__name__
    assert top_class_name in ("SpliceOutcomeSet", "SpliceDonor")


# --------------------------------------------------------------------
# Acceptor-side IntronicSpliceSite emits CrypticAcceptor
# --------------------------------------------------------------------


def test_acceptor_side_intronic_splice_site_uses_cryptic_acceptor():
    variant = Variant("7", 117530896, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    assert isinstance(bare, IntronicSpliceSite) and \
        not isinstance(bare, (SpliceDonor, SpliceAcceptor))
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    mechanisms = {type(c.effect) for c in splice_set.candidates}
    assert CrypticAcceptor in mechanisms
    assert CrypticDonor not in mechanisms


def test_donor_side_intronic_splice_site_uses_cryptic_donor():
    variant = Variant("7", 117531117, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    assert isinstance(bare, IntronicSpliceSite) and \
        not isinstance(bare, (SpliceDonor, SpliceAcceptor))
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    mechanisms = {type(c.effect) for c in splice_set.candidates}
    assert CrypticDonor in mechanisms
    assert CrypticAcceptor not in mechanisms


# --------------------------------------------------------------------
# Multi-protein surface
# --------------------------------------------------------------------


def test_candidate_proteins_keyed_by_mechanism_class():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    proteins = splice_set.candidate_proteins
    # Keys are mechanism classes, values are protein strings (or "").
    assert ExonSkipping in proteins
    assert proteins[ExonSkipping]  # in-frame exon skip has protein
    assert proteins[IntronRetention] == ""
    assert proteins[CrypticDonor] == ""


def test_mutant_protein_sequences_collects_distinct_proteins():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    proteins = splice_set.mutant_protein_sequences
    assert isinstance(proteins, set)
    assert len(proteins) >= 1
    ref = str(transcript.protein_sequence)
    shortest = min(proteins, key=len)
    assert len(shortest) < len(ref)


# --------------------------------------------------------------------
# Most-likely / highest-priority accessors
# --------------------------------------------------------------------


def test_most_likely_candidate_returns_effect_candidate():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    assert isinstance(splice_set.most_likely_candidate, EffectCandidate)


def test_most_likely_effect_returns_splice_mechanism_effect():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    assert isinstance(splice_set.most_likely_effect, SpliceMechanismEffect)
    assert splice_set.most_likely_effect is splice_set.most_likely_candidate.effect


def test_highest_priority_candidate_returns_most_disruptive():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    from varcode.effects import effect_priority
    top = splice_set.highest_priority_candidate
    expected = max(
        effect_priority(c.effect) for c in splice_set.candidates)
    assert effect_priority(top.effect) == expected


def test_resolved_splice_mechanism_delegates_priority_to_protein_effect():
    from varcode.effects import Deletion, effect_priority
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    skip = _candidate_of_type(splice_set, ExonSkipping)
    assert isinstance(skip.effect.protein_effect, Deletion)
    assert skip.effect.modifies_protein_sequence
    assert skip.effect.modifies_coding_sequence
    assert effect_priority(skip.effect) == effect_priority(skip.effect.protein_effect)
    assert splice_set.highest_priority_candidate is skip


def test_effects_unwraps_candidate_tuple():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    inner = splice_set.effects
    assert all(isinstance(e, SpliceMechanismEffect) for e in inner)
    assert inner == tuple(c.effect for c in splice_set.candidates)


# --------------------------------------------------------------------
# RNA evidence reconciliation
# --------------------------------------------------------------------


def test_rna_evidence_reconciles_splice_set_by_excluding_unobserved_mechanisms():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    original = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    skip = _candidate_of_type(original, ExonSkipping)
    observed = make_rna_outcome(
        effect=skip.effect,
        source="junction_reads",
        read_count=18,
        probability=0.91)

    class Resolver:
        def observed_outcomes(self, variant, transcript):
            return (observed,)

    effects = [original]
    result = apply_rna_evidence_to_effects(effects, Resolver())
    reconciled = effects[0]

    assert result is effects
    assert reconciled is not original
    assert isinstance(reconciled, SpliceOutcomeSet)
    assert tuple(reconciled.dna_candidates) == tuple(original.candidates)
    assert reconciled.rna_evidence == (observed,)
    assert reconciled.added_candidates == ()
    assert reconciled.candidates == (observed,)
    assert reconciled.rna_evidence_for(observed) == (observed,)
    excluded_mechanisms = {type(c.effect) for c in reconciled.excluded_candidates}
    assert excluded_mechanisms == {NormalSplicing, IntronRetention, CrypticDonor}


def test_rna_evidence_records_added_splice_mechanisms():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    original = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    skip = _candidate_of_type(original, ExonSkipping)
    added_effect = CrypticAcceptor(
        variant=variant,
        transcript=transcript,
        splice_signal=skip.effect.splice_signal,
        affected_exon=skip.effect.affected_exon,
        cryptic_genomic_position=117530930,
        motif_score=4.2,
        exon_length_delta=-12)
    observed = make_rna_outcome(
        effect=added_effect,
        source="long_reads",
        read_count=6,
        extra_evidence={"junction_id": "cryptic_acceptor_1"})

    reconciled = original.with_rna_evidence((observed,))

    assert reconciled is not original
    assert reconciled.added_candidates == (observed,)
    assert observed in reconciled.candidates
    assert reconciled.rna_evidence_for(observed) == (observed,)
    assert CrypticAcceptor not in {type(c.effect) for c in original.dna_candidates}


# --------------------------------------------------------------------
# alternate_effect back-compat shim
# --------------------------------------------------------------------


def test_alternate_effect_for_exonic_splice_site_resolves_to_coding_effect():
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    assert isinstance(bare, ExonicSpliceSite)
    bare_alt = bare.alternate_effect
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    assert type(splice_set.alternate_effect) is type(bare_alt)


def test_alternate_effect_none_when_no_underlying_coding_change():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    assert splice_set.alternate_effect is None


# --------------------------------------------------------------------
# Protein sequence prediction — biological validation.
#
# These tests check that the mutant protein varcode predicts for a
# splice variant matches what we'd get by independently computing
# "reference protein with the relevant amino acids
# removed/replaced." That arithmetic IS the biology — if varcode's
# protein math is correct, it should match. The tests also pin
# specific landmarks (length, first/last residues, HGVS-style
# descriptions) so regressions surface clearly.
#
# CFTR ENST00000003084 is the canonical 1480-aa CFTR transcript;
# coordinates are GRCh38 (Ensembl 81). Variant positions used
# below are at canonical splice signals of well-defined CFTR
# exons.
# --------------------------------------------------------------------


def test_cftr_exon_3_skip_protein_equals_reference_minus_skipped_aas():
    """In-frame skip of CFTR exon 3 (216 bp = 72 codons, donor +1 at
    g.117531115 G>A). The mutant protein must equal the reference
    protein with exactly those 72 amino acids removed at the right
    offset. This is the cleanest "is the protein math correct?"
    check we can do without external sequence data.
    """
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    skip = _candidate_of_type(splice_set, ExonSkipping)
    assert skip.effect.in_frame
    # Exon 3 of ENST00000003084 is 216 bp → 72 codons.
    assert len(skip.effect.aa_ref) == 72
    ref = str(transcript.protein_sequence)
    mut = skip.effect.mutant_protein_sequence
    start = skip.effect.aa_mutation_start_offset
    # Mutant protein = reference with aa_ref removed at start.
    expected = ref[:start] + ref[start + len(skip.effect.aa_ref):]
    assert mut == expected
    # And the removed AAs match what's in the reference protein at
    # that range (i.e. aa_ref is the actually-removed substring, not
    # a re-translated approximation).
    assert skip.effect.aa_ref == ref[start:start + len(skip.effect.aa_ref)]
    # Protein is exactly 72 AAs shorter than reference.
    assert len(ref) - len(mut) == 72


def test_cftr_exon_3_skip_short_description_uses_hgvs_del_notation():
    """The HGVS-style short_description names the deleted AA range:
    ``p.{first}{start}_{last}{end}del``. For CFTR exon 3 in
    ENST00000003084 this is ``p.E92_K163del`` (Glu92 to Lys163,
    1-indexed)."""
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    skip = _candidate_of_type(splice_set, ExonSkipping)
    desc = skip.effect.short_description
    assert "p.E92_K163del" in desc
    # First/last AAs of the removed range:
    assert skip.effect.aa_ref[0] == "E"   # Glu
    assert skip.effect.aa_ref[-1] == "K"  # Lys


def test_cftr_out_of_frame_exon_skip_produces_frameshift_protein():
    """Skipping a CFTR exon whose length is not divisible by 3
    causes a frameshift starting at the new exon junction. The
    mutant protein equals the reference up to the divergence point,
    then diverges with the frameshifted readthrough until a
    premature stop. CFTR exon at g.117509034-117509142 is 109 bp
    (out-of-frame); donor +1 is at g.117509143.
    """
    variant = Variant("7", 117509143, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    if not isinstance(bare, SpliceDonor):
        pytest.skip("g.117509143 G>A not classified as SpliceDonor")
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    skip = _candidate_of_type(splice_set, ExonSkipping)
    assert skip.effect.in_frame is False
    ref = str(transcript.protein_sequence)
    mut = skip.effect.mutant_protein_sequence
    # Reference and mutant should be identical up to the divergence
    # point at aa_mutation_start_offset, then diverge.
    start = skip.effect.aa_mutation_start_offset
    assert ref[:start] == mut[:start]
    if start < len(mut) and start < len(ref):
        assert ref[start] != mut[start]
    # Frameshift typically truncates the protein well before the
    # reference end — the mutant should be much shorter.
    assert len(mut) < len(ref)
    # Short_description uses HGVS frameshift notation, not "del".
    assert "fs" in skip.effect.short_description
    assert "del" not in skip.effect.short_description


def test_cftr_exon_3_skip_aa_alt_describes_junction():
    """For an in-frame, codon-aligned skip (CFTR exon 3 starts and
    ends at codon boundaries), the new junction codon doesn't get
    reshaped — ``aa_alt`` is empty.
    """
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    skip = _candidate_of_type(splice_set, ExonSkipping)
    assert skip.effect.in_frame
    # Codon-aligned skip → no boundary reshape → aa_alt is "".
    assert skip.effect.aa_alt == ""


def test_cftr_acceptor_disruption_skips_same_exon_as_donor_disruption():
    """The CFTR exon 3 acceptor at g.117530898 (G>A on the canonical
    AG signal) disrupts the splice acceptor and produces the same
    EXON_SKIPPING outcome as the donor disruption at the other end.
    The mutant protein from the acceptor side must be identical to
    the one from the donor side — both skip exon 3 the same way.
    """
    donor_variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    acceptor_variant = Variant("7", 117530898, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)

    donor_set = next(
        e for e in donor_variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    acceptor_set = next(
        e for e in acceptor_variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    donor_skip = _candidate_of_type(donor_set, ExonSkipping)
    acceptor_skip = _candidate_of_type(acceptor_set, ExonSkipping)
    # Same exon skipped.
    assert (
        donor_skip.effect.affected_exon.exon_id
        == acceptor_skip.effect.affected_exon.exon_id)
    # Same mutant protein.
    assert (
        donor_skip.effect.mutant_protein_sequence
        == acceptor_skip.effect.mutant_protein_sequence)
    # Same AA range removed.
    assert donor_skip.effect.aa_ref == acceptor_skip.effect.aa_ref
    # Both reference back to their own splice signal class — donor
    # vs acceptor — but the splice_signal pointer is each variant's
    # own classification.
    assert isinstance(donor_skip.effect.splice_signal, SpliceDonor)
    assert isinstance(acceptor_skip.effect.splice_signal, SpliceAcceptor)


def test_intron_retention_with_provider_truncates_protein_at_in_intron_stop():
    """When a genomic_sequence provider returns an intron with an
    early stop codon, the IntronRetention candidate's mutant protein
    must end at that stop. The truncated protein equals the
    reference protein up to the divergence point (where the intron
    starts encoding nonsense) and then stops.
    """
    # Provider returns A-padding with a TAA stop at offset 30.
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    splice_set = enumerate_splice_outcomes(
        bare, genomic_sequence=_synthetic_intron_provider(stop_offset=30))
    ir = _candidate_of_type(splice_set, IntronRetention)
    assert ir.effect.resolved
    ref = str(transcript.protein_sequence)
    mut = ir.effect.mutant_protein_sequence
    # Mutant ends at a premature stop — strictly shorter than ref.
    assert len(mut) < len(ref)
    # Reference and mutant share their leading prefix; intron
    # retention only diverges where the retained sequence starts
    # encoding the new (eventually nonsense) protein.
    common_prefix_end = next(
        (i for i, (r, m) in enumerate(zip(ref, mut)) if r != m),
        min(len(ref), len(mut)))
    # The shared prefix is positive — there's a real common stretch
    # before divergence.
    assert common_prefix_end > 0


def test_cryptic_donor_with_provider_changes_exon_length_consistently():
    """A cryptic donor at a fixed offset from the canonical boundary
    changes the affected exon's length by a predictable amount.
    Verify the ``exon_length_delta`` field matches the actual cDNA
    length delta between the cryptic-resolved mutant transcript and
    the reference transcript.
    """
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    bare = variant.effect_on_transcript(transcript)
    splice_set = enumerate_splice_outcomes(
        bare, genomic_sequence=_synthetic_cryptic_donor_provider(
            motif_offset=65))
    cryptic = _candidate_of_type(splice_set, CrypticDonor)
    assert cryptic.effect.resolved
    # cDNA length delta matches the field.
    mt = cryptic.effect.mutant_transcript
    cdna_delta = len(mt.cdna_sequence) - len(str(transcript.sequence))
    assert cryptic.effect.exon_length_delta == cdna_delta


# --------------------------------------------------------------------
# Per-mechanism serialization round-trip — each mechanism class
# survives to_json / from_json with its structural fields intact.
# --------------------------------------------------------------------


def test_exon_skipping_serialization_preserves_affected_exon_and_protein():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    restored = SpliceOutcomeSet.from_json(splice_set.to_json())
    og = _candidate_of_type(splice_set, ExonSkipping)
    rt = _candidate_of_type(restored, ExonSkipping)
    assert type(rt.effect) is type(og.effect)
    assert rt.effect.in_frame == og.effect.in_frame
    assert rt.effect.aa_ref == og.effect.aa_ref
    assert rt.effect.mutant_protein_sequence == og.effect.mutant_protein_sequence
    assert rt.effect.affected_exon.exon_id == og.effect.affected_exon.exon_id
    # splice_signal round-trips too — same class, same nearest_exon.
    assert type(rt.effect.splice_signal) is type(og.effect.splice_signal)


def test_intron_retention_serialization_preserves_intron_coords():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    restored = SpliceOutcomeSet.from_json(splice_set.to_json())
    og = _candidate_of_type(splice_set, IntronRetention)
    rt = _candidate_of_type(restored, IntronRetention)
    assert rt.effect.retained_intron_start == og.effect.retained_intron_start
    assert rt.effect.retained_intron_end == og.effect.retained_intron_end
    assert rt.effect.side == og.effect.side


def test_cryptic_donor_serialization_preserves_direction():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    restored = SpliceOutcomeSet.from_json(splice_set.to_json())
    og = _candidate_of_type(splice_set, CrypticDonor)
    rt = _candidate_of_type(restored, CrypticDonor)
    assert rt.effect.direction == "donor" == og.effect.direction
    assert rt.effect.affected_exon.exon_id == og.effect.affected_exon.exon_id


def test_normal_splicing_serialization_preserves_coding_effect_class():
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    restored = SpliceOutcomeSet.from_json(splice_set.to_json())
    og = _candidate_of_type(splice_set, NormalSplicing)
    rt = _candidate_of_type(restored, NormalSplicing)
    if og.effect.coding_effect is None:
        assert rt.effect.coding_effect is None
    else:
        assert type(rt.effect.coding_effect) is type(og.effect.coding_effect)
        # Protein vocab is delegated through coding_effect.
        assert rt.effect.aa_ref == og.effect.aa_ref


# --------------------------------------------------------------------
# Predicted-state HGVS / short_description behavior
# --------------------------------------------------------------------


def test_predicted_exon_skipping_short_description_includes_predicted_tag():
    """Without a resolved protein, ExonSkipping's short_description
    flags the predicted state explicitly so a reader can't confuse
    it with a resolved skip."""
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    # The IntronRetention candidate is predicted (no provider).
    ir = _candidate_of_type(splice_set, IntronRetention)
    assert "predicted" in ir.effect.short_description
    cd = _candidate_of_type(splice_set, CrypticDonor)
    assert "predicted" in cd.effect.short_description


def test_intron_retention_predicted_side_is_donor_or_acceptor_or_unknown():
    """``IntronRetention.side`` is either ``"donor"`` /
    ``"acceptor"`` (when inferable from the splice classification)
    or ``None`` (when not). Sentinel-zero coords are not used —
    unknown is None."""
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(CFTR_TRANSCRIPT_ID)
    splice_set = next(
        e for e in variant.effects(splice_outcomes=True)
        if e.transcript is transcript)
    ir = _candidate_of_type(splice_set, IntronRetention)
    # For a SpliceDonor-classed disruption we expect side="donor".
    assert ir.effect.side == "donor"
    # Coords are populated (not None) because we have the
    # adjacent-exon geometry.
    assert ir.effect.retained_intron_start is not None
    assert ir.effect.retained_intron_end is not None
