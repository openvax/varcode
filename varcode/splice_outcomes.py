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
Multi-mechanism candidate model for splice-disrupting variants — see
openvax/varcode#262, #382.

When a variant disrupts a splice signal, the downstream protein
outcome is not deterministic from DNA alone. This module wraps the
splice classification in a :class:`SpliceOutcomeSet` carrying one
:class:`EffectCandidate` per plausible mechanism. Each candidate's
``effect`` is a :class:`SpliceMechanismEffect` subclass —
:class:`NormalSplicing`, :class:`ExonSkipping`,
:class:`IntronRetention`, :class:`CrypticDonor`, or
:class:`CrypticAcceptor` — describing what the splice machinery
does and (when computable) the resulting protein change.

Class identity = mechanism. Each mechanism Effect carries its own
protein vocab as instance fields (``aa_ref``, ``aa_alt``,
``mutant_protein_sequence``, ``mutant_transcript``), with ``None``
denoting "predicted mechanism, exact protein not computable from
cached cDNA alone." No parallel enum, no placeholder class
hierarchy.

Limitations of this prototype (documented openly so callers know what
the candidate ordering means):

* Varcode-generated splice candidates are ordered by a hard-coded,
  DNA-only mechanism preference. The ordering is intentionally not
  exported as a score or probability.
* :class:`ExonSkipping` and :class:`NormalSplicing` produce concrete
  protein sequences using the cDNA already available from PyEnsembl.
* :class:`IntronRetention` and the cryptic-splice mechanisms
  (:class:`CrypticDonor` / :class:`CrypticAcceptor`) need intron /
  flanking genomic sequence to resolve the protein. Without a
  ``genomic_sequence`` provider they're returned in their
  "predicted" state with protein fields ``None``; consumers can
  still read ``effect.affected_exon`` / ``effect.side`` /
  ``effect.retained_intron_*`` etc. to know what mechanism is
  predicted.
* Candidates whose protein math resolved expose the
  :class:`~varcode.MutantTranscript` they were derived from via
  ``candidate.effect.mutant_transcript``. Unresolved candidates
  carry ``mutant_transcript=None``.
"""
from dataclasses import dataclass, field
from typing import Tuple

from serializable import DataclassSerializable

from .effect_candidates import EffectCandidate
from .effects.effect_classes import (
    CrypticAcceptor,
    CrypticDonor,
    ExonSkipping,
    ExonicSpliceSite,
    IntronRetention,
    MultiOutcomeEffect,
    NormalSplicing,
    SpliceMechanismEffect,
    SpliceSite,
    TranscriptMutationEffect,
)
from .splice_exon_builders import (
    build_exon_skipping_candidate,
    build_normal_splicing_candidate,
    mechanism_order_for,
)
from .splice_intron_builders import (
    build_cryptic_splice_candidate,
    build_intron_retention_candidate,
)


@dataclass(frozen=True)
class SpliceCandidateRNAEvidence(DataclassSerializable):
    """RNA observations supporting one visible splice candidate."""

    candidate: EffectCandidate
    rna_evidence: Tuple[EffectCandidate, ...] = field(default_factory=tuple)


def _exon_signature(exon):
    """Stable value key for exon identity in RNA-evidence audits."""
    if exon is None:
        return None
    return (
        getattr(exon, "exon_id", None),
        getattr(exon, "contig", None),
        getattr(exon, "start", None),
        getattr(exon, "end", None),
        getattr(exon, "strand", None),
    )


def _splice_candidate_signature(candidate):
    """Candidate-specific key for grouping RNA support.

    Mechanism class alone is too coarse: two distinct cryptic donors
    have the same class but different junctions. This key captures the
    mechanism-specific coordinates that make one observed splice
    candidate distinct from another, while ignoring source/evidence so
    multiple observations can support the same visible candidate.
    """
    effect = candidate.effect
    return (
        type(effect),
        _exon_signature(getattr(effect, "affected_exon", None)),
        getattr(effect, "in_frame", None),
        getattr(effect, "retained_intron_start", None),
        getattr(effect, "retained_intron_end", None),
        getattr(effect, "side", None),
        getattr(effect, "cryptic_genomic_position", None),
        getattr(effect, "exon_length_delta", None),
        getattr(effect, "aa_mutation_start_offset", None),
        getattr(effect, "aa_mutation_end_offset", None),
        getattr(effect, "aa_ref", None),
        getattr(effect, "aa_alt", None),
        getattr(effect, "mutant_protein_sequence", None),
    )


class SpliceOutcomeSet(MultiOutcomeEffect, TranscriptMutationEffect):
    """A splice-disrupting variant's effect, expressed as a set of
    plausible mechanisms rather than a single Effect.

    Implements the :class:`MultiOutcomeEffect` protocol: downstream
    consumers can write ``isinstance(effect, MultiOutcomeEffect)`` to
    catch this and any future multi-outcome wrapper (RNA evidence,
    germline-aware, etc.) uniformly.

    This set holds two distinct things — keep them straight:

    * **The disruption** — *where* the splice signal was hit. Stored
      as :attr:`disrupted_signal_class`, which is the
      :class:`~varcode.effects.effect_classes.SpliceSite` **subclass**
      (a *type*, e.g. ``SpliceDonor`` — **not** an Effect instance).
      It's metadata used for priority lookup. The corresponding
      Effect *instance* is reachable on every candidate via
      ``candidate.effect.splice_signal``.
    * **The protein outcomes** — *what* the spliceosome might do.
      Stored as :attr:`candidates`, a ``tuple[EffectCandidate, ...]``.
      Each candidate's ``effect`` is a
      :class:`~varcode.effects.effect_classes.SpliceMechanismEffect`
      subclass (:class:`NormalSplicing`, :class:`ExonSkipping`,
      :class:`IntronRetention`, :class:`CrypticDonor`, or
      :class:`CrypticAcceptor`) — class identity is the mechanism.
      These are the effects that carry a protein consequence:
      ``aa_ref`` / ``aa_alt`` / ``mutant_protein_sequence`` /
      ``mutant_transcript``, populated when the protein math resolved,
      ``None`` when only the mechanism is predicted.

    :attr:`short_description` lists every current candidate in producer
    order, including unresolved predictions. Constructed by
    :func:`enumerate_splice_outcomes` as part of every
    ``Variant.effects()`` / ``Variant.effect_on_transcript(...)`` /
    ``VariantCollection.effects()`` call (always-on as of varcode 6.0).

    RNA evidence reconciliation
    ---------------------------
    :meth:`with_rna_evidence` returns a new set rather than mutating this
    one. The new set keeps ``dna_candidates`` for audit, narrows
    ``candidates`` to RNA-observed mechanisms, records DNA predictions
    absent from RNA as ``excluded_candidates``, records observed
    mechanisms absent from the DNA prior as ``added_candidates``, and
    exposes per-candidate support through ``candidate_rna_evidence`` /
    :meth:`rna_evidence_for`.

    Serialization is inherited from :class:`Serializable` — no
    overrides needed now that there's no enum to stringify.
    """

    # Class-level filter shortcuts. A splice disruption is always
    # *potentially* protein-modifying — any non-NormalSplicing candidate
    # (ExonSkipping, IntronRetention, CrypticDonor/Acceptor) can rewrite
    # the protein. Hardcoded so ``drop_silent_and_noncoding()`` doesn't
    # silently drop splice variants whose NormalSplicing candidate
    # happens to carry a Silent coding_effect, and so the filter path
    # doesn't trigger lazy candidate materialization.
    modifies_protein_sequence = True
    modifies_coding_sequence = True

    def __init__(
            self, variant, transcript, candidates=None,
            disrupted_signal_class=None,
            *,
            splice_effect=None,
            genomic_sequence=None,
            mechanism_order=None,
            dna_candidates=None, rna_evidence=(), added_candidates=(),
            excluded_candidates=(), candidate_rna_evidence=()):
        # SpliceOutcomeSet is both a MultiOutcomeEffect (candidate
        # set) and a TranscriptMutationEffect (carries variant, gene,
        # transcript). Initialising via TranscriptMutationEffect sets
        # all three; MultiOutcomeEffect carries no extra state.
        TranscriptMutationEffect.__init__(self, variant, transcript)
        # The SpliceSite *subclass* (a type, not an instance) of the
        # disruption this set replaced — SpliceDonor, SpliceAcceptor,
        # ExonicSpliceSite, or IntronicSpliceSite. Used for priority
        # lookup; the disruption *instance* lives on each candidate's
        # effect.splice_signal.
        self.disrupted_signal_class = disrupted_signal_class
        if candidates is not None:
            # Materialized path — RNA-reconciliation rebuilds,
            # deserialization, tests with explicit candidate lists.
            self._candidates_cache = tuple(candidates)
            self._normal_splicing_eager = next(
                (c for c in self._candidates_cache
                 if isinstance(c.effect, NormalSplicing)),
                None)
            self._lazy_state = None
        else:
            # Lazy path — enumerate_splice_outcomes uses this. The
            # NormalSplicing candidate is computed eagerly because it's
            # cheap (just wraps the existing coding effect) and is the
            # fast path for ``effect_if_splicing_unchanged`` /
            # ``alternate_effect``. Expensive candidates (ExonSkipping,
            # IntronRetention, cryptic) are deferred until ``.candidates``
            # is read.
            if splice_effect is None or mechanism_order is None:
                raise ValueError(
                    "SpliceOutcomeSet needs either explicit ``candidates`` "
                    "or (``splice_effect`` + ``mechanism_order``) for "
                    "lazy construction.")
            self._normal_splicing_eager = (
                build_normal_splicing_candidate(splice_effect))
            self._candidates_cache = None
            self._lazy_state = (splice_effect, mechanism_order, genomic_sequence)
        self._explicit_dna_candidates = (
            None if dna_candidates is None else tuple(dna_candidates))
        self.rna_evidence = tuple(rna_evidence)
        self.added_candidates = tuple(added_candidates)
        self.excluded_candidates = tuple(excluded_candidates)
        self.candidate_rna_evidence = tuple(candidate_rna_evidence)

    @property
    def _candidates(self):
        if self._candidates_cache is None:
            splice_effect, mechanism_order, genomic_sequence = self._lazy_state
            materialized = []
            for mech_cls in mechanism_order:
                if mech_cls is NormalSplicing:
                    materialized.append(self._normal_splicing_eager)
                else:
                    materialized.append(_build_candidate_for(
                        mech_cls, splice_effect, genomic_sequence))
            self._candidates_cache = tuple(materialized)
            self._lazy_state = None
        return self._candidates_cache

    def _get_candidates(self):
        return self._candidates

    @property
    def dna_candidates(self):
        # When not explicitly set, defaults to the DNA-derived candidate
        # set — same as ``_candidates`` here (additional outcomes appended via
        # ``add_candidates`` live on the public ``candidates`` property).
        # Reading this triggers full materialization if constructed lazily.
        if self._explicit_dna_candidates is not None:
            return self._explicit_dna_candidates
        return self._candidates

    def with_rna_evidence(self, observed_candidates):
        """Return a new set reconciled with RNA-observed splice evidence.

        RNA evidence is treated as a refinement of the DNA-predicted
        mechanism set: observed mechanisms replace their broad DNA
        predictions, unobserved DNA-predicted mechanisms are recorded as
        excluded, and observed mechanisms absent from the DNA prediction
        are recorded as added. The original set is left untouched.
        """
        observed = tuple(observed_candidates)
        if not observed:
            return self

        combined_evidence = self.rna_evidence + observed
        splice_evidence = tuple(
            c for c in combined_evidence
            if isinstance(c.effect, SpliceMechanismEffect))
        if not splice_evidence:
            return SpliceOutcomeSet(
                variant=self.variant,
                transcript=self.transcript,
                candidates=self.candidates + observed,
                disrupted_signal_class=self.disrupted_signal_class,
                dna_candidates=self.dna_candidates,
                rna_evidence=combined_evidence,
                added_candidates=observed,
                excluded_candidates=self.excluded_candidates,
                candidate_rna_evidence=self.candidate_rna_evidence,
            )

        dna_by_mechanism = {
            type(candidate.effect): candidate
            for candidate in self.dna_candidates
        }
        evidence_by_mechanism = {}
        evidence_by_signature = {}
        representative_by_signature = {}
        for candidate in splice_evidence:
            mechanism = type(candidate.effect)
            evidence_by_mechanism.setdefault(mechanism, []).append(candidate)
            signature = _splice_candidate_signature(candidate)
            evidence_by_signature.setdefault(signature, []).append(candidate)
            representative_by_signature.setdefault(signature, candidate)

        current = []
        excluded = []
        added = []
        seen_signatures = set()

        def add_current_observed(candidate):
            signature = _splice_candidate_signature(candidate)
            if signature in seen_signatures:
                return None
            seen_signatures.add(signature)
            representative = representative_by_signature[signature]
            current.append(representative)
            return representative

        for dna_candidate in self.dna_candidates:
            mechanism = type(dna_candidate.effect)
            observed_for_mechanism = tuple(
                evidence_by_mechanism.get(mechanism, ()))
            if observed_for_mechanism:
                for candidate in observed_for_mechanism:
                    add_current_observed(candidate)
            else:
                excluded.append(dna_candidate)

        for candidate in splice_evidence:
            mechanism = type(candidate.effect)
            if mechanism not in dna_by_mechanism:
                representative = add_current_observed(candidate)
                if representative is not None:
                    added.append(representative)

        support = tuple(
            SpliceCandidateRNAEvidence(
                candidate=candidate,
                rna_evidence=tuple(evidence_by_signature.get(
                    _splice_candidate_signature(candidate), ())))
            for candidate in current
        )

        return SpliceOutcomeSet(
            variant=self.variant,
            transcript=self.transcript,
            candidates=current,
            disrupted_signal_class=self.disrupted_signal_class,
            dna_candidates=self.dna_candidates,
            rna_evidence=combined_evidence,
            added_candidates=added,
            excluded_candidates=excluded,
            candidate_rna_evidence=support,
        )

    def rna_evidence_for(self, candidate):
        """RNA-observed candidates supporting ``candidate``."""
        candidate_signature = _splice_candidate_signature(candidate)
        for entry in self.candidate_rna_evidence:
            if entry.candidate is candidate or entry.candidate == candidate:
                return entry.rna_evidence
            if _splice_candidate_signature(entry.candidate) == candidate_signature:
                return entry.rna_evidence
        return ()

    @property
    def priority_class(self):
        """Rank by the signal or its higher-priority exonic coding effect.

        This preserves the previous transcript selection without replacing
        the set with the consequence conditional on normal splicing. Only
        the eager NormalSplicing candidate is needed for this lookup.
        """
        if self.disrupted_signal_class is ExonicSpliceSite:
            from .effects.effect_ordering import effect_priority, transcript_effect_priority_dict

            coding = self.effect_if_splicing_unchanged
            signal_priority = transcript_effect_priority_dict[ExonicSpliceSite]
            if coding is not None and effect_priority(coding) > signal_priority:
                return getattr(coding, "priority_class", None) or type(coding)
        return self.disrupted_signal_class

    @property
    def effect_if_splicing_unchanged(self):
        """The coding consequence that applies *if the spliceosome
        still splices normally* despite the disrupted signal — i.e.
        the protein consequence of the underlying nucleotide change
        on its own (a :class:`MutationEffect` such as
        :class:`Substitution`, not a :class:`SpliceMechanismEffect`).

        Resolves to the :class:`NormalSplicing` candidate's
        :attr:`NormalSplicing.coding_effect`. Returns ``None`` when
        there is no :class:`NormalSplicing` candidate or it carries no
        coding consequence (e.g. a purely intronic disruption, where a
        nucleotide change outside the CDS leaves the protein
        untouched if splicing proceeds).

        This is the canonical "alternative outcome" accessor, alongside
        :attr:`most_likely_effect` (the top predicted splicing
        mechanism) and :attr:`candidates` (every plausible mechanism).
        Unlike the old ``ExonicSpliceSite.alternate_effect`` it works
        for intronic disruptions too.
        """
        # Fast path: NormalSplicing is computed eagerly and stored on
        # the instance, so this read never triggers full materialization
        # of the expensive candidates (ExonSkipping, IntronRetention,
        # cryptic).
        if self._normal_splicing_eager is None:
            return None
        return self._normal_splicing_eager.effect.coding_effect

    # Back-compat alias. Existing callers reading ``effect.alternate_effect``
    # on the default 2-outcome shape (``ExonicSpliceSite``) keep working
    # when the effect is now wrapped in a ``SpliceOutcomeSet``. New code
    # should use ``effect_if_splicing_unchanged``.
    alternate_effect = effect_if_splicing_unchanged

    @property
    def candidate_proteins(self):
        """Mapping from each mechanism class to its computed mutant
        protein sequence (empty string when the protein is not
        computable from cDNA alone — i.e. intron retention and
        cryptic splice without a ``genomic_sequence`` provider).

        Keyed by class object (e.g. ``ExonSkipping``) rather than by
        enum string — match the candidate-iteration shape.
        """
        return {
            type(c.effect): _protein_str(c.effect)
            for c in self.candidates
        }

    @property
    def mutant_protein_sequences(self):
        """Set of distinct non-empty mutant protein sequences across
        all candidates.
        """
        return {
            protein for protein in (
                _protein_str(c.effect) for c in self.candidates)
            if protein
        }

    @property
    def short_description(self) -> str:
        return "splice-set:%s" % " | ".join(
            candidate.effect.short_description for candidate in self.candidates)

    def _format_short_description(self, deletion_notation):
        return "splice-set:%s" % " | ".join(
            candidate.effect.format_short_description(deletion_notation=deletion_notation)
            for candidate in self.candidates)

    def __str__(self) -> str:
        return "SpliceOutcomeSet(variant=%s, transcript=%s, candidates=[%s])" % (
            self.variant,
            getattr(self.transcript, "name", None),
            ", ".join(
                c.effect.short_description
                for c in self.candidates),
        )

    def __repr__(self) -> str:
        return str(self)




def _protein_str(effect):
    """Return ``effect.mutant_protein_sequence`` as a str, or ``""``
    if the effect carries no protein (placeholder stubs, start-loss,
    etc.). Centralizes the None/empty/str coercion that
    :attr:`SpliceOutcomeSet.candidate_proteins` used to inline.
    """
    protein = getattr(effect, "mutant_protein_sequence", None)
    return str(protein) if protein else ""




# ---------------------------------------------------------------------
# Public entry point: enumerate outcomes for a splice effect
# ---------------------------------------------------------------------


def enumerate_splice_outcomes(splice_effect, genomic_sequence=None):
    """Wrap a splice-signal disruption in a SpliceOutcomeSet.

    The recognized inputs are exactly the
    :class:`~varcode.effects.effect_classes.SpliceSite` effects
    (``SpliceDonor`` / ``SpliceAcceptor`` / ``ExonicSpliceSite`` /
    ``IntronicSpliceSite``) — the "where the signal was hit" effects.
    Anything else passes through unchanged.

    Parameters
    ----------
    splice_effect : MutationEffect
        Output of varcode's existing splice classification. Only
        :class:`SpliceSite` instances are wrapped.

    genomic_sequence : Callable[[str, int, int], str], optional
        Caller-supplied provider for genomic sequence, returning
        forward-strand bases for ``(contig, start_1based_inclusive,
        end_1based_inclusive)`` (#296). When present:

        * :class:`IntronRetention` materializes a real
          :class:`MutantTranscript` with the intron inserted and the
          protein truncated at the first in-intron stop.
        * :class:`CrypticDonor` / :class:`CrypticAcceptor` scan a
          ±50 bp window around the canonical splice boundary (both
          exon and intron sides) for the highest-scoring cryptic
          motif using the same consensus tables as
          :mod:`varcode.cryptic_exons`. The chosen cryptic site
          replaces the canonical boundary and the resulting
          truncation / extension is reflected in the
          :class:`MutantTranscript`.

        Providers that raise (missing FASTA, unknown contig, etc.)
        fall back to the predicted state transparently.

    Returns
    -------
    SpliceOutcomeSet or the original effect
        SpliceOutcomeSet wrapping the candidate mechanisms when the
        input is a :class:`SpliceSite` disruption; otherwise the input
        unchanged.
    """
    # Already-wrapped: idempotent for the no-provider case; with a
    # provider, rebuild so the new sequence flows through to lazy
    # IntronRetention / cryptic candidates.
    if isinstance(splice_effect, SpliceOutcomeSet):
        if genomic_sequence is None:
            return splice_effect
        if splice_effect._lazy_state is None:
            raise RuntimeError(
                "enumerate_splice_outcomes(splice_set, genomic_sequence=...) "
                "needs a lazy SpliceOutcomeSet; this one has already "
                "materialized its candidates.")
        raw, mechanism_order, _ = splice_effect._lazy_state
        return SpliceOutcomeSet(
            variant=splice_effect.variant,
            transcript=splice_effect.transcript,
            disrupted_signal_class=splice_effect.disrupted_signal_class,
            splice_effect=raw,
            mechanism_order=mechanism_order,
            genomic_sequence=genomic_sequence,
        )
    # Only splice-signal disruptions (the SpliceSite subclasses) get
    # wrapped; everything else passes through untouched.
    if not isinstance(splice_effect, SpliceSite):
        return splice_effect
    mechanism_order = mechanism_order_for(splice_effect)
    if mechanism_order is None:
        # Defensive: every SpliceSite subclass today is handled by
        # mechanism_order_for, but a future direct SpliceSite subclass
        # without a registered ordering would land here — pass it
        # through rather than crashing on a None mechanism order.
        return splice_effect

    # Lazy construction: the SpliceOutcomeSet eagerly builds only the
    # cheap NormalSplicing candidate (used by the fast
    # effect_if_splicing_unchanged / alternate_effect path) and defers
    # ExonSkipping / IntronRetention / cryptic until ``.candidates`` is
    # accessed.
    return SpliceOutcomeSet(
        variant=splice_effect.variant,
        transcript=getattr(splice_effect, "transcript", None),
        disrupted_signal_class=type(splice_effect),
        splice_effect=splice_effect,
        mechanism_order=mechanism_order,
        genomic_sequence=genomic_sequence,
    )


def _build_candidate_for(mechanism_cls, splice_effect, genomic_sequence):
    """Dispatch from a mechanism class to its builder. Used by the lazy
    materialization path in :class:`SpliceOutcomeSet`.
    """
    if mechanism_cls is NormalSplicing:
        return build_normal_splicing_candidate(splice_effect)
    if mechanism_cls is ExonSkipping:
        return build_exon_skipping_candidate(splice_effect)
    if mechanism_cls is IntronRetention:
        return build_intron_retention_candidate(
            splice_effect, genomic_sequence=genomic_sequence)
    if mechanism_cls in (CrypticDonor, CrypticAcceptor):
        return build_cryptic_splice_candidate(
            splice_effect, mechanism_cls,
            genomic_sequence=genomic_sequence)
    raise ValueError(
        "Unknown splice mechanism class: %r" % (mechanism_cls,))


def wrap_splice_effects_in_collection(effect_collection):
    """Apply :func:`enumerate_splice_outcomes` to every splice-related
    Effect in an EffectCollection. Non-splice effects pass through.
    """
    new_effects = [enumerate_splice_outcomes(e) for e in effect_collection]
    return effect_collection.clone_with_new_elements(new_effects)


# Priority integration happens via the `priority_class` attribute on
# SpliceOutcomeSet — see `effect_priority` in effects.effect_ordering.

# Serialization is inherited from Serializable (the base of
# MutationEffect). to_serializable_repr stamps each nested object's
# {__module__, __name__} and from_serializable_dict resolves them via
# _lookup_value, so the polymorphic EffectCandidate.effect (any
# SpliceMechanismEffect subclass) and SpliceOutcomeSet.disrupted_signal_class
# (a class object) both round-trip without hand-rolled registries.

