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

from .effects.classify import classify_from_protein_diff
from .effects.codon_tables import codon_table_for_transcript, translate_sequence
from .effects.effect_classes import (
    CrypticAcceptor,
    CrypticDonor,
    ExonSkipping,
    ExonicSpliceSite,
    IntronRetention,
    IntronicSpliceSite,
    MultiOutcomeEffect,
    NormalSplicing,
    SpliceAcceptor,
    SpliceDonor,
    SpliceMechanismEffect,
)
from .mutant_transcript import MutantTranscript, TranscriptEdit
from .effect_candidates import EffectCandidate


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


class SpliceOutcomeSet(MultiOutcomeEffect):
    """A splice-disrupting variant's effect, expressed as a set of
    plausible mechanisms rather than a single Effect.

    Implements the :class:`MultiOutcomeEffect` protocol: downstream
    consumers can write ``isinstance(effect, MultiOutcomeEffect)`` to
    catch this and any future multi-outcome wrapper (RNA evidence,
    germline-aware, etc.) uniformly.

    :attr:`candidates` is a ``tuple[EffectCandidate, ...]``. Each
    candidate's ``effect`` is a :class:`SpliceMechanismEffect`
    subclass — :class:`NormalSplicing`, :class:`ExonSkipping`,
    :class:`IntronRetention`, :class:`CrypticDonor`, or
    :class:`CrypticAcceptor`. Class identity = mechanism. Each
    mechanism Effect carries its own protein vocab (``aa_ref`` /
    ``aa_alt`` / ``mutant_protein_sequence`` / ``mutant_transcript``)
    populated when the protein math resolved, ``None`` when not.

    For back-compat, :attr:`short_description` delegates to the most
    likely candidate. Constructed by
    :func:`enumerate_splice_outcomes` when a caller passes
    ``splice_outcomes=True`` to ``Variant.effects()`` or
    ``VariantCollection.effects()``.

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

    def __init__(
            self, variant, transcript, candidates, disrupted_signal_class=None,
            dna_candidates=None, rna_evidence=(), added_candidates=(),
            excluded_candidates=(), candidate_rna_evidence=()):
        MultiOutcomeEffect.__init__(self, variant)
        self.transcript = transcript
        self._candidates = tuple(candidates)
        # Class of the original splice-disrupting effect this set
        # replaced (SpliceDonor, SpliceAcceptor, ExonicSpliceSite, or
        # IntronicSpliceSite). Used for priority lookup.
        self.disrupted_signal_class = disrupted_signal_class
        self.dna_candidates = tuple(
            self._candidates if dna_candidates is None else dna_candidates)
        self.rna_evidence = tuple(rna_evidence)
        self.added_candidates = tuple(added_candidates)
        self.excluded_candidates = tuple(excluded_candidates)
        self.candidate_rna_evidence = tuple(candidate_rna_evidence)

    @property
    def candidates(self):
        return self._combine_with_extra_candidates(self._candidates)

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
        """Delegate priority lookup to the disrupted-signal class so a
        SpliceOutcomeSet sorts as if it were the original splice effect.

        Read by :func:`varcode.effects.effect_priority`.
        """
        return self.disrupted_signal_class

    @property
    def alternate_effect(self):
        """Back-compat shim matching :attr:`ExonicSpliceSite.alternate_effect`.

        Resolves to the underlying coding change carried by the
        :class:`NormalSplicing` candidate (its
        :attr:`NormalSplicing.coding_effect`) — i.e. "the coding
        change that applies when splicing proceeds normally."
        Returns ``None`` when there is no :class:`NormalSplicing`
        candidate or when it carries ``coding_effect=None`` (e.g.
        the variant is purely intronic).
        """
        for candidate in self._candidates:
            if isinstance(candidate.effect, NormalSplicing):
                return candidate.effect.coding_effect
        return None

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
        return "splice-set:%s" % self.most_likely_candidate.effect.short_description

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


# ---------------------------------------------------------------------
# DNA-only mechanism order.
#
# Hand-tuned heuristic ordering. These are not probabilities. The order
# only gives varcode a deterministic leading candidate when DNA alone
# leaves several splice mechanisms possible.
# ---------------------------------------------------------------------


# Mechanism orders are tuples of :class:`SpliceMechanismEffect`
# subclasses — class identity IS the mechanism, no parallel enum.

# Disrupted canonical donor (intronic +1/+2): exon skipping dominates,
# intron retention common, leaky normal splicing rare.
_MECHANISM_ORDER_SPLICE_DONOR = (
    ExonSkipping,
    IntronRetention,
    CrypticDonor,
    NormalSplicing,
)

# Disrupted canonical acceptor (intronic -1/-2): same rough
# distribution as donor.
_MECHANISM_ORDER_SPLICE_ACCEPTOR = (
    ExonSkipping,
    IntronRetention,
    CrypticAcceptor,
    NormalSplicing,
)

# Disrupted exonic splice site (last 3 of exon): often splicing still
# proceeds (the disruption is in the exon, where the spliceosome can
# tolerate more variation), so normal splicing is more competitive.
_MECHANISM_ORDER_EXONIC_SPLICE_SITE = (
    NormalSplicing,
    ExonSkipping,
    CrypticDonor,
    IntronRetention,
)

# Intronic splice site (positions +3 to +6 or -3): less critical
# region; normal splicing dominates. Two variants — donor-side (+3 to +6)
# and acceptor-side (-3) — which differ only in the cryptic direction.
_MECHANISM_ORDER_INTRONIC_SPLICE_SITE_DONOR = (
    NormalSplicing,
    ExonSkipping,
    IntronRetention,
    CrypticDonor,
)

_MECHANISM_ORDER_INTRONIC_SPLICE_SITE_ACCEPTOR = (
    NormalSplicing,
    ExonSkipping,
    IntronRetention,
    CrypticAcceptor,
)

# Order matters for isinstance matching: more specific classes first.
# SpliceDonor and SpliceAcceptor are subclasses of IntronicSpliceSite
# (see effect_classes.py), so they must be checked before the
# IntronicSpliceSite fallback in _mechanism_order_for — otherwise
# a SpliceDonor would incorrectly fall into the intronic-window table.
# IntronicSpliceSite itself is deliberately held out of this tuple and
# dispatched separately so side detection can pick the right cryptic
# direction.
_MECHANISM_ORDER_TABLES = (
    (SpliceDonor, _MECHANISM_ORDER_SPLICE_DONOR),
    (SpliceAcceptor, _MECHANISM_ORDER_SPLICE_ACCEPTOR),
    (ExonicSpliceSite, _MECHANISM_ORDER_EXONIC_SPLICE_SITE),
)


def _intronic_splice_side_is_acceptor(splice_effect):
    """True if an IntronicSpliceSite effect is on the acceptor side of
    its nearest exon.

    The classifier emits IntronicSpliceSite for positions +3–6 (donor
    side, after exon in transcript order) and -3 (acceptor side,
    before exon in transcript order). The side determines whether a
    cryptic donor or cryptic acceptor is the relevant alternative.

    Caller contract: ``splice_effect`` must be an ``IntronicSpliceSite``
    (guaranteed by :func:`_mechanism_order_for`), so ``variant`` and
    ``nearest_exon`` are always present.
    """
    exon = splice_effect.nearest_exon
    variant = splice_effect.variant
    if exon.strand == "+":
        return variant.trimmed_base1_start < exon.start
    # Reverse strand: acceptor-side intronic variants sit past the
    # genomic end of the exon (which is the 5' end in transcript order).
    return variant.trimmed_base1_end > exon.end


def _mechanism_order_for(splice_effect):
    """Return the DNA-only mechanism order for this splice
    effect, or None for effects we don't wrap.

    Uses :func:`isinstance` rather than exact-class dispatch so
    subclasses are handled correctly.
    """
    for cls, order in _MECHANISM_ORDER_TABLES:
        if isinstance(splice_effect, cls):
            return order
    if isinstance(splice_effect, IntronicSpliceSite):
        if _intronic_splice_side_is_acceptor(splice_effect):
            return _MECHANISM_ORDER_INTRONIC_SPLICE_SITE_ACCEPTOR
        return _MECHANISM_ORDER_INTRONIC_SPLICE_SITE_DONOR
    return None


# ---------------------------------------------------------------------
# Splice-mechanism Effect construction
# ---------------------------------------------------------------------


def _wrap_in_candidate(mechanism_effect):
    """Lift a :class:`SpliceMechanismEffect` into an
    unscored :class:`EffectCandidate` with ``source="varcode"``.

    Trivial wrapper kept as a named helper for symmetry with how
    other annotators construct ``EffectCandidate`` objects.
    """
    return EffectCandidate(
        effect=mechanism_effect,
        source="varcode",
    )


def _build_normal_splicing_candidate(splice_effect):
    """Normal splicing: the splice signal was hit but the spliceosome
    handles it; the protein consequence (if any) is whatever the
    underlying coding change normally produces.

    For :class:`ExonicSpliceSite`, ``alternate_effect`` is the
    underlying coding change (a :class:`Substitution` etc.). For
    intronic splice classes there is no underlying coding change —
    the variant is in/near an intron — so ``coding_effect=None``.
    """
    transcript = getattr(splice_effect, "transcript", None)
    coding_effect = getattr(splice_effect, "alternate_effect", None)
    return _wrap_in_candidate(
        NormalSplicing(
            variant=splice_effect.variant,
            transcript=transcript,
            splice_signal=splice_effect,
            coding_effect=coding_effect),
    )


def _aa_fields_from_protein_diff(
        variant, transcript, mutant_transcript, length_delta):
    """Run :func:`classify_from_protein_diff` against the resolved
    mutant transcript and return its ``aa_ref`` / ``aa_alt`` /
    offset fields as a dict, plus the classified protein effect. The
    splice-mechanism Effect builders copy the AA vocab onto the
    mechanism instance and retain the classified effect under
    ``protein_effect`` so severity queries can delegate to it while
    class identity remains the splice mechanism.
    """
    ref_protein = str(transcript.protein_sequence)
    mut_protein = mutant_transcript.mutant_protein_sequence
    classified = classify_from_protein_diff(
        variant=variant,
        transcript=transcript,
        ref_protein=ref_protein,
        mut_protein=mut_protein,
        length_delta=length_delta,
        mutant_transcript=mutant_transcript,
    )
    return {
        "aa_ref": getattr(classified, "aa_ref", None),
        "aa_alt": getattr(classified, "aa_alt", None),
        "aa_mutation_start_offset": getattr(
            classified, "aa_mutation_start_offset", None),
        "aa_mutation_end_offset": getattr(
            classified, "aa_mutation_end_offset", None),
        "mutant_protein_sequence": mut_protein,
        "protein_effect": classified,
    }


def _build_exon_skipping_candidate(splice_effect):
    """Exon skipping: the affected exon is excluded from the mature
    transcript.

    Builds a :class:`MutantTranscript` by deleting the exon's cDNA
    range, translates, runs :func:`classify_from_protein_diff` to
    extract the protein-level AA vocab, and copies those fields onto
    a fresh :class:`ExonSkipping` instance. Returns the
    :class:`ExonSkipping` in its "predicted" state (protein fields
    ``None``) when the protein math couldn't resolve.
    """
    transcript = getattr(splice_effect, "transcript", None)
    exon = _affected_exon(splice_effect)
    if transcript is None or exon is None:
        # We can't even compute frame without the exon.
        return _wrap_in_candidate(
            ExonSkipping(
                variant=splice_effect.variant,
                transcript=transcript,
                splice_signal=splice_effect,
                affected_exon=exon,
                in_frame=False),
        )

    exon_length = exon.end - exon.start + 1
    in_frame = (exon_length % 3 == 0)
    mt = _build_exon_skip_mutant_transcript(
        splice_effect.variant, transcript, exon)

    if mt is None or mt.mutant_protein_sequence is None:
        # Couldn't compute protein (incomplete transcript, exon not
        # found, etc.) — return predicted-state ExonSkipping.
        return _wrap_in_candidate(
            ExonSkipping(
                variant=splice_effect.variant,
                transcript=transcript,
                splice_signal=splice_effect,
                affected_exon=exon,
                in_frame=in_frame,
                mutant_transcript=mt),
        )

    aa_fields = _aa_fields_from_protein_diff(
        splice_effect.variant, transcript, mt, length_delta=-exon_length)
    return _wrap_in_candidate(
        ExonSkipping(
            variant=splice_effect.variant,
            transcript=transcript,
            splice_signal=splice_effect,
            affected_exon=exon,
            in_frame=in_frame,
            mutant_transcript=mt,
            **aa_fields),
    )


def _adjacent_intron_coords(transcript, exon, side):
    """Return ``(intron_start, intron_end)`` (1-based inclusive) of the
    intron on ``side`` of ``exon`` within ``transcript``, or ``None``
    if the exon has no intron on that side (e.g. first or last exon).

    ``side`` is ``"donor"`` (transcript-3' of exon → intron after it
    in transcript order) or ``"acceptor"`` (transcript-5' of exon).
    Handles forward + reverse strand transcripts.
    """
    exons = list(transcript.exons)
    try:
        idx = next(
            i for i, ex in enumerate(exons) if ex.exon_id == exon.exon_id)
    except StopIteration:
        return None
    reverse = transcript.on_backward_strand
    if side == "donor":
        adjacent_idx = idx + 1
    elif side == "acceptor":
        adjacent_idx = idx - 1
    else:
        return None
    if adjacent_idx < 0 or adjacent_idx >= len(exons):
        return None
    adjacent = exons[adjacent_idx]
    # Intron genomic coords are between the two exons on the reference.
    if reverse:
        # Transcript order: higher genomic → lower. Intron on donor
        # side sits below exon; acceptor side above.
        if side == "donor":
            intron_start = adjacent.end + 1
            intron_end = exon.start - 1
        else:
            intron_start = exon.end + 1
            intron_end = adjacent.start - 1
    else:
        if side == "donor":
            intron_start = exon.end + 1
            intron_end = adjacent.start - 1
        else:
            intron_start = adjacent.end + 1
            intron_end = exon.start - 1
    if intron_end < intron_start:
        return None
    return (intron_start, intron_end)


def _build_intron_retention_mutant_transcript(
        variant, transcript, exon, side, genomic_sequence):
    """Construct a :class:`MutantTranscript` for an intron-retention
    outcome using caller-supplied genomic sequence (#296).

    Returns ``None`` if the intron can't be resolved (no adjacent
    exon, provider raised, or transcript is incomplete).
    """
    if not transcript.complete:
        return None
    coords = _adjacent_intron_coords(transcript, exon, side)
    if coords is None:
        return None
    intron_start, intron_end = coords
    try:
        intron_seq = genomic_sequence(
            transcript.contig, intron_start, intron_end)
    except Exception:
        return None
    if not intron_seq:
        return None
    intron_seq = intron_seq.upper()
    if transcript.on_backward_strand:
        # Caller returns forward-strand genomic bases; for a reverse
        # strand transcript we need the reverse complement to get
        # transcript-order cDNA.
        from .nucleotides import reverse_complement
        intron_seq = reverse_complement(intron_seq)
    # Insertion point in the cDNA: immediately after the exon in
    # transcript order when the donor side is retained; immediately
    # before it when the acceptor side is retained.
    try:
        exon_start_offset = _exon_start_offset_in_transcript(transcript, exon)
    except ValueError:
        return None
    exon_length = exon.end - exon.start + 1
    if side == "donor":
        insert_at = exon_start_offset + exon_length
    else:
        insert_at = exon_start_offset
    full_cdna = str(transcript.sequence)
    mutant_cdna = (
        full_cdna[:insert_at] + intron_seq + full_cdna[insert_at:])
    edit = TranscriptEdit(
        cdna_start=insert_at,
        cdna_end=insert_at,
        alt_bases=intron_seq,
        source_variant=variant,
    )
    # Translate from the canonical start codon; premature stop inside
    # the retained intron terminates the ORF.
    cds_start = min(transcript.start_codon_spliced_offsets)
    mut_protein = None
    if cds_start < len(mutant_cdna):
        codon_table = codon_table_for_transcript(transcript)
        coding = mutant_cdna[cds_start:]
        truncated = coding[:(len(coding) // 3) * 3]
        try:
            mut_protein = translate_sequence(
                truncated, codon_table=codon_table, to_stop=True)
        except ValueError:
            mut_protein = None
    return MutantTranscript(
        reference_transcript=transcript,
        edits=(edit,),
        cdna_sequence=mutant_cdna,
        mutant_protein_sequence=mut_protein,
        annotator_name="splice_outcomes",
    )


def _splice_side_for_effect(splice_effect):
    """Return ``"donor"`` / ``"acceptor"`` for the disrupted side, or
    ``None`` if the splice class doesn't distinguish (falls through
    to a single-side default in the caller)."""
    if isinstance(splice_effect, SpliceDonor):
        return "donor"
    if isinstance(splice_effect, SpliceAcceptor):
        return "acceptor"
    if isinstance(splice_effect, IntronicSpliceSite):
        # The existing _intronic_splice_side_is_acceptor helper keys
        # the side decision off the variant's position relative to
        # the exon.
        return "acceptor" if (
            _intronic_splice_side_is_acceptor(splice_effect)) else "donor"
    if isinstance(splice_effect, ExonicSpliceSite):
        # ExonicSpliceSite sits at the end of an exon — default donor.
        return "donor"
    return None


def _build_intron_retention_candidate(splice_effect, genomic_sequence=None):
    """Intron retention: the spliceosome fails and the intron is
    transcribed through to the mature mRNA.

    With a ``genomic_sequence`` provider, materializes the retained-
    intron :class:`MutantTranscript` and fills in the protein vocab
    (#296). Without a provider, the :class:`IntronRetention` is
    returned in its "predicted" state — ``retained_intron_start`` /
    ``retained_intron_end`` / ``side`` are still populated (from the
    adjacent-exon geometry) so consumers know *which* intron, but
    the protein fields are ``None``.
    """
    transcript = getattr(splice_effect, "transcript", None)
    exon = _affected_exon(splice_effect)
    side = _splice_side_for_effect(splice_effect)

    # Compute the retained-intron coords from the adjacent-exon
    # geometry, independent of provider availability. When the
    # geometry isn't resolvable (e.g. first/last exon, transcript
    # missing the adjacent exon), pass None — the IntronRetention
    # instance carries None rather than a sentinel zero so consumers
    # can tell "unknown" from "actually coord 0".
    coords = None
    if transcript is not None and exon is not None and side is not None:
        coords = _adjacent_intron_coords(transcript, exon, side)
    intron_start = coords[0] if coords is not None else None
    intron_end = coords[1] if coords is not None else None

    common = dict(
        variant=splice_effect.variant,
        transcript=transcript,
        splice_signal=splice_effect,
        retained_intron_start=intron_start,
        retained_intron_end=intron_end,
        side=side,
    )

    if (genomic_sequence is not None and transcript is not None
            and exon is not None and side is not None):
        mt = _build_intron_retention_mutant_transcript(
            splice_effect.variant, transcript, exon, side, genomic_sequence)
        if mt is not None and mt.mutant_protein_sequence is not None:
            length_delta = len(mt.cdna_sequence) - len(
                str(transcript.sequence))
            aa_fields = _aa_fields_from_protein_diff(
                splice_effect.variant, transcript, mt, length_delta)
            return _wrap_in_candidate(
                IntronRetention(mutant_transcript=mt, **common, **aa_fields),
            )

    # Predicted state: no protein math.
    return _wrap_in_candidate(
        IntronRetention(**common),
    )


_CRYPTIC_SCAN_FLANK = 50
"""Default number of bases to scan on each side of the canonical
splice boundary when looking for cryptic donor/acceptor sites.
Cryptic activations within ~50 bp are the usual clinical observations;
consumers wanting a wider scan pass ``cryptic_scan_flank`` through
:func:`enumerate_splice_outcomes`."""


def _best_cryptic_site(sequence, canonical_offset, kind):
    """Scan ``sequence`` for the highest-scoring cryptic splice site
    of ``kind`` (``"donor"`` or ``"acceptor"``), ignoring the
    canonical position itself (#296).

    ``sequence`` is in transcript order (5'→3'); the canonical
    boundary sits AFTER position ``canonical_offset`` for donor (the
    exon/intron junction) and BEFORE for acceptor. Returns
    ``(offset, score)`` — ``offset`` is the cryptic equivalent of the
    canonical boundary in ``sequence`` coordinates — or ``None`` when
    no above-threshold site is found.
    """
    from .cryptic_exons import (
        ACCEPTOR_WINDOW,
        DONOR_WINDOW,
        score_acceptor,
        score_donor,
    )
    scorer = score_donor if kind == "donor" else score_acceptor
    window = DONOR_WINDOW if kind == "donor" else ACCEPTOR_WINDOW
    # For donor: offset of returned candidate = position of
    # exon/intron junction in sequence (= end of exonic portion).
    # Donor window is 3 exonic + 6 intronic → junction sits at
    # window_start + 3.
    #
    # For acceptor: window is 3 intronic + 1 exonic → junction
    # (= start of new exon) sits at window_start + 3.
    junction_offset_in_window = 3
    best = None
    for i in range(len(sequence) - window + 1):
        sub = sequence[i:i + window]
        score = scorer(sub)
        if score <= 0:
            continue
        junction = i + junction_offset_in_window
        # Skip the canonical site itself — we're looking for alternatives.
        if junction == canonical_offset:
            continue
        if best is None or score > best[1]:
            best = (junction, score)
    return best


def _build_cryptic_site_mutant_transcript(
        variant, transcript, exon, side, genomic_sequence,
        scan_flank=_CRYPTIC_SCAN_FLANK):
    """Construct a :class:`MutantTranscript` where the canonical
    splice boundary on ``side`` (``"donor"`` or ``"acceptor"``) is
    replaced by the best-scoring cryptic site within ``scan_flank``
    bp on either side (#296).

    Returns ``(mt, score, cryptic_genomic_pos)`` or ``None`` when no
    cryptic site is found or the transcript is incomplete.
    """
    if not transcript.complete:
        return None
    reverse = transcript.on_backward_strand
    # Determine the scan region around the canonical boundary.
    # Donor boundary = exon's transcript-3' end in genomic coords:
    #   forward: exon.end; reverse: exon.start.
    # Acceptor boundary = exon's transcript-5' start:
    #   forward: exon.start; reverse: exon.end.
    if side == "donor":
        boundary_pos = exon.start if reverse else exon.end
        # On forward strand, scan extends into the intron (higher coords)
        # by scan_flank and into the exon (lower coords) by scan_flank.
        # On reverse strand, intron is at LOWER coords.
        if reverse:
            scan_start = boundary_pos - scan_flank
            scan_end = boundary_pos + scan_flank
        else:
            scan_start = boundary_pos - scan_flank
            scan_end = boundary_pos + scan_flank
    else:  # acceptor
        boundary_pos = exon.end if reverse else exon.start
        scan_start = boundary_pos - scan_flank
        scan_end = boundary_pos + scan_flank
    # Fetch forward-strand sequence.
    try:
        forward_seq = genomic_sequence(
            transcript.contig, scan_start, scan_end).upper()
    except Exception:
        return None
    if not forward_seq:
        return None
    # For reverse-strand transcript, transcript-order sequence is the
    # reverse complement of the forward-strand region.
    from .nucleotides import reverse_complement
    if reverse:
        scan_seq = reverse_complement(forward_seq)
    else:
        scan_seq = forward_seq
    # Map genomic positions in the scan window → sequence offsets in
    # transcript order.
    region_length = scan_end - scan_start + 1
    if reverse:
        def genomic_to_seq_offset(g_pos):
            return scan_end - g_pos
    else:
        def genomic_to_seq_offset(g_pos):
            return g_pos - scan_start
    canonical_offset = genomic_to_seq_offset(boundary_pos)
    # Correction for donor/acceptor "junction" convention in
    # _best_cryptic_site: donor's junction = last exonic base → on
    # forward strand that's boundary_pos, offset (boundary_pos - scan_start).
    # For acceptor: junction = first exonic base = boundary_pos, same
    # offset calculation. Both map directly.
    result = _best_cryptic_site(scan_seq, canonical_offset, kind=side)
    if result is None:
        return None
    cryptic_offset_in_seq, score = result
    # Translate back to genomic position.
    if reverse:
        cryptic_genomic_pos = scan_end - cryptic_offset_in_seq
    else:
        cryptic_genomic_pos = scan_start + cryptic_offset_in_seq
    # Build the mutant cDNA by replacing the canonical exon boundary
    # with the cryptic one.
    #
    # Delta in exon length (transcript-order bases):
    #   donor:   cryptic - canonical (positive = extends exon)
    #   acceptor: canonical - cryptic (positive = extends exon at 5' end)
    if side == "donor":
        exon_length_delta = cryptic_offset_in_seq - canonical_offset
    else:
        exon_length_delta = canonical_offset - cryptic_offset_in_seq
    try:
        exon_start_in_tx = _exon_start_offset_in_transcript(transcript, exon)
    except ValueError:
        return None
    exon_len = exon.end - exon.start + 1
    full_cdna = str(transcript.sequence)
    if side == "donor":
        # Donor change shifts the exon's 3' end. Junction in transcript
        # cDNA is at exon_start_in_tx + exon_len.
        junction_cdna = exon_start_in_tx + exon_len
    else:
        # Acceptor change shifts the exon's 5' start.
        junction_cdna = exon_start_in_tx
    if exon_length_delta >= 0:
        # Exon is extended → insert intron bases from the scan window.
        if side == "donor":
            # Added bases sit immediately after the canonical junction.
            # In transcript order they come from scan_seq between
            # canonical_offset and cryptic_offset_in_seq.
            added = scan_seq[canonical_offset:cryptic_offset_in_seq]
            insert_at = junction_cdna
        else:
            # Acceptor extension: intron bases go before the canonical
            # junction.
            added = scan_seq[cryptic_offset_in_seq:canonical_offset]
            insert_at = junction_cdna
        mutant_cdna = (
            full_cdna[:insert_at] + added + full_cdna[insert_at:])
        edit = TranscriptEdit(
            cdna_start=insert_at,
            cdna_end=insert_at,
            alt_bases=added,
            source_variant=variant,
        )
    else:
        # Exon is truncated → remove bases from the exon's end (donor)
        # or start (acceptor).
        removed = -exon_length_delta
        if side == "donor":
            remove_start = junction_cdna - removed
            remove_end = junction_cdna
        else:
            remove_start = junction_cdna
            remove_end = junction_cdna + removed
        mutant_cdna = (
            full_cdna[:remove_start] + full_cdna[remove_end:])
        edit = TranscriptEdit(
            cdna_start=remove_start,
            cdna_end=remove_end,
            alt_bases="",
            source_variant=variant,
        )
    cds_start = min(transcript.start_codon_spliced_offsets)
    mut_protein = None
    if cds_start < len(mutant_cdna):
        codon_table = codon_table_for_transcript(transcript)
        coding = mutant_cdna[cds_start:]
        truncated = coding[:(len(coding) // 3) * 3]
        try:
            mut_protein = translate_sequence(
                truncated, codon_table=codon_table, to_stop=True)
        except ValueError:
            mut_protein = None
    mt = MutantTranscript(
        reference_transcript=transcript,
        edits=(edit,),
        cdna_sequence=mutant_cdna,
        mutant_protein_sequence=mut_protein,
        annotator_name="splice_outcomes",
    )
    return (mt, score, cryptic_genomic_pos)


def _build_cryptic_splice_candidate(
        splice_effect, mechanism_cls, genomic_sequence=None,
        scan_flank=_CRYPTIC_SCAN_FLANK):
    """Cryptic donor/acceptor candidate. When ``genomic_sequence`` is
    provided (#296), scan the canonical splice boundary's flanking
    region on BOTH exon and intron sides for the best-scoring weak
    motif and materialize a real :class:`MutantTranscript` that uses
    the cryptic site in place of the disrupted canonical one.

    Without a provider (or when no motif scores above threshold),
    falls back to the predicted state — :class:`CrypticDonor` /
    :class:`CrypticAcceptor` with protein fields ``None`` and the
    cryptic-site fields (``cryptic_genomic_position``,
    ``motif_score``, ``exon_length_delta``) also ``None``.

    ``mechanism_cls`` is :class:`CrypticDonor` or
    :class:`CrypticAcceptor`; ``direction`` is derived from the
    class.
    """
    direction = mechanism_cls.direction
    transcript = getattr(splice_effect, "transcript", None)
    exon = _affected_exon(splice_effect)

    common = dict(
        variant=splice_effect.variant,
        transcript=transcript,
        splice_signal=splice_effect,
        affected_exon=exon,
    )

    if (genomic_sequence is not None and transcript is not None
            and exon is not None):
        result = _build_cryptic_site_mutant_transcript(
            splice_effect.variant, transcript, exon, direction,
            genomic_sequence, scan_flank=scan_flank)
        if result is not None:
            mt, score, cryptic_pos = result
            length_delta = len(mt.cdna_sequence) - len(
                str(transcript.sequence))
            if mt.mutant_protein_sequence is not None:
                aa_fields = _aa_fields_from_protein_diff(
                    splice_effect.variant, transcript, mt, length_delta)
                return _wrap_in_candidate(
                    mechanism_cls(
                        cryptic_genomic_position=cryptic_pos,
                        motif_score=score,
                        exon_length_delta=length_delta,
                        mutant_transcript=mt,
                        **common,
                        **aa_fields),
                )
            # MT built but protein didn't translate — partial info.
            return _wrap_in_candidate(
                mechanism_cls(
                    cryptic_genomic_position=cryptic_pos,
                    motif_score=score,
                    exon_length_delta=length_delta,
                    mutant_transcript=mt,
                    **common),
            )

    # Predicted state.
    return _wrap_in_candidate(
        mechanism_cls(**common),
    )


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------


def _build_exon_skip_mutant_transcript(variant, transcript, exon):
    """Build a :class:`MutantTranscript` representing an exon skip.

    Deletes the exon's cDNA range from the transcript sequence and
    translates the result. Adjusts the CDS start offset when the
    skipped exon precedes the coding region (5'UTR exon). When the
    skipped exon contains the start codon, returns a MutantTranscript
    with an empty protein so the classifier emits StartLoss.
    """
    if not transcript.complete:
        return None
    try:
        exon_start_in_tx = _exon_start_offset_in_transcript(
            transcript, exon)
    except (StopIteration, ValueError):
        return None

    exon_length = exon.end - exon.start + 1
    full_sequence = str(transcript.sequence)

    post_skip_cdna = (
        full_sequence[:exon_start_in_tx]
        + full_sequence[exon_start_in_tx + exon_length:]
    )

    cds_start = min(transcript.start_codon_spliced_offsets)

    # Adjust CDS start for the deletion.
    if exon_start_in_tx + exon_length <= cds_start:
        # Exon is entirely in the 5' UTR → shift CDS start left.
        new_cds_start = cds_start - exon_length
    elif exon_start_in_tx >= cds_start:
        # Exon is entirely after the CDS start → no shift.
        new_cds_start = cds_start
    else:
        # Exon overlaps the CDS start → start codon lost.
        new_cds_start = None

    mut_protein = ""
    if new_cds_start is not None and new_cds_start < len(post_skip_cdna):
        codon_table = codon_table_for_transcript(transcript)
        coding = post_skip_cdna[new_cds_start:]
        truncated = coding[:(len(coding) // 3) * 3]
        try:
            mut_protein = translate_sequence(
                truncated, codon_table=codon_table, to_stop=True)
        except ValueError:
            mut_protein = ""

    edit = TranscriptEdit(
        cdna_start=exon_start_in_tx,
        cdna_end=exon_start_in_tx + exon_length,
        alt_bases="",
        source_variant=variant,
    )

    return MutantTranscript(
        reference_transcript=transcript,
        edits=(edit,),
        cdna_sequence=post_skip_cdna,
        mutant_protein_sequence=mut_protein,
        annotator_name="splice_outcomes",
    )


def _affected_exon(splice_effect):
    """Best-effort extraction of the affected exon from a splice effect.

    ExonicSpliceSite carries it directly. SpliceDonor/SpliceAcceptor
    carry the nearest exon. IntronicSpliceSite likewise.
    """
    if hasattr(splice_effect, "exon"):
        return splice_effect.exon
    if hasattr(splice_effect, "nearest_exon"):
        return splice_effect.nearest_exon
    return None


def _protein_str(effect):
    """Return ``effect.mutant_protein_sequence`` as a str, or ``""``
    if the effect carries no protein (placeholder stubs, start-loss,
    etc.). Centralizes the None/empty/str coercion that
    :attr:`SpliceOutcomeSet.candidate_proteins` used to inline.
    """
    protein = getattr(effect, "mutant_protein_sequence", None)
    return str(protein) if protein else ""


def _exon_start_offset_in_transcript(transcript, exon):
    """Offset (in transcript coordinates) of the first base of the
    given exon.
    """
    offset = 0
    for ex in transcript.exons:
        if ex.exon_id == exon.exon_id:
            return offset
        offset += ex.end - ex.start + 1
    raise ValueError("Exon %s not found in transcript %s" % (
        exon.exon_id, transcript.transcript_id))


# ---------------------------------------------------------------------
# Public entry point: enumerate outcomes for a splice effect
# ---------------------------------------------------------------------


def enumerate_splice_outcomes(splice_effect, genomic_sequence=None):
    """Wrap a splice-disrupting Effect in a SpliceOutcomeSet.

    Recognized splice effect classes are SpliceDonor, SpliceAcceptor,
    ExonicSpliceSite, and IntronicSpliceSite. Unrecognized classes
    pass through unchanged.

    Parameters
    ----------
    splice_effect : MutationEffect
        Output of varcode's existing splice classification.

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
        input is a splice-disrupting effect; otherwise the input
        unchanged.
    """
    mechanism_order = _mechanism_order_for(splice_effect)
    if mechanism_order is None:
        return splice_effect

    candidates = []
    for mechanism_cls in mechanism_order:
        if mechanism_cls is NormalSplicing:
            candidates.append(_build_normal_splicing_candidate(
                splice_effect))
        elif mechanism_cls is ExonSkipping:
            candidates.append(_build_exon_skipping_candidate(
                splice_effect))
        elif mechanism_cls is IntronRetention:
            candidates.append(_build_intron_retention_candidate(
                splice_effect,
                genomic_sequence=genomic_sequence))
        elif mechanism_cls in (CrypticDonor, CrypticAcceptor):
            candidates.append(_build_cryptic_splice_candidate(
                splice_effect, mechanism_cls,
                genomic_sequence=genomic_sequence))

    return SpliceOutcomeSet(
        variant=splice_effect.variant,
        transcript=getattr(splice_effect, "transcript", None),
        candidates=candidates,
        disrupted_signal_class=type(splice_effect),
    )


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
