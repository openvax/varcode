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


from __future__ import annotations

from collections.abc import Iterable, Sequence
from functools import cached_property
from typing import TYPE_CHECKING, Any, Literal

from serializable import Serializable

from .common import validate_deletion_notation

if TYPE_CHECKING:
    from ..effect_candidates import EffectCandidate
    from ..effect_traversal import Candidate
    from ..mutant_transcript import MutantTranscript
    from ..orf import ORFAnnotation
    from ..protocols import GeneLike, TranscriptLike


class MutationEffect(Serializable):
    """
    Base class for mutation effects.
    """

    def __init__(self, variant):
        self.variant = variant

    def to_dict(self):
        from .structural_serialization import (
            needs_effect_graph_serialization,
            structural_effect_to_dict,
        )

        if self.orf is not None or self.orf_annotation is not None:
            return structural_effect_to_dict(self, schema="orf_effect_schema")
        if needs_effect_graph_serialization(self):
            return structural_effect_to_dict(self, schema="effect_graph_schema")
        if getattr(self.variant, "is_structural", False):
            return structural_effect_to_dict(self)
        return super().to_dict()

    @classmethod
    def from_dict(cls, state_dict):
        if "orf_effect_schema" in state_dict:
            from .structural_serialization import structural_effect_from_dict
            return structural_effect_from_dict(cls, state_dict, schema="orf_effect_schema")
        if "effect_graph_schema" in state_dict:
            from .structural_serialization import structural_effect_from_dict
            return structural_effect_from_dict(cls, state_dict, schema="effect_graph_schema")
        if "structural_effect_schema" in state_dict:
            from .structural_serialization import structural_effect_from_dict
            return structural_effect_from_dict(cls, state_dict)
        return super().from_dict(state_dict)

    @classmethod
    def from_json(cls, json_string):
        """Decode an archived effect using the finite reviewed class allowlist."""
        from .structural_serialization import effect_from_json

        result = effect_from_json(json_string)
        if not isinstance(result, cls):
            raise ValueError(
                f"Expected archive root to decode to {cls.__name__}, "
                f"got {type(result).__name__}")
        return result

    def __hash__(self):
        from .structural_serialization import needs_effect_graph_serialization

        if (getattr(self.variant, "is_structural", False)
                or isinstance(self, MultiOutcomeEffect)
                or "source_evidence" in self.__dict__
                or needs_effect_graph_serialization(self)):
            # Evidence contains dictionaries and candidates can point back
            # to this effect. Equal effects have equal identity fields.
            variants = getattr(self, "variants", None)
            return hash((type(self), tuple(variants) if variants else self.variant, self.transcript))
        return super().__hash__()

    def __str__(self):
        return "%s(variant=%s)" % (self.__class__.__name__, self.variant)

    def __repr__(self):
        return str(self)

    def __lt__(self, other):
        """
        Effects are ordered by their associated variants, which have
        comparison implement in terms of their chromosomal locations.
        """
        return self.variant < other.variant

    @property
    def short_description(self) -> str:
        """
        A short but human-readable description of the effect.
        Defaults to class name for most of the non-coding effects,
        but is more informative for coding ones.
        """
        return self.__class__.__name__.lower()

    def format_short_description(
            self, *, deletion_notation: Literal["sequence", "range"] = "sequence",
    ) -> str:
        """Render the short description with a protein deletion display choice.

        Parameters
        ----------
        deletion_notation : {"sequence", "range"}
            ``"sequence"`` (default) includes every deleted residue and its
            one-based start, e.g. ``p.RSP3del``. ``"range"`` uses the first
            and last residue positions, e.g. ``p.R3_P5del``. Other consequence
            and splice-mechanism labels retain their existing notation.

        Returns
        -------
        str
            A display label. This does not change the effect or its archive.
        """
        validate_deletion_notation(deletion_notation)
        if deletion_notation == "sequence":
            return self.short_description
        return self._format_short_description(deletion_notation)

    def _format_short_description(self, deletion_notation):
        return self.short_description

    transcript: TranscriptLike | None = None
    orf_annotation: ORFAnnotation | None = None
    discovered_orf: Any | None = None

    @property
    def gene(self) -> GeneLike | None:
        """Parent gene for this effect, or ``None`` when unavailable."""
        return self.__dict__.get("gene")

    @property
    def orf(self):
        """The selected reference ORF, or None for ordinary annotation."""
        if self.orf_annotation is not None:
            return self.orf_annotation.reference
        return getattr(self.transcript, "reference_orf", None)

    @property
    def reference_transcript(self) -> TranscriptLike | None:
        """Original annotated RNA, including when an extra coding frame is used."""
        return self.orf.transcript if self.orf is not None else self.transcript

    @property
    def original_protein_sequence(self) -> str | None:
        """Pre-somatic protein: the assembled patient baseline when supplied,
        otherwise the reference protein. None means unavailable.
        """
        if "_original_protein_sequence" in self.__dict__:
            return self.__dict__["_original_protein_sequence"]
        if self.baseline_transcript is not None:
            return self.baseline_transcript.mutant_protein_sequence
        if self.transcript is not None:
            return self.transcript.protein_sequence
        else:
            return None

    @original_protein_sequence.setter
    def original_protein_sequence(self, sequence: str | None) -> None:
        self.__dict__["_original_protein_sequence"] = sequence

    @property
    def gene_name(self) -> str | None:
        if self.gene is not None:
            return self.gene.name
        else:
            return None

    @property
    def gene_id(self) -> str | None:
        if "gene" in self.__dict__:
            return self.gene.id if self.gene is not None else None
        transcript = self.reference_transcript
        if transcript is not None:
            transcript_gene_id = getattr(transcript, "gene_id", None)
            if transcript_gene_id is not None:
                return transcript_gene_id
        if self.gene is not None:
            return self.gene.id
        return None

    @property
    def transcript_name(self) -> str | None:
        if self.reference_transcript is not None:
            return self.reference_transcript.name
        else:
            return None

    @property
    def transcript_id(self) -> str | None:
        if self.reference_transcript is not None:
            return self.reference_transcript.id
        else:
            return None

    # Default tri-state change flags; subclasses override either the class
    # default or the ``_get_*`` hook while keeping instance-level overrides
    # (``_modifies_*``) and archive restoration working cleanly.
    _default_modifies_coding_sequence: bool | None = False
    _default_modifies_protein_sequence: bool | None = False

    def _get_modifies_coding_sequence(self) -> bool | None:
        return self._default_modifies_coding_sequence

    @property
    def modifies_coding_sequence(self) -> bool | None:
        """Tri-state flag: True if coding sequence changes, False if unchanged, None if unknown."""
        if "_modifies_coding_sequence" in self.__dict__:
            return self.__dict__["_modifies_coding_sequence"]
        if "modifies_coding_sequence" in self.__dict__:
            return self.__dict__["modifies_coding_sequence"]
        return self._get_modifies_coding_sequence()

    @modifies_coding_sequence.setter
    def modifies_coding_sequence(self, value: bool | None) -> None:
        self.__dict__["_modifies_coding_sequence"] = value

    def _get_modifies_protein_sequence(self) -> bool | None:
        return self._default_modifies_protein_sequence

    @property
    def modifies_protein_sequence(self) -> bool | None:
        """Tri-state flag: True if protein sequence changes, False if unchanged, None if unknown."""
        if "_modifies_protein_sequence" in self.__dict__:
            return self.__dict__["_modifies_protein_sequence"]
        if "modifies_protein_sequence" in self.__dict__:
            return self.__dict__["modifies_protein_sequence"]
        return self._get_modifies_protein_sequence()

    @modifies_protein_sequence.setter
    def modifies_protein_sequence(self, value: bool | None) -> None:
        self.__dict__["_modifies_protein_sequence"] = value

    def _get_mutant_protein_sequence(self) -> str | None:
        if self.mutant_transcript is not None:
            return self.mutant_transcript.mutant_protein_sequence
        return None

    @property
    def mutant_protein_sequence(self) -> str | None:
        """Protein from the assembled transcript, when one is available."""
        if "_mutant_protein_sequence" in self.__dict__:
            return self._mutant_protein_sequence
        return self._get_mutant_protein_sequence()

    @mutant_protein_sequence.setter
    def mutant_protein_sequence(self, sequence: str | None) -> None:
        self._mutant_protein_sequence = sequence

    # which residues in the mutant protein sequence are different from
    # the original?
    @property
    def aa_mutation_start_offset(self) -> int | None:
        """Zero-based start offset of mutated residues in the protein, or ``None``."""
        return self.__dict__.get("aa_mutation_start_offset")

    @aa_mutation_start_offset.setter
    def aa_mutation_start_offset(self, value: int | None) -> None:
        self.__dict__["aa_mutation_start_offset"] = value

    @property
    def aa_mutation_end_offset(self) -> int | None:
        """Zero-based half-open end offset of mutated residues in the protein, or ``None``."""
        return self.__dict__.get("aa_mutation_end_offset")

    @aa_mutation_end_offset.setter
    def aa_mutation_end_offset(self, value: int | None) -> None:
        self.__dict__["aa_mutation_end_offset"] = value

    # Optional full :class:`~varcode.MutantTranscript` attached to
    # this effect (#269 / #296 / #335). SV effects and splice
    # exon-skip candidates populate this during construction. Other
    # effect kinds (Substitution, FrameShift, etc.) leave it None by
    # default; an Isovar-style :class:`PhaseResolver` can attach an
    # assembled-contig MutantTranscript post-hoc via
    # :func:`varcode.phasing.apply_phase_resolver_to_effects`.
    @property
    def mutant_transcript(self) -> MutantTranscript | None:
        """Assembled mutant transcript model when available, else ``None``."""
        return self.__dict__.get("mutant_transcript")

    @mutant_transcript.setter
    def mutant_transcript(self, value: MutantTranscript | None) -> None:
        self.__dict__["mutant_transcript"] = value

    baseline_transcript: MutantTranscript | None = None
    """The actual pre-somatic transcript used for a composed comparison.

    None means no baseline model was retained. A model whose sequence is None
    instead records an explicitly unavailable product. The model uses the
    same MutantTranscript representation as the post-somatic product.
    """

    def _get_candidates(self) -> tuple[Candidate, ...]:
        return ()

    @property
    def candidates(self) -> tuple[Candidate, ...]:
        """Candidate outcomes in producer order, or ``()`` for bare effects.

        Every effect exposes the same candidate surface: default ``"dna"``
        effects, ``"transcript_model"`` effects, and :class:`MultiOutcomeEffect`
        containers return their candidate tuple here; bare ``"fast"`` effects
        return ``()``. Use :meth:`add_candidates` to append RNA, splice, or
        motif observations without touching private attributes.
        """
        if "candidates" in self.__dict__:
            base: tuple[Candidate, ...] = self.__dict__["candidates"]
        elif "realized_candidates" in self.__dict__:
            base = tuple(self.__dict__["realized_candidates"])
        else:
            base = self._get_candidates()
        extra: tuple[Candidate, ...] = self.__dict__.get("extra_candidates") or ()
        return base + extra if extra else base

    @candidates.setter
    def candidates(self, value: Sequence[Candidate]) -> None:
        self.__dict__["candidates"] = tuple(value)
        self.__dict__.pop("extra_candidates", None)

    def add_candidates(self, candidates: Iterable[Candidate]) -> None:
        """Append additional candidate outcomes in producer order.

        Preserves existing candidate order so :attr:`most_likely_candidate`
        stays unchanged. On effects storing an explicit ``candidates`` tuple,
        appends directly; on dynamic or lazy outcome sets, appends without
        forcing premature materialization of base candidates.
        """
        incoming = tuple(candidates)
        if not incoming:
            return
        if "candidates" in self.__dict__:
            self.__dict__["candidates"] = self.__dict__["candidates"] + incoming
        elif "realized_candidates" in self.__dict__:
            self.__dict__["candidates"] = tuple(self.__dict__["realized_candidates"]) + incoming
        elif type(self)._get_candidates is MutationEffect._get_candidates:
            self.__dict__["candidates"] = incoming
        else:
            existing = self.__dict__.get("extra_candidates") or ()
            self.__dict__["extra_candidates"] = tuple(existing) + incoming

    @property
    def effects(self) -> tuple[MutationEffect, ...]:
        """Tuple of inner :class:`MutationEffect` objects in :attr:`candidates` order."""
        return tuple(c.effect for c in self.candidates)

    @property
    def most_likely_candidate(self) -> Candidate | None:
        """First candidate in producer order, or ``None`` when :attr:`candidates` is empty."""
        candidates = self.candidates
        return candidates[0] if candidates else None

    @property
    def most_likely_effect(self) -> MutationEffect | None:
        """Inner :class:`MutationEffect` of :attr:`most_likely_candidate`, or ``None``."""
        candidate = self.most_likely_candidate
        return candidate.effect if candidate is not None else None

    @property
    def highest_priority_candidate(self) -> Candidate | None:
        """Most protein-disruptive candidate by classification priority, or ``None``."""
        candidates = self.candidates
        if not candidates:
            return None
        from .effect_ordering import classification_priority, transcript_effect_priority_dict

        def _rank(candidate: Candidate) -> int:
            if candidate.effect is self:
                return transcript_effect_priority_dict.get(type(self), -1)
            return classification_priority(candidate.effect)

        return max(candidates, key=_rank)

    @property
    def highest_priority_effect(self) -> MutationEffect | None:
        """Inner :class:`MutationEffect` of :attr:`highest_priority_candidate`, or ``None``."""
        candidate = self.highest_priority_candidate
        return candidate.effect if candidate is not None else None


class MultiOutcomeEffect(MutationEffect):
    """Marker base class for legacy container effects representing a set of
    plausible outcomes rather than a single deterministic effect.

    All candidate accessors (:attr:`candidates`, :attr:`effects`,
    :attr:`most_likely_candidate`, :attr:`most_likely_effect`,
    :attr:`highest_priority_candidate`, :attr:`highest_priority_effect`, and
    :meth:`add_candidates`) live on :class:`MutationEffect` itself so ordinary
    classified effects and legacy outcome containers share one uniform API.

    Subclasses override ``_get_candidates()`` and expose ``priority_class``.
    Downstream consumers can inspect ``effect.candidates`` directly or use
    :func:`varcode.iter_candidate_paths` / :func:`varcode.iter_candidate_effects`.
    """

    def _get_candidates(self) -> tuple[Candidate, ...]:
        raise NotImplementedError(
            "%s must implement _get_candidates" % type(self).__name__)


class Intergenic(MutationEffect):
    """Variant has unknown effect if it occurs between genes"""

    @property
    def short_description(self) -> str:
        return "intergenic"


class Intragenic(MutationEffect):
    """
    Variant within boundaries of a gene but does not overlap
    introns or exons of any transcript. This seems very peculiar but
    apparently does happen sometimes, maybe some genes have two distinct sets
    of exons which are never simultaneously expressed?
    """

    @property
    def short_description(self) -> str:
        return "intragenic"

    def __init__(self, variant, gene):
        MutationEffect.__init__(self, variant)
        self.__dict__["gene"] = gene


class TranscriptMutationEffect(Intragenic):
    _haplotype_fields = (
        "variants", "phase_source", "annotator", "annotator_version",
        "observed_mutant_transcripts")

    def __init__(self, variant, transcript):
        MutationEffect.__init__(self, variant)
        self.transcript = transcript

    @cached_property
    def gene(self) -> GeneLike | None:
        """Resolve the full gene only when requested, not when restoring an effect.

        Transcript identity and retained consequences remain usable without an
        indexed annotation database. Accessing the gene retains the transcript's
        normal lookup and error behavior.
        """
        return self.transcript.gene if self.transcript is not None else None

    def to_dict(self):
        state = super().to_dict()
        if any(key in state for key in ("effect_graph_schema", "structural_effect_schema", "orf_effect_schema")):
            return state
        # Annotators may return any ordinary effect class for a cis group.
        # Preserve the group independently of that class's constructor fields.
        if len(getattr(self, "variants", ())) > 1:
            state["_haplotype_context"] = tuple(
                (key, getattr(self, key)) for key in self._haplotype_fields
                if hasattr(self, key))
        return state

    @classmethod
    def from_dict(cls, state_dict):
        state = dict(state_dict)
        context = dict(state.pop("_haplotype_context", ()))
        result = super().from_dict(state)
        for key in cls._haplotype_fields:
            if key in context:
                setattr(result, key, context[key])
        return result

    def __str__(self):
        orf_label = "" if self.orf is None else ", orf_id=%s" % self.orf.orf_id
        return "%s(variant=%s, transcript_name=%s, transcript_id=%s%s)" % (
            self.__class__.__name__,
            self.variant,
            self.transcript_name,
            self.transcript_id,
            orf_label)


class GermlineAlleleOverlap(TranscriptMutationEffect):
    """This reported allele is already present in the supplied germline.

    No new allele sequence is established. This does not assess allele dosage,
    loss of the other allele, expression, or the allele's pathogenicity.
    """

    @property
    def short_description(self) -> str:
        return "germline-allele-overlap"

    is_germline_overlap = True
    is_loh = None
    loh_status = "not_assessed"


class Failure(TranscriptMutationEffect):
    """
    Special placeholder effect for when we want to suppress errors but still
    need to create a non-empty list of effects for each variant.
    """
    def __init__(self, variant, transcript=None, error=None):
        TranscriptMutationEffect.__init__(self, variant, transcript)
        self.error = error

    def __str__(self):
        return "Failure(variant=%s, transcript_id=%s, error=%s)" % (
            self.variant, self.transcript_id, self.error)


class Unresolved(TranscriptMutationEffect):
    """An effect that could not be resolved with this model or evidence.

    Unlike :class:`Failure`, this is not an annotation error. A hypothesis
    may require missing genomic sequence or an RNA assembly. Alternatively,
    ``mechanism="unsupported_annotation"`` means the selected annotator
    returned ``NotImplemented`` for this input. ``reason`` explains why;
    unresolved does not mean the variant is harmless.

    ``reference_effect``, when supplied, retains a reference-only annotation
    for inspection and severity ranking. It is not the unresolved patient's
    consequence: its sequences never substitute for an unavailable patient
    product. Supplied transcript models may retain known sequences while
    the comparison and change flags remain unresolved.

    Parameters
    ----------
    variant : Variant or StructuralVariant
        The edit whose consequence is unresolved.
    transcript : transcript-like
        The reference transcript used for annotation.
    mechanism : str
        Machine-readable source of uncertainty.
    reason : str, optional
        Human-readable explanation of the missing comparison.
    reference_effect : MutationEffect, optional
        Reference-only diagnostic and ranking anchor, when available.
    """

    _default_modifies_coding_sequence = None
    _default_modifies_protein_sequence = None

    def __init__(self, variant, transcript, mechanism: str, reason: str | None = None,
                 reference_effect: MutationEffect | None = None) -> None:
        TranscriptMutationEffect.__init__(self, variant, transcript)
        self.mechanism = mechanism
        self.reason = reason
        self.reference_effect = reference_effect

    @property
    def priority_class(self) -> type[MutationEffect]:
        """Retain reference severity without asserting a patient consequence."""
        from .effect_ordering import (
            select_between_exonic_splice_site_and_alternate_effect)
        if self.reference_effect is None:
            return Unresolved
        ranked = select_between_exonic_splice_site_and_alternate_effect(
            self.reference_effect)
        return getattr(ranked, "priority_class", None) or type(ranked)

    @property
    def short_description(self):
        return "unresolved-%s" % self.mechanism.replace("_", "-")


class HypothesisLimit(Unresolved):
    """More alternative outcomes than ``max_hypotheses``, left unenumerated.

    Choosing some of the outcomes would present a guess as a prediction, so
    this effect records what is known instead. ``phase`` is the
    :class:`~varcode.PhasePartition` of nearby germline variants (known cis,
    known trans, unphased), and ``reference_effect`` is the variant's effect
    without germline context. It ranks like ``reference_effect``, so germline
    context that cannot be resolved does not change which transcript's
    effect is selected.
    """

    def __init__(self, variant, transcript, max_hypotheses, phase,
                 reference_effect=None):
        reason = "More than %d alternative outcomes" % max_hypotheses
        if phase.unphased:
            reason += " (%d unphased germline variants)" % len(phase.unphased)
        Unresolved.__init__(
            self, variant, transcript, mechanism="hypothesis_limit",
            reason=reason, reference_effect=reference_effect)
        self.max_hypotheses = max_hypotheses
        self.phase = phase


class NoncodingTranscript(TranscriptMutationEffect):
    """
    Any mutation to a transcript with a non-coding biotype
    """

    @property
    def short_description(self) -> str:
        return "non-coding-transcript"


class IncompleteTranscript(TranscriptMutationEffect):
    """
    Any mutation to an incompletely annotated transcript with a coding biotype.

    A stop codon that polyadenylation completes is not incomplete. MT-ND1,
    MT-ND2, MT-CO3, MT-ND3, MT-ND4 and MT-CYB end their reading frame in ``T``
    or ``TA`` and the poly(A) tail adds the rest of a ``TAA`` stop codon
    (NCBI NC_012920.1). Varcode translates them as complete, with the tail
    finishing the last partial codon. This needs PyEnsembl 2.34.1 or later;
    with an earlier PyEnsembl these genes are still reported as incomplete.
    """

    @property
    def short_description(self) -> str:
        return "incomplete"


class FivePrimeUTR(TranscriptMutationEffect):
    """
    Any mutation to the 5' untranslated region (before the start codon) of
    coding transcript.
    """

    @property
    def short_description(self) -> str:
        return "5' UTR"


class ThreePrimeUTR(TranscriptMutationEffect):
    """
    Any mutation to the 3' untranslated region (after the stop codon) of
    coding transcript.
    """

    @property
    def short_description(self) -> str:
        return "3' UTR"


class Intronic(TranscriptMutationEffect):
    """
    Mutation in an intronic region of a coding transcript
    """
    def __init__(self, variant, transcript, nearest_exon, distance_to_exon):
        TranscriptMutationEffect.__init__(self, variant, transcript)
        self.nearest_exon = nearest_exon
        self.distance_to_exon = distance_to_exon

    @property
    def short_description(self) -> str:
        return "intronic"


class SpliceSite(object):
    """Marker base for the four splice-signal *disruption* effects:
    :class:`SpliceDonor`, :class:`SpliceAcceptor`,
    :class:`IntronicSpliceSite`, and :class:`ExonicSpliceSite`.

    A ``SpliceSite`` effect answers **where** a variant hit a splice
    signal — a position fact (which boundary, how close). It does
    *not* describe the protein consequence: from DNA alone you can't
    say how the spliceosome will respond. That downstream protein-level
    outcome is modeled separately by :class:`SpliceMechanismEffect`
    subclasses (exon skipping, intron retention, cryptic-site use,
    normal splicing).

    So the two halves of splice handling are:

    * ``SpliceSite``           → *where* the signal was disrupted.
    * ``SpliceMechanismEffect`` → *what* the spliceosome does, and the
      protein change that results.

    :func:`varcode.splice_outcomes.enumerate_splice_outcomes` bridges
    them: given a ``SpliceSite`` disruption it produces a
    :class:`~varcode.splice_outcomes.SpliceOutcomeSet` of plausible
    mechanisms.
    """
    pass


class IntronicSpliceSite(Intronic, SpliceSite):
    """
    Mutations near exon boundaries, excluding the first two and last two
    nucleotides in an intron, since those are  known to more confidently
    affect splicing and are given their own effect classes below.
    """
    def __init__(self, variant, transcript, nearest_exon, distance_to_exon):
        Intronic.__init__(
            self, variant, transcript, nearest_exon, distance_to_exon)

    @property
    def short_description(self) -> str:
        return "intronic-splice-site"


class SpliceDonor(IntronicSpliceSite):
    """
    Mutation in the first two intron residues.
    """
    def __init__(self, variant, transcript, nearest_exon, distance_to_exon):
        IntronicSpliceSite.__init__(
            self, variant, transcript, nearest_exon, distance_to_exon)

    @property
    def short_description(self) -> str:
        return "splice-donor"


class SpliceAcceptor(IntronicSpliceSite):
    """
    Mutation in the last two intron residues.
    """

    @property
    def short_description(self) -> str:
        return "splice-acceptor"


# =====================================================================
# Splice-mechanism effects (#382).
#
# What the splice machinery does when a disrupted splice signal forces
# a choice. Each subclass represents one mechanism and carries the
# protein-level consequence of that mechanism in its own vocab —
# ``aa_ref`` / ``aa_alt`` / ``mutant_protein_sequence`` /
# ``mutant_transcript``, populated when varcode can resolve the protein
# math, ``None`` otherwise. The "unresolved" state is "protein fields
# are None"; no parallel placeholder class hierarchy.
#
# Every subclass references back to the underlying splice-signal
# disruption via :attr:`splice_signal` — a :class:`SpliceSite` effect
# (SpliceDonor / SpliceAcceptor / IntronicSpliceSite / ExonicSpliceSite)
# describing *where* in the transcript the splice signal was hit. The
# mechanism Effect describes *what* happens downstream of that
# disruption.
#
# Used by :class:`~varcode.splice_outcomes.SpliceOutcomeSet` to
# fill :attr:`EffectCandidate.effect` for each plausible mechanism a
# splice disruption could produce. Class identity = mechanism — no
# parallel enum.
# =====================================================================


class SpliceMechanismEffect(TranscriptMutationEffect):
    """Base for effects describing what the splice machinery does
    downstream of a disrupted splice signal.

    Concrete subclasses: :class:`NormalSplicing`, :class:`ExonSkipping`,
    :class:`IntronRetention`, :class:`CrypticDonor`,
    :class:`CrypticAcceptor`. Class identity is the mechanism;
    instances of this base are not constructed directly.

    Protein vocab fields are optional — ``None`` means "predicted
    mechanism, exact protein not computable from cached cDNA alone."
    When resolved, the fields follow the same conventions as
    :class:`KnownAminoAcidChange` so consumers can read them
    uniformly.

    Parameters
    ----------
    variant : Variant
    transcript : pyensembl.Transcript
    splice_signal : SpliceSite
        The splice-signal-disruption Effect *instance* describing
        *where* the disruption was — a :class:`SpliceSite` (one of
        ``SpliceDonor`` / ``SpliceAcceptor`` / ``IntronicSpliceSite`` /
        ``ExonicSpliceSite``). Lets consumers reach
        ``effect.splice_signal.nearest_exon``,
        ``effect.splice_signal.distance_to_exon``, etc. without
        carrying those fields redundantly on every mechanism subclass.
    protein_effect : MutationEffect or None
        Classified protein-level consequence of this splice mechanism
        (e.g. ``Deletion``, ``FrameShift``, ``PrematureStop``) when the
        protein sequence was resolved. The mechanism object keeps its
        own class identity while delegating severity metadata to this
        classified consequence.
    mutant_transcript : MutantTranscript or None
    aa_ref, aa_alt : str or None
    aa_mutation_start_offset, aa_mutation_end_offset : int or None
    mutant_protein_sequence : str or None
    """

    def __init__(
            self, variant, transcript, splice_signal,
            protein_effect=None,
            mutant_transcript=None,
            aa_ref=None, aa_alt=None,
            aa_mutation_start_offset=None,
            aa_mutation_end_offset=None,
            mutant_protein_sequence=None):
        TranscriptMutationEffect.__init__(self, variant, transcript)
        self.splice_signal = splice_signal
        self.protein_effect = protein_effect
        self.mutant_transcript = mutant_transcript
        self.aa_ref = aa_ref
        self.aa_alt = aa_alt
        self.aa_mutation_start_offset = aa_mutation_start_offset
        self.aa_mutation_end_offset = aa_mutation_end_offset
        self.mutant_protein_sequence = mutant_protein_sequence

    @property
    def resolved(self) -> bool:
        """True when the protein math succeeded for this mechanism —
        i.e. :attr:`mutant_protein_sequence` is populated. Convenience
        for consumers that want to branch on "do we have a concrete
        protein for this mechanism?"
        """
        return self.mutant_protein_sequence is not None

    @property
    def priority_class(self):
        if self.protein_effect is None:
            return None
        return self.protein_effect.__class__

    def _get_modifies_coding_sequence(self) -> bool | None:
        return getattr(self.protein_effect, "modifies_coding_sequence", False)

    def _get_modifies_protein_sequence(self) -> bool | None:
        return getattr(self.protein_effect, "modifies_protein_sequence", False)


class NormalSplicing(SpliceMechanismEffect):
    """Canonical splicing proceeds despite the disruption.

    The splice signal was hit, but the spliceosome handles it; the
    protein consequence (if any) is whatever the underlying
    nucleotide change normally produces — a :class:`Substitution`
    for an :class:`ExonicSpliceSite`, no protein change for a
    purely intronic disruption. The underlying coding change is
    carried as :attr:`coding_effect`; when there isn't one (intronic
    variant), it's ``None``.

    Note: this class is special among :class:`SpliceMechanismEffect`
    subclasses — it doesn't itself transform the protein, so all the
    protein-vocab attributes (``aa_ref``, ``aa_alt``,
    ``mutant_protein_sequence``, ``mutant_transcript``,
    ``aa_mutation_*_offset``) are *delegating descriptors* that read
    through to ``coding_effect``. Construction therefore calls
    :class:`TranscriptMutationEffect`'s ``__init__`` directly rather
    than :class:`SpliceMechanismEffect`'s — the latter would set
    those names as plain instance attributes and shadow the
    descriptors. The other four mechanisms (skipping, retention,
    cryptic donor/acceptor) actually transform the protein and
    carry their own protein vocab as instance state.
    """

    def __init__(self, variant, transcript, splice_signal,
                 coding_effect=None):
        TranscriptMutationEffect.__init__(self, variant, transcript)
        self.splice_signal = splice_signal
        self.coding_effect = coding_effect

    # Delegating descriptors. Read-only by construction (no setters)
    # — NormalSplicing never owns these fields directly; they're a
    # view onto ``coding_effect``. Attempting to assign to one is a
    # programming error and rightly raises AttributeError.

    @property
    def aa_ref(self):
        return getattr(self.coding_effect, "aa_ref", None)

    @property
    def aa_alt(self):
        return getattr(self.coding_effect, "aa_alt", None)

    @property
    def mutant_protein_sequence(self) -> str | None:
        return getattr(self.coding_effect, "mutant_protein_sequence", None)

    @mutant_protein_sequence.setter
    def mutant_protein_sequence(self, value: str | None) -> None:
        raise AttributeError("NormalSplicing.mutant_protein_sequence is read-only")

    @property
    def aa_mutation_start_offset(self) -> int | None:
        return getattr(self.coding_effect, "aa_mutation_start_offset", None)

    @aa_mutation_start_offset.setter
    def aa_mutation_start_offset(self, value: int | None) -> None:
        raise AttributeError("NormalSplicing.aa_mutation_start_offset is read-only")

    @property
    def aa_mutation_end_offset(self) -> int | None:
        return getattr(self.coding_effect, "aa_mutation_end_offset", None)

    @aa_mutation_end_offset.setter
    def aa_mutation_end_offset(self, value: int | None) -> None:
        raise AttributeError("NormalSplicing.aa_mutation_end_offset is read-only")

    @property
    def mutant_transcript(self) -> MutantTranscript | None:
        return getattr(self.coding_effect, "mutant_transcript", None)

    @mutant_transcript.setter
    def mutant_transcript(self, value: MutantTranscript | None) -> None:
        raise AttributeError("NormalSplicing.mutant_transcript is read-only")

    @property
    def priority_class(self):
        if self.coding_effect is None:
            return None
        return self.coding_effect.__class__

    def _get_modifies_coding_sequence(self) -> bool | None:
        return getattr(self.coding_effect, "modifies_coding_sequence", False)

    def _get_modifies_protein_sequence(self) -> bool | None:
        return getattr(self.coding_effect, "modifies_protein_sequence", False)

    @property
    def short_description(self):
        if self.coding_effect is None:
            return "normal-splicing"
        return "normal-splicing (%s)" % self.coding_effect.short_description

    def _format_short_description(self, deletion_notation):
        if self.coding_effect is None:
            return self.short_description
        return "normal-splicing (%s)" % self.coding_effect.format_short_description(
            deletion_notation=deletion_notation)


class ExonSkipping(SpliceMechanismEffect):
    """The affected exon is excluded from the mature transcript.

    When in-frame (the exon's length is divisible by 3), the result
    is a clean amino-acid deletion plus an optional boundary-codon
    reshape. When out-of-frame, a frameshift propagates from the
    new exon junction through the rest of the transcript.

    Parameters beyond the base class
    --------------------------------
    affected_exon : pyensembl.Exon
        The skipped exon.
    in_frame : bool
        ``True`` when ``len(affected_exon) % 3 == 0`` — the
        exon-skip is codon-aligned and the downstream sequence
        keeps its frame.

    Note on ``aa_ref`` / ``aa_alt`` semantics
    -----------------------------------------
    The AA-level fields come from
    :func:`~varcode.effects.classify.classify_from_protein_diff`,
    so what they mean depends on ``in_frame``:

    * **In-frame skip** (``in_frame=True``): ``aa_ref`` is the
      contiguous run of amino acids removed from the reference
      protein (including the boundary-codon reshape when the exon
      boundary lands mid-codon). ``aa_alt`` is the reshaped boundary
      codon's translation (or ``""`` for a clean codon-aligned skip).
    * **Out-of-frame skip** (``in_frame=False``): ``aa_ref`` is the
      reference AA at the divergence point (where the frame breaks
      after the new exon junction); ``aa_alt`` is the first
      frameshifted AA. The downstream protein continues in the new
      frame until a premature stop — read the full sequence from
      :attr:`mutant_protein_sequence`, not from ``aa_alt``.
    """

    def __init__(self, variant, transcript, splice_signal, affected_exon,
                 in_frame, mutant_transcript=None,
                 protein_effect=None,
                 aa_ref=None, aa_alt=None,
                 aa_mutation_start_offset=None,
                 aa_mutation_end_offset=None,
                 mutant_protein_sequence=None):
        SpliceMechanismEffect.__init__(
            self, variant, transcript, splice_signal,
            protein_effect=protein_effect,
            mutant_transcript=mutant_transcript,
            aa_ref=aa_ref, aa_alt=aa_alt,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_mutation_end_offset=aa_mutation_end_offset,
            mutant_protein_sequence=mutant_protein_sequence)
        self.affected_exon = affected_exon
        self.in_frame = in_frame

    @property
    def short_description(self):
        exon_id = getattr(self.affected_exon, "exon_id", "?")
        if not self.resolved:
            suffix = "" if self.in_frame else ", out-of-frame"
            return "exon-skip %s (predicted%s)" % (exon_id, suffix)
        if not self.aa_ref or self.aa_mutation_start_offset is None:
            # Resolved but no usable AA range — fall back to a bare label.
            return "exon-skip %s" % exon_id
        pos = self.aa_mutation_start_offset + 1
        if self.in_frame:
            # HGVS deletion: p.{first}{start}_{last}{end}del
            end = self.aa_mutation_start_offset + len(self.aa_ref)
            return "exon-skip %s (p.%s%d_%s%ddel)" % (
                exon_id, self.aa_ref[0], pos, self.aa_ref[-1], end)
        # HGVS frameshift: p.{aa}{pos}fs — distinct notation from del,
        # not concatenated with it.
        return "exon-skip %s (p.%s%dfs)" % (exon_id, self.aa_ref[0], pos)


class IntronRetention(SpliceMechanismEffect):
    """The intron stays in the mature transcript; translation
    typically hits a premature stop inside the retained intron.

    Parameters beyond the base class
    --------------------------------
    retained_intron_start, retained_intron_end : int or None
        Genomic coordinates of the retained intron (1-based
        inclusive). Populated when the adjacent-exon geometry
        identifies which intron is being retained; ``None`` when
        that can't be determined (e.g. transcript missing the
        adjacent exon, splice boundary ambiguous).
    side : str or None
        ``"donor"`` or ``"acceptor"`` — which splice boundary failed,
        causing the intron to be retained. ``None`` when the side
        can't be inferred from the splice classification.
    """

    def __init__(self, variant, transcript, splice_signal,
                 retained_intron_start=None, retained_intron_end=None,
                 side=None,
                 mutant_transcript=None,
                 protein_effect=None,
                 aa_ref=None, aa_alt=None,
                 aa_mutation_start_offset=None,
                 aa_mutation_end_offset=None,
                 mutant_protein_sequence=None):
        SpliceMechanismEffect.__init__(
            self, variant, transcript, splice_signal,
            protein_effect=protein_effect,
            mutant_transcript=mutant_transcript,
            aa_ref=aa_ref, aa_alt=aa_alt,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_mutation_end_offset=aa_mutation_end_offset,
            mutant_protein_sequence=mutant_protein_sequence)
        self.retained_intron_start = retained_intron_start
        self.retained_intron_end = retained_intron_end
        self.side = side

    @property
    def short_description(self):
        side = self.side or "?"
        if not self.resolved:
            return "intron-retention (%s side, predicted)" % side
        if self.aa_mutation_start_offset is None or not self.aa_ref:
            return "intron-retention (%s side)" % side
        return "intron-retention (%s side, p.%s%dfs*)" % (
            side, self.aa_ref[0], self.aa_mutation_start_offset + 1)


class CrypticSpliceSiteEffect(SpliceMechanismEffect):
    """Base for cryptic donor / acceptor effects. Carries the
    cryptic motif's position and score, plus the resulting exon
    length delta."""

    direction: str | None = None  # "donor" / "acceptor", set by subclass

    def __init__(self, variant, transcript, splice_signal, affected_exon,
                 cryptic_genomic_position=None, motif_score=None,
                 exon_length_delta=None,
                 mutant_transcript=None,
                 protein_effect=None,
                 aa_ref=None, aa_alt=None,
                 aa_mutation_start_offset=None,
                 aa_mutation_end_offset=None,
                 mutant_protein_sequence=None):
        SpliceMechanismEffect.__init__(
            self, variant, transcript, splice_signal,
            protein_effect=protein_effect,
            mutant_transcript=mutant_transcript,
            aa_ref=aa_ref, aa_alt=aa_alt,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_mutation_end_offset=aa_mutation_end_offset,
            mutant_protein_sequence=mutant_protein_sequence)
        self.affected_exon = affected_exon
        self.cryptic_genomic_position = cryptic_genomic_position
        self.motif_score = motif_score
        self.exon_length_delta = exon_length_delta

    @property
    def short_description(self):
        if not self.resolved:
            return "cryptic-%s (predicted)" % self.direction
        delta = self.exon_length_delta or 0
        kind = "ext" if delta > 0 else "trunc" if delta < 0 else "shift"
        return "cryptic-%s %s%dbp" % (
            self.direction, kind, abs(delta))


class CrypticDonor(CrypticSpliceSiteEffect):
    """A cryptic GT donor is used instead of the disrupted canonical
    donor; the exon is extended (cryptic site downstream) or
    truncated (cryptic site upstream)."""
    direction = "donor"


class CrypticAcceptor(CrypticSpliceSiteEffect):
    """A cryptic AG acceptor is used instead of the disrupted
    canonical acceptor; the exon is extended (cryptic site upstream)
    or truncated (cryptic site downstream)."""
    direction = "acceptor"


class Exonic(TranscriptMutationEffect):
    """
    Any mutation which affects the contents of an exon (coding region or UTRs)
    """
    pass


class ExonLoss(Exonic):
    """
    Deletion of one or more exons in a transcript.
    """
    def __init__(self, variant, transcript, exons):
        Exonic.__init__(self, variant, transcript)
        self.exons = exons

    def __str__(self):
        return "ExonLoss(%s, %s, %s)" % (
            self.variant,
            self.transcript,
            "+".join(str(exon) for exon in self.exons))

    @property
    def short_description(self) -> str:
        return "exon-loss"

    _default_modifies_protein_sequence = True
    _default_modifies_coding_sequence = True


class ExonicSpliceSite(Exonic, SpliceSite, MultiOutcomeEffect):
    """
    Mutation in the last three nucleotides before an intron
    or in the first nucleotide after an intron.

    Expresses the two plausible outcomes of a splice-adjacent exonic
    variant — splice-signal disruption vs. the underlying coding
    change if splicing proceeds normally — via the
    :class:`MultiOutcomeEffect` protocol (#299, #382).
    :attr:`candidates` yields a 2-tuple of
    :class:`~varcode.effect_candidates.EffectCandidate` objects;
    ``alternate_effect`` stays on the instance as a first-class field
    for back-compat with callers that depended on it.

    Note: as of varcode 6.0, ``ExonicSpliceSite`` no longer surfaces
    as the top-level effect at the user-facing API — splice-disrupting
    variants are always wrapped in
    :class:`~varcode.splice_outcomes.SpliceOutcomeSet`. ``ExonicSpliceSite``
    survives as ``splice_set.disrupted_signal_class`` (a type) and as
    the ``splice_signal`` reference on each candidate (an instance).
    """
    def __init__(self, variant, transcript, exon, alternate_effect):
        Exonic.__init__(self, variant, transcript)
        self.exon = exon
        self.alternate_effect = alternate_effect

    def __str__(self):
        return "ExonicSpliceSite(exon=%s, alternate_effect=%s)" % (
            self.exon,
            self.alternate_effect)

    @property
    def short_description(self) -> str:
        return "exonic-splice-site"

    def _get_mutant_protein_sequence(self) -> str | None:
        """
        TODO: determine when exonic splice variants cause exon skipping
        vs. translation of the underlying modified coding sequence.

        For now just pretending like there is no effect on splicing.
        """
        return self.alternate_effect.mutant_protein_sequence

    def _get_modifies_protein_sequence(self) -> bool | None:
        return self.alternate_effect.modifies_protein_sequence

    def _get_modifies_coding_sequence(self) -> bool | None:
        return self.alternate_effect.modifies_coding_sequence

    def _get_candidates(self) -> tuple[EffectCandidate, ...]:
        """The two plausible outcomes wrapped as
        :class:`~varcode.effect_candidates.EffectCandidate`:

        Position [0] is the splice-disruption outcome (this effect
        itself); position [1] is the coding change if splicing
        proceeds. Returns a 1-tuple when ``alternate_effect`` is None.
        Additional candidates appended via :meth:`add_candidates`
        (e.g. from RNA evidence) follow automatically.
        """
        from ..effect_candidates import EffectCandidate
        base = [EffectCandidate(effect=self, source="varcode")]
        if self.alternate_effect is not None:
            base.append(EffectCandidate(
                effect=self.alternate_effect, source="varcode"))
        return tuple(base)

    @property
    def priority_class(self):
        """ExonicSpliceSite priority is used directly — no
        delegation, since this class IS the splice-adjacent effect
        rather than a wrapper around one.
        """
        return ExonicSpliceSite


class CodingMutation(Exonic):
    """
    Base class for all mutations which result in a modified coding sequence.
    """
    def __str__(self):
        fields = [
            ("variant", self.variant),
            ("transcript_name", self.transcript_name),
            ("transcript_id", self.transcript_id),
            ("effect_description", self.short_description)
        ]
        if self.orf is not None:
            fields.append(("orf_id", self.orf.orf_id))
        fields_str = ", ".join([
            "%s=%s" % (field_name, field_value)
            for (field_name, field_value)
            in fields
        ])
        return "%s(%s)" % (self.__class__.__name__, fields_str)

    _default_modifies_coding_sequence = True


class Silent(CodingMutation):
    """Mutation to an exon of a coding region which doesn't change the
    amino acid sequence.
    """
    def __init__(
            self,
            variant,
            transcript,
            aa_pos,
            aa_ref):
        """
        Parameters
        ----------
        variant : Variant

        transcript : Transcript

        aa_pos : int
            Offset of first synonymous codon in protein sequence

        aa_ref : str
            Reference amino acid(s) at offset
        """
        CodingMutation.__init__(
            self,
            variant=variant,
            transcript=transcript)
        self.aa_pos = aa_pos
        self.aa_ref = str(aa_ref)

    @property
    def short_description(self):
        # HGVS p.= notation for synonymous changes: "p.{aa_ref}{1-indexed pos}="
        # (http://varnomen.hgvs.org/). If aa_ref is empty (the silent
        # change was fully captured by the shared prefix/suffix and we
        # have nothing to name), fall back to a positional-only form.
        if self.aa_ref:
            return "p.%s%d=" % (self.aa_ref, self.aa_pos + 1)
        return "p.%d=" % (self.aa_pos + 1)


class AlternateStartCodon(Silent):
    """Change to the start codon (e.g. ATG>CTG) but without changing the
    starting amino acid from methionine.
    """
    def __init__(
            self,
            variant,
            transcript,
            ref_codon,
            alt_codon,
            aa_ref: str | None = None):
        Silent.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_pos=0,
            aa_ref=transcript.protein_sequence[0] if aa_ref is None else aa_ref)
        self.ref_codon = str(ref_codon)
        self.alt_codon = str(alt_codon)

    @property
    def short_description(self):
        return "alternate-start-codon (%s>%s)" % (
            self.ref_codon, self.alt_codon)


class NonsilentCodingMutation(CodingMutation):
    """
    All coding mutations other than silent codon substitutions
    """

    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_mutation_end_offset,
            aa_ref):
        """
        variant : Variant

        transcript : Transcript

        aa_mutation_start_offset : int
            Offset of first modified amino acid in protein (starting from 0)

        aa_mutation_end_offset : int
            Offset after last mutated amino acid (half-open coordinates)

        aa_ref : str
            Amino acid string of what used to be at aa_mutation_start_offset
            in the wildtype (unmutated) protein.
        """
        CodingMutation.__init__(
            self,
            variant=variant,
            transcript=transcript)
        self.aa_mutation_start_offset = aa_mutation_start_offset
        self.aa_mutation_end_offset = aa_mutation_end_offset
        self.aa_ref = str(aa_ref)

    _default_modifies_protein_sequence = True


class StartGain(TranscriptMutationEffect):
    """A newly created initiation hypothesis on a paired mutant RNA.

    ``discovered_orf`` describes the mutant frame; ``orf`` remains None because
    no reference ORF is asserted. The conditional peptide is available through
    the usual mutant accessor. Baseline translation and actual initiation are
    unknown, so the protein-change flag and original protein remain None.
    """

    _default_modifies_coding_sequence = True
    _default_modifies_protein_sequence = None

    def __init__(self, variant, transcript, discovered_orf, baseline_transcript,
                 mutant_transcript, discovery_policy, hypothesis=None, source_evidence=None):
        super().__init__(variant, transcript)
        self.discovered_orf = discovered_orf
        self.baseline_transcript = baseline_transcript
        self.mutant_transcript = mutant_transcript
        self.discovery_policy = discovery_policy
        self.hypothesis = hypothesis
        self.source_evidence = source_evidence or {}
        self.annotator = "created_orf"
        from ..version import __version__
        self.annotator_version = __version__

    @property
    def short_description(self) -> str:
        orf = self.discovered_orf
        if orf is None:
            return "start-gain (conditional)"
        return "start-gain %s at RNA offset %d (conditional)" % (
            orf.start_codon, orf.start)


class StartLoss(NonsilentCodingMutation):
    """Loss of the annotated start, with available nucleotide context.

    ``ref_codon`` and ``alt_codon`` are in transcript orientation. Their
    offsets are zero-based in the baseline and mutant cDNA, respectively.
    A removed or unmapped start has no alternate codon/offset. Edits and
    sequence remain on ``mutant_transcript`` when available. Neither these
    details nor a codon's internal amino-acid meaning establish initiation.
    The protein outcome remains unknown; alternative-start discovery is
    outside this effect's scope.
    """
    def __init__(
            self,
            variant,
            transcript,
            mutant_transcript=None,
            ref_codon=None,
            alt_codon=None,
            ref_codon_offset=None,
            alt_codon_offset=None,
            codon_change_status=None,
            aa_ref: str | None = None):
        NonsilentCodingMutation.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=0,
            aa_mutation_end_offset=1,
            aa_ref=transcript.protein_sequence[0] if aa_ref is None else aa_ref)
        if codon_change_status is None:
            from .start_loss import start_loss_details
            details = start_loss_details(variant, transcript, mutant_transcript)
            mutant_transcript = details["mutant_transcript"]
            ref_codon = details["ref_codon"] if ref_codon is None else ref_codon
            alt_codon = details["alt_codon"] if alt_codon is None else alt_codon
            ref_codon_offset = details["ref_codon_offset"] if ref_codon_offset is None else ref_codon_offset
            alt_codon_offset = details["alt_codon_offset"] if alt_codon_offset is None else alt_codon_offset
            codon_change_status = details["codon_change_status"]
        self.mutant_transcript = mutant_transcript
        self.ref_codon = ref_codon
        self.alt_codon = alt_codon
        self.ref_codon_offset = ref_codon_offset
        self.alt_codon_offset = alt_codon_offset
        self.codon_change_status = codon_change_status

    def __hash__(self):
        # Transcript-model evidence can contain dictionaries.
        return hash((type(self), self.variant, self.transcript))

    @property
    def alt_codon_amino_acid(self):
        """Internal decoding of the alternate triplet, not an initiated protein residue."""
        if self.alt_codon is None:
            return None
        from .codon_tables import codon_table_for_transcript
        table = codon_table_for_transcript(self.transcript)
        if self.alt_codon in table.stop_codons:
            return "*"
        return table.forward_table.get(self.alt_codon)

    def _get_mutant_protein_sequence(self) -> str | None:
        return None

    @property
    def short_description(self) -> str:
        return "p.%s1? (start-loss)" % (self.aa_ref,)


class KnownAminoAcidChange(NonsilentCodingMutation):
    """
    Coding mutations in which we can predict what the new/mutant protein
    sequence will be.
    """
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset: int,
            aa_ref,
            aa_alt):
        NonsilentCodingMutation.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_mutation_end_offset=aa_mutation_start_offset + len(aa_alt),
            aa_ref=aa_ref)
        self.aa_alt = str(aa_alt)

    def _required_start_offset(self) -> int:
        offset = self.aa_mutation_start_offset
        if offset is None:
            raise ValueError(
                "%s requires a non-None aa_mutation_start_offset" % type(self).__name__)
        return offset

    def _required_original_protein_sequence(self) -> str:
        original = self.original_protein_sequence
        if original is None:
            raise ValueError(
                "%s requires a non-None original_protein_sequence" % type(self).__name__)
        return original

    @property
    def short_description(self) -> str:
        start = self._required_start_offset()
        if len(self.aa_ref) == 0:
            return "p.%dins%s" % (
                start,
                self.aa_alt)
        elif len(self.aa_alt) == 0:
            # Effect offsets stay zero-based; labels name the one-based start.
            return "p.%s%ddel" % (self.aa_ref, start + 1)
        else:
            return "p.%s%d%s" % (
                self.aa_ref,
                start + 1,
                self.aa_alt)

    @cached_property
    def mutant_protein_sequence(self) -> str | None:
        if self.mutant_transcript is not None:
            return self.mutant_transcript.mutant_protein_sequence
        original = self._required_original_protein_sequence()
        start = self._required_start_offset()
        prefix = original[:start]
        suffix = original[start + len(self.aa_ref):]
        return prefix + self.aa_alt + suffix


class Substitution(KnownAminoAcidChange):
    """
    Single amino acid substitution, e.g. BRAF-001 V600E
    """
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_ref,
            aa_alt):
        if len(aa_ref) != 1:
            raise ValueError(
                "Simple substitution can't have aa_ref='%s'" % (aa_ref,))
        if len(aa_alt) != 1:
            raise ValueError(
                "Simple substitution can't have aa_alt='%s'" % (aa_alt,))
        KnownAminoAcidChange.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt=aa_alt)


class ComplexSubstitution(KnownAminoAcidChange):
    """
    In-frame substitution of multiple amino acids, e.g. TP53-002 p.391FY>QQQ
    Can change the length of the protein sequence but since it has
    non-empty ref and alt strings, is more complicated than an insertion or
    deletion alone.
    """
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_ref,
            aa_alt):
        if len(aa_ref) == 1 and len(aa_alt) == 1:
            raise ValueError(
                "ComplexSubstitution can't have aa_ref='%s' and aa_alt='%s'" % (
                    aa_ref, aa_alt))
        KnownAminoAcidChange.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt=aa_alt)


class Insertion(KnownAminoAcidChange):
    """
    In-frame insertion of one or more amino acids.
    """
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_alt):
        KnownAminoAcidChange.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref="",
            aa_alt=aa_alt)


class Deletion(KnownAminoAcidChange):
    """
    In-frame deletion of one or more amino acids.
    """

    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_ref):
        KnownAminoAcidChange.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt="")

    def _format_short_description(self, deletion_notation):
        if len(self.aa_ref) <= 1:
            return self.short_description
        start = self._required_start_offset() + 1
        return "p.%s%d_%s%ddel" % (
            self.aa_ref[0], start,
            self.aa_ref[-1], start + len(self.aa_ref) - 1)


class PrematureStop(KnownAminoAcidChange):
    """In-frame insertion of codons containing a stop codon. May also involve
    insertion/deletion/substitution of other amino acids preceding the stop."""
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            aa_ref="",
            aa_alt="",
            original_protein_sequence: str | None = None):
        """
        Insertion of premature stop codon, possibly preceded by a substitution
        of `aa_ref` amino acids for `aa_alt` alternative residues.

        ``original_protein_sequence`` supplies the patient baseline when it
        differs from the reference transcript; coordinates refer to that protein.
        """
        if original_protein_sequence is not None:
            self._original_protein_sequence = original_protein_sequence
        if "*" in aa_ref:
            raise ValueError(
                ("Unexpected aa_ref = '%s', should only include amino acids "
                 "before the new stop codon.") % aa_ref)
        if "*" in aa_alt:
            raise ValueError(
                ("Unexpected aa_ref = '%s', should only include amino acids "
                 "before the new stop codon.") % aa_alt)
        KnownAminoAcidChange.__init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt=aa_alt)
        self.stop_codon_offset = aa_mutation_start_offset + len(aa_alt)

        original = self._required_original_protein_sequence()
        if self.stop_codon_offset >= len(original):
            raise ValueError(
                ("Premature stop codon cannot be at position %d"
                 " since the original protein of %s has length %d") % (
                    self.stop_codon_offset,
                    transcript,
                    len(original)))

    @property
    def short_description(self) -> str:
        # When aa_ref is empty, the mutation is an insertion (no residues
        # were replaced) — use HGVS "ins" notation so the description is
        # unambiguous. Otherwise fall through to "{ref}{pos}{alt}*".
        start = self._required_start_offset()
        if len(self.aa_ref) == 0:
            return "p.%dins%s*" % (
                start + 1,
                self.aa_alt)
        return "p.%s%d%s*" % (
            self.aa_ref,
            start + 1,
            self.aa_alt)

    @cached_property
    def mutant_protein_sequence(self) -> str | None:
        if self.mutant_transcript is not None:
            return self.mutant_transcript.mutant_protein_sequence
        prefix = self._required_original_protein_sequence()[:self._required_start_offset()]
        return prefix + self.aa_alt


class StopLoss(KnownAminoAcidChange):
    def __init__(
            self,
            variant,
            transcript,
            aa_ref,
            aa_alt,
            original_protein_sequence: str | None = None):
        """Extend the supplied baseline protein, or the reference when omitted."""
        if original_protein_sequence is not None:
            self._original_protein_sequence = original_protein_sequence
        # StopLoss assumes that we deleted some codons ending with a
        # stop codon
        if "*" in aa_ref:
            raise ValueError(
                "StopLoss aa_ref '%s' should not contain '*'" % (
                    aa_ref,))
        if len(aa_alt) == 0:
            raise ValueError(
                "If no amino acids added by StopLoss then it should be Silent")
        # subtract 1 for the stop codon
        n_ref_amino_acids = len(aa_ref)
        protein_length = len(original_protein_sequence if original_protein_sequence is not None
                             else transcript.protein_sequence)
        aa_mutation_start_offset = protein_length - n_ref_amino_acids
        KnownAminoAcidChange.__init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_alt=aa_alt,
            aa_ref=aa_ref)

    @property
    def short_description(self) -> str:
        return "p.%s*%d%s (stop-loss)" % (
            self.aa_ref,
            self._required_start_offset() + 1,
            self.aa_alt)


class FrameShift(KnownAminoAcidChange):
    def __init__(
            self,
            variant,
            transcript,
            aa_mutation_start_offset,
            shifted_sequence,
            original_protein_sequence: str | None = None):
        """Frameshift mutation preserves all the amino acids up to
        aa_mutation_start_offset and then replaces the rest of the protein with
        new (frameshifted) sequence. Unlike an insertion, where we denote with
        aa_ref as the chracter before the variant sequence, a frameshift starts
        at aa_ref.

        ``original_protein_sequence`` is the patient baseline when supplied.
        """
        if original_protein_sequence is not None:
            self._original_protein_sequence = original_protein_sequence
        protein = (original_protein_sequence if original_protein_sequence is not None
                   else transcript.protein_sequence)
        aa_ref = protein[aa_mutation_start_offset:]
        KnownAminoAcidChange.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=aa_mutation_start_offset,
            aa_ref=aa_ref,
            aa_alt=shifted_sequence)

    @property
    def shifted_sequence(self):
        return self.aa_alt

    @cached_property
    def mutant_protein_sequence(self) -> str | None:
        if self.mutant_transcript is not None:
            return self.mutant_transcript.mutant_protein_sequence
        original_aa_sequence = self._required_original_protein_sequence()
        prefix = original_aa_sequence[:self._required_start_offset()]
        return prefix + self.shifted_sequence

    @property
    def short_description(self) -> str:
        return "p.%s%dfs" % (
            self.aa_ref[0],
            self._required_start_offset() + 1)


class FrameShiftTruncation(PrematureStop, FrameShift):
    """
    A frame-shift mutation which immediately introduces a stop codon.
    """
    def __init__(
            self,
            variant,
            transcript,
            stop_codon_offset,
            aa_ref="",
            original_protein_sequence: str | None = None):
        PrematureStop.__init__(
            self,
            variant=variant,
            transcript=transcript,
            aa_mutation_start_offset=stop_codon_offset,
            aa_ref=aa_ref,
            aa_alt="",
            original_protein_sequence=original_protein_sequence)

    @property
    def shifted_sequence(self):
        return self.aa_alt

    @property
    def short_description(self) -> str:
        return "p.%s%dfs*" % (
            self.aa_ref,
            self._required_start_offset() + 1)


# =====================================================================
# Structural-variant effects (#252 / PR 10)
#
# These describe consequences of SVs at the transcript level. A single
# SV commonly produces multiple plausible consequences (fusion vs exon
# loss vs readthrough), so each of these is a ``MultiOutcomeEffect``
# subclass exposing a :attr:`candidates` tuple per #299 / #382.
#
# The classes are deliberately thin wrappers — the important interface
# is :attr:`candidates`, which carries :class:`EffectCandidate` objects
# with source / evidence provenance so external tools (RNA evidence,
# SpliceAI, long-read assembly) can annotate them without subclassing.
# =====================================================================


class StructuralVariantEffect(TranscriptMutationEffect, MultiOutcomeEffect):
    """Base class for effects of a :class:`StructuralVariant` on a
    specific transcript. Local models contain consequence effects in
    :attr:`candidates`; their description follows the first prediction and
    their priority follows the highest-priority candidate. Gene fusions and
    translocations retain their specialized descriptions.

    :attr:`candidates` always returns a tuple of
    :class:`~varcode.effect_candidates.EffectCandidate` objects.
    Primary candidates carry ``source="varcode"``; cryptic-exon, splice,
    and RNA candidates are appended via :meth:`add_candidates` with their
    respective ``source`` tags (``"varcode_motif"``, ``"varcode_splice"``,
    ``"rna"``, etc.).
    """

    candidate_evidence = None

    @property
    def short_description(self):
        if self.primary_effects and self.primary_effects[0] is not self:
            return self.primary_effects[0].short_description
        return "unresolved-structural-transcript"

    def _format_short_description(self, deletion_notation):
        # Subclasses with dedicated labels (e.g. sv-deletion) retain them.
        if type(self).short_description is not StructuralVariantEffect.short_description:
            return self.short_description
        if self.primary_effects and self.primary_effects[0] is not self:
            return self.primary_effects[0].format_short_description(
                deletion_notation=deletion_notation)
        return self.short_description

    @property
    def priority_class(self):
        own = type(self) if type(self) is not StructuralVariantEffect else Unresolved
        if not self.primary_effects or self.primary_effects[0] is self:
            return own
        from .effect_ordering import effect_priority, transcript_effect_priority_dict

        def rank(candidate):
            # A candidate that is this set stands for its own class; ranking
            # it through this property would never terminate.
            if candidate.effect is self:
                return transcript_effect_priority_dict.get(own, -1)
            return effect_priority(candidate.effect)

        highest = max(self.candidates, key=rank).effect
        if highest is self:
            return own
        return getattr(highest, "priority_class", None) or type(highest)

    def __init__(
            self, variant, transcript,
            primary_effects=None, mutant_transcript=None):
        TranscriptMutationEffect.__init__(self, variant, transcript)
        self.primary_effects: tuple[MutationEffect, ...] = (
            tuple(primary_effects) if primary_effects is not None
            else (self,))
        self.mutant_transcript = mutant_transcript

    def add_primary_effects(self, effects: Iterable[MutationEffect]) -> None:
        """Append additional primary SV classifications (e.g. partner fusions or span effects)."""
        incoming = tuple(effects)
        if incoming:
            self.primary_effects = self.primary_effects + incoming

    def _get_mutant_protein_sequence(self) -> str | None:
        """Protein carried by the realized structural transcript, if any."""
        return getattr(
            self.mutant_transcript, "mutant_protein_sequence", None)

    def _get_modifies_coding_sequence(self) -> bool | None:
        """True if any candidate changes CDS, None if unresolved, else False.

        Protein equality alone does not establish coding-sequence equality.
        """
        from .sequence_change import modification_status
        return modification_status(self, "modifies_coding_sequence")

    def _get_modifies_protein_sequence(self) -> bool | None:
        """True if any candidate changes protein, None if unresolved, else False.

        A set with an unchanged primary and a changed alternative returns True.
        False requires every candidate to be known unchanged; None is not silent.
        """
        from .sequence_change import modification_status
        return modification_status(self, "modifies_protein_sequence")

    def _get_candidates(self) -> tuple[EffectCandidate, ...]:
        """Primary :class:`~varcode.EffectCandidate` entries for :attr:`primary_effects`.

        Additional motif, splice, and RNA candidates appended via
        :meth:`add_candidates` follow automatically.
        """
        from ..effect_candidates import EffectCandidate
        sv_type = getattr(self.variant, "sv_type", None)
        base_evidence = dict(self.candidate_evidence or {})
        if sv_type is not None:
            base_evidence["sv_type"] = sv_type
        return tuple(
            EffectCandidate(
                effect=effect,
                source="varcode",
                evidence=base_evidence)
            for effect in self.primary_effects)


class LargeDeletion(StructuralVariantEffect):
    """A deletion (``<DEL>`` / ``<CN0>``) that removes one or more
    exons — or an entire gene. Carries the list of affected exons
    for downstream analysis; the single ``outcome`` is this effect
    itself (callers that want to add RNA evidence or SpliceAI
    scoring construct additional outcomes and wrap the result)."""

    @property
    def short_description(self) -> str:
        return "sv-deletion"

    def __init__(
            self, variant, transcript, affected_exons,
            primary_effects=None, mutant_transcript=None):
        StructuralVariantEffect.__init__(
            self, variant, transcript,
            primary_effects=primary_effects,
            mutant_transcript=mutant_transcript)
        self.affected_exons = tuple(affected_exons)

    def __str__(self):
        return "LargeDeletion(%s, %s, exons=%d)" % (
            self.variant.short_description,
            self.transcript.id if self.transcript is not None else "?",
            len(self.affected_exons))


class LargeDuplication(StructuralVariantEffect):
    """A tandem duplication (``<DUP>``) overlapping exons.
    Biologically may produce a copy-number increase, a fused
    reading frame if junctions land in-frame, or a regulatory
    effect. Varcode reports the affected exons and leaves scoring
    to downstream tools."""

    @property
    def short_description(self) -> str:
        return "sv-duplication"

    def __init__(
            self, variant, transcript, affected_exons,
            primary_effects=None, mutant_transcript=None):
        StructuralVariantEffect.__init__(
            self, variant, transcript,
            primary_effects=primary_effects,
            mutant_transcript=mutant_transcript)
        self.affected_exons = tuple(affected_exons)


class Inversion(StructuralVariantEffect):
    """An inversion (``<INV>``) that flips a stretch of a transcript.
    Depending on whether breakpoints fall in exons or introns, the
    consequence is very different (exonic: likely disruptive; purely
    intronic: may or may not affect splicing). Reported as a single
    outcome here; subclasses / callers can enrich with cryptic-
    splice candidates per PR 11."""

    @property
    def short_description(self) -> str:
        return "sv-inversion"


class GeneFusion(StructuralVariantEffect):
    """A structural variant joining ``transcript`` to a sense-oriented
    transcript in another gene — the canonical fusion shape. Breakends
    produce it, and so do deletions, duplications and inversions with
    one end in each gene. A 3' partner that starts just past a junction
    end, rather than spanning it, is joined from its second exon and is
    only ever a further candidate behind the primary classification.

    The default annotator retains every compatible annotated partner isoform
    and junction in :attr:`candidates`, with no count cap. The first is the
    historical primary; producer order is not a probability ranking. Distinct
    transcript pairs remain separate even when their proteins match. Missing
    sequence leaves that candidate's transcript/protein unresolved.

    ``transcript`` is the transcript being annotated and
    ``partner_transcript`` the other one; ``five_prime_transcript``
    and ``three_prime_transcript`` say which is which (by default
    ``transcript`` is the 5' partner). When the annotator has enough
    context, :attr:`mutant_transcript` is built from
    :class:`ReferenceSegment` entries describing the fused
    allele. Predicting the exact fused-protein sequence requires
    knowing which exons are retained, which typically needs RNA
    evidence — outcomes beyond "this is a plausible fusion" are
    imported with :func:`varcode.make_fusion_outcome` or supplied by tools
    attaching :class:`EffectCandidate` objects with their own ``source`` tag.
    An RNA-observed gene fusion need not encode a complete protein.
    """

    @property
    def short_description(self) -> str:
        return "sv-gene-fusion"

    def __init__(
            self, variant, transcript, partner_transcript,
            mutant_transcript=None, primary_effects=None,
            five_prime_transcript=None, three_prime_transcript=None):
        StructuralVariantEffect.__init__(
            self, variant, transcript,
            primary_effects=primary_effects,
            mutant_transcript=mutant_transcript)
        self.partner_transcript = partner_transcript
        self.five_prime_transcript = (
            transcript if five_prime_transcript is None
            else five_prime_transcript)
        self.three_prime_transcript = (
            partner_transcript if three_prime_transcript is None
            else three_prime_transcript)


class TranslocationToIntergenic(StructuralVariantEffect):
    """A breakend that doesn't form a gene fusion: its mate lies in
    intergenic space, or the join meets a gene in an orientation that
    can't read sense-to-sense into it (e.g. two genes' 5' ends joined
    head to head). A junction with both ends in the transcript's gene is
    annotated as the local DEL, DUP or INV it describes instead. The downstream consequence depends on whether the
    sequence beyond the breakpoint contains cryptic splice / ORF
    signals. Fusions with a coding transcript that starts just past the
    mate follow as further candidates, and PR 11's cryptic-exon
    enumerator adds candidate outcomes when applicable."""

    @property
    def short_description(self) -> str:
        return "sv-translocation-intergenic"


class TranslocationToUnannotatedContig(TranslocationToIntergenic):
    """A breakend that forms no known gene fusion because its mate lies on
    a contig the annotation doesn't cover: a decoy, EBV, an unplaced or alt
    scaffold, or hg19's chrM under Ensembl's GRCh37 (whose MT is the rCRS,
    at other positions). Whether the mate is intergenic is unknown, so a
    fusion with it is not ruled out. ``partner_contig`` names the mate's
    contig. The local side is annotated as for
    :class:`TranslocationToIntergenic`, which this specializes."""

    @property
    def short_description(self) -> str:
        return "sv-translocation-unannotated-contig"

    def __init__(
            self, variant, transcript, partner_contig,
            mutant_transcript=None, primary_effects=None):
        TranslocationToIntergenic.__init__(
            self, variant, transcript,
            primary_effects=primary_effects,
            mutant_transcript=mutant_transcript)
        self.partner_contig = partner_contig


class CrypticExonCandidate(MutationEffect):
    """A region where an SV has brought novel sequence into range
    of the transcript, and motif scoring flags a plausible new
    splice acceptor / donor pair. Produced by PR 11's cryptic-exon
    enumerator; attached as additional :class:`EffectCandidate` entries on
    SV effects rather than standalone.

    Not a :class:`TranscriptMutationEffect` because the candidate
    region may not overlap any existing transcript — it's a *new*
    exon hypothesis. Carries the contig / interval and the motif
    scores as plain fields; external predictors (SpliceAI,
    Pangolin) attach their own scores via the enclosing
    :class:`EffectCandidate.evidence` dict.
    """

    @property
    def short_description(self) -> str:
        return "sv-cryptic-exon-candidate"

    _default_modifies_coding_sequence = None
    _default_modifies_protein_sequence = None

    def __init__(
            self, variant, contig, interval_start, interval_end,
            donor_score=None, acceptor_score=None):
        MutationEffect.__init__(self, variant)
        self.contig = contig
        self.interval_start = interval_start
        self.interval_end = interval_end
        self.donor_score = donor_score
        self.acceptor_score = acceptor_score


# =====================================================================
# Haplotype effects — joint consequence of multiple cis variants on
# the same transcript (#269).
#
# When a phase resolver (VCF PS, Isovar assembly, long-read haplotype
# caller) confirms two or more variants sit on the same physical copy
# of a transcript, their combined effect on the protein isn't just the
# sum of per-variant effects — they can share a codon, one can rescue
# another's frameshift, etc. HaplotypeEffect carries the joint
# MutantTranscript (via apply_variants_to_transcript) and coexists
# with the per-variant effects on the same EffectCollection so
# downstream consumers can pick the granularity they need.
# =====================================================================


class HaplotypeEffect(TranscriptMutationEffect, MultiOutcomeEffect):
    """Joint effect of two or more cis variants on the same
    transcript (#269).

    Emitted by :meth:`VariantCollection.effects` when a
    ``phase_resolver`` (VCF ``PS`` tag, Isovar assembly, or any
    other :class:`PhaseResolver`) groups cis variants together and
    :func:`varcode.mutant_transcript.apply_variants_to_transcript`
    can build a combined mutant cDNA. The per-variant effects stay
    on the collection — this is additive, not a replacement.

    Downstream consumers that want haplotype-level context (peptide
    prediction, ASE, compound het interpretation) read the
    :attr:`mutant_transcript` here. Consumers that want per-variant
    HGVS or annotation granularity iterate the per-variant effects
    as before.
    """

    def __init__(
            self, variants, transcript, mutant_transcript,
            phase_source=None):
        # The HaplotypeEffect is scoped to a set of variants — the
        # :class:`MutationEffect.variant` slot takes the first one
        # as a back-compat anchor so existing single-variant
        # consumers don't blow up, but real consumers read
        # :attr:`variants`.
        assert len(variants) >= 2, (
            "HaplotypeEffect requires ≥ 2 cis variants; got %d" % len(variants))
        TranscriptMutationEffect.__init__(self, variants[0], transcript)
        self.variants = tuple(variants)
        self.mutant_transcript = mutant_transcript
        # Which resolver produced the phase grouping (e.g. "vcf_ps"
        # for DNA-PS-tag phasing, "read_phasing" for an RNA / long-
        # read source). Useful for downstream filtering and for
        # auditing whether the cis call came from DNA-only phasing
        # or RNA-assembly evidence.
        self.phase_source = phase_source

    def _get_candidates(self) -> tuple[EffectCandidate, ...]:
        """Single-outcome wrapping: the haplotype IS the outcome,
        lifted into the uniform :class:`EffectCandidate` shape so
        consumers iterate :attr:`candidates` the same way across all
        :class:`MutationEffect` instances. Additional candidates appended
        via :meth:`add_candidates` (e.g. from RNA evidence) follow
        automatically.
        """
        from ..effect_candidates import EffectCandidate
        evidence = {}
        if self.phase_source is not None:
            evidence["phase_source"] = self.phase_source
        return (EffectCandidate(
            effect=self, source="varcode", evidence=evidence),)

    def _get_mutant_protein_sequence(self) -> str | None:
        """Joint protein translated from the haplotype's cDNA, or
        ``None`` if the edits don't land after the CDS start."""
        return getattr(
            self.mutant_transcript, "mutant_protein_sequence", None)

    def __str__(self):
        return "HaplotypeEffect(%d variants, transcript=%s, source=%s)" % (
            len(self.variants),
            self.transcript.name if self.transcript is not None else "?",
            self.phase_source or "unknown")

    @property
    def short_description(self):
        """HGVS haplotype notation: ``[c.1A>T;c.5G>C]`` for cis.

        Per-variant short descriptions are synthesized from each
        variant's genomic coords when the per-variant effect isn't
        available here — the transcript-level HGVS would require
        re-running the per-variant classifier, which
        ``HaplotypeEffect`` deliberately doesn't do (that work lives
        on the per-variant effects that coexist in the collection).
        """
        parts = [v.short_description for v in self.variants]
        return "[" + ";".join(parts) + "]"


# =====================================================================
# Phase candidate set — possibility set when somatic + germline
# share a window and phase between them is unknown (#268).
#
# Sibling of HaplotypeEffect: that one captures the *known-cis* case
# (multiple variants composed on the same haplotype), this one
# captures the *unknown-phase* case (multiple hypotheses, each with
# its own resulting effect). Both inherit the same MultiOutcomeEffect
# protocol so consumers iterate ``outcomes`` uniformly. This isn't a
# wrapper around a reference-relative effect — it's the primary
# effect class for the unknown-phase case, just like
# SpliceOutcomeSet is for splice ambiguity.
# =====================================================================


class PhaseCandidateSet(TranscriptMutationEffect, MultiOutcomeEffect):
    """Possibility set across phase hypotheses when a somatic variant
    and one or more germline variants share a window on a transcript
    and phase between them is unknown.

    The somatic effect at this locus depends on which haplotype the
    somatic landed on, and without phase data we can't say which.
    The honest output is the set of plausible effects, one per
    hypothesis. ``most_likely`` returns the highest-priority
    candidate; :attr:`candidates` exposes the full set with
    per-hypothesis evidence keys (``phase_state``, ``haplotype``,
    ``germline_variants``) so consumers — including RNA-evidence
    integrations from #259 — can align across axes.

    Emitted by :func:`varcode.germline.predict_germline_aware_effect`.
    Sibling of :class:`HaplotypeEffect`: that one captures the
    *known-cis* multi-variant case; this one captures the
    *unknown-phase* possibility set.
    """

    def __init__(
            self,
            variant,
            transcript,
            candidates,
            hypotheses,
            germline_variants):
        TranscriptMutationEffect.__init__(self, variant, transcript)
        if len(candidates) != len(hypotheses):
            raise ValueError(
                "PhaseCandidateSet needs one candidate per "
                "hypothesis; got %d candidates and %d hypotheses." % (
                    len(candidates), len(hypotheses)))
        if not candidates:
            raise ValueError(
                "PhaseCandidateSet requires at least one candidate.")
        self.candidate_effects = tuple(candidates)
        self.hypotheses = tuple(hypotheses)
        self.germline_variants = tuple(germline_variants)

    def _get_candidates(self) -> tuple[EffectCandidate, ...]:
        """One :class:`EffectCandidate` per hypothesis, sorted by
        underlying effect priority (most-severe first), carrying the
        phase metadata needed to align with RNA-evidence outcomes
        (#259, #382).

        ``evidence`` keys: ``phase_state`` (``"phased"`` /
        ``"implicit"`` / ``"unknown"``),
        ``haplotype`` (opaque tag), ``germline_variants`` (tuple of
        the cis germline variants on that hypothesis's haplotype).
        """
        from ..effect_candidates import EffectCandidate
        from .effect_ordering import effect_priority
        paired = sorted(
            zip(self.candidate_effects, self.hypotheses),
            key=lambda pair: -effect_priority(pair[0]))
        return tuple(
            EffectCandidate(
                effect=candidate,
                source="varcode_germline",
                evidence={
                    "phase_state": hypothesis.phase_state,
                    "haplotype": hypothesis.haplotype,
                    "germline_variants": tuple(hypothesis.cis),
                })
            for candidate, hypothesis in paired)

    # Note on accessor semantics for this class: the base-class
    # "producer order" contract is replaced here by
    # "sorted highest-impact-first" via
    # :func:`effect_priority`. Consequently
    # :attr:`most_likely_candidate` (== ``candidates[0]``) and
    # :attr:`highest_priority_candidate` coincide for this class —
    # both return the most protein-disruptive hypothesis, which is
    # also the conservative single-effect representation downstream
    # consumers see in ``short_description`` etc.

    @property
    def short_description(self):
        """``"?<most_likely_effect>"`` — the leading ``?`` flags the
        ambiguity. Consumers wanting the full possibility set read
        :attr:`candidates`."""
        top = self.most_likely_effect
        assert top is not None
        return "?" + top.short_description

    def _format_short_description(self, deletion_notation):
        top = self.most_likely_effect
        assert top is not None
        return "?" + top.format_short_description(
            deletion_notation=deletion_notation)

    @property
    def priority_class(self):
        """Rank as the most severe hypothesis, like other outcome sets."""
        highest = self.highest_priority_effect
        assert highest is not None
        return getattr(highest, "priority_class", None) or type(highest)

    def _get_mutant_protein_sequence(self) -> str | None:
        """Most-likely candidate's mutant protein. Consumers iterating
        over hypotheses pull per-candidate sequences from
        :attr:`candidates`."""
        return getattr(self.most_likely_effect, "mutant_protein_sequence", None)

    def _get_modifies_protein_sequence(self) -> bool | None:
        return getattr(
            self.most_likely_effect, "modifies_protein_sequence", False)

    def _get_modifies_coding_sequence(self) -> bool | None:
        return getattr(
            self.most_likely_effect, "modifies_coding_sequence", False)

    def __str__(self):
        return (
            "PhaseCandidateSet(variant=%s, transcript=%s, "
            "%d hypotheses, germline=%s)" % (
                self.variant,
                self.transcript.name if self.transcript is not None else "?",
                len(self.hypotheses),
                ",".join(g.short_description for g in self.germline_variants)))
