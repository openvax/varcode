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

"""Base mutation effect classes and non-coding transcript consequences."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from functools import cached_property
from typing import TYPE_CHECKING, Any, Literal

from serializable import Serializable

from .common import validate_deletion_notation

if TYPE_CHECKING:
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



class Exonic(TranscriptMutationEffect):
    """
    Any mutation which affects the contents of an exon (coding region or UTRs)
    """
    pass



__all__ = [
    "MutationEffect",
    "MultiOutcomeEffect",
    "Intergenic",
    "Intragenic",
    "TranscriptMutationEffect",
    "GermlineAlleleOverlap",
    "Failure",
    "Unresolved",
    "HypothesisLimit",
    "NoncodingTranscript",
    "IncompleteTranscript",
    "FivePrimeUTR",
    "ThreePrimeUTR",
    "Intronic",
    "Exonic",
]
