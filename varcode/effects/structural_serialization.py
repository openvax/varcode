"""Versioned effect graphs, including self-referential realized candidates.

Effects and transcript-model records use graph references. Ordinary effects
retain their existing reference formats. The ORF schema additionally snapshots
reference RNAs and coding models, independently of an annotation database.
"""

from dataclasses import dataclass, fields, replace
from inspect import Parameter, signature

from serializable import DataclassSerializable

from ..effect_candidates import EffectCandidate


@dataclass(frozen=True)
class _EffectReference(DataclassSerializable):
    index: int


# Explicit post-construction annotation state. Do not archive caches, logger
# handles, or arbitrary implementation attributes from __dict__.
_ATTACHED_FIELDS = (
    "mutant_transcript", "baseline_transcript", "affected_exons", "candidate_evidence", "_cryptic_candidates",
    "_splice_candidates", "_extra_candidates", "_original_protein_sequence",
    "_mutant_protein_sequence", "_modifies_coding_sequence",
    "_modifies_protein_sequence", "modifies_coding_sequence",
    "modifies_protein_sequence", "variants", "phase_source", "annotator",
    "annotator_version", "observed_mutant_transcripts",
    "candidates", "realized_candidates", "most_likely_candidate",
    "highest_priority_candidate", "hypotheses", "mechanism_rank",
    "probability", "splice_signal", "germline_unknown",
    "germline_phase_state", "germline_variants_in_window",
    "orf_annotation",
)

# A comparison describes one outcome, independently of candidate containers
# and back-references to an enclosing annotation. Its inert projection retains
# the mechanism, transcript model and evidence but does not archive siblings.
_ALTERNATIVE_FIELDS = {
    "primary_effects", "candidates", "hypotheses", "realized_candidates",
    "_cryptic_candidates", "_splice_candidates", "_extra_candidates",
    "most_likely_candidate", "highest_priority_candidate", "orf_annotation",
}


def _model_types():
    from ..orf_discovery import DiscoveredORF, ORFDiscoveryPolicy
    from ..mutant_transcript import MutantTranscript, ReferenceSegment, TranscriptEdit, _AssembledAllele
    from ..effect_hypotheses import (
        ClassifiedOutcome, EffectHypothesis, RealizedEffectCandidate,
        RealizedTranscriptProduct, SpliceOption, SplicePlan, SpliceSiteKey,
    )
    return (DiscoveredORF, ORFDiscoveryPolicy,
            MutantTranscript, ReferenceSegment, TranscriptEdit, _AssembledAllele,
            ClassifiedOutcome, EffectHypothesis, RealizedEffectCandidate,
            RealizedTranscriptProduct, SpliceOption, SplicePlan, SpliceSiteKey)


def structural_effect_to_dict(root, *, schema="structural_effect_schema"):
    """Archive a native graph, or project one inert ``orf_comparison_schema`` outcome."""
    from .effect_classes import MutationEffect, PhaseCandidateSet, StructuralVariantEffect
    from ..variant import Variant
    from ..structural_variant import StructuralVariant

    indices = {}
    nodes = []
    model_types = _model_types()
    comparison = schema == "orf_comparison_schema"
    orf = schema == "orf_effect_schema" or comparison
    special_types = (Variant, StructuralVariant)
    if orf:
        from orfdata import ReferenceORF
        from pyensembl import Transcript
        from ..orf import _ORFTranscript
        from ..orf_serialization import (
            _ReferenceGenome, _ReferenceTranscript, orf_graph_types, record_fields)
        model_types += orf_graph_types()
        special_types = (ReferenceORF, _ORFTranscript, Variant, StructuralVariant)

    def encode(value):
        if orf and isinstance(value, Transcript):
            return _ReferenceTranscript.capture(value)
        if isinstance(value, MutationEffect) or type(value) in model_types + special_types:
            key = id(value)
            if key not in indices:
                indices[key] = len(nodes)
                nodes.append(None)  # reserve before following self links
                if orf and isinstance(value, ReferenceORF):
                    state = (dict(snapshot=value.to_dict(),
                                  transcript=_ReferenceTranscript.capture(value.transcript))
                             if value.record is None else dict(
                                 record=value.record, transcript=_ReferenceTranscript.capture(value.transcript),
                                 context=value.context, any_transcript=True))
                    nodes[indices[key]] = (type(value), encode(state), {})
                    return _EffectReference(indices[key])
                if orf and isinstance(value, _ORFTranscript):
                    nodes[indices[key]] = (type(value), encode(dict(reference=value.reference_orf)), {})
                    return _EffectReference(indices[key])
                if isinstance(value, Variant):
                    state = value.to_dict()
                    if orf:
                        state = {**state, "genome": _ReferenceGenome.capture(value.genome)}
                    nodes[indices[key]] = (type(value), encode(state), {})
                    return _EffectReference(indices[key])
                if type(value) in model_types:
                    state = record_fields(value) if orf else {f.name: getattr(value, f.name) for f in fields(value)}
                    nodes[indices[key]] = (type(value), encode(state), {})
                    return _EffectReference(indices[key])
                state = {}
                for name, parameter in signature(type(value).__init__).parameters.items():
                    if name == "self" or parameter.kind in (
                            Parameter.VAR_POSITIONAL, Parameter.VAR_KEYWORD):
                        continue
                    if comparison and name in _ALTERNATIVE_FIELDS:
                        continue
                    if name == "primary_effects" and isinstance(value, StructuralVariantEffect):
                        state[name] = value._primary_effects
                    elif isinstance(value, PhaseCandidateSet) and name in ("candidates", "hypotheses"):
                        # Preserve constructor order and native effects, rather
                        # than the public priority-sorted EffectCandidate view.
                        state[name] = (value._candidates_raw if name == "candidates" else value._hypotheses)
                    elif name == "candidates" and hasattr(value, "_candidates"):
                        state[name] = value._candidates
                    elif hasattr(value, name):
                        state[name] = getattr(value, name)
                    elif parameter.default is Parameter.empty:
                        raise ValueError("Missing effect constructor field: %s" % name)
                attached = {name: value.__dict__[name] for name in _ATTACHED_FIELDS
                            if name in value.__dict__ and name not in state
                            and not (comparison and name in _ALTERNATIVE_FIELDS)}
                nodes[indices[key]] = (type(value), encode(state), encode(attached))
            return _EffectReference(indices[key])
        if isinstance(value, EffectCandidate):
            return replace(value, effect=encode(value.effect), evidence=encode(value.evidence))
        if isinstance(value, dict):
            return {key: encode(item) for key, item in value.items()}
        if isinstance(value, (tuple, list)):
            return type(value)(encode(item) for item in value)
        return value

    encode(root)
    return {schema: 1, "nodes": tuple(nodes)}


def structural_effect_from_dict(root_class, state, *, schema="structural_effect_schema"):
    from .effect_classes import MutationEffect
    from ..variant import Variant
    from ..structural_variant import StructuralVariant

    if schema == "orf_comparison_schema":
        raise ValueError("Comparison projections are inert, not reconstructable effect graphs")
    if state.get(schema) != 1:
        raise ValueError("Unsupported effect serialization schema")
    nodes = state["nodes"]
    if not nodes or nodes[0][0] is not root_class:
        raise ValueError("Effect root class does not match archive")
    model_types = _model_types()
    special_types = (Variant, StructuralVariant)
    if schema == "orf_effect_schema":
        from orfdata import ReferenceORF
        from ..orf import _ORFTranscript
        from ..orf_serialization import orf_graph_types
        model_types += orf_graph_types()
        special_types = (ReferenceORF, _ORFTranscript, Variant, StructuralVariant)
    for cls, _, attached in nodes:
        if not isinstance(cls, type) or not (issubclass(cls, MutationEffect) or cls in model_types + special_types):
            raise ValueError("Effect graph contains an unsupported class")
        if cls in model_types and attached:
            raise ValueError("Unexpected transcript-model attachment")
    objects = [object.__new__(cls) for cls, _, _ in nodes]
    started = set()

    def decode(value):
        if isinstance(value, _EffectReference):
            index = value.index
            if type(index) is not int or not 0 <= index < len(nodes):
                raise ValueError("Invalid effect graph reference")
            if index not in started:
                started.add(index)
                cls, fields, attached = nodes[index]
                values = decode(fields)
                if schema == "orf_effect_schema" and cls is ReferenceORF and "snapshot" in values:
                    restored = cls.from_dict(values["snapshot"], values["transcript"])
                    for name, item in restored.__dict__.items():
                        object.__setattr__(objects[index], name, item)
                else:
                    cls.__init__(objects[index], **values)
                for name, item in decode(attached).items():
                    if name not in _ATTACHED_FIELDS:
                        raise ValueError("Unknown effect attachment: %s" % name)
                    setattr(objects[index], name, item)
            return objects[index]
        if isinstance(value, EffectCandidate):
            return replace(value, effect=decode(value.effect), evidence=decode(value.evidence))
        if isinstance(value, dict):
            return {key: decode(item) for key, item in value.items()}
        if isinstance(value, (tuple, list)):
            return type(value)(decode(item) for item in value)
        return value

    return decode(_EffectReference(0))
