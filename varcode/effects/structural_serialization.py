"""Versioned structural-effect graphs, including self-referential candidates.

Only effect/candidate links use graph references. Variants, transcripts, and
MutantTranscript models retain their existing Serializable representations.
This is an object archive, not a portable annotation database.
"""

from dataclasses import dataclass, replace
from inspect import Parameter, signature

from serializable import DataclassSerializable

from ..effect_candidates import EffectCandidate


@dataclass(frozen=True)
class _EffectReference(DataclassSerializable):
    index: int


# Explicit post-construction annotation state. Do not archive caches, logger
# handles, or arbitrary implementation attributes from __dict__.
_ATTACHED_FIELDS = (
    "mutant_transcript", "affected_exons", "candidate_evidence", "_cryptic_candidates",
    "_splice_candidates", "_extra_candidates", "_original_protein_sequence",
    "_mutant_protein_sequence", "_modifies_coding_sequence",
    "_modifies_protein_sequence", "modifies_coding_sequence",
    "modifies_protein_sequence", "variants", "phase_source", "annotator",
    "annotator_version", "observed_mutant_transcripts",
)


def structural_effect_to_dict(root):
    from .effect_classes import MutationEffect, StructuralVariantEffect

    indices = {}
    nodes = []

    def encode(value):
        if isinstance(value, MutationEffect):
            key = id(value)
            if key not in indices:
                indices[key] = len(nodes)
                nodes.append(None)  # reserve before following self links
                fields = {}
                for name, parameter in signature(type(value).__init__).parameters.items():
                    if name == "self" or parameter.kind in (
                            Parameter.VAR_POSITIONAL, Parameter.VAR_KEYWORD):
                        continue
                    if name == "primary_effects" and isinstance(value, StructuralVariantEffect):
                        fields[name] = value._primary_effects
                    elif name == "candidates" and hasattr(value, "_candidates"):
                        fields[name] = value._candidates
                    elif hasattr(value, name):
                        fields[name] = getattr(value, name)
                    elif parameter.default is Parameter.empty:
                        raise ValueError("Missing effect constructor field: %s" % name)
                attached = {name: value.__dict__[name] for name in _ATTACHED_FIELDS
                            if name in value.__dict__ and name not in fields}
                nodes[indices[key]] = (type(value), encode(fields), encode(attached))
            return _EffectReference(indices[key])
        if isinstance(value, EffectCandidate):
            return replace(value, effect=encode(value.effect), evidence=encode(value.evidence))
        if isinstance(value, dict):
            return {key: encode(item) for key, item in value.items()}
        if isinstance(value, (tuple, list)):
            return type(value)(encode(item) for item in value)
        return value

    encode(root)
    return {"structural_effect_schema": 1, "nodes": tuple(nodes)}


def structural_effect_from_dict(root_class, state):
    from .effect_classes import MutationEffect

    if state.get("structural_effect_schema") != 1:
        raise ValueError("Unsupported structural effect serialization schema")
    nodes = state["nodes"]
    if not nodes or nodes[0][0] is not root_class:
        raise ValueError("Structural effect root class does not match archive")
    for cls, _, _ in nodes:
        if not isinstance(cls, type) or not issubclass(cls, MutationEffect):
            raise ValueError("Structural effect graph contains a non-effect class")
    objects = [object.__new__(cls) for cls, _, _ in nodes]
    started = set()

    def decode(value):
        if isinstance(value, _EffectReference):
            index = value.index
            if not isinstance(index, int) or not 0 <= index < len(nodes):
                raise ValueError("Invalid structural effect graph reference")
            if index not in started:
                started.add(index)
                cls, fields, attached = nodes[index]
                cls.__init__(objects[index], **decode(fields))
                for name, item in decode(attached).items():
                    if name not in _ATTACHED_FIELDS:
                        raise ValueError("Unknown structural effect attachment: %s" % name)
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
