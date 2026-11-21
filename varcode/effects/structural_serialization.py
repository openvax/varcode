"""Versioned effect graphs, finite class registry, and safe typed JSON codec.

Effects and transcript-model records use graph references. Ordinary effects
retain their existing reference formats. The ORF schema additionally snapshots
reference RNAs and coding models, independently of an annotation database.

:func:`effect_serialization_classes` exposes the finite reviewed class allowlist
for native Varcode effect and :class:`~varcode.EffectCollection` archives, and
:func:`effect_from_json` / :func:`effect_collection_from_json` decode those
archives without importing arbitrary module or function names from input data.
"""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Mapping
from dataclasses import dataclass, fields, replace
from enum import Enum
from functools import lru_cache
from inspect import Parameter, signature
import json
from types import MappingProxyType
from typing import TYPE_CHECKING, Any

from serializable import DataclassSerializable
from serializable.helpers import (
    SERIALIZED_DICTIONARY_KEYS_FIELD,
    parse_serialized_keys_index,
    to_json as serializable_to_json,
)

from ..effect_candidates import EffectCandidate

if TYPE_CHECKING:
    from .effect_classes import MutationEffect
    from .effect_collection import EffectCollection


@dataclass(frozen=True)
class EffectReference(DataclassSerializable):
    """Zero-based node reference within a serialized native effect graph."""

    index: int


# Backwards-compatible alias for legacy serialized archives.
_EffectReference = EffectReference


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

_PRIMITIVE_TYPES = (type(None), bool, int, float, str)
_COLLECTIONS_BUILT_FROM_VALUE = (set, frozenset, tuple)


def _model_types() -> tuple[type, ...]:
    from ..effect_hypotheses import (
        ClassifiedOutcome,
        EffectHypothesis,
        RealizedEffectCandidate,
        RealizedTranscriptProduct,
        SpliceOption,
        SplicePlan,
        SpliceSiteKey,
    )
    from ..mutant_transcript import (
        AssembledAllele,
        MutantTranscript,
        ReferenceSegment,
        TranscriptEdit,
    )
    from ..orf_discovery import DiscoveredORF, ORFDiscoveryPolicy
    from ..transcript_layout import CrypticSplicePolicy

    return (
        CrypticSplicePolicy,
        DiscoveredORF,
        ORFDiscoveryPolicy,
        MutantTranscript,
        ReferenceSegment,
        TranscriptEdit,
        AssembledAllele,
        ClassifiedOutcome,
        EffectHypothesis,
        RealizedEffectCandidate,
        RealizedTranscriptProduct,
        SpliceOption,
        SplicePlan,
        SpliceSiteKey,
    )


def _all_subclasses(base: type) -> set[type]:
    found: set[type] = {base}
    pending: list[type] = list(base.__subclasses__())
    while pending:
        cls = pending.pop()
        if cls not in found:
            found.add(cls)
            pending.extend(cls.__subclasses__())
    return found


@lru_cache(maxsize=1)
def _registry_bundle() -> tuple[
        Mapping[tuple[str, str], type],
        frozenset[type],
        Mapping[tuple[str, str], Any],
]:
    from pyensembl.ensembl_release import EnsemblRelease
    from pyensembl.exon import Exon
    from pyensembl.gene import Gene
    from pyensembl.genetic_code import GeneticCodeAssignment
    from pyensembl.genome import Genome as PyEnsemblGenome
    from pyensembl.protein import Protein
    from pyensembl.species import Species
    from pyensembl.transcript import Transcript

    from ..effect_hypotheses import SpliceAxis
    from ..genome import Genome
    from ..genotype import Genotype, Zygosity
    from ..germline import (
        Completeness,
        GermlineContext,
        PhaseHypothesis,
        PhasePartition,
    )
    from ..mutant_transcript import AssembledAllele
    from ..orf import (
        ORFAnnotation,
        ORFOutcome,
        ORFProduct,
        ORFProductAnnotation,
        ORFTranscriptView,
    )
    from ..orf_discovery import ORFDiscovery, ORFDiscoveryAssessment
    from ..orf_serialization import (
        ReferenceGenomeSnapshot,
        ReferenceTranscriptSnapshot,
        TranscriptFeatureSnapshot,
        orf_graph_types,
    )
    from ..rna_evidence import RNAEvidence
    from ..splice_outcomes import SpliceCandidateRNAEvidence, SpliceOutcomeSet  # noqa: F401
    from ..structural_variant import Breakend, StructuralVariant
    from ..uorf import ORFTranscript, UORFAnnotation
    from ..variant import Variant
    from ..variant_collection import VariantCollection
    from .effect_classes import MutationEffect
    from .effect_collection import EffectCollection

    reviewed: set[type] = {
        tuple,
        set,
        frozenset,
        OrderedDict,
        EffectCollection,
        Variant,
        VariantCollection,
        StructuralVariant,
        Breakend,
        Genotype,
        Zygosity,
        Genome,
        PyEnsemblGenome,
        EnsemblRelease,
        Transcript,
        Gene,
        Exon,
        Protein,
        Species,
        GeneticCodeAssignment,
        EffectReference,
        EffectCandidate,
        SpliceAxis,
        ORFDiscovery,
        ORFDiscoveryAssessment,
        Completeness,
        GermlineContext,
        PhaseHypothesis,
        PhasePartition,
        RNAEvidence,
        SpliceCandidateRNAEvidence,
        ORFAnnotation,
        ORFOutcome,
        ORFProduct,
        ORFProductAnnotation,
        ORFTranscriptView,
        UORFAnnotation,
        ORFTranscript,
        ReferenceGenomeSnapshot,
        TranscriptFeatureSnapshot,
        ReferenceTranscriptSnapshot,
        *_model_types(),
        *_all_subclasses(MutationEffect),
        *orf_graph_types(),
    }

    try:
        from importlib import import_module

        import orfdata
        from orfdata import ReferenceORF
        from orfdata.catalog import CatalogStatus
        from orfdata.rna import ProteinProduct, RNAFeature
        reviewed.update({ReferenceORF, CatalogStatus, ProteinProduct, RNAFeature})
        # Use public export once openvax/orfdata#210 ships; fall back to
        # _pinned_genome so T2T-CHM13v2.0 and other pinned genomes round-trip.
        pinned_cls = getattr(orfdata, "PinnedEnsemblGenome", None)
        if not isinstance(pinned_cls, type):
            pinned_cls = getattr(
                import_module("orfdata._pinned_genome"),
                "PinnedEnsemblGenome",
                None,
            )
        if isinstance(pinned_cls, type):
            reviewed.add(pinned_cls)
    except ImportError:
        pass

    by_key: dict[tuple[str, str], type] = {
        (cls.__module__, cls.__qualname__): cls for cls in reviewed
    }
    # Explicit legacy aliases for archives written under earlier 11.x names.
    by_key.update({
        ("varcode.effects.structural_serialization", "_EffectReference"): EffectReference,
        ("varcode.mutant_transcript", "_AssembledAllele"): AssembledAllele,
        ("varcode.effects.structural", "_AssembledAllele"): AssembledAllele,
        ("varcode.orf", "_ORFTranscript"): ORFTranscriptView,
        ("varcode.orf_serialization", "_ReferenceGenome"): ReferenceGenomeSnapshot,
        ("varcode.orf_serialization", "_Feature"): TranscriptFeatureSnapshot,
        ("varcode.orf_serialization", "_ReferenceTranscript"): ReferenceTranscriptSnapshot,
    })
    from ..variant import variant_ascending_position_sort_key
    from .effect_collection import default_effect_sort_key
    from .effect_ordering import multi_gene_effect_sort_key

    allowed_callables: dict[tuple[str, str], Any] = {
        ("varcode.effects.effect_collection", "default_effect_sort_key"): default_effect_sort_key,
        ("varcode.effects.effect_collection", "_default_effect_sort_key"): default_effect_sort_key,
        ("varcode.effects.effect_ordering", "multi_gene_effect_sort_key"): multi_gene_effect_sort_key,
        ("varcode.variant", "variant_ascending_position_sort_key"): variant_ascending_position_sort_key,
    }
    return MappingProxyType(by_key), frozenset(by_key.values()), MappingProxyType(allowed_callables)


def effect_serialization_classes() -> Mapping[tuple[str, str], type]:
    """Return the read-only finite allowlist of classes for native effect archives.

    Returns
    -------
    Mapping[tuple[str, str], type]
        Read-only mapping from ``(module_name, qualname)`` to reviewed Python
        classes accepted by :func:`effect_from_json`,
        :func:`effect_from_serializable_repr`,
        :func:`effect_collection_from_json`, and
        :func:`effect_collection_from_serializable_repr` (including legacy
        serialized class aliases such as ``_EffectReference`` and
        ``_AssembledAllele``). No functions or dynamic module imports are
        permitted.
    """
    return _registry_bundle()[0]


def structural_effect_to_dict(
        root: Any, *, schema: str = "structural_effect_schema",
) -> dict[str, Any]:
    """Archive a native graph, or project one inert ``orf_comparison_schema`` outcome."""
    from ..structural_variant import StructuralVariant
    from ..variant import Variant
    from .effect_classes import MutationEffect, PhaseCandidateSet, StructuralVariantEffect

    indices: dict[int, int] = {}
    nodes: list[Any] = []
    model_types = _model_types()
    comparison = schema == "orf_comparison_schema"
    orf = schema == "orf_effect_schema" or comparison
    special_types: tuple[type, ...] = (Variant, StructuralVariant)
    if orf:
        from orfdata import ReferenceORF
        from pyensembl import Transcript

        from ..orf import ORFTranscriptView
        from ..orf_serialization import (
            ReferenceGenomeSnapshot,
            ReferenceTranscriptSnapshot,
            orf_graph_types,
            record_fields,
        )
        model_types += orf_graph_types()
        special_types = (ReferenceORF, ORFTranscriptView, Variant, StructuralVariant)

    def encode(value: Any) -> Any:
        if orf and isinstance(value, Transcript):
            return ReferenceTranscriptSnapshot.capture(value)
        if isinstance(value, MutationEffect) or type(value) in model_types + special_types:
            key = id(value)
            if key not in indices:
                indices[key] = len(nodes)
                nodes.append(None)  # reserve before following self links
                if orf and isinstance(value, ReferenceORF):
                    state = (
                        dict(
                            snapshot=value.to_dict(),
                            transcript=ReferenceTranscriptSnapshot.capture(value.transcript),
                        )
                        if value.record is None
                        else dict(
                            record=value.record,
                            transcript=ReferenceTranscriptSnapshot.capture(value.transcript),
                            context=value.context,
                            any_transcript=True,
                        )
                    )
                    nodes[indices[key]] = (type(value), encode(state), {})
                    return EffectReference(indices[key])
                if orf and isinstance(value, ORFTranscriptView):
                    nodes[indices[key]] = (type(value), encode(dict(reference=value.reference_orf)), {})
                    return EffectReference(indices[key])
                if isinstance(value, Variant):
                    state = value.to_dict()
                    if orf:
                        state = {**state, "genome": ReferenceGenomeSnapshot.capture(value.genome)}
                    nodes[indices[key]] = (type(value), encode(state), {})
                    return EffectReference(indices[key])
                if type(value) in model_types:
                    state = record_fields(value) if orf else {
                        f.name: getattr(value, f.name) for f in fields(value)
                    }
                    nodes[indices[key]] = (type(value), encode(state), {})
                    return EffectReference(indices[key])
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
                        state[name] = (value._candidates_raw if name == "candidates" else value._hypotheses)
                    elif name == "candidates" and hasattr(value, "_candidates"):
                        state[name] = value._candidates
                    elif hasattr(value, name):
                        state[name] = getattr(value, name)
                    elif parameter.default is Parameter.empty:
                        raise ValueError("Missing effect constructor field: %s" % name)
                attached = {
                    name: value.__dict__[name]
                    for name in _ATTACHED_FIELDS
                    if name in value.__dict__
                    and name not in state
                    and not (comparison and name in _ALTERNATIVE_FIELDS)
                }
                nodes[indices[key]] = (type(value), encode(state), encode(attached))
            return EffectReference(indices[key])
        if isinstance(value, EffectCandidate):
            return replace(value, effect=encode(value.effect), evidence=encode(value.evidence))
        if isinstance(value, dict):
            return {key: encode(item) for key, item in value.items()}
        if isinstance(value, (tuple, list)):
            return type(value)(encode(item) for item in value)
        return value

    encode(root)
    return {schema: 1, "nodes": tuple(nodes)}


def structural_effect_from_dict(
        root_class: type, state: Mapping[str, Any], *, schema: str = "structural_effect_schema",
) -> Any:
    """Reconstruct a validated native effect graph from ``state``."""
    from ..structural_variant import StructuralVariant
    from ..variant import Variant
    from .effect_classes import MutationEffect

    if schema == "orf_comparison_schema":
        raise ValueError("Comparison projections are inert, not reconstructable effect graphs")
    if state.get(schema) != 1:
        raise ValueError("Unsupported effect serialization schema")
    nodes = state["nodes"]
    if not nodes or nodes[0][0] is not root_class:
        raise ValueError("Effect root class does not match archive")
    model_types = _model_types()
    special_types: tuple[type, ...] = (Variant, StructuralVariant)
    if schema == "orf_effect_schema":
        from orfdata import ReferenceORF

        from ..orf import ORFTranscriptView
        from ..orf_serialization import orf_graph_types
        model_types += orf_graph_types()
        special_types = (ReferenceORF, ORFTranscriptView, Variant, StructuralVariant)
    for cls, _, attached in nodes:
        if not isinstance(cls, type) or not (
                issubclass(cls, MutationEffect) or cls in model_types + special_types):
            raise ValueError("Effect graph contains an unsupported class")
        if cls in model_types and attached:
            raise ValueError("Unexpected transcript-model attachment")
    objects = [object.__new__(cls) for cls, _, _ in nodes]
    started: set[int] = set()

    def decode(value: Any) -> Any:
        if isinstance(value, EffectReference):
            index = value.index
            if type(index) is not int or not 0 <= index < len(nodes):
                raise ValueError("Invalid effect graph reference")
            if index not in started:
                started.add(index)
                cls, node_fields, attached = nodes[index]
                values = decode(node_fields)
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

    return decode(EffectReference(0))


def _safe_from_reconstructed_dict(
        x: dict[Any, Any],
        allowlist: Mapping[tuple[str, str], type],
        allowed_classes: frozenset[type],
        allowed_callables: Mapping[tuple[str, str], Any],
) -> Any:
    if "__function__" in x:
        raise ValueError("Function metadata is not permitted in effect archives")
    name = x.get("__name__")
    module_name = x.get("__module__")
    if len(x) == 2 and isinstance(name, str) and isinstance(module_name, str):
        key = (module_name, name)
        if key in allowlist:
            return allowlist[key]
        if key in allowed_callables:
            return allowed_callables[key]
        raise ValueError(
            "Unsupported class or function metadata in effect archive: %s.%s"
            % key)

    serialized_keys = x.get(SERIALIZED_DICTIONARY_KEYS_FIELD)
    if (
        isinstance(serialized_keys, list)
        and serialized_keys
        and all(isinstance(serialized_key, str) for serialized_key in serialized_keys)
    ):
        non_string_key_objects = [
            _safe_from_json_any(
                serialized_key, allowlist, allowed_classes, allowed_callables
            )
            for serialized_key in serialized_keys
        ]
        converted_dict: dict[Any, Any] = type(x)()
        for k, v in x.items():
            if k == SERIALIZED_DICTIONARY_KEYS_FIELD:
                continue
            serialized_key_index = parse_serialized_keys_index(k)
            if serialized_key_index is not None and serialized_key_index < len(non_string_key_objects):
                k = non_string_key_objects[serialized_key_index]
            converted_dict[k] = v
        x = converted_dict

    class_object = x.get("__class__")
    if not isinstance(class_object, type):
        if callable(class_object):
            raise ValueError(f"Expected __class__ to be a reviewed class, got {class_object!r}")
        return x
    if class_object not in allowed_classes:
        raise ValueError(
            "Unsupported class in effect archive: %s.%s"
            % (class_object.__module__, class_object.__qualname__))
    state = dict(x)
    del state["__class__"]
    if issubclass(class_object, Enum):
        if "__member__" in state:
            member = state["__member__"]
            if not isinstance(member, str):
                raise ValueError(f"Invalid enum member {member!r} for {class_object.__qualname__}")
            try:
                return class_object[member]
            except KeyError:
                raise ValueError(
                    f"{class_object.__qualname__} has no member named {member!r}"
                ) from None
        if "__value__" in state:
            return class_object(state["__value__"])
    elif "__value__" in state and class_object in _COLLECTIONS_BUILT_FROM_VALUE:
        return class_object(state["__value__"])
    if hasattr(class_object, "from_dict"):
        return class_object.from_dict(state)
    return class_object(**state)


def _safe_from_serializable_repr_any(
        x: Any,
        allowlist: Mapping[tuple[str, str], type],
        allowed_classes: frozenset[type],
        allowed_callables: Mapping[tuple[str, str], Any],
) -> Any:
    if isinstance(x, _PRIMITIVE_TYPES):
        return x
    if isinstance(x, list):
        return [
            _safe_from_serializable_repr_any(
                element, allowlist, allowed_classes, allowed_callables
            )
            for element in x
        ]
    if isinstance(x, dict):
        converted: dict[Any, Any] = type(x)()
        for k, v in x.items():
            if not isinstance(k, str):
                raise ValueError("Serialized dictionary keys must be strings")
            converted[k] = _safe_from_serializable_repr_any(
                v, allowlist, allowed_classes, allowed_callables
            )
        return _safe_from_reconstructed_dict(
            converted, allowlist, allowed_classes, allowed_callables
        )
    raise TypeError(
        f"Cannot convert {x!r} ({type(x).__name__}) from serializable representation")


def _safe_from_json_any(
        json_string: str | bytes,
        allowlist: Mapping[tuple[str, str], type],
        allowed_classes: frozenset[type],
        allowed_callables: Mapping[tuple[str, str], Any],
) -> Any:
    if not isinstance(json_string, (str, bytes)):
        raise TypeError(
            f"Expected JSON str or bytes, got {type(json_string).__name__}")
    return json.loads(
        json_string,
        object_hook=lambda d: _safe_from_reconstructed_dict(
            d, allowlist, allowed_classes, allowed_callables
        ),
    )


def effect_to_json(effect: MutationEffect) -> str:
    """Serialize one :class:`~varcode.MutationEffect` to native JSON.

    Parameters
    ----------
    effect : MutationEffect
        Ordinary, phased, germline, structural, or ORF mutation effect.

    Returns
    -------
    str
        JSON archive string preserving graph cycles, candidate ordering,
        baselines, and evidence.

    Raises
    ------
    TypeError
        If ``effect`` is not a :class:`~varcode.MutationEffect`.
    """
    from .effect_classes import MutationEffect

    if not isinstance(effect, MutationEffect):
        raise TypeError(f"Expected a MutationEffect, got {type(effect).__name__}")
    return serializable_to_json(effect)


def effect_from_serializable_repr(data: object) -> MutationEffect:
    """Decode one :class:`~varcode.MutationEffect` from parsed JSON data using the finite class allowlist.

    Parameters
    ----------
    data : object
        JSON-decoded dictionary produced by ``serializable.to_serializable_repr``
        or ``json.loads(effect_to_json(effect))``.

    Returns
    -------
    MutationEffect
        Reconstructed effect instance.

    Raises
    ------
    TypeError
        If ``data`` contains non-JSON container types.
    ValueError
        If ``data`` names an unreviewed class or function, violates graph or
        schema invariants, or does not decode to a :class:`~varcode.MutationEffect`.
    """
    from .effect_classes import MutationEffect

    allowlist, allowed_classes, allowed_callables = _registry_bundle()
    result = _safe_from_serializable_repr_any(
        data, allowlist, allowed_classes, allowed_callables
    )
    if not isinstance(result, MutationEffect):
        raise ValueError(
            f"Expected archive root to decode to MutationEffect, got {type(result).__name__}")
    return result


def effect_from_json(json_string: str | bytes) -> MutationEffect:
    """Decode one :class:`~varcode.MutationEffect` from JSON using the finite class allowlist.

    Parameters
    ----------
    json_string : str or bytes
        JSON archive produced by :func:`effect_to_json` or ``effect.to_json()``.

    Returns
    -------
    MutationEffect
        Reconstructed effect instance.

    Raises
    ------
    TypeError
        If ``json_string`` is not ``str`` or ``bytes``.
    ValueError
        If the archive is malformed, names an unreviewed class or function, or
        does not decode to a :class:`~varcode.MutationEffect`.
    """
    from .effect_classes import MutationEffect

    allowlist, allowed_classes, allowed_callables = _registry_bundle()
    result = _safe_from_json_any(
        json_string, allowlist, allowed_classes, allowed_callables
    )
    if not isinstance(result, MutationEffect):
        raise ValueError(
            f"Expected archive root to decode to MutationEffect, got {type(result).__name__}")
    return result


def effect_collection_to_json(effects: EffectCollection) -> str:
    """Serialize an :class:`~varcode.EffectCollection` to native JSON.

    Parameters
    ----------
    effects : EffectCollection
        Effect collection to archive.

    Returns
    -------
    str
        JSON archive string preserving collection provenance, ORF annotations,
        and per-effect candidate graphs.

    Raises
    ------
    TypeError
        If ``effects`` is not an :class:`~varcode.EffectCollection`.
    """
    from .effect_collection import EffectCollection

    if not isinstance(effects, EffectCollection):
        raise TypeError(f"Expected an EffectCollection, got {type(effects).__name__}")
    return serializable_to_json(effects)


def effect_collection_from_serializable_repr(data: object) -> EffectCollection:
    """Decode an :class:`~varcode.EffectCollection` from parsed JSON data using the finite class allowlist.

    Parameters
    ----------
    data : object
        JSON-decoded dictionary produced by ``json.loads(effect_collection_to_json(effects))``.

    Returns
    -------
    EffectCollection
        Reconstructed effect collection.

    Raises
    ------
    TypeError
        If ``data`` contains non-JSON container types.
    ValueError
        If ``data`` names an unreviewed class or function, violates graph or
        schema invariants, or does not decode to an :class:`~varcode.EffectCollection`.
    """
    from .effect_collection import EffectCollection

    allowlist, allowed_classes, allowed_callables = _registry_bundle()
    result = _safe_from_serializable_repr_any(
        data, allowlist, allowed_classes, allowed_callables
    )
    if not isinstance(result, EffectCollection):
        raise ValueError(
            f"Expected archive root to decode to EffectCollection, got {type(result).__name__}")
    return result


def effect_collection_from_json(json_string: str | bytes) -> EffectCollection:
    """Decode an :class:`~varcode.EffectCollection` from JSON using the finite class allowlist.

    Parameters
    ----------
    json_string : str or bytes
        JSON archive produced by :func:`effect_collection_to_json` or
        ``effects.to_json()``.

    Returns
    -------
    EffectCollection
        Reconstructed effect collection.

    Raises
    ------
    TypeError
        If ``json_string`` is not ``str`` or ``bytes``.
    ValueError
        If the archive is malformed, names an unreviewed class or function, or
        does not decode to an :class:`~varcode.EffectCollection`.
    """
    from .effect_collection import EffectCollection

    allowlist, allowed_classes, allowed_callables = _registry_bundle()
    result = _safe_from_json_any(
        json_string, allowlist, allowed_classes, allowed_callables
    )
    if not isinstance(result, EffectCollection):
        raise ValueError(
            f"Expected archive root to decode to EffectCollection, got {type(result).__name__}")
    return result


__all__ = [
    "EffectReference",
    "effect_collection_from_json",
    "effect_collection_from_serializable_repr",
    "effect_collection_to_json",
    "effect_from_json",
    "effect_from_serializable_repr",
    "effect_serialization_classes",
    "effect_to_json",
    "structural_effect_from_dict",
    "structural_effect_to_dict",
]
