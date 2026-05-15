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

"""RNA-evidence resolver protocol for refining DNA-predicted variant
effects with observed isoform / read-level evidence (#259).

Varcode's effect classifiers predict possibilities from genomic
coordinates alone — for many variants, especially structural ones,
the actual protein consequence depends on what's transcribed and
spliced. This module defines the contract that lets RNA-aware tools
(Isovar, Exacto, custom long-read pipelines) attach observed
:class:`~varcode.effect_candidates.EffectCandidate` objects to existing effects
without subclassing or churning the core classes.

The shape mirrors :mod:`varcode.phasing`: a runtime-checkable
:class:`Protocol`, an :func:`apply_rna_evidence_to_effects` walk that
refines multi-candidate effects, and a small :func:`make_rna_outcome`
factory for the common provenance fields. Like :mod:`varcode.phasing`,
nothing here imports Isovar — implementers duck-type the protocol.

Usage::

    resolver = MyIsovarResolver(reads=...)
    # Per-call:
    effects = vc.effects(rna_resolver=resolver)
    # Or apply post-hoc to an existing collection:
    apply_rna_evidence_to_effects(effects, resolver)

    # Each affected effect gains observed candidates alongside its
    # DNA-predicted ones. Filter by source:
    for effect in effects:
        for candidate in getattr(effect, "candidates", ()):
            if candidate.source != "varcode":
                # An RNA-observed isoform.
                ...
"""
from __future__ import annotations

from typing import Any, Iterable, Mapping, Optional, Protocol, Sequence, runtime_checkable

from .effect_candidates import EffectCandidate


@runtime_checkable
class RNAEvidenceResolver(Protocol):
    """Source of RNA-observed outcomes for a ``(variant, transcript)``
    pair.

    Implementers return zero or more :class:`~varcode.effect_candidates.EffectCandidate`
    objects describing isoforms, fusions, or RNA-level events that were
    actually observed in reads. An empty sequence means "no evidence
    for this pair" — the existing DNA-predicted outcomes are left
    alone.

    Returned outcomes should set ``source`` to a producer-specific
    string (the name of the RNA assembler, long-read caller, fusion
    detector, etc.) and populate ``evidence`` with whatever shape that
    producer natively emits (transcript model IDs, junction read
    counts, etc.). See :func:`make_rna_outcome` for a convenience
    factory that fills the common fields.
    """

    def observed_outcomes(self, variant, transcript) -> Sequence[EffectCandidate]:
        """Return RNA-observed outcomes for ``variant`` on
        ``transcript``, or an empty sequence when no evidence is
        available. Must not raise on unknown ``(variant, transcript)``
        pairs — return an empty sequence instead."""
        ...


class NullRNAEvidenceResolver:
    """No-op resolver that always reports "no evidence".

    Useful as a default in pipelines where an RNA resolver is optional
    and as a baseline in tests. ``apply_rna_evidence_to_effects`` is
    safe to call with this resolver — it's a no-op walk."""

    def observed_outcomes(self, variant, transcript) -> Sequence[EffectCandidate]:
        return ()


def make_rna_outcome(
        effect,
        *,
        source: str = "rna",
        transcript_model_id: Optional[str] = None,
        read_count: Optional[int] = None,
        extra_evidence: Optional[Mapping[str, Any]] = None) -> EffectCandidate:
    """Construct an :class:`~varcode.effect_candidates.EffectCandidate`
    carrying RNA-derived provenance.

    Convenience factory for the common fields a reads-based or
    long-read assembly tool wants on each observed outcome — keeps
    consumers from hand-rolling the ``evidence`` dict shape and lets
    downstream code rely on a small set of well-known keys.

    Parameters
    ----------
    effect : MutationEffect
        The effect this RNA-observed outcome represents.
    source : str
        Producer name; defaults to ``"rna"``. Set to a tool-specific
        string (RNA assembler, long-read caller, etc.) for downstream
        filtering. Opaque to varcode.
    transcript_model_id : str or None
        Stable ID of the observed transcript model from the producer.
        Stored under ``evidence["transcript_model_id"]``.
    read_count : int or None
        Supporting read count. Stored under ``evidence["read_count"]``.
    extra_evidence : Mapping or None
        Producer-specific extra fields, merged into the evidence dict
        on top of the well-known keys above. Allows tool-native fields
        (e.g. ``"tpm"``, ``"junction_id"``) to ride along without
        forcing a schema here.
    """
    evidence: dict = {}
    if transcript_model_id is not None:
        evidence["transcript_model_id"] = transcript_model_id
    if read_count is not None:
        evidence["read_count"] = read_count
    if extra_evidence:
        evidence.update(extra_evidence)
    return EffectCandidate(
        effect=effect,
        source=source,
        evidence=evidence,
    )


def _replace_effects(effects, replacements):
    if hasattr(effects, "elements") and hasattr(effects, "effects"):
        effects.elements = replacements
        effects.effects = replacements
        return effects
    try:
        effects[:] = replacements
        return effects
    except TypeError:
        return tuple(replacements)


def apply_rna_evidence_to_effects(effects: Iterable, resolver) -> Iterable:
    """Attach RNA-observed candidates from ``resolver`` to each effect.

    Walks ``effects`` and, for any effect with a resolvable
    ``(variant, transcript)``, asks ``resolver.observed_outcomes``
    for RNA-observed candidates.

    Splice mechanism sets use RNA evidence as a reconciliation signal:
    a new set replaces the old one, retaining an audit trail of raw RNA
    evidence, added candidates, excluded DNA-predicted candidates, and
    per-current-candidate RNA support. Other multi-outcome effects keep
    the additive side-channel behavior: observed candidates are stashed
    on ``_extra_candidates`` and exposed through ``.candidates``.

    Single-outcome effects (Missense, FrameShift, etc.) are left
    untouched even when the resolver has evidence — those classes
    don't expose a multi-candidate view, and replacing them with a
    multi-outcome wrapper would break downstream ``isinstance`` checks.
    Producers that need to surface RNA observations on point variants
    should report them as a separate :class:`MultiOutcomeEffect` rather
    than mutating an existing single-outcome one. (The point-variant
    diff is generally already correct from DNA, so this is rarely an
    issue in practice.)

    Safe to call on a mixed collection where only some variants have
    RNA evidence; no-op when ``resolver`` is None or doesn't implement
    the protocol.

    Returns ``effects`` for chaining convenience.
    """
    if resolver is None:
        return effects
    if not hasattr(resolver, "observed_outcomes"):
        return effects

    # Lazy import to avoid an import cycle (effect_classes imports
    # from varcode.effect_candidates, which sits below us).
    from .effects.effect_classes import MultiOutcomeEffect

    replacements = []
    changed = False
    for effect in effects:
        if not isinstance(effect, MultiOutcomeEffect):
            replacements.append(effect)
            continue
        variant = getattr(effect, "variant", None)
        transcript = getattr(effect, "transcript", None)
        if variant is None or transcript is None:
            replacements.append(effect)
            continue
        observed = resolver.observed_outcomes(variant, transcript)
        if not observed:
            replacements.append(effect)
            continue
        # Coerce to tuple eagerly so we don't keep a generator that
        # would silently exhaust on a second candidates() read.
        observed_tuple = tuple(observed)
        if not observed_tuple:
            replacements.append(effect)
            continue
        if hasattr(effect, "with_rna_evidence"):
            refined = effect.with_rna_evidence(observed_tuple)
            replacements.append(refined)
            changed = changed or refined is not effect
        else:
            existing = getattr(effect, "_extra_candidates", ())
            effect._extra_candidates = tuple(existing) + observed_tuple
            replacements.append(effect)
    if changed:
        return _replace_effects(effects, replacements)
    return effects
