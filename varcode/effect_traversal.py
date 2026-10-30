"""Traverse native effect alternatives without discarding their provenance."""

from __future__ import annotations

from collections.abc import Iterator
from typing import TypeAlias

from .effect_candidates import EffectCandidate
from .effect_hypotheses import RealizedEffectCandidate
from .effects.effect_classes import MutationEffect


Candidate: TypeAlias = EffectCandidate | RealizedEffectCandidate
CandidatePath: TypeAlias = tuple[Candidate, ...]


def candidate_children(effect: MutationEffect) -> tuple[Candidate, ...]:
    """Return direct alternatives in producer order, including attached evidence.

    Parameters
    ----------
    effect : MutationEffect
        An ordinary classification, shared transcript comparison, or legacy
        multi-outcome wrapper. ORF effects use the same interface.

    Returns
    -------
    tuple of EffectCandidate or RealizedEffectCandidate
        Existing candidate objects, without copying or merging them. Shared
        transcript candidates retain their complete ``outcomes`` (including
        each patient baseline, mutant product and hypothesis). A plain effect
        returns an empty tuple. A candidate may point back to ``effect``.

    Notes
    -----
    ``MultiOutcomeEffect`` is not a prerequisite for having alternatives.
    Enumeration can materialize lazily computed structural/splice candidates;
    prediction errors propagate rather than silently dropping a candidate.
    """
    realized = getattr(effect, "realized_candidates", None)
    if realized is not None:
        return tuple(realized) + tuple(getattr(effect, "_extra_candidates", ()))
    return tuple(getattr(effect, "candidates", ()))


def iter_effects(effect: MutationEffect) -> Iterator[MutationEffect]:
    """Yield the selected effect and all reachable effects once by identity.

    Parameters
    ----------
    effect : MutationEffect
        Root result to inspect, including ordinary candidate-bearing effects.

    Yields
    ------
    MutationEffect
        Root first, followed by depth-first producer order. Wrappers remain
        present. Shared objects and cycles are visited once; equal-but-distinct
        objects are not merged. Use this iterator for filtering, not for
        counting hypotheses or reporting provenance.

    Notes
    -----
    Each effect's own sequence/change fields describe its selected comparison.
    No protein, baseline, or probability is inferred by traversal. To retain
    distinct routes to an effect, use :func:`iter_effect_paths` instead.
    """
    pending = [effect]
    seen: set[int] = set()
    while pending:
        current = pending.pop()
        if id(current) in seen:
            continue
        seen.add(id(current))
        yield current
        pending.extend(c.effect for c in reversed(candidate_children(current)))


def iter_effect_paths(effect: MutationEffect) -> Iterator[tuple[MutationEffect, CandidatePath]]:
    """Yield every effect occurrence with its original candidate path.

    Parameters
    ----------
    effect : MutationEffect
        Root result to inspect.

    Yields
    ------
    tuple of MutationEffect and tuple of candidates
        Root with an empty path, then depth-first occurrences in producer order.
        Shared effects under different evidence remain separate occurrences.
        A self/cycle-closing candidate is yielded with its evidence, but is not
        expanded again. This is a path traversal, not a unique-effect iterator
        or an enumeration of translated proteins.

    Notes
    -----
    Read ``source``/``evidence`` from EffectCandidate path entries; read
    ``outcomes`` from RealizedEffectCandidate entries for complete baseline,
    mutant and hypothesis records. Neither producer order nor ordinal rank
    implies a calibrated probability. No global deduplication erases paths.
    """
    pending: list[tuple[MutationEffect, CandidatePath, frozenset[int]]] = [
        (effect, (), frozenset())]
    while pending:
        current, path, ancestors = pending.pop()
        yield current, path
        if id(current) in ancestors:
            continue
        ancestors = ancestors | {id(current)}
        pending.extend((candidate.effect, (*path, candidate), ancestors)
                       for candidate in reversed(candidate_children(current)))
