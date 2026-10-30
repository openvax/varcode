"""Traverse classified comparisons without discarding candidate provenance."""

from __future__ import annotations

from collections.abc import Iterator
from typing import TYPE_CHECKING, TypeAlias

from .effect_candidates import EffectCandidate
from .effect_hypotheses import RealizedEffectCandidate

if TYPE_CHECKING:
    from .effects.effect_classes import MutationEffect

Candidate: TypeAlias = EffectCandidate | RealizedEffectCandidate
CandidatePath: TypeAlias = tuple[Candidate, ...]


def get_effect_candidates(effect: MutationEffect | None) -> tuple[Candidate, ...]:
    """Return immediate candidates in producer order, retaining original wrappers.

    Realized candidates are authoritative on transcript-model results; their
    externally attached candidates follow them. Legacy outcome sets include
    those extras in their own ``candidates`` property. A bare effect has none.
    This does not change probabilities, merge outcomes or infer missing products.
    """
    realized = getattr(effect, "realized_candidates", None)
    if realized is not None:
        return tuple(realized) + tuple(getattr(effect, "_extra_candidates", ()))
    return tuple(getattr(effect, "candidates", ()) or ())


def iter_candidate_paths(
        effect: MutationEffect | None,
) -> Iterator[tuple[MutationEffect | None, CandidatePath]]:
    """Yield terminal comparisons and their original enclosing candidate wrappers.

    Parameters
    ----------
    effect : MutationEffect or None
        An ordinary classified effect, a legacy outcome set, or no prediction.

    Yields
    ------
    (MutationEffect or None, tuple of candidates)
        A terminal effect and its path in producer order. A bare effect has an
        empty path. An ordinary selected comparison is also emitted when its
        candidates contain only alternatives. A primary self-candidate terminates at its own effect; a
        longer cycle or empty legacy outcome set yields ``None``. That unknown
        terminal is not an unchanged prediction or an invented comparison.

    Notes
    -----
    Distinct paths are never deduplicated, even when they share an effect.
    Realized wrappers retain every outcome's baseline, mutant and hypothesis.
    Traversal does not split or merge those outcomes. Use
    :func:`iter_candidate_effects` only when provenance is unnecessary.
    """
    from .effects.effect_classes import MultiOutcomeEffect

    pending: list[tuple[MutationEffect | None, CandidatePath, frozenset[int], bool]] = [
        (effect, (), frozenset(), False)]
    while pending:
        current, path, ancestors, primary = pending.pop()
        if current is None or primary:
            yield current, path
            continue
        if id(current) in ancestors:
            yield None, path
            continue
        children = get_effect_candidates(current)
        if not children:
            yield (None if isinstance(current, MultiOutcomeEffect) else current), path
            continue
        if not isinstance(current, MultiOutcomeEffect) and not any(
                child.effect is current for child in children):
            yield current, path
        ancestors = ancestors | {id(current)}
        pending.extend(
            (child.effect, path + (child,), ancestors, child.effect is current)
            for child in reversed(children))


def iter_candidate_effects(effect: MutationEffect | None) -> Iterator[MutationEffect | None]:
    """Yield terminal effects once by identity, including at most one unknown.

    This convenience view is suitable for sequence availability and filtering.
    It deliberately drops path provenance; use :func:`iter_candidate_paths`
    when different patient baselines, phases or evidence must remain distinct.
    """
    seen: set[int] = set()
    for terminal, _path in iter_candidate_paths(effect):
        if id(terminal) not in seen:
            seen.add(id(terminal))
            yield terminal
