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

"""Annotator registry with one shared default and explicit compatibility choices.

``transcript_model`` composes transcript consequences. Other annotators
can decline individual inputs by returning ``NotImplemented``; no static
capability declaration or implicit fallback is used.
"""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager
from typing import TYPE_CHECKING, TypeVar, cast

if TYPE_CHECKING:
    from . import EffectAnnotator

from .fast import FastEffectAnnotator
from .protein_diff import ProteinDiffEffectAnnotator
from ..transcript_model import TranscriptModelEffectAnnotator


_AnnotatorT = TypeVar("_AnnotatorT", bound="EffectAnnotator")
_REGISTRY: dict[str, EffectAnnotator] = {}
_DEFAULT_NAME = "transcript_model"


def register_annotator(annotator: _AnnotatorT) -> _AnnotatorT:
    """Add an annotator to the process-global registry, keyed by its
    ``.name``. Re-registering under the same name overrides the
    previous entry — this is deliberate so callers can swap
    implementations in tests.

    Parameters
    ----------
    annotator : EffectAnnotator
        Implementation with a nonempty name and an annotation method.

    Returns
    -------
    EffectAnnotator
        The same instance, preserving its concrete type.

    Raises
    ------
    ValueError
        The instance has no nonempty name.
    """
    name = getattr(annotator, "name", None)
    if not name:
        raise ValueError(
            "Annotator %r has no .name attribute; cannot register." % annotator)
    _REGISTRY[name] = annotator
    return annotator


def get_annotator(name: str) -> EffectAnnotator:
    """Look up a registered annotator by name. Raises ``KeyError``
    if no annotator is registered under that name.
    """
    return _REGISTRY[name]


def get_default_annotator() -> EffectAnnotator:
    """Return the annotator currently configured as the default
    (``"transcript_model"`` unless changed via :func:`set_default_annotator` or
    :func:`use_annotator`).
    """
    return _REGISTRY[_DEFAULT_NAME]


def set_default_annotator(name: str) -> None:
    """Swap the process-wide default annotator. ``name`` must refer
    to a registered annotator.
    """
    global _DEFAULT_NAME
    if name not in _REGISTRY:
        raise KeyError(
            "No annotator registered under %r — call register_annotator() "
            "first or pick from %r." % (name, sorted(_REGISTRY)))
    _DEFAULT_NAME = name


def resolve_annotator(
        annotator_or_name: str | EffectAnnotator | None,
) -> EffectAnnotator:
    """Normalize a per-call ``annotator=`` argument to an instance.

    Accepts ``None`` (use the current default), a string (look up in
    the registry), or an object that already implements the
    :class:`EffectAnnotator` protocol (used directly without
    validation — we trust duck-typing).

    Raises :class:`KeyError` for unknown names so callers that
    mistype a string get an early, legible error rather than an
    obscure attribute failure downstream.
    """
    if annotator_or_name is None:
        return get_default_annotator()
    if isinstance(annotator_or_name, str):
        if annotator_or_name not in _REGISTRY:
            raise KeyError(
                "No annotator registered under %r — known annotators: %r. "
                "Register one with varcode.register_annotator() or pick "
                "from the list." % (
                    annotator_or_name, sorted(_REGISTRY)))
        return _REGISTRY[annotator_or_name]
    return annotator_or_name


@contextmanager
def use_annotator(name_or_instance: str | EffectAnnotator) -> Iterator[None]:
    """Context manager that temporarily swaps the default annotator.

    Useful for sequential comparisons and scoped overrides::

        with varcode.use_annotator("protein_diff"):
            effects = variant_collection.effects()

    Accepts the same argument shape as the ``annotator=`` kwarg:
    a registered-name string, or an annotator instance. Passing an
    instance registers it temporarily under its ``.name`` so that
    name-based lookups inside the block find it; on exit the
    previous default and any previously-registered annotator under
    that name are restored.

    Parameters
    ----------
    name_or_instance : str or EffectAnnotator
        Registered name or implementation to select for the block.

    Yields
    ------
    None
        The context manager does not bind an annotator with ``as``.

    Notes
    -----
    The override is process-global, not thread- or task-local. Concurrent
    callers should pass ``annotator=`` explicitly. Restoration also occurs
    when the block raises an exception.
    """
    global _DEFAULT_NAME
    prior_default = _DEFAULT_NAME

    if isinstance(name_or_instance, str):
        if name_or_instance not in _REGISTRY:
            raise KeyError(
                "No annotator registered under %r — register one first or "
                "pass an instance." % name_or_instance)
        _DEFAULT_NAME = name_or_instance
        prior_registration = None
    else:
        # Keep the runtime check for callers without static checking.
        name = cast(str, getattr(name_or_instance, "name", None))
        if not name:
            raise ValueError(
                "Annotator instance has no .name attribute; cannot scope.")
        prior_registration = _REGISTRY.get(name)
        _REGISTRY[name] = name_or_instance
        _DEFAULT_NAME = name

    try:
        yield
    finally:
        _DEFAULT_NAME = prior_default
        if not isinstance(name_or_instance, str):
            if prior_registration is None:
                _REGISTRY.pop(name, None)
            else:
                _REGISTRY[name] = prior_registration


# Register the default and opt-in experiments. No input-dependent routing.
register_annotator(FastEffectAnnotator())
register_annotator(ProteinDiffEffectAnnotator())
register_annotator(TranscriptModelEffectAnnotator())
# Backward-compatible selection; new results carry the canonical name.
_REGISTRY["realized"] = _REGISTRY["transcript_model"]
