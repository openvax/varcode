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

"""Process-global registry for :class:`EffectAnnotator` instances.

Kept as a module-level dict (not a class) to match the flat registry
pattern in the rest of varcode. Default selection is
``"protein_diff"`` — the protein-diff annotator is self-consistent by
construction and produces HGVS-canonical output everywhere. ``"fast"``
remains available as an opt-in for performance-sensitive pipelines or
byte-for-byte compatibility with 2.x output. See #271.
"""

from contextlib import contextmanager

from .fast import FastEffectAnnotator
from .protein_diff import ProteinDiffEffectAnnotator


class UnsupportedVariantError(ValueError):
    """Raised when an :class:`EffectAnnotator` is asked to handle a
    variant kind outside its declared ``supports`` set.

    Prefer this over silent mis-annotation — the whole point of the
    pluggable-annotator design is that callers can see exactly which
    annotator handles which variant kinds.
    """
    pass


_REGISTRY = {}
_DEFAULT_NAME = "protein_diff"


def register_annotator(annotator):
    """Add an annotator to the process-global registry, keyed by its
    ``.name``. Re-registering under the same name overrides the
    previous entry — this is deliberate so callers can swap
    implementations in tests.
    """
    name = getattr(annotator, "name", None)
    if not name:
        raise ValueError(
            "Annotator %r has no .name attribute; cannot register." % annotator)
    _REGISTRY[name] = annotator
    return annotator


def get_annotator(name):
    """Look up a registered annotator by name. Raises ``KeyError``
    if no annotator is registered under that name.
    """
    return _REGISTRY[name]


def get_default_annotator():
    """Return the annotator currently configured as the default.

    Current default is ``"protein_diff"`` (#322–#327 closed the last
    known correctness bugs between the two). ``"fast"`` stays
    available as an opt-in.
    """
    return _REGISTRY[_DEFAULT_NAME]


def set_default_annotator(name):
    """Swap the process-wide default annotator. ``name`` must refer
    to a registered annotator.
    """
    global _DEFAULT_NAME
    if name not in _REGISTRY:
        raise KeyError(
            "No annotator registered under %r — call register_annotator() "
            "first or pick from %r." % (name, sorted(_REGISTRY)))
    _DEFAULT_NAME = name


def resolve_annotator(annotator_or_name):
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
def use_annotator(name_or_instance):
    """Context manager that temporarily swaps the default annotator.

    Useful for A/B comparisons and scoped overrides without mutating
    global state across the codebase::

        with varcode.use_annotator("protein_diff"):
            effects = variant_collection.effects()

    Accepts the same argument shape as the ``annotator=`` kwarg:
    a registered-name string, or an annotator instance. Passing an
    instance registers it temporarily under its ``.name`` so that
    name-based lookups inside the block find it; on exit the
    previous default and any previously-registered annotator under
    that name are restored.
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
        name = getattr(name_or_instance, "name", None)
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


# Register built-in annotators at import time.
register_annotator(FastEffectAnnotator())
register_annotator(ProteinDiffEffectAnnotator())
