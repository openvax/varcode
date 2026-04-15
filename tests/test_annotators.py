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

"""Tests for the pluggable EffectAnnotator interface and registry
(openvax/varcode#271, stage 1)."""

import pytest
from pyensembl import cached_release

import varcode
from varcode import (
    EffectAnnotator,
    LegacyEffectAnnotator,
    UnsupportedVariantError,
    Variant,
    get_annotator,
    get_default_annotator,
    register_annotator,
    set_default_annotator,
)


ensembl_grch38 = cached_release(81)


# ====================================================================
# Protocol shape
# ====================================================================


def test_legacy_annotator_satisfies_protocol():
    annotator = LegacyEffectAnnotator()
    # @runtime_checkable Protocol — isinstance works structurally.
    assert isinstance(annotator, EffectAnnotator)
    assert annotator.name == "legacy"
    assert {"snv", "indel", "mnv"}.issubset(annotator.supports)


def test_duck_typed_annotator_satisfies_protocol():
    # Third parties don't need to inherit from anything; matching the
    # shape is enough.
    class DuckAnnotator:
        name = "duck"
        supports = frozenset({"snv"})
        def annotate_on_transcript(self, variant, transcript):
            return None
    assert isinstance(DuckAnnotator(), EffectAnnotator)


# ====================================================================
# Registry
# ====================================================================


def test_legacy_annotator_is_registered_by_default():
    assert get_annotator("legacy").__class__ is LegacyEffectAnnotator
    assert get_default_annotator().__class__ is LegacyEffectAnnotator


def test_register_and_retrieve_custom_annotator():
    class CustomAnnotator:
        name = "test_custom"
        supports = frozenset({"snv"})
        def annotate_on_transcript(self, variant, transcript):
            return None
    register_annotator(CustomAnnotator())
    try:
        assert get_annotator("test_custom").__class__ is CustomAnnotator
    finally:
        # Clean up so we don't leak state into other tests.
        from varcode.annotators.registry import _REGISTRY
        _REGISTRY.pop("test_custom", None)


def test_register_rejects_nameless_annotators():
    class Unnamed:
        supports = frozenset()
        def annotate_on_transcript(self, variant, transcript):
            return None
    with pytest.raises(ValueError):
        register_annotator(Unnamed())


def test_set_default_annotator_swaps_registry_default():
    class Swap:
        name = "test_swap_default"
        supports = frozenset({"snv"})
        def annotate_on_transcript(self, variant, transcript):
            return None
    register_annotator(Swap())
    try:
        set_default_annotator("test_swap_default")
        assert get_default_annotator().__class__ is Swap
    finally:
        set_default_annotator("legacy")
        from varcode.annotators.registry import _REGISTRY
        _REGISTRY.pop("test_swap_default", None)


def test_set_default_annotator_rejects_unknown_name():
    with pytest.raises(KeyError):
        set_default_annotator("nonexistent_annotator")


# ====================================================================
# LegacyEffectAnnotator end-to-end — byte-for-byte match with the
# existing Variant.effect_on_transcript API.
# ====================================================================


def test_legacy_annotator_matches_effect_on_transcript():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id("ENST00000003084")
    direct = variant.effect_on_transcript(transcript)
    annotated = LegacyEffectAnnotator().annotate_on_transcript(
        variant, transcript)
    assert type(annotated) is type(direct)
    assert annotated.short_description == direct.short_description


# ====================================================================
# UnsupportedVariantError is available as an exception class for the
# sequence-diff annotator to raise. No code throws it yet (no
# annotator currently checks `.supports` at runtime), but downstream
# code can already catch it.
# ====================================================================


def test_unsupported_variant_error_is_a_value_error():
    # Users can `except ValueError` and catch unsupported-variant
    # errors alongside other validation failures if they want.
    assert issubclass(UnsupportedVariantError, ValueError)
    with pytest.raises(UnsupportedVariantError):
        raise UnsupportedVariantError("test")


# ====================================================================
# Package-level exports
# ====================================================================


def test_annotator_types_exported_at_package_level():
    assert varcode.EffectAnnotator is EffectAnnotator
    assert varcode.LegacyEffectAnnotator is LegacyEffectAnnotator
    assert varcode.UnsupportedVariantError is UnsupportedVariantError
    # Registry functions too
    assert varcode.get_annotator is get_annotator
    assert varcode.get_default_annotator is get_default_annotator
    assert varcode.register_annotator is register_annotator
    assert varcode.set_default_annotator is set_default_annotator
