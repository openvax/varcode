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
    FastEffectAnnotator,
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


def test_fast_annotator_annotator_satisfies_protocol():
    annotator = FastEffectAnnotator()
    # @runtime_checkable Protocol — isinstance works structurally.
    assert isinstance(annotator, EffectAnnotator)
    assert annotator.name == "fast"
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


def test_fast_annotator_annotator_is_registered():
    assert get_annotator("fast").__class__ is FastEffectAnnotator


def test_default_annotator_matches_active_configuration(request):
    configured = request.config.getoption("--annotator") or "fast"
    default = get_default_annotator()
    assert default.name == configured
    assert default is get_annotator(configured)


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
        set_default_annotator("fast")
        from varcode.annotators.registry import _REGISTRY
        _REGISTRY.pop("test_swap_default", None)


def test_set_default_annotator_rejects_unknown_name():
    with pytest.raises(KeyError):
        set_default_annotator("nonexistent_annotator")


# ====================================================================
# FastEffectAnnotator end-to-end — byte-for-byte match with the
# existing Variant.effect_on_transcript API.
# ====================================================================


def test_fast_annotator_annotator_matches_effect_on_transcript():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id("ENST00000003084")
    direct = variant.effect_on_transcript(transcript)
    annotated = FastEffectAnnotator().annotate_on_transcript(
        variant, transcript)
    assert type(annotated) is type(direct)
    assert annotated.short_description == direct.short_description


# ====================================================================
# UnsupportedVariantError is available as an exception class for the
# protein-diff annotator to raise. No code throws it yet (no
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
    assert varcode.FastEffectAnnotator is FastEffectAnnotator
    assert varcode.UnsupportedVariantError is UnsupportedVariantError
    # Registry functions too
    assert varcode.get_annotator is get_annotator
    assert varcode.get_default_annotator is get_default_annotator
    assert varcode.register_annotator is register_annotator
    assert varcode.set_default_annotator is set_default_annotator


# ====================================================================
# annotator= kwarg on Variant.effects() and VariantCollection.effects()
# (see #271 stage 3a)
# ====================================================================


def test_variant_effects_default_uses_fast_annotator():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    # No annotator kwarg — should still produce the same output the
    # fast annotator produces, since fast is the default.
    default_effects = variant.effects()
    explicit_fast_annotator = variant.effects(annotator="fast")
    assert [type(e).__name__ for e in default_effects] == \
           [type(e).__name__ for e in explicit_fast_annotator]


def test_variant_effects_accepts_annotator_name_string():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    effects = variant.effects(annotator="fast")
    assert len(effects) > 0


def test_variant_effects_accepts_annotator_instance():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    effects = variant.effects(annotator=FastEffectAnnotator())
    assert len(effects) > 0


def test_variant_effects_rejects_unknown_annotator_name():
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    with pytest.raises(KeyError) as exc_info:
        variant.effects(annotator="nonexistent_annotator")
    # Error message lists known annotators so the caller can fix the typo.
    assert "fast" in str(exc_info.value)


def test_variant_collection_effects_accepts_annotator():
    from varcode import VariantCollection
    vc = VariantCollection([
        Variant("7", 117531115, "G", "A", ensembl_grch38),
        Variant("7", 117531114, "G", "T", ensembl_grch38),
    ])
    default_effects = vc.effects()
    with_annotator = vc.effects(annotator="fast")
    assert len(default_effects) == len(with_annotator)


def test_variant_collection_effects_preserves_annotator_provenance():
    from varcode import VariantCollection

    vc = VariantCollection([
        Variant("7", 117531115, "G", "A", ensembl_grch38),
    ])

    effects = vc.effects(annotator="protein_diff")
    assert effects.annotator == "protein_diff"
    assert effects.annotator_version is not None
    assert effects.annotated_at is not None


# ====================================================================
# use_annotator context manager
# ====================================================================


def test_use_annotator_swaps_default_within_scope():
    from varcode import use_annotator
    class Scoped:
        name = "test_scoped"
        supports = frozenset({"snv"})
        def annotate_on_transcript(self, variant, transcript):
            return None
    scoped_instance = Scoped()
    prior_default = get_default_annotator()
    with use_annotator(scoped_instance):
        assert get_default_annotator() is scoped_instance
    # On exit, default is restored to whatever it was before.
    assert get_default_annotator() is prior_default
    # And the scoped-instance registration is cleaned up.
    from varcode.annotators.registry import _REGISTRY
    assert "test_scoped" not in _REGISTRY


def test_use_annotator_with_registered_name_preserves_registration():
    from varcode import use_annotator, register_annotator
    class Persisted:
        name = "test_persisted"
        supports = frozenset({"snv"})
        def annotate_on_transcript(self, variant, transcript):
            return None
    persisted = Persisted()
    register_annotator(persisted)
    prior_default = get_default_annotator()
    try:
        with use_annotator("test_persisted"):
            assert get_default_annotator() is persisted
        assert get_default_annotator() is prior_default
        # Registration persists after the context exit since it
        # wasn't created by the context manager.
        from varcode.annotators.registry import _REGISTRY
        assert "test_persisted" in _REGISTRY
    finally:
        from varcode.annotators.registry import _REGISTRY
        _REGISTRY.pop("test_persisted", None)


def test_use_annotator_affects_variant_effects_default():
    from varcode import use_annotator
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)

    # Capture fast output inside the context via a wrapping annotator
    # that records the calls.
    calls = []
    class Recording:
        name = "test_recording"
        supports = frozenset({"snv", "indel", "mnv"})
        def annotate_on_transcript(self, variant, transcript):
            calls.append((variant, transcript))
            return FastEffectAnnotator().annotate_on_transcript(
                variant, transcript)
    with use_annotator(Recording()):
        variant.effects()  # no annotator= kwarg → uses scoped default
    assert len(calls) > 0, (
        "Variant.effects() inside use_annotator(...) should have "
        "dispatched to the scoped annotator")


def test_use_annotator_rejects_instance_without_name():
    from varcode import use_annotator
    class Unnamed:
        supports = frozenset()
        def annotate_on_transcript(self, variant, transcript):
            return None
    with pytest.raises(ValueError):
        with use_annotator(Unnamed()):
            pass


def test_resolve_annotator_passes_through_instance():
    from varcode import resolve_annotator
    instance = FastEffectAnnotator()
    assert resolve_annotator(instance) is instance


def test_resolve_annotator_resolves_none_to_default():
    from varcode import resolve_annotator
    resolved = resolve_annotator(None)
    assert resolved is get_default_annotator()


# ====================================================================
# ProteinDiffEffectAnnotator (#309, stage 3d)
# ====================================================================


def test_protein_diff_annotator_is_registered():
    from varcode.annotators import ProteinDiffEffectAnnotator
    annotator = get_annotator("protein_diff")
    assert isinstance(annotator, ProteinDiffEffectAnnotator)


def test_protein_diff_annotator_satisfies_protocol():
    from varcode.annotators import ProteinDiffEffectAnnotator
    assert isinstance(ProteinDiffEffectAnnotator(), EffectAnnotator)


def test_protein_diff_usable_via_kwarg():
    variant = Variant("7", 117531095, "T", "A", ensembl_grch38)
    effects = variant.effects(annotator="protein_diff")
    assert len(effects) > 0


def test_protein_diff_parity_on_coding_snv():
    variant = Variant("7", 117531095, "T", "A", ensembl_grch38)
    fast_effects = list(variant.effects(annotator="fast"))
    sdiff = list(variant.effects(annotator="protein_diff"))
    assert len(fast_effects) == len(sdiff)
    for le, se in zip(fast_effects, sdiff):
        assert type(le).__name__ == type(se).__name__, (
            "Class mismatch: %s vs %s" % (type(le).__name__, type(se).__name__))
        assert le.short_description == se.short_description, (
            "short_description mismatch: %r vs %r" % (
                le.short_description, se.short_description))


def test_protein_diff_parity_on_splice_donor():
    # SpliceDonor: pure-intronic → fast-only path. Both annotators
    # should agree because protein_diff delegates to fast.
    variant = Variant("7", 117531115, "G", "A", ensembl_grch38)
    fast_effects = list(variant.effects(annotator="fast"))
    sdiff = list(variant.effects(annotator="protein_diff"))
    for le, se in zip(fast_effects, sdiff):
        assert type(le).__name__ == type(se).__name__
        assert le.short_description == se.short_description


def test_protein_diff_parity_on_exonic_splice_site():
    # ExonicSpliceSite: dual-dispatch. Splice class from fast,
    # alternate_effect from protein-diff's protein diff.
    variant = Variant("7", 117531114, "G", "T", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id("ENST00000003084")
    legacy = variant.effect_on_transcript(transcript)
    from varcode.annotators import ProteinDiffEffectAnnotator
    sdiff = ProteinDiffEffectAnnotator().annotate_on_transcript(
        variant, transcript)
    from varcode.effects import ExonicSpliceSite
    assert isinstance(legacy, ExonicSpliceSite)
    assert isinstance(sdiff, ExonicSpliceSite)
    # alternate_effect comes from protein-diff in the new annotator.
    assert sdiff.alternate_effect is not None
    assert type(sdiff.alternate_effect).__name__ == type(legacy.alternate_effect).__name__


def test_protein_diff_parity_on_reverse_strand_mnv():
    # BRCA1 reverse-strand MNV — exercises the strand-flip path.
    variant = Variant("17", 43082570, "CCT", "GGG", ensembl_grch38)
    fast_effects = list(variant.effects(annotator="fast"))
    sdiff = list(variant.effects(annotator="protein_diff"))
    for le, se in zip(fast_effects, sdiff):
        assert type(le).__name__ == type(se).__name__
        assert le.short_description == se.short_description


def test_protein_diff_parity_on_mt_codon_table():
    # MT-CO1 TCA→TGA: Substitution to Trp under mt table.
    variant = Variant("MT", 6739, "C", "G", ensembl_grch38)
    fast_effects = list(variant.effects(annotator="fast"))
    sdiff = list(variant.effects(annotator="protein_diff"))
    for le, se in zip(fast_effects, sdiff):
        assert type(le).__name__ == type(se).__name__
        assert le.short_description == se.short_description


def test_protein_diff_parity_on_mt_premature_stop():
    # MT-CO1 CGA→AGA: PrematureStop under mt table.
    variant = Variant("MT", 6015, "C", "A", ensembl_grch38)
    fast_effects = list(variant.effects(annotator="fast"))
    sdiff = list(variant.effects(annotator="protein_diff"))
    for le, se in zip(fast_effects, sdiff):
        assert type(le).__name__ == type(se).__name__
        assert le.short_description == se.short_description


@pytest.mark.parametrize(
    ("variant", "expected_aa_ref", "expected_aa_alt", "expected_start", "expected_end"),
    [
        (Variant("chr1", 99772782, "A", "ATGA", "mm10"), "", "", 6, 6),
        (Variant("chr1", 99772782, "A", "ACCCTGA", "mm10"), "", "P", 6, 7),
    ],
)
def test_protein_diff_matches_fast_on_premature_stop_insertions(
        variant,
        expected_aa_ref,
        expected_aa_alt,
        expected_start,
        expected_end):
    fast_effect = variant.effects(annotator="fast").top_priority_effect()
    protein_diff_effect = variant.effects(
        annotator="protein_diff").top_priority_effect()

    assert type(fast_effect).__name__ == "PrematureStop"
    assert type(protein_diff_effect) is type(fast_effect)
    assert fast_effect.aa_ref == protein_diff_effect.aa_ref == expected_aa_ref
    assert fast_effect.aa_alt == protein_diff_effect.aa_alt == expected_aa_alt
    assert (
        fast_effect.aa_mutation_start_offset,
        fast_effect.aa_mutation_end_offset,
    ) == (
        protein_diff_effect.aa_mutation_start_offset,
        protein_diff_effect.aa_mutation_end_offset,
    ) == (expected_start, expected_end)
