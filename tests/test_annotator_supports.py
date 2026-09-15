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

"""An annotator is only given variants whose kind is in its ``supports``
set. Before this was enforced, ``annotator="fast"`` on a structural
deletion returned ``ExonLoss`` instead of ``LargeDeletion`` (#412)."""

import pytest
from pyensembl import cached_release

from varcode import (
    FastEffectAnnotator,
    StructuralVariant,
    UnsupportedVariantError,
    Variant,
    VariantCollection,
    use_annotator,
)
from varcode.annotators.registry import variant_kind
from varcode.effects.effect_classes import LargeDeletion


ensembl_grch38 = cached_release(81)

CFTR_ID = "ENST00000003084"


def _cftr_deletion():
    cftr = ensembl_grch38.transcript_by_id(CFTR_ID)
    return StructuralVariant(
        contig="7",
        start=cftr.start + 100,
        end=cftr.start + 50_000,
        sv_type="DEL",
        genome=ensembl_grch38)


def _cftr_snv():
    return Variant("7", 117531115, "G", "A", ensembl_grch38)


@pytest.mark.parametrize("variant, kind", [
    (Variant("7", 117531115, "G", "A", ensembl_grch38), "snv"),
    (Variant("7", 117531115, "G", "GA", ensembl_grch38), "indel"),
    (Variant("7", 117531115, "GA", "G", ensembl_grch38), "indel"),
    (Variant("7", 117531115, "GA", "TC", ensembl_grch38), "mnv"),
    (Variant("7", 117531115, "GA", "C", ensembl_grch38), "mnv"),
])
def test_variant_kind_of_point_variants(variant, kind):
    assert variant_kind(variant) == kind


def test_variant_kind_of_structural_variant_is_its_sv_type():
    sv = _cftr_deletion()
    # The placeholder ref/alt looks like an SNV; the kind must not.
    assert sv.is_snv
    assert variant_kind(sv) == "DEL"


@pytest.mark.parametrize("annotator", ["fast", "protein_diff"])
def test_point_annotator_refuses_structural_variant(annotator):
    with pytest.raises(UnsupportedVariantError, match="DEL"):
        _cftr_deletion().effects(annotator=annotator)


@pytest.mark.parametrize("annotator", ["fast", "protein_diff"])
def test_point_annotator_refuses_structural_variant_in_collection(annotator):
    collection = VariantCollection([_cftr_snv(), _cftr_deletion()])
    with pytest.raises(UnsupportedVariantError):
        collection.effects(annotator=annotator)


def test_refusal_is_not_swallowed_by_raise_on_error_false():
    with pytest.raises(UnsupportedVariantError):
        _cftr_deletion().effects(annotator="fast", raise_on_error=False)


def test_annotator_instance_is_checked_too():
    with pytest.raises(UnsupportedVariantError):
        _cftr_deletion().effects(annotator=FastEffectAnnotator())


def test_structural_annotator_refuses_point_variant():
    with pytest.raises(UnsupportedVariantError, match="snv"):
        _cftr_snv().effects(annotator="structural_variant")


def test_default_annotator_routes_structural_variant():
    effects = _cftr_deletion().effects()
    assert effects.annotator == "structural_variant"
    assert any(isinstance(e, LargeDeletion) for e in effects)


def test_scoped_point_default_still_routes_structural_variant():
    with use_annotator("protein_diff"):
        effects = _cftr_deletion().effects()
    assert effects.annotator == "structural_variant"


def test_annotator_without_supports_is_not_checked():
    class NoSupports:
        name = "no_supports"

        def annotate_on_transcript(self, variant, transcript):
            return FastEffectAnnotator().annotate_on_transcript(
                variant, transcript)

    effects = _cftr_snv().effects(annotator=NoSupports())
    assert effects.annotator == "no_supports"
    assert len(effects) > 0
