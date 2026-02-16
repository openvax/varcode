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

import pytest

from varcode.effects import EffectCollection
from varcode.effects.effect_classes import (
    Intronic,
    IntronicSpliceSite,
    SpliceAcceptor,
    SpliceDonor,
)
from varcode.effects.effect_prediction import (
    choose_intronic_effect_class,
    predict_variant_effects,
)


class VariantWithGeneLookupError(object):
    @property
    def gene_ids(self):
        raise RuntimeError("simulated lookup failure")

    @property
    def transcripts(self):
        return []


class DummyVariant(object):
    def __init__(self, start, end, is_insertion=False):
        self.trimmed_base1_start = start
        self.trimmed_base1_end = end
        self.is_insertion = is_insertion


class DummyExon(object):
    def __init__(self, strand, start, end):
        self.strand = strand
        self.start = start
        self.end = end


def test_predict_variant_effects_returns_collection_on_lookup_error():
    effects = predict_variant_effects(
        variant=VariantWithGeneLookupError(),
        raise_on_error=False)
    assert isinstance(effects, EffectCollection)
    assert len(effects) == 0


def test_predict_variant_effects_reraises_on_lookup_error():
    with pytest.raises(RuntimeError):
        predict_variant_effects(
            variant=VariantWithGeneLookupError(),
            raise_on_error=True)


@pytest.mark.parametrize(
    "variant,exon,distance_to_exon,expected_effect_class",
    [
        # + strand, before exon, distance 1-2 => splice acceptor
        (DummyVariant(start=9, end=9), DummyExon(strand="+", start=10, end=20), 1, SpliceAcceptor),
        # + strand, after exon, distance 1-2 => splice donor
        (DummyVariant(start=21, end=21), DummyExon(strand="+", start=10, end=20), 2, SpliceDonor),
        # + strand, after exon, distance 3-6 => intronic splice site
        (DummyVariant(start=21, end=21), DummyExon(strand="+", start=10, end=20), 6, IntronicSpliceSite),
        # + strand, before exon, distance 3 => intronic splice site
        (DummyVariant(start=7, end=7), DummyExon(strand="+", start=10, end=20), 3, IntronicSpliceSite),
        # + strand, far from exon => intronic
        (DummyVariant(start=6, end=6), DummyExon(strand="+", start=10, end=20), 4, Intronic),
        # + strand insertion exactly at exon start is before exon
        (DummyVariant(start=10, end=10, is_insertion=True), DummyExon(strand="+", start=10, end=20), 1, SpliceAcceptor),
        # - strand, genomic position after exon is before exon in transcript direction
        (DummyVariant(start=21, end=21), DummyExon(strand="-", start=10, end=20), 1, SpliceAcceptor),
        # - strand, genomic position before exon is after exon in transcript direction
        (DummyVariant(start=9, end=9), DummyExon(strand="-", start=10, end=20), 2, SpliceDonor),
        # - strand insertion exactly at exon end counts as before exon
        (DummyVariant(start=20, end=20, is_insertion=True), DummyExon(strand="-", start=10, end=20), 2, SpliceAcceptor),
        # - strand, after exon (transcript direction), distance 3-6 => intronic splice site
        (DummyVariant(start=8, end=8), DummyExon(strand="-", start=10, end=20), 6, IntronicSpliceSite),
        # - strand, far from exon => intronic
        (DummyVariant(start=8, end=8), DummyExon(strand="-", start=10, end=20), 7, Intronic),
    ],
)
def test_choose_intronic_effect_class_paths(
        variant,
        exon,
        distance_to_exon,
        expected_effect_class):
    effect_class = choose_intronic_effect_class(
        variant=variant,
        nearest_exon=exon,
        distance_to_exon=distance_to_exon)
    assert effect_class is expected_effect_class


def test_choose_intronic_effect_class_rejects_nonpositive_distance():
    with pytest.raises(AssertionError):
        choose_intronic_effect_class(
            variant=DummyVariant(start=9, end=9),
            nearest_exon=DummyExon(strand="+", start=10, end=20),
            distance_to_exon=0)
