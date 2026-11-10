"""Transcript names break ties between effects, and unnamed transcripts
don't crash effect selection (#727).

``name_of_associated_transcript`` returned ``""`` for every named transcript
and ``None`` for an unnamed one, so ``top_priority_effect`` raised
``AttributeError`` on annotations without transcript names, and the
documented preference for ``TP53-201`` over ``TP53-206`` never applied.
"""

from types import SimpleNamespace

import pytest

from varcode import Variant
from varcode.effects.effect_ordering import (
    name_of_associated_transcript,
    tie_breaking_sort_key_for_single_gene_effects,
)

from . import native_reference_fixture
from .native_reference_fixture import CASES

native_genomes = native_reference_fixture.native_genomes


def _effect(name, exons=3):
    return SimpleNamespace(transcript=SimpleNamespace(name=name, exons=[None] * exons))


def test_name_of_associated_transcript():
    assert name_of_associated_transcript(_effect("TP53-201")) == "TP53-201"
    assert name_of_associated_transcript(_effect(None)) == ""
    assert name_of_associated_transcript(SimpleNamespace(transcript=None)) == ""


@pytest.mark.parametrize("names, best", [
    (["TP53-206", "TP53-201", "TP53-202"], "TP53-201"),
    (["TP53-203", "TP53-202"], "TP53-202"),
    ([None, "TP53-201"], "TP53-201"),
])
def test_transcript_names_break_ties(names, best):
    effects = [_effect(name) for name in names]
    chosen = max(effects, key=tie_breaking_sort_key_for_single_gene_effects)
    assert chosen.transcript.name == best


def test_more_exons_still_win_over_the_name():
    effects = [_effect("TP53-201", exons=3), _effect("TP53-206", exons=4)]
    chosen = max(effects, key=tie_breaking_sort_key_for_single_gene_effects)
    assert chosen.transcript.name == "TP53-206"


@pytest.mark.parametrize("case", CASES, ids=lambda case: case["id"])
def test_top_priority_effect_on_native_annotations(native_genomes, case):
    # Ensembl's rapid-release CHM13 annotation in the fixture names no
    # transcripts; top_priority_effect raised AttributeError on it.
    variant = Variant(case["contig"], case["position"], case["ref"], case["alt"],
                      genome=native_genomes[case["genome"]])
    effects = variant.effects()
    assert effects.top_priority_effect() in effects
    assert effects.top_priority_effect_per_variant()[variant] in effects
