"""Varcode annotates every small variant the OpenVax libraries share.

openvax-v2, the OpenVax libraries' shared Sid test data published by
osteosarc (iskandr/osteosarc#56), lists in its recipe each target the
libraries test. Its 179 current catalogue alleles are the ready alleles of
``tests/data/osteosarc_variants.json``. The rest are alleles Isovar, Topiary
and Vaxrank depend on: the historical MAP2 deletion, NTF3's compound
substitution, three count-export indels, NR2F2 at its GRCh37 position and
MT_ND5 in two mitochondrial conventions. Both annotators must annotate each
one and agree on every transcript. The first run downloads the bundle
(28 MB) into the osteosarc cache; later runs reuse it offline.
"""

from collections import Counter
import json
from pathlib import Path

import pytest
from pyensembl import cached_release

from varcode import Variant
from varcode.errors import ReferenceMismatchError

from .osteosarc_variants import read_fixture


ENSEMBL_RELEASES = {"GRCh38": 81, "GRCh37": 75}
# hg19's chrM is NC_001807, not the rCRS that Ensembl names MT: hg19
# chrM:12995 is rCRS MT:12994, the target MT_ND5-chrM-12994-GRCh37. Varcode
# refuses the hg19 coordinate rather than annotate a different base.
HG19_MITOCHONDRIAL = "MT_ND5-chrM-12994-hg19"


@pytest.fixture(scope="module")
def small_variant_targets():
    osteosarc = pytest.importorskip("osteosarc", reason="install .[test-data]")
    recipe = json.loads(
        (Path(osteosarc.fetch_bundle("openvax-v2")) / "recipe.json").read_text())
    return {name: target for name, target in recipe["targets"].items()
            if target["kind"] == "small_variant"}


def shared_variant(target):
    return Variant(
        target["contig"], target["position"], target["ref"], target["alt"],
        genome=cached_release(ENSEMBL_RELEASES[target["assembly"]]),
        convert_ucsc_contig_names=True)


def effects_by_transcript(variant, annotator):
    effects = variant.effects(raise_on_error=True, annotator=annotator)
    return {getattr(getattr(effect, "transcript", None), "transcript_id", None):
            (type(effect).__name__, effect.short_description)
            for effect in effects}


def test_current_targets_are_the_fixture_alleles(small_variant_targets):
    current = {(t["contig"], t["position"], t["ref"], t["alt"])
               for t in small_variant_targets.values() if t["label"] == "current"}
    ready = {tuple(allele) for entry in read_fixture()["entries"]
             if entry["status"] == "ready" for allele in entry["alleles"]}
    assert current == ready
    assert Counter(t["label"] for t in small_variant_targets.values()) == {
        "current": 179, "count_export": 3, "mitochondrial": 2,
        "historical": 1, "native_grch37": 1, "library_allele": 1}


def test_every_shared_small_variant_annotates_and_the_annotators_agree(
        small_variant_targets):
    for name, target in sorted(small_variant_targets.items()):
        variant = shared_variant(target)
        if name == HG19_MITOCHONDRIAL:
            with pytest.raises(ReferenceMismatchError):
                variant.effects(annotator="fast")
            continue
        fast = effects_by_transcript(variant, "fast")
        assert fast, name
        assert effects_by_transcript(variant, "protein_diff") == fast, name
