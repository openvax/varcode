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

"""Runtime behavioral tests for typed collection traversal and effect contracts (#653)."""

from __future__ import annotations

import pytest
from pyensembl import cached_release

from varcode import (
    EffectCollection,
    GeneLike,
    NormalSplicing,
    TranscriptLike,
    Unresolved,
    Variant,
    VariantCollection,
    effect_from_json,
    effect_to_json,
)


ENSEMBL81 = cached_release(81)


class _TaggedVariantCollection(VariantCollection):
    pass


class _TaggedEffectCollection(EffectCollection):
    pass


def test_variant_collection_subclass_and_metadata_preservation():
    v1 = Variant("7", 117559590, "A", "G", genome=ENSEMBL81)
    v2 = Variant("17", 7675088, "C", "T", genome=ENSEMBL81)
    meta = {"test.vcf": {v1: {"id": "rs1"}, v2: {"id": "rs2"}}}
    vc = _TaggedVariantCollection([v1, v2], source_to_metadata_dict=meta)

    assert isinstance(vc.filter(lambda v: v.contig == "7"), _TaggedVariantCollection)
    assert isinstance(vc.clone_with_new_elements([v1]), _TaggedVariantCollection)
    assert isinstance(vc.union(vc), _TaggedVariantCollection)
    assert isinstance(vc.intersection(vc), _TaggedVariantCollection)
    assert isinstance(vc.difference(VariantCollection([v2])), _TaggedVariantCollection)
    assert isinstance(vc.clone_without_ucsc_data(), _TaggedVariantCollection)

    # Slicing preserves Sercol's underlying tuple/sequence return rather than
    # silently wrapping into a collection.
    sliced = vc[:1]
    assert not isinstance(sliced, VariantCollection)
    assert list(sliced) == [vc[0]]

    by_gene = vc.groupby_gene_name()
    assert "CFTR" in by_gene
    assert isinstance(by_gene["CFTR"], _TaggedVariantCollection)
    assert by_gene["CFTR"].metadata[v1]["id"] == "rs1"


def test_effect_collection_subclass_and_provenance_preservation():
    v1 = Variant("7", 117559590, "A", "G", genome=ENSEMBL81)
    effects = v1.effects()
    tagged = _TaggedEffectCollection(
        list(effects),
        annotator="custom_annotator",
        annotator_version="1.0",
        annotated_at="2026-11-20T00:00:00+00:00",
    )
    filtered = tagged.drop_silent_and_noncoding()
    assert isinstance(filtered, _TaggedEffectCollection)
    assert filtered.annotator == "custom_annotator"
    assert filtered.annotator_version == "1.0"

    by_gene = tagged.groupby_gene_id()
    for group in by_gene.values():
        assert isinstance(group, _TaggedEffectCollection)
        assert group.annotator == "custom_annotator"

    top_by_var = tagged.top_priority_effect_per_variant()
    assert v1 in top_by_var
    assert top_by_var[v1].variant == v1


def test_normal_splicing_delegated_properties_remain_read_only():
    tx = ENSEMBL81.transcript_by_id("ENST00000003084")
    v = Variant("7", 117559590, "A", "G", genome=ENSEMBL81)
    coding = v.effect_on_transcript(tx)
    normal = NormalSplicing(v, tx, splice_signal=None, coding_effect=coding)
    assert normal.mutant_protein_sequence == coding.mutant_protein_sequence
    assert normal.aa_mutation_start_offset == coding.aa_mutation_start_offset
    with pytest.raises(AttributeError):
        normal.mutant_protein_sequence = "M"
    with pytest.raises(AttributeError):
        normal.aa_mutation_start_offset = 0
    with pytest.raises(AttributeError):
        normal.aa_mutation_end_offset = 1
    with pytest.raises(AttributeError):
        normal.mutant_transcript = None


def test_unresolved_mutable_change_flags_round_trip_through_json():
    tx = ENSEMBL81.transcript_by_id("ENST00000003084")
    v = Variant("7", 117559590, "A", "G", genome=ENSEMBL81)
    eff = Unresolved(v, tx, mechanism="missing_sequence", reason="no FASTA")
    assert eff.modifies_coding_sequence is None
    assert eff.modifies_protein_sequence is None
    eff.modifies_coding_sequence = True
    eff.modifies_protein_sequence = False
    assert eff.modifies_coding_sequence is True
    assert eff.modifies_protein_sequence is False

    restored = effect_from_json(effect_to_json(eff))
    assert restored.modifies_coding_sequence is True
    assert restored.modifies_protein_sequence is False


def test_protocols_exported_and_match_concrete_and_custom_views():
    # Protocols are exported at package root and describe the actual attributes
    # present on PyEnsembl transcripts/genes without requiring inheritance.
    assert GeneLike is not None
    assert TranscriptLike is not None
    tx = ENSEMBL81.transcript_by_id("ENST00000003084")
    gene = tx.gene
    for attr in ("id", "name", "contig", "start", "end", "strand"):
        assert getattr(gene, attr) is not None
    for attr in (
        "id", "name", "contig", "start", "end", "strand",
        "gene_id", "protein_sequence", "complete", "is_protein_coding", "gene",
    ):
        assert hasattr(tx, attr)
