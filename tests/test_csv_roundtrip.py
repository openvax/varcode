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

"""
Round-trip tests for VariantCollection.from_csv and
EffectCollection.from_csv.

These are the lightweight inverse of the existing .to_csv() path. For
VariantCollection the round-trip is exact (ref/alt/position are
preserved). For EffectCollection the round-trip is semantic — it
records (variant, transcript_id) pairs and re-runs annotation on read,
so the resulting effects match whenever annotation is deterministic.
"""

import os
import tempfile

from varcode import Variant, VariantCollection
from varcode.effects import EffectCollection


def _tmp_csv():
    fd, path = tempfile.mkstemp(suffix=".csv")
    os.close(fd)
    return path


# -----------------------------------------------------------------------
# VariantCollection round-trip
# -----------------------------------------------------------------------


def test_variant_collection_csv_roundtrip_preserves_variants():
    variants = [
        Variant("17", 43082575 - 5, "CCT", "GGG", "GRCh38"),
        Variant("7", 117531114, "G", "T", "GRCh38"),
        Variant("7", 117530898, "G", "A", "GRCh38"),
    ]
    original = VariantCollection(variants=variants)

    path = _tmp_csv()
    try:
        original.to_csv(path)
        loaded = VariantCollection.from_csv(path, genome="GRCh38")
    finally:
        os.unlink(path)

    # Sort key on both collections produces the same ordering.
    assert len(loaded) == len(original)
    for original_variant, loaded_variant in zip(original, loaded):
        assert loaded_variant.contig == original_variant.contig
        assert loaded_variant.start == original_variant.start
        assert loaded_variant.ref == original_variant.ref
        assert loaded_variant.alt == original_variant.alt


def test_variant_collection_csv_roundtrip_handles_indels():
    # Insertion (empty ref) and deletion (empty alt) should survive the
    # CSV round-trip even though pandas may read '' as NaN.
    variants = [
        Variant("17", 43082404, "C", "", "GRCh38"),          # deletion
        Variant("17", 43082575 - 6, "", "AAA", "GRCh38"),    # insertion
    ]
    original = VariantCollection(variants=variants)

    path = _tmp_csv()
    try:
        original.to_csv(path)
        loaded = VariantCollection.from_csv(path, genome="GRCh38")
    finally:
        os.unlink(path)

    assert len(loaded) == 2
    deletions = [v for v in loaded if v.alt == ""]
    insertions = [v for v in loaded if v.ref == ""]
    assert len(deletions) == 1
    assert len(insertions) == 1
    assert deletions[0].ref == "C"
    assert insertions[0].alt == "AAA"


def test_variant_collection_from_csv_rejects_missing_columns():
    # If the CSV is missing required columns, from_csv should raise
    # rather than silently produce empty or wrong variants.
    import pandas as pd
    path = _tmp_csv()
    try:
        pd.DataFrame({"chr": ["17"], "start": [100]}).to_csv(path, index=False)
        try:
            VariantCollection.from_csv(path, genome="GRCh38")
        except ValueError as e:
            assert "ref" in str(e) and "alt" in str(e), \
                "Error should name the missing columns, got %r" % str(e)
        else:
            raise AssertionError("Expected ValueError for missing columns")
    finally:
        os.unlink(path)


# -----------------------------------------------------------------------
# EffectCollection round-trip
# -----------------------------------------------------------------------


def test_effect_collection_csv_roundtrip_preserves_effect_types():
    # Use a variant with a well-known effect signature across several
    # transcripts so we can compare the class distribution pre and post.
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", "GRCh38")
    original = variant.effects()

    path = _tmp_csv()
    try:
        original.to_csv(path)
        loaded = EffectCollection.from_csv(path, genome="GRCh38")
    finally:
        os.unlink(path)

    assert len(loaded) == len(original), (
        "Effect counts differ: original=%d loaded=%d" % (
            len(original), len(loaded))
    )

    # Same class distribution and same short descriptions.
    original_types = sorted(e.__class__.__name__ for e in original)
    loaded_types = sorted(e.__class__.__name__ for e in loaded)
    assert original_types == loaded_types

    original_descs = sorted(e.short_description for e in original)
    loaded_descs = sorted(e.short_description for e in loaded)
    assert original_descs == loaded_descs


def test_effect_collection_csv_roundtrip_preserves_transcripts():
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", "GRCh38")
    original = variant.effects()

    path = _tmp_csv()
    try:
        original.to_csv(path)
        loaded = EffectCollection.from_csv(path, genome="GRCh38")
    finally:
        os.unlink(path)

    # Transcript IDs should be preserved (with None appearing exactly
    # once per original None).
    def transcript_ids(ec):
        ids = [
            str(e.transcript_id) if e.transcript_id else None
            for e in ec
        ]
        return sorted(ids, key=lambda x: (x is None, x or ""))

    assert transcript_ids(loaded) == transcript_ids(original)


def test_effect_collection_from_csv_skips_comment_lines():
    # from_csv should treat extra '#'-prefixed lines as comments that
    # don't interfere with body parsing.
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", "GRCh38")
    original = variant.effects()

    path = _tmp_csv()
    try:
        original.to_csv(path)
        # Prepend a few extra comment lines in addition to the header
        # that to_csv already wrote.
        with open(path, "r") as f:
            body = f.read()
        with open(path, "w") as f:
            f.write("# annotator=legacy\n")
            f.write("# timestamp=2026-04-12T14:30:00Z\n")
            f.write(body)
        loaded = EffectCollection.from_csv(path)
    finally:
        os.unlink(path)

    assert len(loaded) == len(original)


# -----------------------------------------------------------------------
# Metadata header: genome is optional when the header carries it,
# required otherwise.
# -----------------------------------------------------------------------


def test_variant_collection_to_csv_writes_metadata_header():
    variants = [
        Variant("17", 43082404, "C", "T", "GRCh38"),
    ]
    vc = VariantCollection(variants=variants)
    path = _tmp_csv()
    try:
        vc.to_csv(path)
        with open(path) as f:
            head = "".join(f.readline() for _ in range(4))
    finally:
        os.unlink(path)
    assert "# varcode_version=" in head
    assert "# reference_name=GRCh38" in head


def test_variant_collection_from_csv_reads_genome_from_header():
    variants = [
        Variant("17", 43082404, "C", "T", "GRCh38"),
        Variant("7", 117531114, "G", "T", "GRCh38"),
    ]
    original = VariantCollection(variants=variants)
    path = _tmp_csv()
    try:
        original.to_csv(path)
        # Don't pass genome — should be read from the header.
        loaded = VariantCollection.from_csv(path)
    finally:
        os.unlink(path)
    assert len(loaded) == len(original)
    # The loaded variants should have resolved the same reference.
    assert all(v.reference_name == "GRCh38" for v in loaded)


def test_variant_collection_from_csv_explicit_genome_overrides_header():
    # Explicit genome argument takes precedence over header.
    variants = [Variant("17", 43082404, "C", "T", "GRCh38")]
    path = _tmp_csv()
    try:
        VariantCollection(variants=variants).to_csv(path)
        loaded = VariantCollection.from_csv(path, genome="GRCh38")
    finally:
        os.unlink(path)
    assert len(loaded) == 1


def test_variant_collection_from_csv_raises_without_genome_or_header():
    # If the CSV has no header metadata and the caller doesn't pass
    # genome, we should raise a clear error rather than guessing.
    variants = [Variant("17", 43082404, "C", "T", "GRCh38")]
    path = _tmp_csv()
    try:
        VariantCollection(variants=variants).to_csv(path, include_header=False)
        try:
            VariantCollection.from_csv(path)
        except ValueError as e:
            assert "genome" in str(e) and "reference_name" in str(e), \
                "Expected error to mention genome and reference_name, got %r" % str(e)
        else:
            raise AssertionError("Expected ValueError when no genome is available")
    finally:
        os.unlink(path)


def test_effect_collection_to_csv_writes_metadata_header():
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", "GRCh38")
    ec = variant.effects()
    path = _tmp_csv()
    try:
        ec.to_csv(path)
        with open(path) as f:
            head = "".join(f.readline() for _ in range(4))
    finally:
        os.unlink(path)
    assert "# varcode_version=" in head
    assert "# reference_name=GRCh38" in head


def test_effect_collection_from_csv_reads_genome_from_header():
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", "GRCh38")
    original = variant.effects()
    path = _tmp_csv()
    try:
        original.to_csv(path)
        # Don't pass genome — should be read from the header.
        loaded = EffectCollection.from_csv(path)
    finally:
        os.unlink(path)
    assert len(loaded) == len(original)


def test_effect_collection_from_csv_raises_without_genome_or_header():
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", "GRCh38")
    ec = variant.effects()
    path = _tmp_csv()
    try:
        ec.to_csv(path, include_header=False)
        try:
            EffectCollection.from_csv(path)
        except ValueError as e:
            assert "genome" in str(e) and "reference_name" in str(e)
        else:
            raise AssertionError("Expected ValueError when no genome is available")
    finally:
        os.unlink(path)


def test_csv_to_csv_without_header_is_opt_out_legacy_format():
    # to_csv(include_header=False) produces a plain CSV for legacy
    # consumers that don't tolerate comment lines.
    variants = [Variant("17", 43082404, "C", "T", "GRCh38")]
    path = _tmp_csv()
    try:
        VariantCollection(variants=variants).to_csv(path, include_header=False)
        with open(path) as f:
            first_line = f.readline()
    finally:
        os.unlink(path)
    # First line should be the column header, not a '#' comment.
    assert not first_line.startswith("#"), \
        "include_header=False should not produce leading comment lines"
