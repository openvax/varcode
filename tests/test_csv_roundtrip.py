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


# -----------------------------------------------------------------------
# Edge cases: empty collections, blank-line tolerance in header parser,
# intergenic fallback path in EffectCollection.from_csv.
# -----------------------------------------------------------------------


def test_variant_collection_empty_csv_roundtrip():
    # Empty collection round-trips without crashing. No reference is
    # available to write in the header, so from_csv needs an explicit
    # genome.
    empty = VariantCollection(variants=[])
    path = _tmp_csv()
    try:
        empty.to_csv(path)
        loaded = VariantCollection.from_csv(path, genome="GRCh38")
    finally:
        os.unlink(path)
    assert len(loaded) == 0


def test_effect_collection_empty_csv_roundtrip():
    empty = EffectCollection(effects=[])
    path = _tmp_csv()
    try:
        empty.to_csv(path)
        loaded = EffectCollection.from_csv(path, genome="GRCh38")
    finally:
        os.unlink(path)
    assert len(loaded) == 0


def test_header_parser_tolerates_blank_lines():
    # Hand-edited CSVs sometimes have a blank line between the comment
    # block and the body; the header parser should skip it rather than
    # treat it as end-of-header.
    variants = [Variant("17", 43082404, "C", "T", "GRCh38")]
    path = _tmp_csv()
    try:
        VariantCollection(variants=variants).to_csv(path)
        # Inject a blank line in the middle of the header.
        with open(path, "r") as f:
            content = f.read()
        lines = content.splitlines(keepends=True)
        # Insert a blank line after the first header line.
        new_content = lines[0] + "\n" + "".join(lines[1:])
        with open(path, "w") as f:
            f.write(new_content)
        # Should still read the genome from the header and work fine.
        loaded = VariantCollection.from_csv(path)
    finally:
        os.unlink(path)
    assert len(loaded) == 1
    assert loaded[0].reference_name == "GRCh38"


def test_effect_collection_roundtrip_preserves_intergenic_effect():
    # An intergenic variant produces an Intergenic effect with no
    # transcript. The CSV row will have an empty transcript_id, and
    # from_csv should reconstruct it via the intergenic fallback path
    # without warning.
    import warnings

    # chr1:100 is far from any gene on GRCh38.
    variant = Variant("1", 100, "A", "T", "GRCh38")
    original = variant.effects()
    intergenic_effects = [
        e for e in original if e.__class__.__name__ == "Intergenic"
    ]
    assert len(intergenic_effects) >= 1, (
        "Test variant should produce at least one Intergenic effect, got %r"
        % [e.__class__.__name__ for e in original]
    )

    path = _tmp_csv()
    try:
        original.to_csv(path)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            loaded = EffectCollection.from_csv(path)
    finally:
        os.unlink(path)

    assert len(loaded) == len(original)
    assert any(e.__class__.__name__ == "Intergenic" for e in loaded)
    # No spurious "could not recover" warnings on the normal path.
    unexpected = [
        w for w in caught if "Could not recover effect" in str(w.message)
    ]
    assert len(unexpected) == 0, (
        "Unexpected fallback warnings during normal intergenic round-trip: %r"
        % [str(w.message) for w in unexpected]
    )


def test_effect_collection_from_csv_warns_when_fallback_cant_match():
    # Hand-craft a CSV with an empty transcript_id and a made-up
    # effect_type that the annotation won't produce. from_csv should
    # warn and drop the row rather than silently omitting it.
    import warnings

    import pandas as pd

    path = _tmp_csv()
    try:
        df = pd.DataFrame([{
            "variant": "chr1 g.100A>T",
            "contig": "1",
            "start": 100,
            "ref": "A",
            "alt": "T",
            "is_snv": True,
            "is_transversion": True,
            "is_transition": False,
            "gene_id": None,
            "gene_name": None,
            "transcript_id": None,
            "transcript_name": None,
            "effect_type": "DefinitelyNotARealEffectClass",
            "effect": "p.fake",
        }])
        # Write a metadata header by hand so from_csv can resolve the
        # genome, then append the body.
        with open(path, "w") as f:
            f.write("# varcode_version=test\n")
            f.write("# reference_name=GRCh38\n")
            df.to_csv(f, index=False)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            loaded = EffectCollection.from_csv(path)
    finally:
        os.unlink(path)

    # The unmatched row should be dropped with a warning.
    assert len(loaded) == 0
    matching = [
        w for w in caught if "Could not recover effect" in str(w.message)
    ]
    assert len(matching) >= 1, \
        "Expected a warning when fallback could not match the effect_type"
    assert "DefinitelyNotARealEffectClass" in str(matching[0].message)
