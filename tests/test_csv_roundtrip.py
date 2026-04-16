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


# -----------------------------------------------------------------------
# Annotator provenance round-trip (#271 stage 3b)
# -----------------------------------------------------------------------


def test_effect_collection_carries_annotator_provenance():
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", "GRCh38")
    ec = variant.effects()
    # predict_variant_effects should have populated the annotator name
    # from whichever annotator is currently the default.
    assert ec.annotator is not None
    assert ec.annotator_version is not None
    assert ec.annotated_at is not None


def test_effect_collection_to_csv_emits_annotator_header():
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", "GRCh38")
    ec = variant.effects()
    path = _tmp_csv()
    try:
        ec.to_csv(path)
        with open(path) as f:
            head = "".join(f.readline() for _ in range(8))
    finally:
        os.unlink(path)
    assert "# annotator=" in head
    assert "# annotator_version=" in head
    assert "# annotated_at=" in head


def test_effect_collection_from_csv_recovers_annotator_metadata():
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", "GRCh38")
    original = variant.effects()
    path = _tmp_csv()
    try:
        original.to_csv(path)
        loaded = EffectCollection.from_csv(path)
    finally:
        os.unlink(path)
    assert loaded.annotator == original.annotator
    assert loaded.annotator_version == original.annotator_version
    # annotated_at is preserved verbatim from the header, not
    # refreshed — so the restored collection remembers when it was
    # originally produced.
    assert loaded.annotated_at == original.annotated_at


def test_effect_collection_clone_preserves_annotator_metadata():
    # filter / groupby / clone_with_new_elements flow through to_dict;
    # the provenance fields must survive so derived collections don't
    # lose their origin.
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", "GRCh38")
    ec = variant.effects()
    cloned = ec.clone_with_new_elements(list(ec.elements))
    assert cloned.annotator == ec.annotator
    assert cloned.annotator_version == ec.annotator_version
    assert cloned.annotated_at == ec.annotated_at


def test_effect_collection_from_csv_warns_on_annotator_mismatch():
    import warnings
    from varcode import get_default_annotator
    variant = Variant("17", 43082575 - 5, "CCT", "GGG", "GRCh38")
    ec = variant.effects()
    path = _tmp_csv()
    # Use a fake annotator name that's guaranteed to differ from
    # whatever the current default is.
    fake_name = "fake_annotator_for_test"
    current_default = get_default_annotator().name
    try:
        ec.to_csv(path)
        with open(path) as f:
            lines = f.readlines()
        with open(path, "w") as f:
            for line in lines:
                if line.startswith("# annotator="):
                    f.write("# annotator=%s\n" % fake_name)
                else:
                    f.write(line)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            EffectCollection.from_csv(path)
    finally:
        os.unlink(path)

    messages = [str(w.message) for w in caught]
    assert any(
        fake_name in m and current_default in m for m in messages), (
        "Expected a warning about annotator mismatch, got: %r" % messages)


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


# -----------------------------------------------------------------------
# #274: VariantCollection and EffectCollection use different column
# names for the contig field ("chr" vs "contig"). Readers should accept
# either so CSVs are interchangeable across the two types.
# -----------------------------------------------------------------------


def test_variant_collection_from_csv_accepts_contig_column_name():
    # A CSV written by EffectCollection-style code uses "contig" but
    # should still load as a VariantCollection.
    import pandas as pd
    path = _tmp_csv()
    try:
        df = pd.DataFrame([
            {"contig": "17", "start": 43082404, "ref": "C", "alt": "T"},
            {"contig": "7", "start": 117531114, "ref": "G", "alt": "T"},
        ])
        with open(path, "w") as f:
            f.write("# reference_name=GRCh38\n")
            df.to_csv(f, index=False)
        loaded = VariantCollection.from_csv(path)
    finally:
        os.unlink(path)
    assert len(loaded) == 2
    starts = sorted(v.start for v in loaded)
    assert starts == [43082404, 117531114]


def test_effect_collection_from_csv_accepts_chr_column_name():
    # An EffectCollection CSV with the VariantCollection-style "chr"
    # column should still load — we only need one contig column name
    # present.
    import pandas as pd
    path = _tmp_csv()
    try:
        df = pd.DataFrame([
            {
                "variant": "chr1 g.100A>T",
                "chr": "1",  # the alternate name
                "start": 100,
                "ref": "A",
                "alt": "T",
                "transcript_id": None,
                "effect_type": "Intergenic",
                "effect": "intergenic",
            }
        ])
        with open(path, "w") as f:
            f.write("# reference_name=GRCh38\n")
            df.to_csv(f, index=False)
        loaded = EffectCollection.from_csv(path)
    finally:
        os.unlink(path)
    # Intergenic variant at chr1:100 produces at least one Intergenic
    # effect after re-annotation.
    assert any(e.__class__.__name__ == "Intergenic" for e in loaded)


def test_variant_collection_from_csv_errors_without_contig_column():
    import pandas as pd
    path = _tmp_csv()
    try:
        df = pd.DataFrame([
            # No 'chr' or 'contig' column at all.
            {"start": 100, "ref": "A", "alt": "T"},
        ])
        with open(path, "w") as f:
            f.write("# reference_name=GRCh38\n")
            df.to_csv(f, index=False)
        try:
            VariantCollection.from_csv(path)
        except ValueError as e:
            assert "contig" in str(e) and "chr" in str(e), \
                "Error should list both aliases, got %r" % str(e)
        else:
            raise AssertionError("Expected ValueError for missing contig column")
    finally:
        os.unlink(path)


def test_effect_collection_from_csv_errors_without_contig_column():
    import pandas as pd
    path = _tmp_csv()
    try:
        df = pd.DataFrame([
            {
                "variant": "chr1 g.100A>T",
                # No 'chr' or 'contig'.
                "start": 100,
                "ref": "A",
                "alt": "T",
                "transcript_id": None,
                "effect_type": "Intergenic",
                "effect": "intergenic",
            }
        ])
        with open(path, "w") as f:
            f.write("# reference_name=GRCh38\n")
            df.to_csv(f, index=False)
        try:
            EffectCollection.from_csv(path)
        except ValueError as e:
            assert "contig" in str(e) and "chr" in str(e)
        else:
            raise AssertionError("Expected ValueError for missing contig column")
    finally:
        os.unlink(path)


# -----------------------------------------------------------------------
# #275: warn on major-version drift when loading serialized collections.
# Because from_csv re-runs annotation, mismatch between the varcode
# version that wrote the file and the one reading it can produce
# different effects. We only warn on major-version drift since minor
# and patch are semver-compatible.
# -----------------------------------------------------------------------


def _write_csv_with_header_version(path, version, body_df):
    with open(path, "w") as f:
        f.write("# varcode_version=%s\n" % version)
        f.write("# reference_name=GRCh38\n")
        body_df.to_csv(f, index=False)


def test_from_csv_warns_on_major_version_drift():
    import warnings
    import pandas as pd

    path = _tmp_csv()
    try:
        body = pd.DataFrame([
            {"chr": "17", "start": 43082404, "ref": "C", "alt": "T"},
        ])
        # Pretend this CSV was written by an ancient 1.0.0.
        _write_csv_with_header_version(path, "1.0.0", body)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            VariantCollection.from_csv(path)
    finally:
        os.unlink(path)

    drift = [w for w in caught if "major versions" in str(w.message)]
    assert len(drift) == 1, \
        "Expected a single major-version drift warning, got %d: %r" % (
            len(drift), [str(w.message) for w in caught])
    assert "1.0.0" in str(drift[0].message)


def test_from_csv_does_not_warn_on_patch_or_minor_drift():
    import warnings
    import pandas as pd
    import varcode

    path = _tmp_csv()
    try:
        body = pd.DataFrame([
            {"chr": "17", "start": 43082404, "ref": "C", "alt": "T"},
        ])
        # Same major (2.x) as current install; different minor/patch.
        # Compute the major from varcode.__version__.
        current_major = varcode.__version__.split(".")[0]
        fake_version = "%s.0.1" % current_major
        _write_csv_with_header_version(path, fake_version, body)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            VariantCollection.from_csv(path)
    finally:
        os.unlink(path)

    drift = [w for w in caught if "major versions" in str(w.message)]
    assert len(drift) == 0, \
        "Should not warn on within-major drift, got: %r" % (
            [str(w.message) for w in drift])


def test_from_csv_ignores_malformed_version_header():
    import warnings
    import pandas as pd

    path = _tmp_csv()
    try:
        body = pd.DataFrame([
            {"chr": "17", "start": 43082404, "ref": "C", "alt": "T"},
        ])
        # Malformed version string — the drift check should not
        # raise, just silently skip the comparison.
        _write_csv_with_header_version(path, "not-a-version", body)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            VariantCollection.from_csv(path)
    finally:
        os.unlink(path)

    drift = [w for w in caught if "major versions" in str(w.message)]
    assert len(drift) == 0


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
