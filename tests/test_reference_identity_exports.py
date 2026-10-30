"""Exports pin assembly and annotation, including heterogeneous and custom data."""

import gzip
import json
import pickle

import pytest
from pyensembl import cached_release

from varcode import EffectCollection, Genome, Variant, VariantCollection, StructuralVariant, to_proteome
from varcode.reference_metadata import REFERENCE_COLUMNS, reference_metadata
from .native_reference_fixture import CASES
from . import native_reference_fixture
from .test_native_reference_consequences import variant_for

native_genomes = native_reference_fixture.native_genomes


@pytest.mark.parametrize("structural", [False, True])
def test_archives_preserve_input_alias_without_re_resolving_it(structural):
    original = (StructuralVariant("chr7", 100, "DEL", end=105, genome="hg19:2013_09")
                if structural else Variant("chr7", 100, "A", "T", genome="hg19:2013_09"))
    restored = type(original).from_json(original.to_json())
    assert restored.genome.release == "2013_09"
    assert restored.original_reference_name == "hg19"
    assert restored.original_genome == "hg19:2013_09"
    assert restored.original_genome_was_ucsc
    assert restored.to_json() == original.to_json()


@pytest.mark.parametrize("header", [True, False])
def test_mixed_ensembl_releases_survive_variant_csv(tmp_path, header):
    variants = VariantCollection([
        Variant("7", 117531115, "G", "A", cached_release(release))
        for release in (81, 93)])
    path = tmp_path / "mixed.csv"
    if header:
        variants.to_csv(path)
    else:
        variants.to_dataframe().to_csv(path, index=False)
    restored = VariantCollection.from_csv(path)
    assert [v.genome.release for v in restored] == [81, 93]
    with pytest.raises(ValueError, match="annotation_version"):
        VariantCollection.from_csv(path, genome=cached_release(81))


@pytest.mark.parametrize("collection", [VariantCollection, EffectCollection])
@pytest.mark.parametrize("name", ["mouse", "dog", "chm13"])
def test_native_csv_records_and_validates_every_identity_field(native_genomes, tmp_path, collection, name):
    case = next(c for c in CASES if c["genome"] == name and "protein_expectation" not in c)
    genome = native_genomes[name]
    variant = variant_for(case, native_genomes)
    original = (collection([variant]) if collection is VariantCollection else
                collection([variant.effect_on_transcript(genome.transcript_by_id(case["transcript"]))]))
    path = tmp_path / "native.csv"
    original.to_csv(path)
    data = original.to_dataframe()
    assert all(column in data for column in REFERENCE_COLUMNS)
    identity = reference_metadata(genome)
    assert all(str(data.iloc[0][key]) == str(value) for key, value in identity.items())
    restored = collection.from_csv(path, genome=genome)
    restored_genome = restored[0].genome if collection is VariantCollection else restored[0].variant.genome
    assert reference_metadata(restored_genome) == identity
    # Supplying the wrong annotation of the same assembly must not reannotate silently.
    metadata = genome.to_dict()
    wrong = type(genome)(**dict(metadata, annotation_version="wrong-release"))
    with pytest.raises(ValueError, match="annotation_version"):
        collection.from_csv(path, genome=wrong)
    if name == "chm13" and collection is VariantCollection:
        with pytest.raises(ValueError, match="custom annotation"):
            collection.from_csv(path)
    contents = path.read_text()
    marker = "# annotation_version=" + identity["annotation_version"]
    assert marker in contents
    path.write_text(contents.replace(marker, "# annotation_version=conflicting-header"))
    with pytest.raises(ValueError, match="disagrees with header"):
        collection.from_csv(path, genome=genome)


@pytest.mark.parametrize("structural", [False, True])
def test_native_json_and_pickle_keep_exact_genome(native_genomes, structural):
    genome = native_genomes["chm13"]
    original = (StructuralVariant("22", 21798486, "DEL", end=21798490, genome=genome)
                if structural else Variant("22", 21798486, "T", "C", genome))
    for restored in (type(original).from_json(original.to_json()), pickle.loads(pickle.dumps(original))):
        assert reference_metadata(restored.genome) == reference_metadata(genome)
        assert restored.genome.reference_metadata == genome.reference_metadata


def test_wrapper_archive_retains_explicit_chromosome_dna(native_genomes, tmp_path):
    path = tmp_path / "dna.fa"
    path.write_text(">22\nACGT\n")
    original = Genome(native_genomes["chm13"], fasta=path, verify=False)
    try:
        restored = Genome.from_json(original.to_json())
        try:
            assert restored.sequence("22", 1, 4) == "ACGT"
            assert reference_metadata(restored) == reference_metadata(original)
        finally:
            restored.close()
    finally:
        original.close()
    with pytest.raises(TypeError, match="FASTA"):
        native_genomes["chm13_dna"].to_json()


def test_proteome_json_retains_native_reference_context(native_genomes, tmp_path):
    import orfdata
    case = next(c for c in CASES if c["id"].startswith("1KGP-"))
    exported = to_proteome(variant_for(case, native_genomes).effects())
    path = tmp_path / "native-proteome.json.gz"
    exported.save(path)
    encoded = gzip.decompress(path.read_bytes()).decode()
    assert "GCA_009914755.4" in encoded
    assert "2022_07" in encoded
    assert "ensembl-rapid" in encoded
    assert json.loads(encoded)
    assert orfdata.Proteome.load(path).identity == exported.identity
