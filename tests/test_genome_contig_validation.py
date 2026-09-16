"""Contig validation belongs to the annotation dataset, not its assembly name."""

import pytest
from pyensembl import Genome as AnnotationGenome, cached_release

from varcode import EffectCollection, Genome, Variant


@pytest.fixture
def make_genome(tmp_path):
    def make(name, contigs, version=1):
        gtf = tmp_path / (name + ".gtf")
        rows = []
        for contig in contigs:
            attributes = (
                'gene_id "g%s"; gene_name "G%s"; gene_biotype "lncRNA"; '
                'transcript_id "t%s"; transcript_name "T%s"; '
                'transcript_biotype "lncRNA"; exon_id "e%s"; exon_number "1";'
            ) % ((contig,) * 5)
            for feature in ("gene", "transcript", "exon"):
                rows.append(
                    "%s\ttest\t%s\t100\t200\t.\t+\t.\t%s\n"
                    % (contig, feature, attributes))
        gtf.write_text("".join(rows))
        genome = AnnotationGenome(
            reference_name="GRCh38",
            annotation_name="test",
            annotation_version=version,
            gtf_path_or_url=str(gtf),
            cache_directory_path=str(tmp_path / name))
        genome.index()
        return genome

    return make


@pytest.mark.parametrize("subset_first", [True, False])
@pytest.mark.parametrize("full_version", [1, 2])
@pytest.mark.parametrize("wrapped", [False, True])
@pytest.mark.parametrize("attribute", ["genes", "gene_ids", "gene_names", "transcripts"])
def test_contigs_are_dataset_specific(
        make_genome, subset_first, full_version, wrapped, attribute):
    subset = make_genome("subset", ["1"])
    full = make_genome("full", ["1", "17"], version=full_version)
    if wrapped:
        subset, full = Genome(subset), Genome(full)
    datasets = [(subset, "1"), (full, "17")]
    if not subset_first:
        datasets.reverse()

    # Fresh variants on each pass also exercise repeated dataset queries.
    for _ in range(2):
        for genome, contig in datasets:
            variant = Variant(contig, 150, "A", "T", genome=genome)
            assert variant.genome is genome
            assert len(getattr(variant, attribute)) == 1
        for genome in (subset, full):
            variant = Variant("99", 150, "A", "T", genome=genome)
            with pytest.raises(ValueError, match="Invalid contig name '99'"):
                getattr(variant, attribute)
        # A full-dataset query must not make chr17 valid in the subset.
        variant = Variant("17", 150, "A", "T", genome=subset)
        with pytest.raises(ValueError, match="Invalid contig name '17'"):
            getattr(variant, attribute)


def test_invalid_contig_effects_preserve_collection_contract(make_genome, caplog):
    variant = Variant("99", 150, "A", "T", genome=make_genome("subset", ["1"]))
    with pytest.raises(ValueError, match="Invalid contig name '99'"):
        variant.effects(raise_on_error=True)

    effects = variant.effects(raise_on_error=False, annotator="fast")
    assert isinstance(effects, EffectCollection)
    assert len(effects) == 0
    assert effects.annotator == "fast"
    assert effects.annotator_version
    assert effects.annotated_at
    assert isinstance(effects.clone_with_new_elements([]), EffectCollection)
    assert "Invalid contig name '99'" in caplog.text


def test_subset_then_full_grch38_tp53_annotation(make_genome):
    subset = make_genome("subset", ["1"])
    assert Variant("1", 150, "A", "T", genome=subset).gene_ids == ["g1"]

    # The same full Ensembl release used elsewhere in the test suite.
    genome = cached_release(81)
    # The coding deletion from isovar's original downstream regression.
    variant = Variant("17", 7676589, "CTC", "", genome=genome)
    effects = variant.effects()
    tp53 = next(effect for effect in effects if effect.transcript_id == "ENST00000269305")
    assert tp53.gene_name == "TP53"
    assert tp53.short_description == "p.E2del"
    assert tp53.mutant_protein_sequence == (
        tp53.original_protein_sequence[:1] + tp53.original_protein_sequence[2:])
