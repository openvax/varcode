"""Contig validation belongs to the annotation dataset, not its assembly name."""

import pytest
from pyensembl import Genome as AnnotationGenome, cached_release

from varcode import EffectCollection, Genome, Variant, VariantCollection, StructuralVariant, load_vcf
from varcode.effects import Failure


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
    assert len(effects) == 1
    assert isinstance(effects[0], Failure)
    assert effects[0].transcript is None
    assert "Invalid contig name '99'" in effects[0].error
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
    assert tp53.short_description == "p.E3del"
    assert tp53.mutant_protein_sequence == (
        tp53.original_protein_sequence[:2] + tp53.original_protein_sequence[3:])


@pytest.mark.parametrize("name,genome,expected", [
    ("chrUn_gl000216", 75, "GL000216.1"),
    ("chr1_gl000191_random", 75, "GL000191.1"),
    ("chrUn_KI270438v1", 81, "KI270438.1"),
    ("chrUn_GL000216v2", 81, "GL000216.2"),
    ("chr12_GL383550v2_alt", 81, "GL383550.2"),
    ("chr2_KN538362v1_fix", 81, "KN538362.1"),
    ("chrUn_gl000216", 81, "chrUn_gl000216"),
    ("chr6_apd_hap1", 75, "chr6_apd_hap1"),
    ("chrUn_unknown", 81, "chrUn_unknown"),
])
def test_ucsc_scaffold_accessions_preserve_versions(name, genome, expected):
    variant = Variant(name, 150, "A", "T", genome=genome, convert_ucsc_contig_names=True)
    assert variant.contig == expected
    assert variant.original_contig == name
    assert Variant(name, 150, "A", "T", genome=genome, convert_ucsc_contig_names=False).contig == name
    assert Variant.from_json(variant.to_json()).contig == expected


def test_mapping_does_not_substitute_a_different_accession_version(make_genome):
    genome = make_genome("version-two-only", ["KI270438.2"])
    variant = Variant("chrUn_KI270438v1", 150, "A", "T", genome=genome, convert_ucsc_contig_names=True)
    with pytest.raises(ValueError, match="KI270438.1"):
        variant.effects()
    sv = StructuralVariant("1", 150, "BND", genome=genome, mate_contig="chrUn_KI270438v1",
                           mate_start=150, convert_ucsc_contig_names=True)
    assert sv.mate_contig == "KI270438.1"


@pytest.mark.parametrize("raise_on_error", [None, False])
def test_batch_missing_contigs_warn_once_and_keep_failures(make_genome, raise_on_error):
    genome = make_genome("batch", ["1"])
    variants = VariantCollection([Variant(c, 150, "A", "T", genome=genome) for c in ("1", "98", "99")])
    with pytest.warns(UserWarning, match="2 variants on contigs absent") as warnings:
        effects = variants.effects(raise_on_error=raise_on_error)
    assert len(warnings) == 1
    assert len(effects) == 3
    assert sum(isinstance(e, Failure) for e in effects) == 2
    with pytest.warns(UserWarning, match="2 variants on contigs absent"):
        table = variants.to_dataframe(raise_on_error=raise_on_error)
    assert list(table.annotation_status) == ["ok", "Failure", "Failure"]
    assert table.iloc[1:].gene_id.isna().all()
    for operation in (variants.effects, variants.to_dataframe):
        with pytest.raises(ValueError, match="Invalid contig"):
            operation(raise_on_error=True)


def test_strelka_unplaced_contigs_keep_every_record():
    from pathlib import Path
    variants = load_vcf(str(Path(__file__).parent / "data/strelka-example.vcf"), genome=75)
    assert len(variants) == 20
    assert "GL000216.1" in {v.contig for v in variants}
    with pytest.warns(UserWarning, match="GL000226.1"):
        effects = variants.effects()
    assert {e.variant for e in effects} == set(variants)
    failures = [e for e in effects if isinstance(e, Failure)]
    assert len(failures) == 1
    assert failures[0].variant.contig == "GL000226.1"
    with pytest.warns(UserWarning, match="GL000226.1"):
        table = variants.to_dataframe()
    assert len(table) == 20
    assert (table.annotation_status == "Failure").sum() == 1


def test_missing_contig_policy_does_not_suppress_other_reference_errors(make_genome, monkeypatch):
    genome = make_genome("unavailable", ["1"])
    variants = VariantCollection([Variant("1", 150, "A", "T", genome=genome)])

    def unavailable():
        raise OSError("annotation database unavailable")

    monkeypatch.setattr(genome, "contigs", unavailable)
    for operation in (variants.effects, variants.to_dataframe):
        with pytest.raises(OSError, match="database unavailable"):
            operation()
    effects = variants.effects(raise_on_error=False)
    assert isinstance(effects[0], Failure)
    assert "database unavailable" in effects[0].error
    table = variants.to_dataframe(raise_on_error=False)
    assert table.iloc[0].annotation_status == "Failure"


def test_contig_normalization_does_not_open_reference_data(monkeypatch):
    genome = cached_release(75)

    def unexpected_access():
        raise AssertionError("construction must not read the annotation")

    monkeypatch.setattr(genome, "contigs", unexpected_access)
    variant = Variant("chrUn_gl000216", 150, "A", "T", genome=genome, convert_ucsc_contig_names=True)
    assert variant.contig == "GL000216.1"
