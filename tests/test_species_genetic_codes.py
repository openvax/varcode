"""Independent codon examples and NCBI species/compartment assignments."""

from types import SimpleNamespace

import pytest

from varcode.effects.codon_tables import codon_table_by_id, codon_table_for_transcript, translate_sequence


@pytest.mark.parametrize("species,contig,code", [
    ("human", "MT", 2), ("mus_musculus", "chrM", 2),
    ("canis_familiaris", "MT", 2), ("danio_rerio", "MT", 2),
    ("saccharomyces_cerevisiae", "Mito", 3),
    ("caenorhabditis_elegans", "MtDNA", 5),
    ("drosophila_melanogaster", "MT", 5),
    ("drosophila_melanogaster", "mitochondrion_genome", 5),
    ("mouse", "MT", 2),
    ("schizosaccharomyces_pombe", "MT", 4),
    ("candida_albicans", "1", 12),
    ("arabidopsis_thaliana", "Mt", 1),
    ("arabidopsis_thaliana", "Pt", 11),
    ("toxoplasma_gondii", "apicoplast", 4),
])
def test_ncbi_assignment(species, contig, code):
    transcript = SimpleNamespace(contig=contig, genome=SimpleNamespace(species=species))
    assert codon_table_for_transcript(transcript).ncbi_id == code


@pytest.mark.parametrize("code,sequence,protein", [
    (2, "ATGAGATGA", "M"),
    (3, "ATGCTATGA", "MTW"),
    (5, "ATGAGATGA", "MSW"),
    (4, "ATGATATGA", "MIW"),
    (12, "ATGCTGTAA", "MS"),
])
def test_meaning_of_distinguishing_codons(code, sequence, protein):
    assert translate_sequence(sequence, codon_table=codon_table_by_id(code)) == protein


def test_unknown_species_requires_explicit_assignment():
    genome = SimpleNamespace(species="unregistered_species")
    transcript = SimpleNamespace(genome=genome, contig="MT")
    with pytest.raises(ValueError, match="No mitochondrial genetic-code assignment"):
        codon_table_for_transcript(transcript)
    genome.mitochondrial_translation_table = 5
    assert codon_table_for_transcript(transcript).ncbi_id == 5
    transcript.translation_table = 3
    assert codon_table_for_transcript(transcript).ncbi_id == 3


def test_source_compartment_metadata_handles_accession_contigs():
    genome = SimpleNamespace(species="human", contig_compartments={"NC_012920.1": "mitochondrial"})
    assert codon_table_for_transcript(SimpleNamespace(genome=genome, contig="NC_012920.1")).ncbi_id == 2


@pytest.mark.parametrize("species,code", [("human", 2), ("yeast", 3)])
def test_mirrored_reference_retains_assembly_species_assignment(tmp_path, species, code):
    from pyensembl import EnsemblRelease
    genome = EnsemblRelease(115, species=species).with_mirror(tmp_path)
    assert codon_table_for_transcript(SimpleNamespace(genome=genome, contig="MT")).ncbi_id == code


def test_custom_assembly_without_species_does_not_guess_mitochondrial_code():
    genome = SimpleNamespace(reference_name="custom-GRCh38")
    with pytest.raises(ValueError, match="No mitochondrial genetic-code assignment"):
        codon_table_for_transcript(SimpleNamespace(genome=genome, contig="MT"))


@pytest.mark.parametrize("code", [27, 28, 31])
def test_context_dependent_stop_codes_are_explicitly_unsupported(code):
    with pytest.raises(ValueError, match="context-dependent stop"):
        codon_table_by_id(code)


def test_taxonomy_conflicts_and_unknown_codes_raise():
    transcript = SimpleNamespace(contig="MT", genome=SimpleNamespace(species="human", taxon_id=4932))
    with pytest.raises(ValueError, match="species and taxon_id disagree"):
        codon_table_for_transcript(transcript)
    for code in (999, True, "2"):
        with pytest.raises(ValueError, match="Unsupported explicit"):
            codon_table_by_id(code)


def test_yeast_orf_modification_retains_declared_code():
    from dataclasses import replace
    from pyensembl import cached_release
    from varcode import Variant, annotate_orf_variants
    from .test_uorf import Transcript
    pytest.importorskip("orfdata")
    transcript = Transcript()
    transcript.contig = "MT"
    transcript.sequence = "CCCATGCTATGGGGGTAAATGAAATAA"
    transcript.genome = cached_release(115, species="yeast")
    record = replace(
        transcript.record(), protein_sequence="MTWG", translation_table=3,
        contig="MT", taxon_id=4932, species="saccharomyces_cerevisiae",
        genome_build=transcript.genome.reference_name)
    reference = record.bind(transcript)
    variant = Variant("MT", 106, "C", "T", genome=transcript.genome)
    annotation = annotate_orf_variants([variant], reference)
    assert annotation.effect.aa_ref == "T"
    assert annotation.effect.aa_alt == "L"
    assert annotation.outcomes[0].mutant.protein_sequence == "MLWG"
