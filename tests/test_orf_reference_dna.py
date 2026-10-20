"""Reference DNA reaches ORF consequences through the ordinary public APIs."""

from dataclasses import replace

import pytest
from pyensembl import Genome as AnnotationGenome

from varcode import Genome, Variant, VariantCollection, annotate_orf_variants
from varcode.genome_sequence import reference_range
from varcode.nucleotides import reverse_complement
from .test_uorf import Transcript

pytest.importorskip("orfdata")


@pytest.fixture(params=[False, True], ids=["plus", "minus"])
def reference_dna(request, tmp_path):
    model = Transcript(reverse=request.param, split=True)
    dna = list("A" * 500)
    for position, base in zip(model.positions, model.sequence):
        dna[position - 1] = reverse_complement(base) if request.param else base
    first, second = sorted(model.exons, key=lambda e: e.start)
    dna[first.end:first.end + 2] = "CT" if request.param else "GT"
    dna[second.start - 3:second.start - 1] = "AC" if request.param else "AG"
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">1\n" + "".join(dna) + "\n")
    cdna = tmp_path / "transcripts.fa"
    cdna.write_text(">test\n" + model.sequence + "\n")
    gtf = tmp_path / "reference.gtf"
    attrs = ('gene_id "gene"; gene_name "GENE"; transcript_id "test"; '
             'transcript_name "test"; transcript_biotype "lncRNA";')
    rows = [("gene", model.start, model.end, ""),
            ("transcript", model.start, model.end, "")]
    rows += [("exon", e.start, e.end, ' exon_id "e%d"; exon_number "%d";' % (i, i))
             for i, e in enumerate(model.exons, 1)]
    gtf.write_text("".join('1\ttest\t%s\t%d\t%d\t.\t%s\t.\t%s%s\n' %
                         (kind, start, end, model.strand, attrs, extra)
                         for kind, start, end, extra in rows))
    genomes = []

    def build(source):
        native = AnnotationGenome(
            reference_name="GRCh38", annotation_name="synthetic", annotation_version="1",
            gtf_path_or_url=str(gtf), transcript_fasta_paths_or_urls=[str(cdna)],
            genome_fasta_path_or_url=str(fasta) if source == "native" else None,
            cache_directory_path=str(tmp_path / source))
        native.species = "homo_sapiens"
        native.index()
        genomes.append(native)
        genome = Genome(native, fasta=fasta, verify=False) if source == "override" else native
        if genome is not native:
            genomes.append(genome)
        transcript = native.transcript_by_id("test")
        record = replace(model.record(), taxon_id=9606, species="homo_sapiens", source="test")
        return genome, record.bind(transcript)

    yield build, 150 if not request.param else 250
    for genome in reversed(genomes):
        genome.close()


@pytest.mark.parametrize("source", ["native", "override"])
def test_ordinary_orf_and_collection_calls_use_available_intronic_dna(reference_dna, source):
    build, position = reference_dna
    genome, reference = build(source)
    assert reference_range(genome, "1", position, position) == "A"
    variant = Variant("1", position, "A", "G", genome=genome)
    direct = annotate_orf_variants([variant], reference)
    automatic = variant.effects(extra_orfs=[reference]).orf_annotations[0]
    collection = VariantCollection([variant]).effects(extra_orfs=[reference]).orf_annotations[0]
    for annotation in (direct, automatic, collection):
        assert type(annotation.effect).__name__ == "Intronic"
        assert len(annotation.outcomes) == 1
        assert annotation.outcomes[0].mutant.protein_sequence == reference.protein_sequence


def test_missing_dna_stays_unresolved_and_explicit_provider_is_honored(reference_dna):
    build, position = reference_dna
    genome, reference = build("missing")
    variant = Variant("1", position, "A", "G", genome=genome)
    result = variant.effects(extra_orfs=[reference]).orf_annotations[0]
    assert type(result.effect).__name__ == "Unresolved"
    calls = []

    def provider(contig, start, end):
        calls.append((contig, start, end))
        return "A" * (end - start + 1)

    explicit = annotate_orf_variants([variant], reference, sequence_provider=provider)
    assert type(explicit.effect).__name__ == "Intronic"
    assert ("1", position, position) in calls


@pytest.mark.parametrize("source", ["native", "override"])
@pytest.mark.parametrize("kind", ["splice_snv", "local_sv", "cis_group"])
def test_splice_structural_and_cis_products_use_intronic_dna(reference_dna, source, kind):
    from varcode import StructuralVariant
    build, _ = reference_dna
    genome, reference = build(source)
    transcript = reference.transcript
    exons = sorted(transcript.exons, key=lambda e: e.start)
    reverse = transcript.on_backward_strand
    donor = exons[-1].start - 1 if reverse else exons[0].end + 1
    if kind == "local_sv":
        start, end = (donor - 5, donor) if reverse else (donor, donor + 5)
        variant = StructuralVariant("1", start, "DEL", end=end, genome=genome)
    else:
        ref = reference_range(genome, "1", donor, donor)
        variant = Variant("1", donor, ref, "T" if reverse else "A", genome=genome)
    if kind == "cis_group":
        partner_position = donor - 15 if reverse else donor + 15
        partner = Variant("1", partner_position, "A", "G", genome=genome)

        class Phase:
            def in_cis(self, first, second, transcript=None):
                return True

        annotations = VariantCollection([variant, partner]).effects(
            extra_orfs=[reference], phase_resolver=Phase()).orf_annotations
        annotation = next(a for a in annotations if len(a.variants) == 2)
    else:
        annotation = variant.effects(extra_orfs=[reference]).orf_annotations[0]
    retained = [
        outcome for outcome in annotation.outcomes
        if outcome.hypothesis is not None and outcome.hypothesis.splice_plan is not None
        and any(option.mechanism == "intron_retention"
                for _, option in outcome.hypothesis.splice_plan.choices)]
    assert retained
    assert any(outcome.mutant is not None and outcome.mutant.cdna_sequence is not None for outcome in retained)
