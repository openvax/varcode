"""Public mammalian deletions and a separately pinned dog mitochondrial genome."""

import json

import pytest

from varcode import StructuralVariant, Variant, iter_candidate_effects
from varcode.effects import StartLoss
from .native_reference_fixture import DATA
from . import native_reference_fixture

native_genomes = native_reference_fixture.native_genomes
MITO = json.loads((DATA / "dog-mito.json").read_text())
STRUCTURAL = json.loads((DATA / "structural-cases.json").read_text())


@pytest.fixture(scope="module")
def dog_mito(tmp_path_factory):
    from orfdata._pinned_genome import PinnedEnsemblGenome

    path = tmp_path_factory.mktemp("dog-mito")
    rows, cdna, pep = [], [], []
    for record in MITO["records"]:
        a, b, strand = record["start"], record["end"], record["strand"]
        attributes = (
            'gene_id "{gene}"; gene_name "{gene}"; transcript_id "{id}"; '
            'transcript_name "{gene}"; transcript_biotype "protein_coding"; '
            'protein_id "{protein_id}";'
        ).format(**record)
        spans = [("gene", a, b), ("transcript", a, b), ("exon", a, b)]
        spans += (
            [("CDS", a, b - 3), ("start_codon", a, a + 2), ("stop_codon", b - 2, b)]
            if strand == "+"
            else [
                ("CDS", a + 3, b),
                ("start_codon", b - 2, b),
                ("stop_codon", a, a + 2),
            ]
        )
        for feature, start, end in spans:
            extra = (
                ' exon_id "%s-exon"; exon_number "1";' % record["id"]
                if feature == "exon"
                else ""
            )
            rows.append(
                "MT\trefseq\t%s\t%d\t%d\t.\t%s\t0\t%s%s\n"
                % (feature, start, end, strand, attributes, extra)
            )
        cdna.append(">%s\n%s\n" % (record["id"], record["sequence"]))
        pep.append(">%s\n%s\n" % (record["protein_id"], record["protein"]))
    for filename, text in [
        ("reference.gtf", "".join(rows)),
        ("rna.fa", "".join(cdna)),
        ("protein.fa", "".join(pep)),
        ("dna.fa", ">MT\n" + MITO["sequence"] + "\n"),
    ]:
        (path / filename).write_text(text)
    genome = PinnedEnsemblGenome(
        **{
            key: MITO[key]
            for key in (
                "reference_name",
                "annotation_name",
                "annotation_version",
                "species",
            )
        },
        gtf_path_or_url=str(path / "reference.gtf"),
        transcript_fasta_paths_or_urls=[str(path / "rna.fa")],
        protein_fasta_paths_or_urls=[str(path / "protein.fa")],
        genome_fasta_path_or_url=str(path / "dna.fa"),
        cache_directory_path=str(path),
    )
    genome.index()
    yield genome
    genome.close()


@pytest.mark.parametrize("case", MITO["cases"], ids=lambda c: c["id"])
@pytest.mark.parametrize("annotator", ["fast", "transcript_model"])
def test_public_dog_mitochondrial_variants(dog_mito, case, annotator):
    assert dog_mito.reference_name == "NC_002008.4"  # never relabeled ROS_Cfam
    assert dog_mito.reference_metadata["translation_tables"]["mitochondrial"] == 2
    assert dog_mito.sequence("MT", case["position"], case["position"]) == case["ref"]
    transcript = dog_mito.transcript_by_id(case["transcript"])
    effect = Variant(
        "MT", case["position"], case["ref"], case["alt"], dog_mito
    ).effect_on_transcript(transcript, annotator=annotator)
    assert effect.mutant_protein_sequence == case["expected_protein"]
    if case["id"] == "dog-mito-ND6":
        # The paper labels this I106V; NC_002008.4 places it at I107V.
        assert transcript.strand == "-"
        assert transcript.protein_sequence[106] == "I"
        assert case["expected_protein"][106] == "V"
    if case["id"] == "dog-mito-ND5":
        assert transcript.sequence[:3] == "ATA"
        assert effect.mutant_protein_sequence.startswith("M")


@pytest.mark.parametrize("case", STRUCTURAL, ids=lambda c: c["id"])
def test_public_native_mouse_and_dog_deletions(native_genomes, case):
    genome = native_genomes[case["genome"]]
    transcript = genome.transcript_by_id(case["transcript"])
    assert transcript.protein_sequence == case["reference_protein"]
    variant = StructuralVariant(
        case["contig"], case["start"], "DEL", end=case["end"], genome=genome
    )
    candidates = list(iter_candidate_effects(variant.effect_on_transcript(transcript)))
    matching = [
        e
        for e in candidates
        if e.mutant_transcript is not None
        and e.mutant_transcript.cdna_sequence == case["expected_cdna"]
    ]
    assert matching
    if case["expectation"] == "start_loss":
        assert all(
            case["start"] <= p <= case["end"] for p in transcript.start_codon_positions
        )
        assert any(
            isinstance(e, StartLoss) and e.mutant_protein_sequence is None
            for e in matching
        )
    else:
        assert len(case["expected_protein"]) == 237
        assert any(
            e.mutant_protein_sequence == case["expected_protein"] for e in matching
        )
