"""New starts are conditional products, with exact ancestry and explicit bounds."""

from dataclasses import replace
import json
from pathlib import Path
import pickle

from Bio.Seq import Seq
import pytest
from pyensembl import cached_release

from varcode import (
    EffectCollection, GermlineContext, ORFDiscovery, ORFDiscoveryPolicy, Variant,
    discover_created_orfs, discover_created_orfs_in_products, iter_candidate_effects, to_proteome,
)
from varcode.effect_hypotheses import RealizedTranscriptProduct
from varcode.effects import StartGain
from varcode.orf_serialization import _Feature, _ReferenceGenome, _ReferenceTranscript


def host(sequence, reverse=False, table=1):
    genome = _ReferenceGenome(dict(reference_name="fixture", species="homo_sapiens", taxon_id=9606,
                                  annotation_name="synthetic", annotation_version="1"))
    return _ReferenceTranscript(
        dict(id="noncoding-host", name="NCRNA", contig="1", strand="-" if reverse else "+",
             start=1, end=len(sequence), sequence=sequence, protein_sequence=None,
             biotype="lncRNA", is_protein_coding=False, complete=False,
             contains_start_codon=False, contains_stop_codon=False,
             on_backward_strand=reverse, on_forward_strand=not reverse, translation_table=table),
        genome, (_Feature(1, len(sequence)),), _Feature(1, len(sequence), "gene", "GENE"))


def product(tx, sequence, origins=None, **evidence):
    if origins is None:
        positions = range(1, len(sequence) + 1)
        if tx.on_backward_strand:
            positions = reversed(tuple(positions))
        origins = tuple(("1", p, tx.strand) for p in positions)
    return RealizedTranscriptProduct(tx, sequence, None, evidence=dict(base_origins=tuple(origins), **evidence))


def discover(before="AGGAAATAA", after="ATGAAATAA", *, reverse=False, table=1, policy=None):
    tx = host(before, reverse, table)
    return discover_created_orfs_in_products(product(tx, before), product(tx, after), policy=policy)


@pytest.fixture(scope="module")
def cdkn2a():
    fixture = json.loads((Path(__file__).parent / "data/created_orfs/CDKN2A.json").read_text())
    native = fixture["native"]
    genome = cached_release(native["annotation_version"])
    transcript = genome.transcript_by_id(native["transcript_id"])
    allele = fixture["genomic"]
    variant = Variant(allele["contig"], allele["start"], allele["ref"], allele["alt"], genome=genome)
    return fixture, transcript, variant


def test_published_cdkn2a_start_gain_has_independent_refseq_peptide(cdkn2a):
    fixture, transcript, variant = cdkn2a
    expected, native = fixture["refseq"], fixture["native"]
    sequence = expected["sequence"]
    offset = expected["variant_offset"]
    assert sequence[offset] == "G" and offset == expected["cds_start"] - 34
    modified = sequence[:offset] + "T" + sequence[offset + 1:]
    coding = modified[expected["new_start"]:expected["stop_offset"] + 3]
    assert str(Seq(coding).translate(to_stop=True)) == expected["expected_protein"]
    # Different 5' transcript extents retain the exact same new coding path.
    start, end = native["new_start"], native["stop_offset"] + 3
    native_modified = transcript.sequence[:start + 1] + "T" + transcript.sequence[start + 2:]
    assert native_modified[start:end] == coding
    ordinary = variant.effect_on_transcript(transcript)
    report = discover_created_orfs(ordinary)
    assert report.complete
    gained, = report.effects
    assert type(ordinary).__name__ == "FivePrimeUTR"
    assert gained.discovered_orf.start == start and gained.discovered_orf.stop_offset == native["stop_offset"]
    assert gained.discovered_orf.baseline_codons == ("AGG",)
    assert gained.mutant_protein_sequence == expected["expected_protein"]
    assert gained.original_protein_sequence is None and gained.modifies_protein_sequence is None
    assert gained.orf is None and gained.reference_transcript is transcript
    assert gained.transcript_id == transcript.id and gained.gene_id == transcript.gene_id
    assert gained.source_evidence["mutant_translation_evidence"] == "not_assessed"
    assert gained.source_evidence["downstream_cds_regulation"] == "not_assessed"
    assert list(iter_candidate_effects(gained)) == [gained]


def test_cdkn2a_archives_csv_and_proteome_preserve_discovery(cdkn2a, tmp_path):
    fixture, transcript, variant = cdkn2a
    report = discover_created_orfs(variant.effect_on_transcript(transcript))
    for restored in (ORFDiscovery.from_json(report.to_json()), pickle.loads(pickle.dumps(report))):
        gained, = restored.effects
        assert restored.complete
        assert gained.discovery_policy == report.policy
        assert gained.discovered_orf == report.effects[0].discovered_orf
        assert gained.baseline_transcript.cdna_sequence == transcript.sequence
        assert gained.hypothesis == report.effects[0].hypothesis
    effects = EffectCollection(report.effects)
    path = tmp_path / "created.csv"
    effects.to_csv(path)
    restored, = EffectCollection.from_csv(path)
    assert isinstance(restored, StartGain)
    assert restored.discovered_orf == report.effects[0].discovered_orf
    text = path.read_text()
    path.write_text(text.replace("created-orf:", "edited-orf:", 1))
    with pytest.raises(ValueError, match="discovered_orf_id disagrees"):
        EffectCollection.from_csv(path)
    proteome = to_proteome(effects)
    entry, = proteome
    assert entry.sequence == fixture["refseq"]["expected_protein"]
    assert entry.producer.origin == "nc_orf" and entry.producer.orf_kind == "variant_created"
    assert entry.producer.start_complete and entry.producer.stop_complete
    assert entry.producer.reference.annotation_version == "95"
    assert entry.producer.tool == "varcode"
    proteome.save(tmp_path / "products.json.gz")
    from orfdata import Proteome
    assert Proteome.load(tmp_path / "products.json.gz").identity == proteome.identity


@pytest.mark.parametrize("reverse", [False, True])
def test_noncoding_host_and_junction_spanning_start(reverse):
    tx = host("AGGAAATAA", reverse)
    # A mapped start need not occupy three contiguous genomic coordinates.
    positions = (3, 2, 100, 99, 98, 97, 96, 95, 94) if reverse else (1, 2, 100, 101, 102, 103, 104, 105, 106)
    origins = tuple(("1", p, tx.strand) for p in positions)
    result = discover_created_orfs_in_products(product(tx, "AGGAAATAA", origins), product(tx, "ATGAAATAA", origins))
    effect, = result.effects
    assert result.status == "complete"
    assert effect.mutant_protein_sequence == "MK"
    assert effect.discovered_orf.base_origins[:3] == origins[:3]
    assert effect.discovered_orf.initiation_context["upstream_bases_available"] == 0
    assert effect.discovered_orf.initiation_context["biological_five_prime_completeness"] == "unknown"


def test_insertion_coordinate_shift_does_not_create_existing_start():
    tx = host("CCATGAAATAA")
    before = product(tx, tx.sequence)
    after = product(tx, "CCCCATGAAATAA", (None, None, *before.evidence["base_origins"]), inserted_base_offsets=(0, 1))
    result = discover_created_orfs_in_products(before, after)
    assert result.status == "complete" and not result.effects


def test_identical_rna_with_different_source_placement_has_no_new_start():
    tx = host("ATGAAATAA")
    before = product(tx, tx.sequence)
    after = product(tx, tx.sequence, tuple(("1", p, "+") for p in range(101, 110)))
    result = discover_created_orfs_in_products(before, after)
    assert result.status == "complete" and not result.effects


@pytest.mark.parametrize("insertion", [False, True])
def test_indels_create_new_start_adjacency(insertion):
    tx = host("AGAAATAA" if insertion else "ACCTGAAATAA")
    before = product(tx, tx.sequence)
    origins = before.evidence["base_origins"]
    new = (origins[0], None, *origins[1:]) if insertion else (origins[0], *origins[3:])
    after = product(tx, "ATGAAATAA", new, inserted_base_offsets=(1,) if insertion else ())
    result = discover_created_orfs_in_products(before, after)
    gained, = result.effects
    assert gained.mutant_protein_sequence == "MK"
    assert gained.discovered_orf.creation == "new_rna_adjacency_or_sequence"


def test_unchanged_duplicated_start_is_not_new_and_new_junction_can_be():
    tx = host("ATGAAATAA")
    before = product(tx, tx.sequence)
    after = product(tx, tx.sequence * 2, before.evidence["base_origins"] * 2)
    assert not discover_created_orfs_in_products(before, after).effects


def test_duplicated_baseline_origins_with_different_codons_are_ambiguous():
    tx = host("ATGAAATAAAGGAAATAA")
    origins = tuple(("1", p, "+") for p in range(1, 10))
    before = product(tx, tx.sequence, origins * 2)
    after = product(tx, "ATGAAATAAATGAAATAA", origins * 2)
    result = discover_created_orfs_in_products(before, after)
    assert result.status == "partial" and "ambiguous_baseline_ancestry" in result.reasons
    assert not result.effects


def test_policy_controls_non_aug_candidates():
    assert not discover("TTGAAATAA", "CTGAAATAA").effects
    result = discover("TTGAAATAA", "CTGAAATAA", policy=ORFDiscoveryPolicy(start_codons=("ATG", "CTG")))
    gained, = result.effects
    assert gained.mutant_protein_sequence == "MK" and gained.discovered_orf.start_codon == "CTG"


@pytest.mark.parametrize("table,expected", [(1, "M"), (2, "MW")])
def test_reference_genetic_code_controls_new_orf_translation(table, expected):
    gained, = discover("AGGTGATAA", "ATGTGATAA", table=table).effects
    assert gained.mutant_protein_sequence == expected
    assert gained.discovered_orf.translation_table == table


@pytest.mark.parametrize("after,termination,protein", [
    ("ATGAAA", "sequence_boundary", "MK"),
    ("ATGNNNTAA", "ambiguous_sequence", "M"),
])
def test_partial_peptides_keep_their_limit(after, termination, protein):
    gained, = discover("A" * len(after), after).effects
    assert gained.mutant_protein_sequence == protein
    assert gained.discovered_orf.termination == termination
    assert gained.discovered_orf.stop_offset is None
    entry, = to_proteome([gained])
    assert entry.producer.stop_complete is False


def test_bounds_do_not_turn_incomplete_search_into_no_orfs():
    result = discover(policy=ORFDiscoveryPolicy(max_scan_bases=5))
    assert result.status == "bounded" and result.reasons == ("baseline_scan_limit",)
    result = discover("AGGAAATAAAGGAAATAA", "ATGAAATAAATGAAATAA", policy=ORFDiscoveryPolicy(max_candidates=1))
    assert result.status == "bounded" and result.reasons == ("candidate_limit",)
    assert len(result.effects) == 1
    gained, = discover(policy=ORFDiscoveryPolicy(max_protein_length=1)).effects
    assert gained.mutant_protein_sequence == "M" and gained.discovered_orf.termination == "protein_length_limit"
    tx = host("AGGAAATAA")
    before = product(tx, tx.sequence)
    after = product(tx, "ATGAAATAACCC", (*before.evidence["base_origins"], None, None, None),
                    inserted_base_offsets=(9, 10, 11))
    result = discover_created_orfs_in_products(before, after, policy=ORFDiscoveryPolicy(max_scan_bases=9))
    assert result.status == "bounded" and result.reasons == ("mutant_scan_limit",)
    assert len(result.effects) == 1


@pytest.mark.parametrize("kind,reason", [
    ("missing_rna", "paired_rna_sequence_unavailable"),
    ("missing_origins", "base_origins_unavailable"),
    ("partial", "partial_rna_comparison"),
])
def test_missing_and_partial_products_are_unavailable(kind, reason):
    tx = host("AGGAAATAA")
    before, after = product(tx, tx.sequence), product(tx, "ATGAAATAA")
    if kind == "missing_rna":
        before = replace(before, cdna_sequence=None)
    elif kind == "missing_origins":
        before = replace(before, evidence={})
    else:
        before = replace(before, evidence={**before.evidence, "partial": True})
    result = discover_created_orfs_in_products(before, after)
    assert result.status == "unavailable" and result.reasons == (reason,) and not result.effects


def test_unknown_insertion_ancestry_and_ambiguous_baseline_do_not_prove_gain():
    tx = host("AGGAAATAA")
    before = product(tx, tx.sequence)
    origins = list(before.evidence["base_origins"])
    origins[1] = None
    result = discover_created_orfs_in_products(before, product(tx, "ATGAAATAA", origins))
    assert not result.effects and "ambiguous_start_ancestry" in result.reasons
    result = discover("ANGAAATAA", "ATGAAATAA")
    assert not result.effects and "ambiguous_baseline_codon" in result.reasons


def test_malformed_origins_and_reference_mismatch_raise():
    tx = host("AGGAAATAA")
    before, after = product(tx, tx.sequence), product(tx, "ATGAAATAA")
    with pytest.raises(ValueError, match="every RNA nucleotide"):
        discover_created_orfs_in_products(before, replace(after, evidence={"base_origins": ()}))
    other = host("AGGAAATAA")
    other.genome.metadata["annotation_version"] = "2"
    with pytest.raises(ValueError, match="same host transcript and reference"):
        discover_created_orfs_in_products(before, product(other, "ATGAAATAA"))


def test_phase_alternatives_use_actual_baseline_and_preserve_hypotheses(cdkn2a):
    _, tx, variant = cdkn2a
    # G>T in the same start triplet abolishes the gained AUG if inherited in cis.
    germline = Variant("9", 21974860, "C", "A", genome=variant.genome)
    root = variant.effect_on_transcript(tx, germline=GermlineContext.from_variants([germline]))
    report = discover_created_orfs(root)
    assert len(report.assessments) == 2 and report.complete
    assert sorted(len(item.effects) for item in report.assessments) == [0, 1]
    gained, = report.effects
    assert any(state == "trans" for _, state in gained.hypothesis.phase)
    limited = discover_created_orfs(root, policy=ORFDiscoveryPolicy(max_products=1))
    assert not limited.complete and limited.assessments[-1].reasons == ("product_limit",)


def test_distinct_phase_baselines_with_same_peptide_are_not_deduplicated(cdkn2a):
    _, tx, variant = cdkn2a
    positions = tuple(p for exon in sorted(tx.exons, key=lambda exon: exon.start, reverse=True)
                      for p in range(exon.end, exon.start - 1, -1))
    ref = str(Seq(tx.sequence[500]).complement())
    alt = next(base for base in "ACGT" if base != ref)
    germline = Variant("9", positions[500], ref, alt, genome=variant.genome)
    report = discover_created_orfs(variant.effect_on_transcript(
        tx, germline=GermlineContext.from_variants([germline])))
    assert len(report.effects) == 2
    assert len({e.baseline_transcript.cdna_sequence for e in report.effects}) == 2
    assert len({e.mutant_protein_sequence for e in report.effects}) == 1
    assert len(EffectCollection(report.effects)) == 2


def test_selected_reference_orf_keeps_original_host_for_discovery(cdkn2a):
    from orfdata import TranscriptORF
    from varcode import annotate_orf_variants
    _, tx, variant = cdkn2a
    reference = TranscriptORF.from_annotation(tx)
    report = discover_created_orfs(annotate_orf_variants([variant], reference).effect)
    gained, = report.effects
    assert gained.transcript is tx and gained.reference_transcript is tx
    assert gained.discovered_orf.transcript_id == tx.id
    assert gained.orf is None


def test_legacy_annotation_does_not_silently_claim_no_new_orfs(cdkn2a):
    _, transcript, variant = cdkn2a
    report = discover_created_orfs(variant.effect_on_transcript(transcript, annotator="fast"))
    assert not report.complete and not report.effects
    assert report.assessments[0].reasons == ("paired_products_unavailable",)


def test_mapped_rna_without_focal_variant_retains_reference_in_csv(tmp_path):
    gained, = discover().effects
    collection = EffectCollection([gained])
    path = tmp_path / "rna.csv"
    collection.to_csv(path)
    restored, = EffectCollection.from_csv(path)
    assert restored.variant is None
    assert restored.discovered_orf == gained.discovered_orf
    assert restored.reference_transcript.genome.reference_name == "fixture"


@pytest.mark.parametrize("kwargs", [dict(max_candidates=0), dict(max_products=True),
                                     dict(start_codons=()), dict(start_codons=("AUG",))])
def test_invalid_policy_rejected(kwargs):
    with pytest.raises(ValueError):
        ORFDiscoveryPolicy(**kwargs)
