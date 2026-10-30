"""Public splice alleles, observed junctions and chromosome-dependent hypotheses."""

from types import SimpleNamespace

import pytest
from Bio.Seq import Seq

from varcode import (
    MutantTranscript, RNAEvidence, StructuralVariant, Variant, iter_candidate_effects, make_rna_outcome,
)
from varcode.effects import Unresolved
from varcode.genome_sequence import reference_range
from .native_reference_fixture import CASES
from . import native_reference_fixture

native_genomes = native_reference_fixture.native_genomes


SPLICE_CASES = [c for c in CASES if c["id"].startswith("ROS-splice-") or c["id"] == "SweGen-DHRS4L2"]


@pytest.mark.parametrize("case", SPLICE_CASES, ids=lambda c: c["id"])
def test_native_intronic_dna_yields_exact_retained_intron_rna(native_genomes, case):
    genome = native_genomes[case["genome"] + "_dna"]
    transcript = genome.transcript_by_id(case["transcript"])
    variant = Variant(case["contig"], case["position"], case["ref"], case["alt"], genome)
    assert reference_range(genome, variant.contig, case["position"],
                           case["position"] + len(case["ref"]) - 1) == case["ref"]
    edit_position, ref, alt = case["position"], case["ref"], case["alt"]
    while ref and alt and ref[0] == alt[0]:
        edit_position, ref, alt = edit_position + 1, ref[1:], alt[1:]
    exons = sorted(transcript.exons, key=lambda e: e.start, reverse=transcript.strand == "-")
    rna_offset = 0
    for left, right in zip(exons, exons[1:]):
        rna_offset += len(left)
        start, end = min(left.end, right.end) + 1, max(left.start, right.start) - 1
        if start <= edit_position <= end:
            break
    else:
        pytest.fail("Public splice fixture does not lie in an annotated intron")
    dna = reference_range(genome, transcript.contig, start, end)
    assert len(dna) == end - start + 1
    offset = edit_position - start
    altered = dna[:offset] + alt + dna[offset + len(ref):]
    intron = str(Seq(altered).reverse_complement()) if transcript.strand == "-" else altered
    expected_rna = transcript.sequence[:rna_offset] + intron + transcript.sequence[rna_offset:]
    cds_start = min(transcript.start_codon_spliced_offsets)
    if rna_offset <= cds_start:
        cds_start += len(intron)
    coding = expected_rna[cds_start:]
    expected_protein = str(Seq(coding[:len(coding) // 3 * 3]).translate(to_stop=True))
    effect = variant.effect_on_transcript(transcript)
    candidates = list(iter_candidate_effects(effect))
    retained = [e for e in candidates if e.mutant_transcript is not None
                and e.mutant_transcript.cdna_sequence == expected_rna]
    assert retained
    assert any(e.mutant_protein_sequence == expected_protein for e in retained)
    assert any(e.mutant_protein_sequence == case["reference_protein"] for e in candidates)


def test_ap2s1_observed_del17_junction_survives_rna_evidence(native_genomes):
    case = next(c for c in CASES if c["id"] == "Ap2s1-del17")
    genome = native_genomes["mouse_dna"]
    transcript = genome.transcript_by_id(case["transcript"])
    variant = Variant(case["contig"], case["position"], case["ref"], case["alt"], genome)
    start = min(transcript.start_codon_spliced_offsets)
    # Paper RT-PCR/sequencing establishes the 51-nt deletion, not full-length RNA.
    before, after = start + 34 * 3, start + 51 * 3
    rna = transcript.sequence[:before] + transcript.sequence[after:]
    assert len(transcript.sequence) - len(rna) == 51
    coding = rna[start:]
    assert str(Seq(coding[:len(coding)//3*3]).translate(to_stop=True)) == case["expected_protein"]
    observed = Unresolved(variant, transcript, "published_splice_junction")
    observed.mutant_transcript = MutantTranscript.from_sequence(
        rna, reference_transcript=transcript, mutant_protein_sequence=case["expected_protein"],
        annotator_name="published_junction_plus_reference",
        evidence={"reference_completion": True, "source": case["source"],
                  "observed_deleted_cdna_interval": (before, after)})
    results = variant.effects(
        rna_resolver=RNAEvidence([make_rna_outcome(observed, source="PMID:29479578")]),
        transcript_filter=lambda t: t.id == transcript.id)
    candidates = [e for root in results for e in iter_candidate_effects(root)
                  if e.mutant_protein_sequence == case["expected_protein"]]
    assert candidates
    # Reference completion does not manufacture a fully observed patient baseline.
    assert any(isinstance(e, Unresolved) and e.mutant_transcript.cdna_sequence == rna for e in candidates)


def test_chm13_structural_deletion_matches_independent_real_sequence(native_genomes):
    # Controlled rearrangement on real duplicated/newly resolved DNA, not a patient call.
    genome = native_genomes["chm13_dna"]
    transcript = genome.transcript_by_id("ENST05220004619")
    positions = [p for e in sorted(transcript.exons, key=lambda e: e.start)
                 for p in range(e.start, e.end + 1)]
    cds_start = min(transcript.start_codon_spliced_offsets)
    offset = cds_start + 30
    a, b = positions[offset], positions[offset + 8]
    assert b - a == 8
    variant = StructuralVariant(transcript.contig, a, "DEL", end=b, genome=genome)
    rna = transcript.sequence[:offset] + transcript.sequence[offset + 9:]
    protein = transcript.protein_sequence[:10] + transcript.protein_sequence[13:]
    assert reference_range(genome, transcript.contig, a, b) == transcript.sequence[offset:offset+9]
    effects = variant.effects(transcript_filter=lambda t: t.id == transcript.id)
    assert any(e.mutant_protein_sequence == protein and e.mutant_transcript.cdna_sequence == rna
               for root in effects for e in iter_candidate_effects(root))
