"""Public variants on native assemblies, checked against independent sequences.

See data/native_references/README.md for evidence levels and annotation differences.
"""

from types import SimpleNamespace

import pytest

from varcode import Variant, VariantCollection, iter_candidate_effects
from .native_reference_fixture import CASES, DATA
from . import native_reference_fixture

native_genomes = native_reference_fixture.native_genomes


def variant_for(case, genomes):
    return Variant(case["contig"], case["position"], case["ref"], case["alt"],
                   genome=genomes[case["genome"]])


CODING_CASES = [case for case in CASES if "expected_cdna" in case
                and "protein_expectation" not in case]


@pytest.mark.parametrize("case", CODING_CASES, ids=lambda c: c["id"])
@pytest.mark.parametrize("mode", ["default", "fast"])
def test_public_coding_variant_full_protein(native_genomes, case, mode):
    variant = variant_for(case, native_genomes)
    transcript = variant.genome.transcript_by_id(case["transcript"])
    assert transcript.strand == case["strand"]
    assert transcript.protein_sequence == case["reference_protein"]
    effects = variant.effects(**({"annotator": "fast"} if mode == "fast" else {}))
    candidates = [candidate for root in effects for candidate in iter_candidate_effects(root)
                  if candidate.transcript_id == transcript.id]
    assert any(effect.mutant_protein_sequence == case["expected_protein"]
               for effect in candidates), [(type(e).__name__, e.short_description)
                                          for e in candidates]


def test_brindled_mouse_allele_combines_two_published_cis_changes(native_genomes):
    cases = [c for c in CASES if c["id"].startswith("Atp7a-")]
    variants = [variant_for(c, native_genomes) for c in cases]
    protein = cases[0]["reference_protein"]
    expected = protein[:513] + "T" + protein[514:800] + protein[802:]
    assert protein[513] == "A" and protein[800:802] == "LA"
    effects = VariantCollection(variants).effects(
        phase_resolver=SimpleNamespace(in_cis=lambda *args, **kwargs: True))
    assert any(e.mutant_protein_sequence == expected
               for root in effects for e in iter_candidate_effects(root))


def test_mouse_mitochondrial_minus_strand_uses_vertebrate_code(native_genomes):
    case = next(c for c in CASES if c["id"] == "mt-Nd6-m1Dwa")
    assert case["strand"] == "-"
    # The pinned GRCm39 sequence is not the original LA9/C3H mutant background.
    assert case["published_stop_position"] == 79
    assert len(case["expected_protein"]) == case["reference_stop_position"] - 1 == 79
    assert case["expected_protein"][62:] == "DVGCVWMYDCYSYWGMS"
    genome = native_genomes["mouse"]
    assert genome.reference_metadata["translation_tables"]["mitochondrial"] == 2
    effect = variant_for(case, native_genomes).effect_on_transcript(
        genome.transcript_by_id(case["transcript"]))
    assert effect.mutant_protein_sequence == case["expected_protein"]


def test_dog_foxi3_isoform_difference_is_explicit(native_genomes):
    case = next(c for c in CASES if c["id"] == "FOXI3-dup7")
    refseq = "".join((DATA / "FOXI3-NP_001129118.1.fa").read_text().splitlines()[1:])
    extension = case["ensembl_n_terminal_extension"]
    assert case["reference_protein"][extension:] == refseq
    assert case["reference_protein"][extension + 22] == "A"
    assert case["expected_protein"][extension + 22] == "R"
    # RefSeq p.A23Rfs*219: 22 unchanged residues, then 218 translated residues.
    assert len(case["expected_protein"]) - extension == 22 + 218


def test_chm13_new_sequence_keeps_both_duplicated_transcript_identities(native_genomes):
    cases = [c for c in CASES if c["id"].startswith("1KGP-") and c["position"] == 21798486]
    assert len(cases) == 2 and all(c["newly_resolved"] for c in cases)
    assert len({c["transcript"] for c in cases}) == 2
    effects = variant_for(cases[0], native_genomes).effects()
    for case in cases:
        assert any(e.transcript_id == case["transcript"]
                   and e.mutant_protein_sequence == case["expected_protein"]
                   for root in effects for e in iter_candidate_effects(root))
