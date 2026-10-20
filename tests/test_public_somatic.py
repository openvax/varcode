"""Complete public call sets, with independent canine consequence checks.

See data/public_somatic/README.md for provenance, interpretation and limits.
References are installed explicitly before pytest (see CONTRIBUTING.md).
"""

from collections import Counter
import hashlib
import json

import pytest
from pyensembl import cached_release

from varcode import load_vcf
from varcode.effects import (
    AlternateStartCodon, IncompleteTranscript, NoncodingTranscript, Silent,
    SpliceAcceptor, StopLoss,
)
from varcode.splice_outcomes import SpliceOutcomeSet
from varcode.transforms import pair_breakends

from .public_somatic import (
    DATA, MANIFEST, assert_published_protein_change, canine_variant,
    hcc1395_structural_variant, materialize_effects, table_rows, vcf_rows,
)


VCFS = [name for name in MANIFEST["fixtures"] if ".vcf" in name]
CODING_CLASSES = {
    "Missense_Mutation": {"Substitution", "ComplexSubstitution"},
    "Nonsense_Mutation": {"PrematureStop"},
    "Frame_Shift_Ins": {"FrameShift", "FrameShiftTruncation"},
    "Frame_Shift_Del": {"FrameShift", "FrameShiftTruncation"},
    "In_Frame_Ins": {"Insertion", "ComplexSubstitution"},
    "In_Frame_Del": {"Deletion", "ComplexSubstitution"},
    "Silent": {"Silent"},
}


@pytest.mark.parametrize("name", VCFS)
def test_every_public_vcf_allele_loads_and_annotates(name):
    spec = MANIFEST["fixtures"][name]
    assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == spec["sha256"]
    rows = vcf_rows(name)
    assert len(rows) == spec["records"]
    assert Counter(row[6] for row in rows) == spec["filters"]
    # Every pinned VCF is biallelic. Fail if a future fixture needs expansion.
    assert all("," not in row[4] for row in rows)
    genome = cached_release(spec["ensembl_release"])
    variants = load_vcf(
        str(DATA / name), genome=genome, only_passing=False,
        parse_structural_variants=True)
    assert len(variants) == len(rows)
    expected = Counter((row[0].removeprefix("chr"), int(row[1]), row[3], row[4])
                       for row in rows)
    actual = Counter((v.contig, v.original_start, v.original_ref, v.original_alt)
                     for v in variants)
    assert actual == expected, "Source alleles must not be dropped or rewritten"
    for variant in variants:
        effects = variant.effects(raise_on_error=True)
        assert effects, str(variant)
        if getattr(variant, "is_structural", False):
            materialize_effects(effects)
    if name == "colo829-sv-truth.vcf":
        events = pair_breakends(variants)
        assert len(events) == 68  # 66 reciprocal pairs + two unpaired INS records
        assert Counter(v for event in events for v in (event.source_variants or (event,))) == Counter(variants)
        for event in events:
            effects = event.effects(raise_on_error=True)
            assert effects
            materialize_effects(effects)


def test_every_hcc1395_consensus_sv_loads_and_annotates():
    rows = table_rows("hcc1395-consensus-sv.tsv")
    assert len(rows) == 1788
    assert len({row["SV_id"] for row in rows}) == 1788
    assert Counter(row["SV_type"] for row in rows) == {
        "DEL": 717, "INS": 551, "DUP": 230, "TRA": 146, "INV": 133, "BND": 11}
    genome = cached_release(93)
    for row in rows:
        variant = hcc1395_structural_variant(row, genome)
        assert variant.info["published_consensus"] == row
        assert variant.alt_assembly is None
        if row["SV_type"] in ("BND", "TRA"):
            assert variant.mate_orientation is None
            assert variant.mate_start == int(row["Pos2"])
        effects = variant.effects(raise_on_error=True)
        assert effects, row["SV_id"]
        materialize_effects(effects)


def test_published_canine_consequences():
    rows = table_rows("canine-osa-vep99.tsv")
    assert len(rows) == 2141
    genome = cached_release(99, species="dog")
    limitations = json.loads((DATA / "canine-unresolved.json").read_text())
    observed_unresolved, unique_variants, counts = {}, set(), Counter()
    for source_row, row in enumerate(rows, 4):  # original spreadsheet row
        variant = canine_variant(row, genome)
        unique_variants.add(variant)
        transcript = genome.transcript_by_id(row["Transcript_ID"])
        effect = variant.effect_on_transcript(transcript)
        # VEP's exonic coding prediction assumes normal splicing. Preserve the
        # original outcome set and compare that explicitly conditional product.
        conditional = getattr(effect, "effect_if_splicing_unchanged", None)
        if conditional is not None:
            effect = conditional
        if isinstance(effect, (IncompleteTranscript, NoncodingTranscript)):
            reason = {
                "transcript": transcript.id,
                "biotype": transcript.biotype,
                "missing_start": not transcript.contains_start_codon,
                "missing_stop": not transcript.contains_stop_codon,
            }
            assert reason["missing_start"] or reason["missing_stop"]
            assert reason == limitations[str(source_row)], (source_row, row)
            observed_unresolved[str(source_row)] = reason
            counts["unresolved"] += 1
            continue

        hgvs = row["HGVSp_Short"]
        original = transcript.protein_sequence
        if row["Variant_Classification"] == "Splice_Site":
            assert row["Consequence"] == "splice_acceptor_variant"
            assert isinstance(effect, SpliceOutcomeSet)
            assert effect.disrupted_signal_class is SpliceAcceptor
            materialize_effects([effect])
            counts["splice"] += 1
        elif hgvs == "p.M1?":
            # The paper leaves the protein outcome unknown. Varcode recognizes
            # TTG as an alternative initiator; this is NOT protein concordance.
            assert isinstance(effect, AlternateStartCodon)
            assert (effect.ref_codon, effect.alt_codon) == ("ATG", "TTG")
            assert original[0] == "M"
            counts["start_uncertain"] += 1
        elif "ext" in hgvs:
            assert hgvs == "p.*136Qext*?"
            assert isinstance(effect, StopLoss)
            assert len(original) == 135
            assert effect.mutant_protein_sequence.startswith(original + "Q")
            counts["stop_loss"] += 1
        else:
            assert type(effect).__name__ in CODING_CLASSES[row["Variant_Classification"]], (
                source_row, row["Variant_Classification"], type(effect).__name__)
            mutant = original if isinstance(effect, Silent) else effect.mutant_protein_sequence
            try:
                assert_published_protein_change(hgvs, original, mutant)
            except AssertionError as error:
                raise AssertionError(
                    "Published spreadsheet row %s, %s: %s versus %s"
                    % (source_row, transcript.id, hgvs, effect.short_description)) from error
            counts["protein_match"] += 1
    assert len(unique_variants) == 1936
    assert observed_unresolved == limitations
    assert counts == {
        "protein_match": 1954, "stop_loss": 1, "splice": 1,
        "start_uncertain": 1, "unresolved": 184,
    }
