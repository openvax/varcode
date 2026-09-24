"""Observed RNA imports never invent sequence, phase, or a coding fusion."""

import csv
import json
from io import StringIO
from pathlib import Path
from types import SimpleNamespace

import pytest

from varcode import (
    MutantTranscript, RNAEvidence, RNAEvidenceResolver, StructuralVariant,
    apply_rna_evidence_to_effects, load_exacto_fusions, make_fusion_outcome,
)
from varcode.effects.effect_classes import (
    GeneFusion, Intronic, StructuralVariantEffect, TranslocationToIntergenic,
)


def tx(identifier, strand, gene, contig="1"):
    return SimpleNamespace(id=identifier, name=identifier, strand=strand,
                           gene_id=gene, gene=SimpleNamespace(id=gene), contig=contig)


@pytest.fixture
def context():
    first, second = tx("T1", "-", "G1"), tx("T2", "+", "G2", "2")

    class Genome:
        def transcript_by_id(self, identifier):
            for transcript in [first, second]:
                if transcript.id == identifier:
                    return transcript
            raise ValueError("Unknown transcript " + identifier)

    variant = SimpleNamespace(is_structural=True, genome=Genome())
    return variant, first, second


def table(rows, fields=None):
    result = StringIO()
    writer = csv.DictWriter(result, fieldnames=fields or list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
    result.seek(0)
    return result


def structures(model="1", refs="T1,T2", sequence="ATGAAAGGGTAA"):
    rows = []
    for i, (start, end, tid, strand, chrom) in enumerate([
            (0, 5, "T1", "-", "1"), (6, 11, "T2", "+", "2")]):
        rows.append(dict(transcript_model_id=model, reference_transcript_ids=refs,
                         index=str(i), read_start=str(start), read_end=str(end),
                         sequence=sequence[start:end+1], type="base", kind="match",
                         context="exonic", chromosome_1=chrom, position_1="100",
                         strand_1=strand, chromosome_2=chrom, position_2="105",
                         strand_2=strand, transcript_id_1=tid, transcript_id_2=tid))
    return rows


def link(model="1", dna="D1", refs="T1,T2", rna="R1"):
    return dict(transcript_model_id=model, reference_transcript_ids=refs,
                rna_variant_call_id=rna, dna_variant_call_id=dna, distance="20")


def load(context, rows=None, links=None, **kwargs):
    return load_exacto_fusions(
        table(structures() if rows is None else rows),
        table([link()] if links is None else links),
        variants_by_id={"D1": context[0]}, **kwargs)


def test_sequence_constructor_reuses_structural_shape():
    mt = MutantTranscript.from_sequence("AACGT", annotator_name="rna")
    assert mt.is_structural and not mt.is_identical_to_reference
    assert mt.reference_segments[0].source.sequence == mt.cdna_sequence == "AACGT"
    assert mt.mutant_protein_sequence is None
    assert mt.reference_segments[0].length == 5


def test_existing_assembled_allele_pickles_remain_readable(monkeypatch):
    import pickle
    from varcode.effects import structural
    original = structural._AssembledAllele
    legacy_type = type("_AssembledAllele", (), {
        "__slots__": ("sequence",), "__module__": structural.__name__})
    monkeypatch.setattr(structural, "_AssembledAllele", legacy_type)
    legacy = legacy_type()
    legacy.sequence = "ACGT"
    serialized = pickle.dumps(legacy)
    monkeypatch.setattr(structural, "_AssembledAllele", original)
    assert pickle.loads(serialized).sequence == "ACGT"


def test_import_uses_existing_effects_without_double_reverse_complement(context):
    evidence = load(context)
    assert isinstance(evidence, RNAEvidenceResolver)
    (candidate,) = evidence.candidates
    assert isinstance(candidate.effect, GeneFusion)
    assert candidate.effect.five_prime_transcript is context[1]
    assert candidate.effect.three_prime_transcript is context[2]
    assert candidate.effect.mutant_transcript.cdna_sequence == "ATGAAAGGGTAA"
    assert candidate.effect.mutant_protein_sequence is None
    assert candidate.source == "exacto"
    assert candidate.evidence["dna_variant_call_id"] == "D1"
    assert candidate.evidence["rna_variant_call_ids"] == ["R1"]
    assert "read_count" not in candidate.evidence
    assert evidence.observed_outcomes(context[0], context[1]) == (candidate,)
    assert evidence.observed_outcomes(context[0], context[2]) == (candidate,)
    assert evidence.observed_outcomes(context[0], tx("OTHER", "+", "G3")) == ()


def primary_rows(sequence="ATGAAAGGGTAA", peptide="1", start=0):
    from varcode.effects.codon_tables import STANDARD
    rows = []
    for i, nt in enumerate(sequence):
        codon = sequence[i // 3 * 3:i // 3 * 3 + 3]
        aa = "*" if codon in STANDARD.stop_codons else STANDARD.forward_table.get(codon, "")
        rows.append(dict(peptide_id=peptide, primary_structure_index=str(i),
                         type="base", amino_acid=aa, amino_acid_index=str(i // 3),
                         codon_index=str(i % 3), nucleotide=nt.lower(),
                         transcript_model_id="1", reference_transcript_ids="T1,T2",
                         transcript_structure_index=str((start + i) // 6),
                         read_start=str(start + i), read_end=str(start + i),
                         codon_dna_variant_call_ids="D1"))
    return rows


def test_native_primary_structure_import(context):
    rows = primary_rows()
    (candidate,) = load(context, primary_structures_path=table(rows[::-1])).candidates
    assert candidate.effect.mutant_protein_sequence == "MKG"
    assert candidate.evidence["protein_completeness"] == "start_to_stop"
    assert candidate.evidence["translation_evidence"] == "sequence_prediction_only"
    assert candidate.evidence["exacto_primary_structure"] == rows
    assert candidate.effect.mutant_transcript.evidence == candidate.evidence


@pytest.mark.parametrize("reference_side", ["five_prime", "three_prime"])
@pytest.mark.parametrize("protein_input", ["native", "explicit_start"])
@pytest.mark.parametrize("reference_coding,reference_protein,changed", [
    ("ATGCCAATGAAAGGGTAA", "MPMKG", None),
    ("ATGCCACAATAA", "MPQ", True),
    ("ATGAAAGGGTAA", "MKG", False),
    ("ATGAAAGGGTGGTAA", "MKGW", True),
])
def test_start_to_stop_import_does_not_infer_missing_five_prime_coverage(
        context, reference_side, protein_input, reference_coding, reference_protein, changed):
    from varcode.effects import EffectCollection

    for transcript in context[1:]:
        transcript.protein_sequence = reference_protein
        transcript.coding_sequence = transcript.sequence = reference_coding
        transcript.complete = True
        transcript.start_codon_spliced_offsets = [0, 1, 2]
    options = ({"primary_structures_path": table(primary_rows())} if protein_input == "native"
               else {"cds_starts": {("1", ("T1", "T2")): 0}})
    candidate, = load(context, **options).candidates
    effect = candidate.effect
    if reference_side == "three_prime":
        effect = GeneFusion(context[0], context[2], context[1],
                            mutant_transcript=effect.mutant_transcript,
                            five_prime_transcript=context[1], three_prime_transcript=context[2])
    before = dict(candidate.evidence)
    assert effect.modifies_coding_sequence is changed
    assert effect.modifies_protein_sequence is changed
    assert effect.mutant_protein_sequence == "MKG"
    assert candidate.evidence == before == effect.mutant_transcript.evidence
    if protein_input == "native":
        assert candidate.evidence["protein_completeness"] == "start_to_stop"
    assert candidate.evidence["sequence_status"] == "observed_model_completeness_unknown"
    assert len(EffectCollection([effect]).drop_silent_and_noncoding()) == (changed is not False)
    assert len(EffectCollection([effect]).drop_silent_and_noncoding(False)) == (changed is True)


@pytest.mark.parametrize("partner_coding,partner_protein,changed", [
    ("ATGCCAATGAAAGGGTAA", "MPMKG", None),
    ("ATGCCACAATAA", "MPQ", True),
])
def test_observed_fusion_comparison_uses_each_partners_own_reference(
        context, partner_coding, partner_protein, changed):
    for transcript, coding, protein in [
            (context[1], "ATGAAAGGGTAA", "MKG"),
            (context[2], partner_coding, partner_protein)]:
        transcript.protein_sequence = protein
        transcript.coding_sequence = transcript.sequence = coding
        transcript.complete = True
        transcript.start_codon_spliced_offsets = [0, 1, 2]
    candidate, = load(context, primary_structures_path=table(primary_rows())).candidates
    five_prime = candidate.effect
    three_prime = GeneFusion(context[0], context[2], context[1],
                             mutant_transcript=five_prime.mutant_transcript,
                             five_prime_transcript=context[1], three_prime_transcript=context[2])
    assert five_prime.modifies_coding_sequence is False
    assert five_prime.modifies_protein_sequence is False
    assert three_prime.modifies_coding_sequence is changed
    assert three_prime.modifies_protein_sequence is changed


def test_primary_structure_keeps_multiple_orfs(context):
    rows = primary_rows() + primary_rows("AAAGGGTAA", peptide="2", start=3)
    candidates = load(context, primary_structures_path=table(rows)).candidates
    assert [c.effect.mutant_protein_sequence for c in candidates] == ["MKG", "KG"]
    assert candidates[1].evidence["protein_completeness"] == "partial_start"
    assert candidates[1].evidence["cds_start"] == 3


@pytest.mark.parametrize("length,protein,trailing", [(9, "MKG", 0), (10, "MKG", 1), (11, "MKG", 2)])
def test_incomplete_primary_protein_remains_visible(context, length, protein, trailing):
    rows = primary_rows("ATGAAAGGGTAA"[:length])
    (candidate,) = load(context, primary_structures_path=table(rows)).candidates
    assert candidate.effect.mutant_protein_sequence == protein
    assert candidate.evidence["protein_completeness"] == "partial_end"
    assert candidate.evidence["trailing_partial_codon_bases"] == trailing


@pytest.mark.parametrize("field,value", [
    ("primary_structure_index", "1"), ("codon_index", "2"),
    ("amino_acid_index", "1"), ("nucleotide", "C"), ("amino_acid", "K"),
    ("read_end", "2"), ("transcript_structure_index", "99"), ("peptide_id", "")])
def test_invalid_primary_records_fail_closed(context, field, value):
    rows = primary_rows()
    rows[0][field] = value
    with pytest.raises(ValueError):
        load(context, primary_structures_path=table(rows))


def test_primary_and_manual_start_are_mutually_exclusive(context):
    with pytest.raises(ValueError, match="not both"):
        load(context, primary_structures_path=table(primary_rows()), cds_starts={})


def test_model_without_primary_rows_stays_untranslated(context):
    rows = primary_rows()
    for row in rows:
        row["transcript_model_id"] = "unselected"
    (candidate,) = load(context, primary_structures_path=table(rows)).candidates
    assert candidate.effect.mutant_protein_sequence is None


def test_primary_events_preserve_codon_order_and_provenance(context):
    structure = structures()
    structure[1]["index"] = "2"
    event = dict(structure[0], index="1", type="event", kind="splicing",
                 read_start="5", read_end="6", sequence="")
    structure.insert(1, event)
    primary = primary_rows()
    for row in primary[6:]:
        row["primary_structure_index"] = str(int(row["primary_structure_index"]) + 1)
        row["transcript_structure_index"] = "2"
    primary.insert(6, dict(primary[0], primary_structure_index="6", type="event",
                           transcript_structure_index="1", read_start="5", read_end="6",
                           amino_acid="", nucleotide="", codon_index="-1", amino_acid_index="-1"))
    (candidate,) = load(context, rows=structure, primary_structures_path=table(primary)).candidates
    assert candidate.effect.mutant_protein_sequence == "MKG"
    assert candidate.evidence["exacto_primary_structure"][6]["type"] == "event"


def test_primary_structure_cannot_translate_past_stop(context):
    sequence = "ATGTAAGGGTAA"
    with pytest.raises(ValueError, match="after a stop"):
        load(context, rows=structures(sequence=sequence),
             primary_structures_path=table(primary_rows(sequence)))


def test_one_dna_event_keeps_multiple_models_and_reference_groups(context):
    rows = structures() + structures(model="2") + structures(refs="T1.1,T2.2")
    links = [link(), link(model="2"), link(refs="T1.1,T2.2")]
    evidence = load(context, rows=rows, links=links)
    assert len(evidence.candidates) == 3
    assert {tuple(c.evidence["reference_transcript_ids"]) for c in evidence.candidates} == {
        ("T1", "T2"), ("T1.1", "T2.2")}


def test_duplicate_integration_rows_do_not_duplicate_models(context):
    evidence = load(context, links=[link(), link(), link(rna="R2")])
    assert len(evidence.candidates) == 1
    assert evidence.candidates[0].evidence["rna_variant_call_ids"] == ["R1", "R2"]


def test_explicit_complete_orf_is_predicted_not_translation_evidence(context):
    (candidate,) = load(context, cds_starts={("1", ("T1", "T2")): 0}).candidates
    assert candidate.effect.mutant_protein_sequence == "MKG"
    assert candidate.evidence["protein_status"] == "predicted_from_observed_rna"
    assert candidate.evidence["cds_end"] == 12


@pytest.mark.parametrize("sequence,start", [
    ("ATGAAAGGG", 0), ("ATGNNNTAA", 0), ("CCCAAATAA", 0),
    ("ATGAAATAA", -1), ("ATGAAATAA", True), ("ATGAAATAA", 100)])
def test_invalid_or_incomplete_orf_is_not_translated(context, sequence, start):
    with pytest.raises(ValueError):
        make_fusion_outcome(context[0], context[1], sequence=sequence,
                            transcript_model_id="x", cds_start=start)


@pytest.mark.parametrize("sequence", ["", "ATGZ", "ATG U", None])
def test_invalid_sequence_rejected(context, sequence):
    with pytest.raises(ValueError, match="sequence"):
        make_fusion_outcome(context[0], context[1], sequence=sequence, transcript_model_id="x")


@pytest.mark.parametrize("count", [-1, 2.5, True])
def test_invalid_read_count_rejected(context, count):
    with pytest.raises(ValueError, match="read_count"):
        make_fusion_outcome(context[0], context[1], sequence="ACGT",
                            transcript_model_id="x", read_count=count)


def test_antisense_partner_does_not_become_a_coding_fusion(context):
    rows = structures()
    rows[-1].update(strand_1="-", strand_2="-")
    (candidate,) = load(context, rows=rows).candidates
    assert isinstance(candidate.effect, TranslocationToIntergenic)
    assert candidate.evidence["partner_status"] == "antisense"
    assert candidate.effect.mutant_protein_sequence is None


def test_intergenic_partner_and_splice_event_preserved(context):
    rows = structures()
    rows[-1].update(transcript_id_1="", transcript_id_2="", context="intergenic", index="2")
    event = dict(rows[0], index="1", type="event", kind="splicing", sequence="",
                 read_start="5", read_end="6", context="noncanonical")
    (candidate,) = load(context, rows=[rows[1], event, rows[0]]).candidates
    assert isinstance(candidate.effect, TranslocationToIntergenic)
    assert candidate.effect.mutant_transcript.cdna_sequence == "ATGAAAGGGTAA"
    assert candidate.evidence["exacto_structure"][1]["kind"] == "splicing"


@pytest.mark.parametrize("change", [
    {"index": "0"}, {"index": "3"}, {"read_start": "7"},
    {"read_start": "5"}, {"sequence": ""}, {"type": "mystery"},
    {"type": "event"}, {"kind": "circular"}, {"transcript_id_2": "UNKNOWN"}])
def test_malformed_or_unsupported_structures_fail_explicitly(context, change):
    rows = structures()
    rows[-1].update(change)
    with pytest.raises(ValueError):
        load(context, rows=rows)


def test_unknown_model_does_not_silently_disappear(context):
    with pytest.raises(ValueError, match="Missing structures"):
        load(context, links=[link(model="missing")])


def test_same_gene_and_antisense_anchor_rejected(context):
    rows = structures()
    rows[-1].update(transcript_id_1="T1", transcript_id_2="T1", strand_1="-", strand_2="-",
                    chromosome_1="1", chromosome_2="1", position_1="90", position_2="95")
    with pytest.raises(ValueError, match="two-gene"):
        load(context, rows=rows)
    rows = structures()
    rows[0]["strand_1"] = "+"
    with pytest.raises(ValueError, match="sense 5-prime"):
        rows[0]["strand_2"] = "+"
        load(context, rows=rows)


def test_backsplice_is_not_flattened_into_a_linear_fusion(context):
    rows = structures()
    repeat = dict(rows[0], index="1", read_start="6", read_end="11")
    rows[-1].update(index="2", read_start="12", read_end="17")
    with pytest.raises(ValueError, match="back-spliced"):
        load(context, rows=[rows[0], repeat, rows[1]])


def test_return_to_first_locus_is_not_a_two_locus_model(context):
    rows = structures()
    repeat = dict(rows[0], index="2", read_start="12", read_end="17")
    with pytest.raises(ValueError, match="More than two loci"):
        load(context, rows=[*rows, repeat])


def test_selected_variant_ids_only(context):
    assert load(context, links=[link(dna="unselected")]).candidates == ()


def test_intronic_dna_prediction_does_not_hide_observed_sv(context):
    evidence = load(context)
    original = Intronic(context[0], context[1], nearest_exon=None, distance_to_exon=10)
    effects = [original]
    result = apply_rna_evidence_to_effects(effects, evidence)
    assert result is effects
    assert isinstance(effects[0], StructuralVariantEffect)
    assert effects[0].candidates[0].effect is original
    assert effects[0].candidates[-1] is evidence.candidates[0]


def test_no_observation_leaves_intronic_sv_unchanged(context):
    original = Intronic(context[0], context[1], nearest_exon=None, distance_to_exon=10)
    effects = [original]
    assert apply_rna_evidence_to_effects(effects, RNAEvidence()) == [original]


def test_point_variant_still_rejected_by_fusion_factory(context):
    with pytest.raises(ValueError, match="StructuralVariant"):
        make_fusion_outcome(SimpleNamespace(is_structural=False), context[1],
                            sequence="ACGT", transcript_model_id="x")


def test_tsv_headers_and_ragged_rows_fail(context):
    for stream in [StringIO("bad\nvalue\n"), StringIO("transcript_model_id\ttranscript_model_id\n1\t1\n")]:
        with pytest.raises(ValueError, match="columns"):
            load_exacto_fusions(stream, table([link()]), variants_by_id={"D1": context[0]})


def test_gzip_paths(context, tmp_path):
    import gzip
    path = tmp_path / "structures.tsv.gz"
    with gzip.open(path, "wt") as h:
        h.write(table(structures()).getvalue())
    assert len(load_exacto_fusions(path, table([link()]),
                                  variants_by_id={"D1": context[0]}).candidates) == 1


def test_audited_osteosarc_junction_fragments_remain_untranslated():
    path = Path(__file__).parent / "data/osteosarc_observed_junctions.json"
    records = json.loads(path.read_text())["records"]
    assert len(records) == 3
    for row in records:
        variant = StructuralVariant(
            contig=row["contig"], start=row["start"], sv_type="BND",
            mate_contig=row["mate_contig"], mate_start=row["mate_start"],
            alt="[%s:%s[N" % (row["mate_contig"], row["mate_start"]),
            genome="GRCh38")
        # Annotation identity/strand is pinned in the fixture, not fetched from
        # a potentially different Ensembl release during a regression test.
        transcript = tx(row["transcript_id"], row["strand"], row["gene"], row["contig"])
        candidate = make_fusion_outcome(
            variant, transcript, sequence=row["sequence"],
            transcript_model_id=row["read_name"], source=row["source"],
            read_count=1, extra_evidence=row)
        assert isinstance(candidate.effect, TranslocationToIntergenic)
        assert candidate.effect.mutant_transcript.cdna_sequence == row["sequence"]
        assert candidate.effect.mutant_protein_sequence is None
        assert candidate.evidence["sequence_status"] == "junction_fragment"
        assert len(variant.junctions) == 1
        if row["label"] == "OTUD7A-FMN1":
            assert row["sequence"][row["left_end"]:row["right_start"]] == "AG"
            assert row["partner_status"] == "antisense"


def test_point_variant_rna_evidence_contract_unchanged(context):
    context[0].is_structural = False
    original = Intronic(context[0], context[1], nearest_exon=None, distance_to_exon=10)
    class UnexpectedResolver:
        def observed_outcomes(self, variant, transcript):
            raise AssertionError("Do not change deterministic point-variant behavior")
    assert apply_rna_evidence_to_effects([original], UnexpectedResolver()) == [original]
