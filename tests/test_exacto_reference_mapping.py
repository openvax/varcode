"""Exacto coordinates anchor observed RNA; missing coverage is not a deletion.

The short references below have independently specified codons and exon maps.
Native mismatch/insertion bounds follow Exacto's identify_records at revision
307c08670d5e706734bddf393bcebc84db497f9f (see docs/rna_structures.md).
"""

from dataclasses import replace
from types import SimpleNamespace

import pytest

from varcode import load_exacto_fusions
from varcode.effects import EffectCollection, GeneFusion

from .test_exacto_fusions import link, table


def transcript(identifier, strand, contig, sequence):
    # Split a codon at the first exon boundary on either strand.
    intervals = [(100, 106), (200, 216)]
    if strand == "-":
        intervals = [(100, 116), (200, 206)]
    positions = [p for lo, hi in intervals for p in range(lo, hi + 1)]
    if strand == "-":
        positions.reverse()
    return SimpleNamespace(
        id=identifier, name=identifier, gene_id=identifier,
        gene=SimpleNamespace(id=identifier), strand=strand, contig=contig,
        sequence=sequence, coding_sequence=sequence,
        protein_sequence={"T1": "MKPGFPK", "T2": "MPKGKFP"}[identifier],
        exon_intervals=intervals, spliced_offset=positions.index,
        start_codon_spliced_offsets=[0, 1, 2],
        stop_codon_spliced_offsets=[21, 22, 23], complete=True,
        positions=positions)


@pytest.fixture(params=["+", "-"])
def context(request):
    first = transcript("T1", request.param, "1", "ATGAAACCCGGGTTTCCCAAATAA")
    second = transcript("T2", "+" if request.param == "-" else "-", "2",
                        "ATGCCCAAAGGGAAATTTCCCTAA")

    class Genome:
        def transcript_by_id(self, identifier):
            for t in (first, second):
                if identifier == t.id:
                    return t
            raise ValueError("Unknown transcript " + identifier)

    return SimpleNamespace(is_structural=True, genome=Genome()), first, second


def rows_for(parts):
    """Write native runs, using flanks for non-reference bases.

    Each part is (transcript, start, end, kind, observed_or_None). Offsets
    refer to the short references above; insertions have start == end.
    """
    rows = []
    read_offset = 0
    for t, start, end, kind, observed in parts:
        sequence = t.sequence[start:end] if observed is None else observed
        if kind == "match":
            # A native match row stays within one contiguous genomic run.
            cuts = [0] + [i for i in range(1, end - start)
                          if abs(t.positions[start + i] - t.positions[start + i - 1]) != 1]
            cuts.append(len(sequence))
        else:
            cuts = [0, len(sequence)]
        for a, b in zip(cuts, cuts[1:]):
            positions = (t.positions[start + a], t.positions[start + b - 1]) if kind == "match" else (
                t.positions[start - 1], t.positions[end])
            lo, hi = sorted(positions)
            rows.append(dict(
                transcript_model_id="1", reference_transcript_ids="T1,T2",
                index=str(len(rows)), read_start=str(read_offset), read_end=str(read_offset + b - a - 1),
                sequence=sequence[a:b].lower(), type="base", kind=kind, context="exonic",
                chromosome_1=t.contig, chromosome_2=t.contig,
                position_1=str(lo), position_2=str(hi), strand_1=t.strand, strand_2=t.strand,
                transcript_id_1=t.id + ".1", transcript_id_2=t.id + ".1"))
            read_offset += b - a
    return rows


def primary(rows, start=0):
    # Explicit codon oracle for these synthetic sequences, independent of
    # Varcode's translator. A partial terminal codon carries no amino acid.
    code = {"ATG": "M", "AAA": "K", "CCC": "P", "GGG": "G", "TTT": "F",
            "CCA": "P", "CGG": "R", "AGG": "R", "GTT": "V", "TCC": "S",
            "ATT": "I", "TTC": "F", "TAA": "*"}
    sequence = "".join(r["sequence"] for r in rows).upper()
    sources = [r["index"] for r in rows for _ in r["sequence"]]
    result = []
    for i in range(start, len(sequence)):
        j = i - start
        codon = sequence[start + j // 3 * 3:start + j // 3 * 3 + 3]
        result.append(dict(
            peptide_id="P1", primary_structure_index=str(j), type="base",
            amino_acid=code[codon] if len(codon) == 3 else "",
            amino_acid_index=str(j // 3), codon_index=str(j % 3), nucleotide=sequence[i],
            transcript_model_id="1", reference_transcript_ids="T1,T2",
            transcript_structure_index=sources[i], read_start=str(i), read_end=str(i)))
    return result


def load(context, rows, translate=False):
    options = {"primary_structures_path": table(primary(rows))} if translate else {}
    candidate, = load_exacto_fusions(
        table(rows), table([link()]), variants_by_id={"D1": context[0]}, **options).candidates
    return candidate


def test_native_matches_map_both_transcripts_without_changing_observed_rna(context):
    _, first, second = context
    rows = rows_for([(first, 0, 9, "match", None), (second, 9, 18, "match", None)])
    candidate = load(context, rows)
    model = candidate.effect.mutant_transcript
    assert model.cdna_sequence == "ATGAAACCCGGGAAATTT"
    assert [(s.source.id, s.start, s.end, s.strand) for s in model.reference_segments] == [
        ("T1", 0, 7, "+"), ("T1", 7, 9, "+"), ("T2", 9, 18, "+")]
    assert "".join(s.source.sequence[s.start:s.end] for s in model.reference_segments) == model.cdna_sequence
    assert candidate.evidence["exacto_structure"] == rows
    assert model.evidence == candidate.evidence
    assert model.mutant_protein_sequence is None


@pytest.mark.parametrize("kind", ["mismatch", "insertion"])
def test_native_nonreference_bounds_do_not_reject_a_linear_path(context, kind):
    _, first, second = context
    end = 10 if kind == "mismatch" else 9
    rows = rows_for([(first, 0, 9, "match", None), (first, 9, end, kind, "A"),
                     (first, end, 15, "match", None), (second, 15, 18, "match", None)])
    candidate = load(context, rows)
    model = candidate.effect.mutant_transcript
    assert model.cdna_sequence == "".join(r["sequence"] for r in rows).upper()
    assert "".join(s.source.sequence[s.start:s.end] for s in model.reference_segments) == model.cdna_sequence
    assert any(s.source is first for s in model.reference_segments)
    assert any(s.source is second for s in model.reference_segments)
    assert any(s.source is not first and s.source is not second for s in model.reference_segments)


def test_partial_import_establishes_changed_junction_sequence(context):
    _, first, second = context
    rows = rows_for([(first, 0, 9, "match", None), (second, 9, 18, "match", None)])
    candidate = load(context, rows, translate=True)
    assert candidate.effect.mutant_protein_sequence == "MKPGKF"
    assert candidate.evidence["protein_completeness"] == "partial_end"
    assert candidate.effect.modifies_coding_sequence is True
    assert candidate.effect.modifies_protein_sequence is True
    assert len(EffectCollection([candidate.effect]).drop_silent_and_noncoding(False)) == 1
    # The same observation must be interpretable relative to its 3' partner.
    other = GeneFusion(context[0], second, first,
                       mutant_transcript=candidate.effect.mutant_transcript,
                       five_prime_transcript=first, three_prime_transcript=second)
    assert other.modifies_protein_sequence is True


def test_partial_import_with_identical_observed_codons_stays_unresolved(context):
    _, first, second = context
    # T2 codons 3-4 are GGG AAA, also present at T1's matching coordinates
    # after choosing a reference with the same observed continuation.
    first.sequence = first.coding_sequence = "ATGAAACCCGGGAAATTTCCCTAA"
    first.protein_sequence = "MKPGKFP"
    rows = rows_for([(first, 0, 9, "match", None), (second, 9, 18, "match", None)])
    effect = load(context, rows, translate=True).effect
    assert effect.mutant_protein_sequence == "MKPGKF"
    assert effect.modifies_coding_sequence is None
    assert effect.modifies_protein_sequence is None


@pytest.mark.parametrize("change,expected", [("insertion", "MKPRVSF"), ("deletion", "MKPIFP")])
def test_partial_import_recognizes_an_observed_frameshift(context, change, expected):
    _, first, second = context
    parts = [(first, 0, 9, "match", None)]
    if change == "insertion":
        parts += [(first, 9, 9, "insertion", "A"), (first, 9, 17, "match", None),
                  (second, 15, 18, "match", None)]
    else:
        first.sequence = first.coding_sequence = "ATGAAACCCGATTTTCCCAAATAA"
        first.protein_sequence = "MKPDFPK"
        parts += [(first, 10, 16, "match", None), (second, 18, 21, "match", None)]
    effect = load(context, rows_for(parts), translate=True).effect
    assert effect.mutant_protein_sequence == expected
    assert effect.modifies_coding_sequence is True
    assert effect.modifies_protein_sequence is True


def test_native_mismatch_changes_a_partial_protein_without_replacing_the_observation(context):
    _, first, second = context
    rows = rows_for([(first, 0, 9, "match", None), (first, 9, 10, "mismatch", "A"),
                     (first, 10, 18, "match", None), (second, 18, 21, "match", None)])
    candidate = load(context, rows, translate=True)
    assert candidate.effect.mutant_protein_sequence == "MKPRFPP"
    assert candidate.effect.modifies_protein_sequence is True
    assert candidate.evidence["exacto_structure"] == rows


def test_inframe_deletion_compares_observed_continuation_not_each_codons_origin(context):
    _, first, second = context
    # Remove the P codon; all retained bases individually match reference.
    rows = rows_for([(first, 0, 6, "match", None), (first, 9, 18, "match", None),
                     (second, 18, 21, "match", None)])
    effect = load(context, rows, translate=True).effect
    assert effect.mutant_protein_sequence == "MKGFPP"
    assert effect.modifies_protein_sequence is True


def test_mapped_bases_do_not_establish_an_out_of_frame_orf(context):
    _, first, second = context
    rows = rows_for([(first, 0, 9, "match", None), (second, 9, 18, "match", None)])
    effect = load(context, rows).effect
    model = effect.mutant_transcript
    effect.mutant_transcript = replace(model, evidence=dict(
        model.evidence, cds_start=1, cds_end=18, protein_completeness="partial_both"))
    assert effect.modifies_coding_sequence is None
    assert effect.modifies_protein_sequence is None


def test_bases_after_an_observed_stop_cannot_supply_a_frame_anchor(context):
    _, first, second = context
    rows = rows_for([(first, 0, 6, "match", "ccctaa"), (second, 9, 18, "match", None)])
    effect = load(context, rows).effect
    model = effect.mutant_transcript
    effect.mutant_transcript = replace(model, evidence=dict(
        model.evidence, cds_start=0, cds_end=18, protein_completeness="partial_both"))
    assert effect.modifies_coding_sequence is None
    assert effect.modifies_protein_sequence is None


def test_native_validation_still_rejects_backtracking_around_an_insertion(context):
    _, first, second = context
    rows = rows_for([(first, 0, 12, "match", None), (first, 9, 9, "insertion", "A"),
                     (second, 9, 18, "match", None)])
    with pytest.raises(ValueError, match="Overlapping or back-spliced"):
        load(context, rows)


@pytest.mark.parametrize("unmapped", ["wrong_base", "ambiguous_base", "intronic", "antisense", "insertion", "length"])
def test_annotations_do_not_substitute_for_reference_sequence(context, unmapped):
    _, first, second = context
    rows = rows_for([(first, 0, 6, "match", None), (second, 9, 15, "match", None)])
    target = rows[-1]
    if unmapped == "wrong_base":
        target["sequence"] = "accaaa"
    elif unmapped == "ambiguous_base":
        target["sequence"] = "nnnaaa"
    elif unmapped == "intronic":
        target.update(position_1="150", position_2="155")
    elif unmapped == "antisense":
        target.update(strand_1=first.strand, strand_2=first.strand)
    elif unmapped == "insertion":
        target.update(kind="insertion", position_1="150", position_2="151")
    else:
        target["position_2"] = str(int(target["position_2"]) + 1)
    candidate = load(context, rows)
    model = candidate.effect.mutant_transcript
    assert model.cdna_sequence == "".join(r["sequence"] for r in rows).upper()
    assert "".join(s.source.sequence[s.start:s.end] for s in model.reference_segments) == model.cdna_sequence
    mapped = sum(s.length for s in model.reference_segments if s.source is second)
    assert mapped == (3 if unmapped in {"wrong_base", "ambiguous_base"} else 0)
    assert model.mutant_protein_sequence is None


def test_real_transcript_mapping_survives_effect_json_roundtrip():
    from pyensembl import cached_release
    from varcode import StructuralVariant

    genome = cached_release(81)
    first = genome.transcript_by_id("ENST00000003084")
    second = genome.transcript_by_id("ENST00000357654")
    # Use two literal exon pieces; expected source offsets come independently
    # from the transcript's public genomic-to-cDNA coordinate API.
    rows = []
    expected = []
    for t in (first, second):
        lo, hi = t.exon_intervals[1]
        positions = [lo, lo + 8] if t.strand == "+" else [hi - 8, hi]
        left = min(t.spliced_offset(p) for p in positions)
        expected.append((t.id, left, left + 9))
        rows.append(dict(
            transcript_model_id="1", reference_transcript_ids="T1,T2", index=str(len(rows)),
            read_start=str(9 * len(rows)), read_end=str(9 * len(rows) + 8),
            sequence=t.sequence[left:left + 9].lower(), type="base", kind="match", context="exonic",
            chromosome_1=t.contig, chromosome_2=t.contig,
            position_1=str(positions[0]), position_2=str(positions[1]),
            strand_1=t.strand, strand_2=t.strand, transcript_id_1=t.id, transcript_id_2=t.id))
    variant = StructuralVariant(first.contig, first.start, "BND", genome=genome,
                                mate_contig=second.contig, mate_start=second.start)
    effect = load((variant, first, second), rows).effect
    restored = GeneFusion.from_json(effect.to_json())
    model = restored.mutant_transcript
    assert [(s.source.id, s.start, s.end) for s in model.reference_segments] == expected
    assert model.cdna_sequence == effect.mutant_transcript.cdna_sequence
    assert model.evidence == effect.mutant_transcript.evidence
