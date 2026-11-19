"""Effect classification under non-standard NCBI genetic codes (#701).

Expected protein sequences and codon meanings come independently from Biotite's
NCBI tables (with the gc.prt 4.6 corrections from tests/test_genetic_codes_independent.py),
never from Varcode's own codon tables or translation helpers.
"""

from itertools import product

import pytest
from biotite.sequence import CodonTable as BiotiteCodonTable
from pyensembl import Genome, cached_release

from .common import reverse_complement
from varcode import MutantTranscript, ReferenceSegment, StructuralVariant, Variant, VariantCollection
from varcode.effects import (
    AlternateStartCodon,
    Deletion,
    FrameShift,
    FrameShiftTruncation,
    GeneFusion,
    PrematureStop,
    Silent,
    StopLoss,
    StructuralVariantEffect,
    Substitution,
)
from varcode.effects.effect_classes import HaplotypeEffect
from varcode.phasing import MolecularPhaseResolver

ANNOTATORS = ("transcript_model", "fast", "protein_diff")

NCBI_CODES = (1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16,
              21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33)
CONTEXT_DEPENDENT = (27, 28, 31)
NON_STANDARD_CODES = tuple(c for c in NCBI_CODES if c != 1 and c not in CONTEXT_DEPENDENT)
FROM_GC_PRT = {
    32: (11, {"TAG": "W"}, ()),
    33: (24, {"TAA": "Y"}, ()),
    3: (3, {}, ("GTG",)),
    29: (29, {"CTG": "L"}, ()),
    30: (30, {"CTG": "L"}, ()),
}
CODONS = tuple("".join(b) for b in product("TCAG", repeat=3))


def ncbi_table(code: int):
    """Independent (codon -> residue_or_star, frozenset_of_starts) from Biotite + gc.prt."""
    base, changes, added_starts = FROM_GC_PRT.get(code, (code, {}, ()))
    table = BiotiteCodonTable.load(base)
    return {**table.codon_dict(), **changes}, frozenset(table.start_codons()) | set(added_starts)


def biotite_translate_cds(dna: str, code: int, *, initiate_as_met: bool = True) -> str:
    """Translate ``dna`` from codon 0 up to the first stop in NCBI ``code`` via Biotite."""
    residues, starts = ncbi_table(code)
    seq = dna.upper()
    out = []
    for idx in range(0, len(seq) // 3 * 3, 3):
        codon = seq[idx:idx + 3]
        aa = residues[codon]
        if aa == "*":
            break
        if idx == 0 and initiate_as_met and codon in starts:
            aa = "M"
        out.append(aa)
    return "".join(out)


def _single_nt_neighbors(codon: str):
    for i in range(3):
        for base in "ACGT":
            if base != codon[i]:
                yield i, codon[i], base, codon[:i] + base + codon[i + 1:]


def _distinguishing_case(code: int):
    """Pick a deterministic 1-nt nonsynonymous mutation whose consequence under ``code`` differs from code 1."""
    std_res, std_starts = ncbi_table(1)
    code_res, code_starts = ncbi_table(code)

    # Prefer a non-silent internal codon mutation where code vs standard code disagree.
    for ref_codon in CODONS:
        if code_res[ref_codon] == "*" or std_res[ref_codon] == "*":
            continue
        for offset, ref_nt, alt_nt, alt_codon in _single_nt_neighbors(ref_codon):
            if code_res[alt_codon] != std_res[alt_codon] and code_res[alt_codon] != code_res[ref_codon]:
                return ("internal", ref_codon, offset, ref_nt, alt_nt, alt_codon)

    # Codes like 11 share internal sense/stop assignments with code 1 and differ in start codons.
    extra_starts = sorted(code_starts - std_starts)
    for alt_start in extra_starts:
        for offset, ref_nt, alt_nt, cand in _single_nt_neighbors("ATG"):
            if cand == alt_start:
                return ("start", "ATG", offset, ref_nt, alt_nt, alt_start)
    raise AssertionError(f"No single-nt distinguishing mutation found for NCBI code {code}")


@pytest.fixture(scope="module")
def all_codes_genome(tmp_path_factory):
    """Custom indexed Genome with one annotated transcript per non-standard NCBI code."""
    work = tmp_path_factory.mktemp("ncbi-codes-genome")
    gtf_lines = []
    cdna_lines = []
    prot_lines = []

    for idx, code in enumerate(NON_STANDARD_CODES):
        code_res, _ = ncbi_table(code)
        kind, ref_codon, _, _, _, _ = _distinguishing_case(code)
        stop_codon = next(c for c in ("TGA", "TAG", "TAA", "AGA", "TCA", "TTA") if code_res[c] == "*")
        utr5 = "CCCCCC"
        body_codon = "GCT" if kind == "start" else ref_codon
        cds_sense = "ATG" + "GCT" + body_codon + "GGA" + "TTC"
        utr3 = "GGT" + "CCA" + stop_codon + "CCCCCC"
        tx_seq = utr5 + cds_sense + stop_codon + utr3
        protein = biotite_translate_cds(cds_sense + stop_codon, code)
        assert len(protein) == 5, (code, cds_sense, stop_codon, protein)

        start = idx * 200 + 1
        end = start + len(tx_seq) - 1
        cds_start = start + len(utr5)
        cds_end = cds_start + len(cds_sense) - 1
        stop_start = cds_end + 1
        stop_end = stop_start + 2
        gid, tid, pid = f"G_{code}", f"T_{code}", f"P_{code}"
        attrs = (
            f'gene_id "{gid}"; gene_name "{gid}"; gene_biotype "protein_coding"; '
            f'transcript_id "{tid}"; transcript_name "{tid}"; '
            f'transcript_biotype "protein_coding"; transl_table "{code}";'
        )
        gtf_lines.extend([
            f'1\ttest\tgene\t{start}\t{end}\t.\t+\t.\tgene_id "{gid}"; gene_name "{gid}"; gene_biotype "protein_coding";',
            f'1\ttest\ttranscript\t{start}\t{end}\t.\t+\t.\t{attrs}',
            f'1\ttest\texon\t{start}\t{end}\t.\t+\t.\t{attrs} exon_number "1"; exon_id "E_{code}";',
            f'1\ttest\tCDS\t{cds_start}\t{cds_end}\t.\t+\t0\t{attrs} protein_id "{pid}";',
            f'1\ttest\tstart_codon\t{cds_start}\t{cds_start + 2}\t.\t+\t0\t{attrs}',
            f'1\ttest\tstop_codon\t{stop_start}\t{stop_end}\t.\t+\t0\t{attrs}',
        ])
        cdna_lines.append(f">{tid}\n{tx_seq}")
        prot_lines.append(f">{pid}\n{protein}")

    gtf_path = work / "codes.gtf"
    cdna_path = work / "codes.cdna.fa"
    prot_path = work / "codes.pep.fa"
    gtf_path.write_text("\n".join(gtf_lines) + "\n")
    cdna_path.write_text("\n".join(cdna_lines) + "\n")
    prot_path.write_text("\n".join(prot_lines) + "\n")

    genome = Genome(
        reference_name="custom_codes",
        annotation_name="ncbi_codes",
        annotation_version="1",
        gtf_path_or_url=str(gtf_path),
        transcript_fasta_paths_or_urls=[str(cdna_path)],
        protein_fasta_paths_or_urls=[str(prot_path)],
        cache_directory_path=str(work / "cache"),
    )
    genome.index()
    return genome


@pytest.mark.parametrize("annotator", ANNOTATORS)
@pytest.mark.parametrize("code", NON_STANDARD_CODES)
def test_every_non_standard_ncbi_code_drives_transcript_effects(
        all_codes_genome, code, annotator):
    transcript = all_codes_genome.transcript_by_id(f"T_{code}")
    assert transcript.translation_table == code

    std_res, _ = ncbi_table(1)
    code_res, _ = ncbi_table(code)
    kind, ref_codon, codon_nt_offset, ref_nt, alt_nt, alt_codon = _distinguishing_case(code)

    cds_genomic_start = min(transcript.start_codon_positions)
    if kind == "start":
        pos = cds_genomic_start + codon_nt_offset
        variant = Variant("1", pos, ref_nt, alt_nt, all_codes_genome)
        effect = variant.effect_on_transcript(transcript, annotator=annotator)
        assert isinstance(effect, AlternateStartCodon)
        assert effect.ref_codon == "ATG" and effect.alt_codon == alt_codon
        return

    # Codon index 2 (0-based) in CDS ("ATG", "GCT", target, "GGA", "TTC")
    pos = cds_genomic_start + 6 + codon_nt_offset
    variant = Variant("1", pos, ref_nt, alt_nt, all_codes_genome)
    effect = variant.effect_on_transcript(transcript, annotator=annotator)

    mutated_cdna = (
        transcript.coding_sequence[:6 + codon_nt_offset]
        + alt_nt
        + transcript.coding_sequence[6 + codon_nt_offset + 1:]
        + transcript.sequence[transcript.last_stop_codon_spliced_offset - 2:]
    )
    expected_protein = biotite_translate_cds(mutated_cdna, code)
    standard_protein = biotite_translate_cds(mutated_cdna, 1)

    expected_aa = code_res[alt_codon]
    ref_aa = code_res[ref_codon]
    if expected_aa == "*":
        assert isinstance(effect, PrematureStop)
        assert effect.mutant_protein_sequence == expected_protein == "MA"
        assert standard_protein != expected_protein
    elif expected_aa == ref_aa:
        assert isinstance(effect, Silent)
        assert effect.mutant_protein_sequence == expected_protein
        assert std_res[alt_codon] != std_res[ref_codon]
    else:
        assert isinstance(effect, Substitution)
        assert effect.aa_ref == ref_aa
        assert effect.aa_alt == expected_aa
        assert effect.mutant_protein_sequence == expected_protein
        assert expected_aa != std_res[alt_codon] or ref_aa != std_res[ref_codon]


@pytest.mark.parametrize("annotator", ANNOTATORS)
def test_mt_co1_frameshifts_use_vertebrate_mitochondrial_code(annotator):
    """Frameshifts on human MT-CO1 read TGA as Trp and terminate at AGA/AGG (#701)."""
    genome = cached_release(81)
    transcript = genome.genes_by_name("MT-CO1")[0].transcripts[0]
    cds = str(transcript.sequence[transcript.first_start_codon_spliced_offset:])

    readthrough_case = None
    early_aga_case = None
    for cds_idx in range(6, len(cds) - 30):
        anchor = cds[cds_idx - 1]
        deleted = cds[cds_idx]
        if anchor != deleted and readthrough_case is None:
            mut_cds = cds[:cds_idx] + cds[cds_idx + 1:]
            prot_code2 = biotite_translate_cds(mut_cds, 2)
            prot_code1 = biotite_translate_cds(mut_cds, 1)
            if len(prot_code2) > len(prot_code1) + 2:
                readthrough_case = (cds_idx - 1, cds[cds_idx - 1:cds_idx + 1], anchor, prot_code2, prot_code1)
        for ins_base in "ACGT":
            if ins_base != anchor and ins_base != deleted and early_aga_case is None:
                mut_cds = cds[:cds_idx] + ins_base + cds[cds_idx:]
                prot_code2 = biotite_translate_cds(mut_cds, 2)
                prot_code1 = biotite_translate_cds(mut_cds, 1)
                if len(prot_code2) + 2 < len(prot_code1):
                    early_aga_case = (cds_idx - 1, anchor, anchor + ins_base, prot_code2, prot_code1)
        if readthrough_case and early_aga_case:
            break

    assert readthrough_case is not None and early_aga_case is not None

    for anchor_idx, ref, alt, expected_code2, wrong_code1 in (readthrough_case, early_aga_case):
        genomic_pos = transcript.start_codon_positions[0] + anchor_idx
        variant = Variant("MT", genomic_pos, ref, alt, genome)
        effect = variant.effect_on_transcript(transcript, annotator=annotator)
        assert isinstance(effect, (FrameShift, FrameShiftTruncation))
        assert effect.mutant_protein_sequence == expected_code2
        assert effect.mutant_protein_sequence != wrong_code1


@pytest.fixture(scope="module")
def custom_organelle_genome(tmp_path_factory):
    """Indexed Genome with chromosome FASTA for stop-loss, frameshift, and SV tests."""
    work = tmp_path_factory.mktemp("organelle-sv")

    t_mt1_utr5 = "CCCCCC"
    t_mt1_cds = "ATGGCTGAACCAGGATAA"  # M A E P G *
    t_mt1_utr3 = "TGAGCTAGACCCCCC"     # code 2: W A * ; code 1: * A R P P
    t_mt1_seq = t_mt1_utr5 + t_mt1_cds + t_mt1_utr3  # 39 bp

    t_mt2_utr5 = "TTTTTT"
    t_mt2_cds = "ATGTGAAAGATAGGAAGG"  # code 2: M W K M G *
    t_mt2_utr3 = "CCCCCC"
    t_mt2_seq = t_mt2_utr5 + t_mt2_cds + t_mt2_utr3  # 30 bp

    mt_chrom = "A" * 20 + t_mt1_seq + "A" * (150 - 20 - len(t_mt1_seq)) + t_mt2_seq + "A" * 40
    assert mt_chrom[20:20 + len(t_mt1_seq)] == t_mt1_seq
    assert mt_chrom[150:150 + len(t_mt2_seq)] == t_mt2_seq

    t_c6_seq = "CCCCCC" + "ATGGCTTGA" + "TAATAGTGACCC"
    t_c14_seq = "CCCCCC" + "ATGGCTTAG" + "TAATGAAAATAGCCC"
    t_c22_seq = "CCCCCC" + "ATGGCTCCAGGATAA" + "TTAGTCAGGGCCCCC"

    chr1_chrom = "A" * 20 + t_c6_seq + "A" * 30 + t_c14_seq + "A" * 30 + t_c22_seq + "A" * 40

    fasta_path = work / "genome.fa"
    fasta_path.write_text(f">MT\n{mt_chrom}\n>1\n{chr1_chrom}\n")

    specs = [
        ("MT", 21, t_mt1_seq, 6, 15, 2, "G_mt1", "T_mt1", "P_mt1"),
        ("MT", 151, t_mt2_seq, 6, 15, 2, "G_mt2", "T_mt2", "P_mt2"),
        ("1", 21, t_c6_seq, 6, 6, 6, "G_c6", "T_c6", "P_c6"),
        ("1", 21 + len(t_c6_seq) + 30, t_c14_seq, 6, 6, 14, "G_c14", "T_c14", "P_c14"),
        ("1", 21 + len(t_c6_seq) + 30 + len(t_c14_seq) + 30, t_c22_seq, 6, 12, 22, "G_c22", "T_c22", "P_c22"),
    ]
    gtf_lines, cdna_lines, prot_lines = [], [], []
    for contig, start, seq, utr5_len, cds_sense_len, table_id, gid, tid, pid in specs:
        end = start + len(seq) - 1
        cds_start = start + utr5_len
        cds_end = cds_start + cds_sense_len - 1
        stop_start = cds_end + 1
        stop_end = stop_start + 2
        prot = biotite_translate_cds(seq[utr5_len:utr5_len + cds_sense_len + 3], table_id)
        attrs = (
            f'gene_id "{gid}"; gene_name "{gid}"; gene_biotype "protein_coding"; '
            f'transcript_id "{tid}"; transcript_name "{tid}"; '
            f'transcript_biotype "protein_coding"; transl_table "{table_id}";'
        )
        gtf_lines.extend([
            f'{contig}\ttest\tgene\t{start}\t{end}\t.\t+\t.\tgene_id "{gid}"; gene_name "{gid}"; gene_biotype "protein_coding";',
            f'{contig}\ttest\ttranscript\t{start}\t{end}\t.\t+\t.\t{attrs}',
            f'{contig}\ttest\texon\t{start}\t{end}\t.\t+\t.\t{attrs} exon_number "1"; exon_id "E_{tid}";',
            f'{contig}\ttest\tCDS\t{cds_start}\t{cds_end}\t.\t+\t0\t{attrs} protein_id "{pid}";',
            f'{contig}\ttest\tstart_codon\t{cds_start}\t{cds_start + 2}\t.\t+\t0\t{attrs}',
            f'{contig}\ttest\tstop_codon\t{stop_start}\t{stop_end}\t.\t+\t0\t{attrs}',
        ])
        cdna_lines.append(f">{tid}\n{seq}")
        prot_lines.append(f">{pid}\n{prot}")

    gtf_path = work / "organelle.gtf"
    cdna_path = work / "organelle.cdna.fa"
    prot_path = work / "organelle.pep.fa"
    gtf_path.write_text("\n".join(gtf_lines) + "\n")
    cdna_path.write_text("\n".join(cdna_lines) + "\n")
    prot_path.write_text("\n".join(prot_lines) + "\n")

    py_genome = Genome(
        reference_name="custom_organelle",
        annotation_name="organelle_test",
        annotation_version="1",
        gtf_path_or_url=str(gtf_path),
        transcript_fasta_paths_or_urls=[str(cdna_path)],
        protein_fasta_paths_or_urls=[str(prot_path)],
        cache_directory_path=str(work / "cache"),
    )
    py_genome.index()
    from varcode import Genome as VarcodeGenome
    return VarcodeGenome(py_genome, fasta=str(fasta_path))


@pytest.mark.parametrize("annotator", ANNOTATORS)
@pytest.mark.parametrize("tid,code,alt_nt,expected_extension", [
    ("T_mt1", 2, "C", "QWA"),     # TAA->CAA (Q), reads TGA (W), GCT (A), stops at AGA
    ("T_c6", 6, "C", "CQQ"),      # TGA->TGC (C), reads TAA (Q), TAG (Q), stops at TGA
    ("T_c14", 14, "C", "YYWN"),   # TAG->TAC (Y), reads TAA (Y), TGA (W), AAA (N), stops at TAG
])
def test_stop_loss_readthrough_uses_transcript_genetic_code(
        custom_organelle_genome, annotator, tid, code, alt_nt, expected_extension):
    transcript = custom_organelle_genome.transcript_by_id(tid)
    stop_positions = transcript.stop_codon_positions
    idx = 0 if tid == "T_mt1" else 2
    pos = stop_positions[idx]
    ref_nt = str(transcript.sequence[transcript.spliced_offset(pos)])
    variant = Variant(transcript.contig, pos, ref_nt, alt_nt, custom_organelle_genome)
    effect = variant.effect_on_transcript(transcript, annotator=annotator)

    mut_cdna = (
        str(transcript.sequence[transcript.first_start_codon_spliced_offset:transcript.spliced_offset(pos)])
        + alt_nt
        + str(transcript.sequence[transcript.spliced_offset(pos) + 1:])
    )
    expected_protein = biotite_translate_cds(mut_cdna, code)
    standard_protein = biotite_translate_cds(mut_cdna, 1)

    assert isinstance(effect, StopLoss)
    assert effect.aa_alt == expected_extension
    assert effect.mutant_protein_sequence == expected_protein
    assert expected_protein != standard_protein


@pytest.mark.parametrize("annotator", ANNOTATORS)
def test_code_22_frameshift_reads_tag_as_leucine_and_stops_at_tca(
        custom_organelle_genome, annotator):
    transcript = custom_organelle_genome.transcript_by_id("T_c22")
    cds_start = min(transcript.start_codon_positions)
    variant = Variant("1", cds_start + 5, "TC", "T", custom_organelle_genome)
    effect = variant.effect_on_transcript(transcript, annotator=annotator)

    cds_and_utr3 = str(transcript.sequence[transcript.first_start_codon_spliced_offset:])
    mut_cds = cds_and_utr3[:6] + cds_and_utr3[7:]
    expected_code22 = biotite_translate_cds(mut_cds, 22)
    wrong_code1 = biotite_translate_cds(mut_cds, 1)

    assert isinstance(effect, FrameShift)
    assert effect.mutant_protein_sequence == expected_code22 == "MAQDNL"
    assert wrong_code1 == "MAQDN"


def test_non_identity_mitochondrial_structural_variants_and_composed_haplotypes(
        custom_organelle_genome):
    t1 = custom_organelle_genome.transcript_by_id("T_mt1")
    t2 = custom_organelle_genome.transcript_by_id("T_mt2")

    # 1. Local in-frame deletion on T_mt2 removing Lys3 (AAG, genomic 163..165),
    # retaining TGA (Trp2) and ATA (Met4) and terminating at AGG: M W M G.
    sv_del = StructuralVariant("MT", 163, "DEL", end=165, genome=custom_organelle_genome)
    del_effect = sv_del.effect_on_transcript(t2)
    assert isinstance(del_effect, Deletion)
    expected_del = biotite_translate_cds("ATGTGAATAGGAAGG", 2)
    assert expected_del == "MWMG"
    assert biotite_translate_cds("ATGTGAATAGGAAGG", 1) == "M"
    assert del_effect.mutant_protein_sequence == expected_del

    # 2. Local in-frame tandem duplication on T_mt2 duplicating TGA AAG ATA (genomic 160..168).
    sv_dup = StructuralVariant("MT", 160, "DUP", end=168, genome=custom_organelle_genome)
    dup_effect = sv_dup.effect_on_transcript(t2)
    expected_dup = biotite_translate_cds("ATGTGAAAGATATGAAAGATAGGAAGG", 2)
    assert expected_dup == "MWKMWKMG"
    assert biotite_translate_cds("ATGTGAAAGATATGAAAGATAGGAAGG", 1) == "M"
    assert dup_effect.mutant_protein_sequence == expected_dup

    # 3. Local frameshifting SV deletion on T_mt1 deleting 2 bp (37..38) so the
    # shifted frame translates ATA as Met under code 2 (MAERMMS) vs Ile under code 1 (MAERIMS).
    sv_fs = StructuralVariant("MT", 37, "DEL", end=38, genome=custom_organelle_genome)
    fs_eff = sv_fs.effect_on_transcript(t1)
    expected_fs = biotite_translate_cds("ATGGCTGAACGGATAATGAGCTAGACCCCCC", 2)
    assert expected_fs == "MAERMMS"
    assert biotite_translate_cds("ATGGCTGAACGGATAATGAGCTAGACCCCCC", 1) == "MAERIMS"
    assert fs_eff.mutant_protein_sequence == expected_fs

    # 4. Local in-frame stop-loss SV deletion on T_mt1 deleting TAA (42..44)
    # so translation reads through 3' UTR TGA (Trp) GCT (Ala) and stops at AGA (*).
    sv_stoploss = StructuralVariant("MT", 42, "DEL", end=44, genome=custom_organelle_genome)
    stoploss_eff = sv_stoploss.effect_on_transcript(t1)
    expected_stoploss = biotite_translate_cds("ATGGCTGAACCAGGATGAGCTAGACCCCCC", 2)
    assert expected_stoploss == "MAEPGWA"
    assert biotite_translate_cds("ATGGCTGAACCAGGATGAGCTAGACCCCCC", 1) == "MAEPG"
    assert stoploss_eff.mutant_protein_sequence == expected_stoploss

    # 5. Large mitochondrial deletion fusing 5' of T_mt1 (ATG GCT, 27..32; VCF anchor at 32)
    # in-frame to T_mt2 at TGA (160..174: TGA AAG ATA GGA AGG) -> M A W K M G.
    sv_fusion = StructuralVariant("MT", 32, "DEL", end=159, genome=custom_organelle_genome)
    fusion_effects = [e for e in sv_fusion.effects() if isinstance(e, GeneFusion)]
    assert len(fusion_effects) >= 1
    expected_fusion = biotite_translate_cds("ATGGCTTGAAAGATAGGAAGG", 2)
    assert expected_fusion == "MAWKMG"
    assert biotite_translate_cds("ATGGCTTGAAAGATAGGAAGG", 1) == "MA"
    for fusion_eff in fusion_effects:
        assert fusion_eff.mutant_protein_sequence == expected_fusion

    # 6. Composed cis HaplotypeEffect on T_mt2: mutate AAG (Lys3, 163..165) via two cis SNVs
    # 163 A->T and 164 A->G -> TGG (Trp), preserving TGA (Trp2) and ATA (Met4) -> M W W M G.
    v1 = Variant("MT", 163, "A", "T", custom_organelle_genome)
    v2 = Variant("MT", 164, "A", "G", custom_organelle_genome)

    class _CisSource:
        source = "test_cis"

        def has_evidence(self, variant):
            return variant in (v1, v2)

        def partners_in_cis(self, variant):
            return (v2,) if variant == v1 else (v1,)

    resolver = MolecularPhaseResolver(_CisSource())
    expected_hap = biotite_translate_cds("ATGTGGTGGATAGGAAGG", 2)
    assert expected_hap == "MWWMG"
    for hap_annotator, expected_cls in (("transcript_model", Substitution), ("fast", HaplotypeEffect)):
        coll_effects = VariantCollection([v1, v2]).effects(
            annotator=hap_annotator, phase_resolver=resolver)
        joint = [e for e in coll_effects if getattr(e, "variants", None) == (v1, v2)]
        assert len(joint) == 1
        assert isinstance(joint[0], expected_cls)
        assert joint[0].mutant_protein_sequence == expected_hap


def test_human_mitochondrial_poly_a_and_att_initiator_structural_models():
    """Poly(A)-terminated MT genes and ATT initiators preserve code-2 flags and fusions (#701)."""
    genome = cached_release(81)
    nd1 = genome.genes_by_name("MT-ND1")[0].transcripts[0]
    nd2 = genome.genes_by_name("MT-ND2")[0].transcripts[0]

    # Identity segment model on polyadenylated MT-ND1 (ends in TA, complete=False) is unchanged.
    sv_bnd = StructuralVariant("MT", nd1.start, "BND", genome=genome)
    identity = StructuralVariantEffect(
        sv_bnd,
        nd1,
        mutant_transcript=MutantTranscript(
            reference_transcript=nd1,
            reference_segments=(ReferenceSegment(nd1, 0, len(nd1.sequence)),),
            cdna_sequence=nd1.sequence,
        ),
    )
    assert identity.modifies_coding_sequence is False
    assert identity.modifies_protein_sequence is False

    # Partial fragment not reaching the 3' end of MT-ND1 does not invent a poly(A) stop.
    partial = StructuralVariantEffect(
        sv_bnd,
        nd1,
        mutant_transcript=MutantTranscript(
            reference_transcript=nd1,
            reference_segments=(ReferenceSegment(nd1, 0, 90),),
            cdna_sequence=nd1.sequence[:90],
        ),
    )
    assert partial.modifies_coding_sequence is None
    assert partial.modifies_protein_sequence is None

    # Local in-frame duplication inside MT-ND1 translates with code 2 + poly(A) stop.
    dup = StructuralVariant("MT", nd1.start + 6, "DUP", end=nd1.start + 17, genome=genome)
    dup_effect = dup.effect_on_transcript(nd1)
    expected_dup_cdna = nd1.sequence[:18] + nd1.sequence[6:18] + nd1.sequence[18:] + "A"
    assert dup_effect.mutant_protein_sequence == biotite_translate_cds(expected_dup_cdna, 2)
    assert dup_effect.modifies_coding_sequence is True
    assert dup_effect.modifies_protein_sequence is True

    # Fusion from MT-ND2 (ATT start codon -> Met) into MT-ND1 (polyadenylated TA -> TAA stop).
    nd2_keep_pos = nd2.start + 11  # first 4 codons (12 bp) of MT-ND2
    nd1_join_pos = nd1.end - 13    # last 14 bp (4 codons + TA) of MT-ND1
    fusion_sv = StructuralVariant(
        "MT",
        nd2_keep_pos,
        "BND",
        alt=f"N[MT:{nd1_join_pos}[",
        mate_contig="MT",
        mate_start=nd1_join_pos,
        genome=genome,
    )
    fusion_eff = fusion_sv.effect_on_transcript(nd2)
    fusions = [c.effect for c in getattr(fusion_eff, "candidates", ()) if isinstance(c.effect, GeneFusion)]
    if isinstance(fusion_eff, GeneFusion):
        fusions.insert(0, fusion_eff)
    assert fusions
    expected_fused = biotite_translate_cds(nd2.sequence[:12] + nd1.sequence[-14:] + "A", 2)
    assert expected_fused.startswith("M")
    assert not biotite_translate_cds(nd2.sequence[:12] + nd1.sequence[-14:] + "A", 2, initiate_as_met=False).startswith("M")
    assert fusions[0].mutant_protein_sequence == expected_fused
    assert fusions[0].modifies_coding_sequence is True
    assert fusions[0].modifies_protein_sequence is True


@pytest.mark.parametrize("annotator", ANNOTATORS)
def test_human_mt_co1_and_mt_nd6_terminal_aga_agg_follow_ncbi_table_2(annotator):
    """Documented scope: MT-CO1 (AGA) and MT-ND6 (AGG) use NCBI table 2's static stop set (#701)."""
    genome = cached_release(81)
    co1 = genome.genes_by_name("MT-CO1")[0].transcripts[0]
    nd6 = genome.genes_by_name("MT-ND6")[0].transcripts[0]

    co1_stop_syn = Variant("MT", 7445, "A", "G", genome).effect_on_transcript(co1, annotator=annotator)
    co1_stop_loss = Variant("MT", 7444, "G", "A", genome).effect_on_transcript(co1, annotator=annotator)
    co1_pre_stop = Variant("MT", 7442, "T", "C", genome).effect_on_transcript(co1, annotator=annotator)
    assert isinstance(co1_stop_syn, Silent)
    assert isinstance(co1_stop_loss, StopLoss) and co1_stop_loss.short_description == "p.*514K (stop-loss)"
    assert isinstance(co1_pre_stop, Silent) and co1_pre_stop.short_description == "p.S513="

    assert reverse_complement("CCTATT") == "AATAGG"
    nd6_stop_syn = Variant("MT", 14149, "C", "T", genome).effect_on_transcript(nd6, annotator=annotator)
    nd6_stop_loss = Variant("MT", 14151, "T", "C", genome).effect_on_transcript(nd6, annotator=annotator)
    nd6_pre_stop = Variant("MT", 14152, "A", "G", genome).effect_on_transcript(nd6, annotator=annotator)
    assert isinstance(nd6_stop_syn, Silent)
    assert isinstance(nd6_stop_loss, StopLoss) and nd6_stop_loss.short_description == "p.*175G (stop-loss)"
    assert isinstance(nd6_pre_stop, Silent) and nd6_pre_stop.short_description == "p.N174="
