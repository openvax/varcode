# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Gene-fusion regression tests from a real osteosarcoma genome.

Three somatic rearrangements from Sid Sijbrandij's public osteosarcoma
dataset (https://osteosarc.com). Each is a single DNA adjacency (one
LINX SV, not a chained event) with split reads in tumor DNA and none in
matched blood. Breakpoints, orientations and transcripts come from the
nf-core/oncoanalyser LINX output (``*.linx.neoepitope.tsv``,
``*.linx.fusion.tsv``, ``*.linx.breakend.tsv``) in the dataset's S3
bucket. Coordinates are hg38 with ``chr`` dropped to match Ensembl
contig names.

Expected proteins were checked three ways: varcode's output, LINX's
independent exon and phase call for the same junction, and a
from-scratch translation (plain exon walk plus a standard codon table,
no varcode helpers). Assertions are on proteins and exon boundaries
rather than raw cDNA offsets, which shift with UTR edits between
Ensembl releases.

None of these junctions has more than one RNA split read, so the
proteins are predicted products of DNA rearrangements, not observed
transcripts.

The tests at the end load the same junctions from
``tests/data/osteosarc_esvee_somatic.vcf``: the esvee records as
published (UCSC ``chr`` contig names) plus two single breakends, with
only the FORMAT and sample columns removed. esvee labels the
CPEB2::FAM193A pair a duplication and OTX1::KIF3C an inversion;
``pair_breakends`` turns each pair into that event, and both the
records and the paired event must still come out as fusions.
"""

import pytest
from pyensembl import cached_release

from varcode import StructuralVariant, StructuralVariantAnnotator, load_vcf
from varcode.annotators.structural_variant import (
    _build_fusion_mutant_transcript,
)
from varcode.effects import GeneFusion, Intronic
from varcode.transforms import pair_breakends

from .data import data_path

# Release 95 rather than the suite's usual 81: LINX's FAM193A
# transcript (ENST00000637812) isn't in 81. CI installs both.
ensembl_grch38 = cached_release(95)
_ANNOTATOR = StructuralVariantAnnotator()


def _sum_exon_lengths(transcript, exon_numbers):
    return sum(
        exon.end - exon.start + 1
        for i, exon in enumerate(transcript.exons, start=1)
        if i in exon_numbers)


def test_otx1_kif3c_out_of_frame_fusion():
    """OTX1 (chr2, forward) intron 3 joined to KIF3C (chr2, reverse)
    intron 2, from the T2 biopsy. The junction carries a 2 bp
    untemplated insert (TC) that splices out with the intron.

    LINX: OTX1 exon 3 -> KIF3C exon 3, OUT_OF_FRAME (phases 1 / 0).
    OTX1 keeps 32 residues, then the shifted frame reads three KIF3C
    residues before a stop.
    """
    otx1 = ensembl_grch38.transcript_by_id("ENST00000282549")
    kif3c = ensembl_grch38.transcript_by_id("ENST00000264712")
    sv = StructuralVariant(
        contig="2",
        start=63_053_516,
        sv_type="BND",
        alt="N]2:25955847]",
        mate_contig="2",
        mate_start=25_955_847,
        mate_orientation="]]",
        genome=ensembl_grch38)
    effect = _ANNOTATOR.annotate_on_transcript(sv, otx1)
    assert isinstance(effect, GeneFusion)
    # The annotator picks the same 3' transcript LINX reported.
    assert effect.partner_transcript.id == kif3c.id

    mt = effect.mutant_transcript
    assert mt.reference_segments[0].end == _sum_exon_lengths(
        otx1, {1, 2, 3})
    assert mt.reference_segments[1].start == _sum_exon_lengths(
        kif3c, {1, 2})

    protein = mt.mutant_protein_sequence
    assert protein == "MMSYLKQPPYGMNGLGLAGPAMDLLHPSVGYPETS"
    assert protein[:32] == otx1.protein_sequence[:32]


def test_sema6a_parm1_reverse_strand_5p_partner():
    """SEMA6A (chr5, reverse) intron 16 joined to PARM1 (chr4,
    forward) intron 3, from the T0 resection.

    LINX's breakend record puts the SEMA6A donor at exon 16 and the
    PARM1 acceptor at exon 4, phases 1 / 2: out of frame.
    ``linx.fusion.tsv`` lists exon 15 instead because LINX skips exon
    16 to find an in-frame alternative (SKIPPED_EXONS); the direct
    junction is what's modeled here.
    """
    sema6a = ensembl_grch38.transcript_by_id("ENST00000343348")
    parm1 = ensembl_grch38.transcript_by_id("ENST00000307428")
    sv = StructuralVariant(
        contig="5",
        start=116_474_281,
        sv_type="BND",
        alt="[4:75037126[N",
        mate_contig="4",
        mate_start=75_037_126,
        mate_orientation="[[",
        genome=ensembl_grch38)
    effect = _ANNOTATOR.annotate_on_transcript(sv, sema6a)
    assert isinstance(effect, GeneFusion)
    assert effect.partner_transcript.id == parm1.id

    mt = effect.mutant_transcript
    assert mt.reference_segments[0].end == _sum_exon_lengths(
        sema6a, set(range(1, 17)))
    assert mt.reference_segments[1].start == _sum_exon_lengths(
        parm1, {1, 2, 3})

    protein = mt.mutant_protein_sequence
    assert len(protein) == 623
    assert protein[:569] == sema6a.protein_sequence[:569]
    # PARM1 read in the shifted frame until the first stop.
    assert protein[569:] == (
        "SIPPMEDFWTTMTTGPGETTTTLCTMTPNNGIWPGMRINCSLFISAYPVELIST")


def test_cpeb2_fam193a_in_frame_fusion():
    """CPEB2 (chr4, forward) intron 3 joined to FAM193A (chr4,
    forward) intron 12, from the T1 organoid.

    LINX: CPEB2 exon 3 -> FAM193A exon 13 on ENST00000637812, INFRAME
    (phases 0 / 0), so the product is CPEB2's first 678 residues
    followed by FAM193A from residue 694 on.

    The annotator's partner lookup takes the first coding transcript
    at the mate locus, a different FAM193A isoform here, so the LINX
    transcript is passed to the builder directly.
    """
    cpeb2 = ensembl_grch38.transcript_by_id("ENST00000538197")
    fam193a = ensembl_grch38.transcript_by_id("ENST00000637812")
    sv = StructuralVariant(
        contig="4",
        start=15_012_987,
        sv_type="BND",
        alt="N[4:2664478[",
        mate_contig="4",
        mate_start=2_664_478,
        mate_orientation="[[",
        genome=ensembl_grch38)

    effect = _ANNOTATOR.annotate_on_transcript(sv, cpeb2)
    assert isinstance(effect, GeneFusion)
    assert effect.partner_transcript.gene_name == "FAM193A"

    mt = _build_fusion_mutant_transcript(
        cpeb2, cpeb2, 15_012_987, fam193a, 2_664_478)
    assert mt.reference_segments[0].end == _sum_exon_lengths(
        cpeb2, {1, 2, 3})
    assert mt.reference_segments[1].start == _sum_exon_lengths(
        fam193a, set(range(1, 13)))

    protein = mt.mutant_protein_sequence
    assert protein == (
        cpeb2.protein_sequence[:678] + fam193a.protein_sequence[693:])
    assert len(protein) == 1500
    assert protein[668:688] == "TAAGTSRIDQAEQAPNTCEC"


def _esvee_svs():
    return load_vcf(
        data_path("osteosarc_esvee_somatic.vcf"),
        genome=ensembl_grch38,
        parse_structural_variants=True)


def _esvee_sv(record_id):
    vc = _esvee_svs()
    (sv,) = [v for v in vc if vc.metadata[v]["id"] == record_id]
    return sv


def _fusion_on(sv, transcript_id):
    (fusion,) = [
        effect for effect in sv.effects()
        if isinstance(effect, GeneFusion)
        and effect.transcript.id == transcript_id]
    return fusion


def test_esvee_vcf_loads_every_record():
    """All eight records load, single breakends included, with Ensembl
    contig names and the genome given to load_vcf. Each record loads as
    the breakend it is, at its own position, whatever SVTYPE its caller
    labeled the pair with."""
    vc = _esvee_svs()
    assert len(vc) == 8
    assert all(isinstance(v, StructuralVariant) for v in vc)
    assert {v.genome for v in vc} == {ensembl_grch38}
    parsed = {
        vc.metadata[v]["id"]: (
            v.sv_type, v.contig, v.start, v.end, v.mate_contig)
        for v in vc}
    assert parsed == {
        "4055": ("BND", "2", 25_955_847, 25_955_847, "2"),
        "4571": ("BND", "2", 63_053_516, 63_053_516, "2"),
        "11585": ("BND", "4", 2_664_478, 2_664_478, "4"),
        "11828": ("BND", "4", 15_012_987, 15_012_987, "4"),
        "11309": ("BND", "4", 75_037_126, 75_037_126, "5"),
        "15093": ("BND", "5", 116_474_281, 116_474_281, "4"),
        "30830": ("BND", "12", 72_342_468, 72_342_468, None),
        "32579": ("BND", "13", 56_435_058, 56_435_058, None),
    }
    assert {v.info.get("svtype") for v in vc} == {"INV", "DUP", "BND", "SGL"}
    singles = [v for v in vc if v.mate_contig is None]
    assert all(v.info["single_breakend"] for v in singles)


def test_esvee_vcf_breakend_records_annotate():
    """Every record, single breakends included, goes through
    ``effects()``. (The paired DUP and INV span tens of megabases, so
    they're checked transcript by transcript below.)"""
    for sv in _esvee_svs():
        assert len(sv.effects()) > 0


def test_esvee_vcf_without_sv_parsing_skips_records():
    """Without ``parse_structural_variants`` the single breakends are
    skipped along with the other SV rows instead of failing nucleotide
    validation."""
    with pytest.warns(UserWarning, match="Skipped 8 symbolic/breakend"):
        vc = load_vcf(
            data_path("osteosarc_esvee_somatic.vcf"), genome=ensembl_grch38)
    assert len(vc) == 0


def test_esvee_vcf_otx1_kif3c_inversion_is_fusion_on_both_partners():
    """Either record of the INV-labeled OTX1::KIF3C junction, and the
    inversion they pair into, give the same fusion on OTX1 and on
    KIF3C: the 35-residue out-of-frame protein from the direct test
    above."""
    otx1 = ensembl_grch38.transcript_by_id("ENST00000282549")
    kif3c = ensembl_grch38.transcript_by_id("ENST00000264712")
    (inversion,) = [
        v for v in pair_breakends(_esvee_svs()) if v.sv_type == "INV"]
    for sv in (_esvee_sv("4055"), _esvee_sv("4571"), inversion):
        for transcript in (otx1, kif3c):
            fusion = _ANNOTATOR.annotate_on_transcript(sv, transcript)
            assert isinstance(fusion, GeneFusion)
            assert fusion.five_prime_transcript.id == otx1.id
            assert fusion.three_prime_transcript.id == kif3c.id
            assert fusion.mutant_transcript.mutant_protein_sequence == (
                "MMSYLKQPPYGMNGLGLAGPAMDLLHPSVGYPETS")


def test_esvee_vcf_cpeb2_fam193a_duplication_is_fusion():
    """esvee labels CPEB2::FAM193A a 12 Mb tandem duplication. CPEB2
    holds its right end and FAM193A its left, so once paired each gets a
    fusion with CPEB2 as 5' partner rather than an internal duplication
    of its own exons."""
    cpeb2 = ensembl_grch38.transcript_by_id("ENST00000538197")
    fam193a = ensembl_grch38.transcript_by_id("ENST00000637812")
    (sv,) = [v for v in pair_breakends(_esvee_svs()) if v.sv_type == "DUP"]
    assert (sv.start, sv.end) == (2_664_477, 15_012_987)

    on_cpeb2 = _ANNOTATOR.annotate_on_transcript(sv, cpeb2)
    assert isinstance(on_cpeb2, GeneFusion)
    assert on_cpeb2.five_prime_transcript.id == cpeb2.id
    assert on_cpeb2.three_prime_transcript.gene_name == "FAM193A"
    protein = on_cpeb2.mutant_transcript.mutant_protein_sequence
    assert protein[:678] == cpeb2.protein_sequence[:678]

    on_fam193a = _ANNOTATOR.annotate_on_transcript(sv, fam193a)
    assert isinstance(on_fam193a, GeneFusion)
    assert on_fam193a.three_prime_transcript.id == fam193a.id
    assert on_fam193a.five_prime_transcript.gene_name == "CPEB2"


def test_esvee_vcf_sema6a_parm1_fusion_from_either_record():
    """Each breakend record annotates its own gene: the SEMA6A record
    (through ``effects()``) as 5' partner, the PARM1 record as 3'
    partner, with the same fused protein."""
    parm1 = ensembl_grch38.transcript_by_id("ENST00000307428")
    on_sema6a = _fusion_on(_esvee_sv("15093"), "ENST00000343348")
    on_parm1 = _ANNOTATOR.annotate_on_transcript(_esvee_sv("11309"), parm1)
    for fusion in (on_sema6a, on_parm1):
        assert isinstance(fusion, GeneFusion)
        assert fusion.five_prime_transcript.gene_name == "SEMA6A"
        assert fusion.three_prime_transcript.id == parm1.id
        protein = fusion.mutant_transcript.mutant_protein_sequence
        assert len(protein) == 623
        assert protein.endswith(
            "SIPPMEDFWTTMTTGPGETTTTLCTMTPNNGIWPGMRINCSLFISAYPVELIST")


def test_esvee_vcf_pairs_into_one_variant_per_junction():
    """pair_breakends collapses each breakend pair, keeping the DUP and
    INV types and spans; the single breakends pass through."""
    paired = pair_breakends(_esvee_svs())
    assert sorted((v.sv_type, v.contig, v.start, v.end) for v in paired) == [
        ("BND", "12", 72_342_468, 72_342_468),
        ("BND", "13", 56_435_058, 56_435_058),
        ("BND", "4", 75_037_126, 75_037_126),
        ("DUP", "4", 2_664_477, 15_012_987),
        ("INV", "2", 25_955_847, 63_053_516),
    ]


def _load_esvee_records(tmp_path, records):
    """Load a few esvee records (sites only) from the dataset."""
    path = tmp_path / "esvee_records.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"type\">\n"
        "##INFO=<ID=MATEID,Number=1,Type=String,Description=\"mate\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join("\t".join(record) + "\n" for record in records))
    return load_vcf(
        str(path), genome=ensembl_grch38, parse_structural_variants=True)


def test_esvee_deletion_within_cbx5_is_not_a_fusion(tmp_path):
    """A 3.3 kb deletion inside CBX5 (T1 organoid, records 34402/34403).
    CBX5's long isoform holds both ends of the paired deletion; two
    short isoforms hold one end each, and what lies at the other end is
    CBX5 itself, so no isoform reports a fusion."""
    vc = _load_esvee_records(tmp_path, [
        ("chr12", "54259010", "34402", "T", "T[chr12:54262307[", "57",
         "PASS", "SVTYPE=DEL;MATEID=34403"),
        ("chr12", "54262307", "34403", "C", "]chr12:54259010]C", "57",
         "PASS", "SVTYPE=DEL;MATEID=34402"),
    ])
    (sv,) = pair_breakends(vc)
    assert (sv.sv_type, sv.start, sv.end) == ("DEL", 54_259_010, 54_262_306)
    long_isoform = ensembl_grch38.transcript_by_id("ENST00000209875")
    assert isinstance(
        _ANNOTATOR.annotate_on_transcript(sv, long_isoform), Intronic)
    for transcript_id in ("ENST00000439541", "ENST00000550411"):
        effect = _ANNOTATOR.annotate_on_transcript(
            sv, ensembl_grch38.transcript_by_id(transcript_id))
        assert not isinstance(effect, GeneFusion), transcript_id


def test_esvee_insertion_breakends_are_not_a_fusion(tmp_path):
    """esvee writes an insertion in COPB2 at chr3:139358615 as two
    breakends a base apart (T1 organoid, records 10440/10441,
    SVTYPE=INS), so they load as BNDs. The mate lies in COPB2 too, so
    no COPB2 transcript reports a fusion with itself."""
    insert = (
        "AATTAGCCTCAGCACCGAGAACGAATTGTATATGGGCGGGCAGGGGAGAAGATGGGCAGAGAAATG"
        "AGAATGAAGGGGAGTTCTTGGGATTCTGGTC")
    vc = _load_esvee_records(tmp_path, [
        ("chr3", "139358615", "10440", "A", "]chr3:139358616]" + insert + "A",
         "80", "PASS", "SVTYPE=INS;MATEID=10441"),
        ("chr3", "139358616", "10441", "G", "G" + insert + "[chr3:139358615[",
         "80", "PASS", "SVTYPE=INS;MATEID=10440"),
    ])
    assert [v.sv_type for v in vc] == ["BND", "BND"]
    sv = vc[0]
    transcripts = [
        t for t in ensembl_grch38.transcripts_at_locus("3", sv.start, sv.start)
        if t.is_protein_coding]
    assert {t.gene_name for t in transcripts} == {"COPB2"}
    for transcript in transcripts:
        effect = _ANNOTATOR.annotate_on_transcript(sv, transcript)
        assert not isinstance(effect, GeneFusion), transcript.id
