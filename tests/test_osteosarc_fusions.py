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
CPEB2::FAM193A pair a duplication and OTX1::KIF3C an inversion, so
those load as typed SVs and must still come out as fusions.
"""

import pytest
from pyensembl import cached_release

from varcode import StructuralVariant, StructuralVariantAnnotator, load_vcf
from varcode.annotators.structural_variant import (
    _build_fusion_mutant_transcript,
)
from varcode.effects import GeneFusion
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

    mt = _build_fusion_mutant_transcript(sv, cpeb2, fam193a)
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
    contig names and the genome given to load_vcf. The DUP- and
    INV-labeled breakend pairs load as typed SVs spanning the event."""
    vc = _esvee_svs()
    assert len(vc) == 8
    assert all(isinstance(v, StructuralVariant) for v in vc)
    assert {v.genome for v in vc} == {ensembl_grch38}
    parsed = {
        vc.metadata[v]["id"]: (
            v.sv_type, v.contig, v.start, v.end, v.mate_contig)
        for v in vc}
    assert parsed == {
        "4055": ("INV", "2", 25_955_847, 63_053_516, "2"),
        "4571": ("INV", "2", 25_955_847, 63_053_516, "2"),
        "11585": ("DUP", "4", 2_664_477, 15_012_987, "4"),
        "11828": ("DUP", "4", 2_664_477, 15_012_987, "4"),
        "11309": ("BND", "4", 75_037_126, 75_037_126, "5"),
        "15093": ("BND", "5", 116_474_281, 116_474_281, "4"),
        "30830": ("BND", "12", 72_342_468, 72_342_468, None),
        "32579": ("BND", "13", 56_435_058, 56_435_058, None),
    }
    singles = [v for v in vc if v.mate_contig is None]
    assert all(v.info["single_breakend"] for v in singles)


def test_esvee_vcf_breakend_records_annotate():
    """Breakend records, single breakends included, go through
    ``effects()``. (The typed DUP and INV span tens of megabases, so
    they're checked transcript by transcript below.)"""
    for sv in _esvee_svs():
        if sv.sv_type == "BND":
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
    """Either record of the INV-labeled OTX1::KIF3C junction gives the
    same fusion on OTX1 and on KIF3C: the 35-residue out-of-frame
    protein from the direct test above."""
    otx1 = ensembl_grch38.transcript_by_id("ENST00000282549")
    kif3c = ensembl_grch38.transcript_by_id("ENST00000264712")
    for record_id in ("4055", "4571"):
        sv = _esvee_sv(record_id)
        for transcript in (otx1, kif3c):
            fusion = _ANNOTATOR.annotate_on_transcript(sv, transcript)
            assert isinstance(fusion, GeneFusion)
            assert fusion.five_prime_transcript.id == otx1.id
            assert fusion.three_prime_transcript.id == kif3c.id
            assert fusion.mutant_transcript.mutant_protein_sequence == (
                "MMSYLKQPPYGMNGLGLAGPAMDLLHPSVGYPETS")


def test_esvee_vcf_cpeb2_fam193a_duplication_is_fusion():
    """esvee labels CPEB2::FAM193A a 12 Mb tandem duplication. CPEB2
    holds its right end and FAM193A its left, so each gets a fusion with
    CPEB2 as 5' partner rather than an internal duplication of its own
    exons."""
    cpeb2 = ensembl_grch38.transcript_by_id("ENST00000538197")
    fam193a = ensembl_grch38.transcript_by_id("ENST00000637812")
    sv = _esvee_sv("11828")

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
