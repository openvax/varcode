"""Varcode annotates every fusion the OpenVax libraries share.

openvax-v2, the OpenVax libraries' shared Sid test data (iskandr/osteosarc#56),
lists seven fusion targets. Each is a junction between two breakends, named
5' partner first. Four give the side of each break the junction keeps
(``retained_side``). The three from Isovar's sv-regressions-v1 panel give
``orientation`` instead, the strand the fusion follows at each breakend: +
keeps the left side at the first breakend and the right side at the second,
and - the reverse (iskandr/osteosarc#96). Built this way, the five junctions
esvee called in IPISRC044_tumor_T2_ucla match its records.

ATP5MG--KMT2A and TPST1--CRCP join exon 1 of the 5' gene to exon 2 of the
3' gene, and Varcode fuses every coding isoform of the 5' gene with the 3'
gene. Isovar translates ATP5MG--KMT2A's junction reads to MAQFVRNLVEKTPALVNG
and a stop, and Varcode builds that protein from each ATP5MG isoform. No
TPST1 isoform keeps an annotated start codon, so TPST1--CRCP has no fusion
protein.

At the far end of the other five, no protein-coding transcript is read in
the direction the junction continues. The 5' gene then keeps only its 5'
fragment. The first run downloads the bundle (28 MB) into the osteosarc
cache; later runs reuse it offline.
"""

import json
from pathlib import Path

import pytest
from pyensembl import cached_release

from varcode.effects import GeneFusion, TranslocationToIntergenic
from varcode.sv_allele_parser import breakend_sides, parse_symbolic_alt

# Release 95, as for the other esvee tests. CI installs it.
ensembl_grch38 = cached_release(95)

# target -> (5' gene, 3' gene, fusion protein from each 5' isoform)
FUSIONS = {
    "ATP5MG--KMT2A": ("ATP5MG", "KMT2A", "MAQFVRNLVEKTPALVNG"),
    "TPST1--CRCP": ("TPST1", "CRCP", None),
}
# target -> 5' gene, where the far end has no sense-oriented coding partner
# in Ensembl 95.
NO_PARTNER = {
    # The junction reads chr17's + strand; CCDC47 and AC046185.1 are on -.
    "FOXO3--STRADA-CCDC47": "FOXO3",
    # The junction reads chr6's + strand from 597 bases before SLC29A1
    # starts. LINX would fuse SLC29A1 from exon 2 (#550).
    "GABBR1--SLC29A1": "GABBR1",
    # Intergenic; PIK3CB, the next gene to the right, is on -.
    "KLF15--PPIAP72": "KLF15",
    # The junction reads chr15's + strand; FMN1 is on -.
    "OTUD7A--FMN1": "OTUD7A",
    # CDKN2B is on -, CDKN2B-AS1 is noncoding and AL359922.1's only
    # transcript is NMD.
    "PARD3B--CDKN2B-AS1-CDKN2B": "PARD3B",
}
# target -> the esvee record (CHROM, POS, REF, ALT) at its 5' breakend
ESVEE_RECORDS = {
    "FOXO3--STRADA-CCDC47": ("chr6", 108_611_245, "G", "G[chr17:63745980["),  # 16765
    "GABBR1--SLC29A1": ("chr6", 29_612_906, "C", "[chr6:44218908[C"),  # 15486
    "KLF15--PPIAP72": ("chr3", 126_343_949, "C", "[chr3:138649012[C"),  # 8618
    "OTUD7A--FMN1": ("chr15", 31_743_670, "G", "[chr15:33043479[CTG"),  # 35324
    "PARD3B--CDKN2B-AS1-CDKN2B": (
        "chr2", 205_310_103, "G", "GTGTGCTGATATG[chr9:22007648["),  # 6316
}


@pytest.fixture(scope="module")
def fusion_targets():
    osteosarc = pytest.importorskip("osteosarc", reason="install .[test-data]")
    recipe = json.loads(
        (Path(osteosarc.fetch_bundle("openvax-v2")) / "recipe.json").read_text())
    return {name: target for name, target in recipe["targets"].items()
            if target["kind"] == "sv" and target["label"] == "fusion"}


def breakend(contig, position, ref, alt):
    return parse_symbolic_alt(
        contig, position, ref, alt, genome=ensembl_grch38,
        convert_ucsc_contig_names=True)


def kept_sides(target):
    """Each breakend's retained side, given or read from its orientation."""
    sides = []
    for index, end in enumerate(target["breakends"]):
        side = end.get("retained_side")
        if side is None:
            side = "left" if (end["orientation"] == "+") == (index == 0) else "right"
        sides.append(side)
    return sides


def fusion_breakend(target):
    """The target's junction as a BND record at its 5' breakend. A kept
    left side ends at the zero-based interbase position; a kept right side
    starts at the next base."""
    (five, three), (five_side, three_side) = target["breakends"], kept_sides(target)
    position = five["position"] + int(five_side == "right")
    bracket = "[" if three_side == "right" else "]"
    mate = "%s%s:%d%s" % (
        bracket, three["contig"],
        three["position"] + int(three_side == "right"), bracket)
    alt = "N" + mate if five_side == "left" else mate + "N"
    return breakend(five["contig"], position, "N", alt)


def geometry(variant):
    return (variant.contig, variant.start, variant.mate_contig,
            variant.mate_start, breakend_sides(variant.symbolic_alt))


def coding_effects(variant, gene_name):
    effects = [effect for effect in variant.effects()
               if effect.transcript.gene_name == gene_name
               and effect.transcript.is_protein_coding]
    assert effects
    return effects


def test_every_fusion_target_is_checked(fusion_targets):
    assert set(fusion_targets) == set(FUSIONS) | set(NO_PARTNER)


@pytest.mark.parametrize("name", sorted(ESVEE_RECORDS))
def test_junction_matches_esvee_record(fusion_targets, name):
    assert geometry(fusion_breakend(fusion_targets[name])) == geometry(
        breakend(*ESVEE_RECORDS[name]))


@pytest.mark.parametrize("name", sorted(FUSIONS))
def test_exon_junction_fuses_every_coding_isoform(fusion_targets, name):
    five_gene, three_gene, protein = FUSIONS[name]
    for effect in coding_effects(fusion_breakend(fusion_targets[name]), five_gene):
        assert isinstance(effect, GeneFusion)
        assert effect.five_prime_transcript.gene_name == five_gene
        assert effect.three_prime_transcript.gene_name == three_gene
        assert effect.mutant_transcript.mutant_protein_sequence == protein


@pytest.mark.parametrize("name", sorted(NO_PARTNER))
def test_junction_without_sense_partner_keeps_five_prime_fragment(fusion_targets, name):
    variant = fusion_breakend(fusion_targets[name])
    assert not any(isinstance(effect, GeneFusion) for effect in variant.effects())
    for effect in coding_effects(variant, NO_PARTNER[name]):
        assert type(effect) is TranslocationToIntergenic
        segment, = effect.mutant_transcript.reference_segments
        assert segment.label == "translocation_5p"
