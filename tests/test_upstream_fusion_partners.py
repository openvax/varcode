"""3' fusion partners that start just past a junction end (#550).

LINX fuses a 5' partner with a coding transcript whose 5' end lies up to
10 kb past the other breakend, when the junction reads into it sense to
sense, joining the partner's second exon because its first has no splice
acceptor. Varcode adds these fusions as further candidates. The examples
sit around openvax-v2's GABBR1--SLC29A1 junction (esvee 15486,
chr6:29612906 ``[chr6:44218908[C``) in Ensembl 95. GABBR1 (-) keeps its
5' end, SLC29A1 (+) starts 597 bases past the mate, and MYMX (+) ends 672
bases before it.
"""

from pyensembl import cached_release

from varcode import StructuralVariant
from varcode.effects import GeneFusion, StartLoss, TranslocationToIntergenic
from varcode.effects.structural import MAX_UPSTREAM_PARTNER_DISTANCE

# Release 95, as for the other esvee tests. CI installs it.
genome = cached_release(95)
GABBR1_201 = genome.transcript_by_id("ENST00000355973")
MYMX_201 = genome.transcript_by_id("ENST00000573382")
SLC29A1_208 = genome.transcript_by_id("ENST00000393844")
BREAKPOINT = 29_612_906  # intronic in GABBR1-201; the right side is its 5' end
MATE = 44_218_908


def gabbr1_breakend(alt, position=BREAKPOINT):
    return StructuralVariant("chr6", position, "BND", ref="C", alt=alt, genome=genome,
                             convert_ucsc_contig_names=True)


def fusion_candidates(effect):
    return {candidate.effect.three_prime_transcript.transcript_id: candidate.effect
            for candidate in effect.candidates
            if isinstance(candidate.effect, GeneFusion)}


def partners_starting_past(position):
    """Coding chr6 + transcripts with a second exon whose start is 1 to
    10,000 bases past ``position``, read straight from Ensembl."""
    return {t.transcript_id: t.start - position
            for t in genome.transcripts(contig="6", strand="+")
            if t.biotype == "protein_coding" and len(t.exons) > 1
            and 0 < t.start - position <= MAX_UPSTREAM_PARTNER_DISTANCE}


def exon_length(exon):
    return exon.end - exon.start + 1


def test_breakend_stays_translocation_with_a_candidate_per_partner_isoform():
    effect = gabbr1_breakend("[chr6:%d[C" % MATE).effect_on_transcript(GABBR1_201)
    assert type(effect) is TranslocationToIntergenic
    assert effect.most_likely_effect is effect
    assert effect.priority_class is TranslocationToIntergenic
    fusions = fusion_candidates(effect)
    expected = partners_starting_past(MATE)
    assert set(fusions) == set(expected)
    assert {f.three_prime_transcript.gene_name for f in fusions.values()} == {"SLC29A1"}
    assert min(expected.values()) == 597
    five_prime = GABBR1_201.sequence[:sum(
        exon_length(e) for e in GABBR1_201.exons if e.start > BREAKPOINT)]
    for transcript_id, fusion in fusions.items():
        partner = fusion.three_prime_transcript
        model = fusion.mutant_transcript
        assert fusion.five_prime_transcript.transcript_id == GABBR1_201.transcript_id
        # The first exon has no splice acceptor: the join reaches the second.
        assert model.cdna_sequence == five_prime + partner.sequence[
            exon_length(partner.exons[0]):]
        assert model.evidence["three_prime_breakend"] == "upstream"
        assert model.evidence["three_prime_upstream_distance"] == expected[transcript_id]
        assert model.mutant_protein_sequence.startswith(GABBR1_201.protein_sequence[:20])


def test_partner_starts_at_most_the_distance_limit_past_the_mate():
    at_limit = SLC29A1_208.start - MAX_UPSTREAM_PARTNER_DISTANCE
    for mate, included in ((at_limit, True), (at_limit - 1, False)):
        effect = gabbr1_breakend("[chr6:%d[C" % mate).effect_on_transcript(GABBR1_201)
        fusions = fusion_candidates(effect)
        assert set(fusions) == set(partners_starting_past(mate))
        assert (SLC29A1_208.transcript_id in fusions) is included
        # MYMX (+) starts 7.4 kb past either mate.
        assert "MYMX" in {f.three_prime_transcript.gene_name for f in fusions.values()}


def test_partner_must_be_read_away_from_the_mate():
    # Keeping the mate's left side reads the - strand leftward, away from
    # SLC29A1 and against MYMX (+).
    effect = gabbr1_breakend("]chr6:%d]C" % MATE).effect_on_transcript(GABBR1_201)
    assert type(effect) is TranslocationToIntergenic
    assert not fusion_candidates(effect)


def test_only_a_retained_five_prime_end_joins_a_partner_past_the_mate():
    # Keeping GABBR1's left side keeps its 3' end.
    effect = gabbr1_breakend(
        "C[chr6:%d[" % MATE, position=BREAKPOINT - 1).effect_on_transcript(GABBR1_201)
    assert type(effect) is TranslocationToIntergenic
    assert not fusion_candidates(effect)


def test_deletion_keeps_its_consequence_first_and_adds_fusion_candidates():
    # MYMX-201 intron 1 to 505 bases before SLC29A1: the deletion removes
    # MYMX's second exon, with its start codon.
    deletion = StructuralVariant("6", 44_217_000, "DEL", end=44_218_999, genome=genome)
    effect = deletion.effect_on_transcript(MYMX_201)
    assert type(effect.most_likely_effect) is StartLoss
    fusions = fusion_candidates(effect)
    assert set(fusions) == set(partners_starting_past(44_219_000))
    assert {f.five_prime_transcript.gene_name for f in fusions.values()} == {"MYMX"}
    # A structural effect set is ranked by its most disruptive candidate.
    assert isinstance(effect.highest_priority_effect, GeneFusion)
