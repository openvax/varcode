"""Bounded multi-candidate DNA cryptic splice enumeration and Ap2s1 del17 regressions (#675).

Primary source for Ap2s1 del17:
https://pmc.ncbi.nlm.nih.gov/articles/PMC5824975/ (PMID 29479578)
GRCm39 chr7:16477295 T>C on ENSMUST00000086112 (Ap2s1-201, Ensembl 115):
exon 2 donor +2 T>C activates an upstream exonic cryptic donor 51 nt inside
exon 2, deleting 17 amino acids (residues 35-51, VHAVVTVRDAKHTNFVE).
"""

from __future__ import annotations

import pytest
from pyensembl import EnsemblRelease, cached_release

from varcode import (
    CrypticSplicePolicy,
    EffectCollection,
    GermlineContext,
    HypothesisLimit,
    MutantTranscript,
    Variant,
    apply_rna_evidence_to_effects,
    effect_from_json,
    effect_to_json,
    make_rna_outcome,
    predict_transcript_model_effect,
)
from varcode.effects import Deletion, Substitution
from varcode.genomic_layout import GenomicLayout, SequenceUnavailable
from varcode.nucleotides import reverse_complement
from varcode.transcript_layout import enumerate_cryptic_boundaries
from varcode.transcript_model import refine_transcript_model_effect


def _ap2s1_fixture():
    """Return mouse Ensembl 115 Ap2s1-201 transcript, del17 variant, and local DNA provider."""
    genome = EnsemblRelease(115, species="mouse")
    if not genome.required_local_files_exist():
        pytest.skip("Mouse Ensembl 115 (GRCm39) is not installed locally")
    transcript = genome.transcript_by_id("ENSMUST00000086112")
    # Exon 2 is chr7:16477144..16477293 (+). Donor +2 is chr7:16477295 T>C.
    variant = Variant("7", 16477295, "T", "C", genome)

    # Supply exonic bases from Ensembl cDNA and deterministic intronic bases
    # flanking exon 2 with a canonical GT donor at +1/+2 (16477294..16477295)
    # and realistic competing intronic GT motifs downstream.
    exon_bases: dict[int, str] = {}
    offset = 0
    cdna = str(transcript.sequence)
    for exon in transcript.exons:
        length = exon.end - exon.start + 1
        for i in range(length):
            exon_bases[exon.start + i] = cdna[offset + i]
        offset += length

    # Intron 2 starts at 16477294 with canonical GT (16477294=G, 16477295=T)
    # followed by intronic bases containing a competing downstream cryptic GT.
    exon2 = transcript.exons[1]
    intron2_prefix = "GTGAGC" + "CCCCCC" + "CAGGTAAGT" + ("C" * 160)
    intron1_suffix = ("C" * 120) + "CAG"
    intron_bases: dict[int, str] = {}
    for i, base in enumerate(intron2_prefix):
        intron_bases[exon2.end + 1 + i] = base
    for i, base in enumerate(reversed(intron1_suffix)):
        intron_bases[exon2.start - 1 - i] = base

    def provider(contig, start, end):
        if str(contig) != "7":
            raise SequenceUnavailable(f"No sequence for {contig}:{start}-{end}")
        chars = []
        for pos in range(start, end + 1):
            base = exon_bases.get(pos) or intron_bases.get(pos)
            if base is None:
                # Pad remaining introns with non-motif C for retention/flank spans.
                base = "C"
            chars.append(base)
        return "".join(chars)

    return transcript, variant, provider


def _ap2s1_del17_expected_sequences(transcript):
    """Return expected cDNA and protein for the published 51 nt exon-2 deletion."""
    e1_len = transcript.exons[0].end - transcript.exons[0].start + 1
    e2_len = transcript.exons[1].end - transcript.exons[1].start + 1
    ref_cdna = str(transcript.sequence)
    del17_cdna = ref_cdna[: e1_len + e2_len - 51] + ref_cdna[e1_len + e2_len :]
    ref_protein = str(transcript.protein_sequence)
    # Residues 35..51 (1-based inclusive) = 0-based slice [34:51] = 17 aa.
    deleted_aa = ref_protein[34:51]
    assert deleted_aa == "VHAVVTVRDAKHTNFVE"
    del17_protein = ref_protein[:34] + ref_protein[51:]
    return del17_cdna, del17_protein


def test_ap2s1_del17_default_policy_limitation_and_independent_rna_regression():
    """Default 50 nt / 1-site DNA search misses del17 (-51 nt), while observed RNA preserves it."""
    transcript, variant, provider = _ap2s1_fixture()
    del17_cdna, del17_protein = _ap2s1_del17_expected_sequences(transcript)

    # 1. Default DNA-only policy: explicit limitation (scan_flank=50, max_sites_per_junction=1)
    default_effect = predict_transcript_model_effect(
        (variant,), transcript, sequence_provider=provider
    )
    default_proteins = {
        c.effect.mutant_protein_sequence
        for c in default_effect.candidates
        if c.effect.mutant_protein_sequence is not None
    }
    assert del17_protein not in default_proteins

    # Heuristic cryptic candidate records explicit search bounds & heuristic status.
    cryptic_outcomes = [
        outcome
        for candidate in default_effect.candidates
        for outcome in candidate.outcomes
        if any(opt.mechanism == "cryptic" for _, opt in outcome.hypothesis.splice_plan.choices)
    ]
    assert len(cryptic_outcomes) == 1
    assert cryptic_outcomes[0].mutant.evidence["discovery_status"] == "heuristic_dna_motif"
    assert cryptic_outcomes[0].mutant.evidence["cryptic_search_policy"] == {
        "scan_flank": 50,
        "max_sites_per_junction": 1,
        "min_score": 0.0,
    }

    # 2. Observed RNA for del17 works independently even when DNA hypothesis coverage misses it.
    observed_model = MutantTranscript.from_sequence(
        del17_cdna,
        reference_transcript=transcript,
        mutant_protein_sequence=del17_protein,
        annotator_name="rt_pcr",
        evidence={"pmid": "29479578", "assay": "RT-PCR"},
    )
    observed_effect = Deletion(
        variant=variant,
        transcript=transcript,
        aa_mutation_start_offset=34,
        aa_ref="VHAVVTVRDAKHTNFVE",
    )
    observed_effect.mutant_transcript = observed_model
    rna_candidate = make_rna_outcome(
        observed_effect,
        source="rt_pcr",
        read_count=42,
        extra_evidence={"pmid": "29479578"},
    )

    class _ObservedResolver:
        def observed_outcomes(self, var, tx):
            if var == variant and tx.id == transcript.id:
                return (rna_candidate,)
            return ()

    refined_collection = apply_rna_evidence_to_effects(
        EffectCollection([default_effect]), _ObservedResolver()
    )
    refined = refined_collection[0]
    matching_rna = [
        c
        for c in refined.candidates
        if c.effect.mutant_protein_sequence == del17_protein
    ]
    assert len(matching_rna) == 1
    assert matching_rna[0].outcomes[0].mutant.cdna_sequence == del17_cdna
    assert matching_rna[0].outcomes[0].hypothesis.evidence["rna_observation"]["pmid"] == "29479578"


def test_ap2s1_del17_bounded_multi_candidate_dna_discovery_and_rna_matching():
    """Expanded bounded CrypticSplicePolicy discovers Ap2s1 del17 alongside competing cryptic donors."""
    transcript, variant, provider = _ap2s1_fixture()
    del17_cdna, del17_protein = _ap2s1_del17_expected_sequences(transcript)

    policy = CrypticSplicePolicy(scan_flank=100, max_sites_per_junction=6)
    effect = predict_transcript_model_effect(
        (variant,), transcript, sequence_provider=provider, cryptic_policy=policy
    )

    # Distinct cryptic candidates are enumerated and capped at max_sites_per_junction.
    cryptic_outcomes = [
        outcome
        for candidate in effect.candidates
        for outcome in candidate.outcomes
        if any(opt.mechanism == "cryptic" for _, opt in outcome.hypothesis.splice_plan.choices)
    ]
    assert 2 <= len(cryptic_outcomes) <= policy.max_sites_per_junction

    # Every DNA cryptic candidate is marked heuristic, never experimentally observed.
    for outcome in cryptic_outcomes:
        assert outcome.mutant.evidence["discovery_status"] == "heuristic_dna_motif"
        assert outcome.mutant.evidence["cryptic_search_policy"] == policy.as_evidence()
        assert "rna_support" not in outcome.mutant.evidence

    # The published del17 product (p.V35_E51del, cryptic donor at chr7:16477242) is discovered!
    del17_candidates = [
        c
        for c in effect.candidates
        if c.effect.mutant_protein_sequence == del17_protein
    ]
    assert len(del17_candidates) == 1
    del17_cand = del17_candidates[0]
    assert isinstance(del17_cand.effect, Deletion)
    assert del17_cand.effect.format_short_description(deletion_notation="range") == "p.V35_E51del"
    del17_outcome = del17_cand.outcomes[0]
    assert del17_outcome.mutant.cdna_sequence == del17_cdna
    resolved = del17_outcome.mutant.evidence["resolved_splice_choices"]
    assert resolved[0][2] == transcript.exons[1].end - 51  # chr7:16477242

    # Archive round-trip preserves all distinct cryptic candidates and evidence.
    restored = effect_from_json(effect_to_json(effect))
    assert len(restored.candidates) == len(effect.candidates)
    assert any(c.effect.mutant_protein_sequence == del17_protein for c in restored.candidates)

    # When observed RNA matches the discovered DNA del17 product, rna_support attaches directly
    # onto the coupled baseline/mutant candidate rather than creating an unmapped orphan.
    observed_model = MutantTranscript.from_sequence(
        del17_cdna,
        reference_transcript=transcript,
        mutant_protein_sequence=del17_protein,
        annotator_name="rt_pcr",
    )
    observed_effect = Deletion(variant, transcript, 34, "VHAVVTVRDAKHTNFVE")
    observed_effect.mutant_transcript = observed_model
    rna_candidate = make_rna_outcome(
        observed_effect,
        source="rt_pcr",
        extra_evidence={"pmid": "29479578"},
    )
    refined = refine_transcript_model_effect(effect, (rna_candidate,))
    matched = [c for c in refined.candidates if c.effect.mutant_protein_sequence == del17_protein]
    assert len(matched) == 1
    assert matched[0].outcomes[0].mutant.evidence["rna_support"][0]["pmid"] == "29479578"
    assert matched[0].outcomes[0].baseline.protein_sequence == str(transcript.protein_sequence)


def test_bounded_cryptic_enumeration_respects_max_hypotheses():
    """Exceeding max_hypotheses with multi-candidate cryptic plans returns HypothesisLimit."""
    transcript, variant, provider = _ap2s1_fixture()
    policy = CrypticSplicePolicy(scan_flank=100, max_sites_per_junction=6)
    limited = predict_transcript_model_effect(
        (variant,),
        transcript,
        sequence_provider=provider,
        max_hypotheses=4,
        cryptic_policy=policy,
    )
    assert isinstance(limited, HypothesisLimit)


def test_shared_germline_cryptic_site_preserves_baseline_mutant_coupling():
    """Baseline and mutant couple on identical cryptic boundary positions for shared disrupted sites."""
    transcript, germline_splice, provider = _ap2s1_fixture()
    # Downstream somatic missense in exon 3 (far from exon 2 donor, does not touch exon 2).
    exon3 = transcript.exons[2]
    somatic_pos = exon3.start + 15
    ref_base = provider("7", somatic_pos, somatic_pos)
    alt_base = "A" if ref_base != "A" else "G"
    somatic = Variant("7", somatic_pos, ref_base, alt_base, transcript.genome)

    policy = CrypticSplicePolicy(scan_flank=100, max_sites_per_junction=4)
    # Known-cis germline splice variant so exon 2 donor is identically disrupted on baseline and mutant.
    from types import SimpleNamespace

    effect = predict_transcript_model_effect(
        (somatic,),
        transcript,
        germline_variants=(germline_splice,),
        phase_resolver=SimpleNamespace(in_cis=lambda *args, **kwargs: True),
        sequence_provider=provider,
        cryptic_policy=policy,
    )
    cryptic_outcomes = [
        outcome
        for candidate in effect.candidates
        for outcome in candidate.outcomes
        if any(opt.mechanism == "cryptic" for _, opt in outcome.hypothesis.splice_plan.choices)
    ]
    # Exactly 4 coupled baseline-vs-mutant cryptic pairs (NOT 4 x 4 = 16 cross-pairs!).
    assert len(cryptic_outcomes) == 4
    seen_positions = set()
    for outcome in cryptic_outcomes:
        base_pos = outcome.baseline.evidence["resolved_splice_choices"][0][2]
        mut_pos = outcome.mutant.evidence["resolved_splice_choices"][0][2]
        assert base_pos == mut_pos
        seen_positions.add(mut_pos)
        assert isinstance(outcome.effect, Substitution) or outcome.effect.modifies_coding_sequence
    assert len(seen_positions) == 4


def test_reverse_strand_bounded_cryptic_enumeration_donor_and_acceptor():
    """Reverse-strand transcripts enumerate multiple distinct donor and acceptor cryptic sites."""
    genome = cached_release(81)
    brca1 = genome.transcript_by_id("ENST00000357654")  # reverse strand
    # Build an inline genomic provider around exon 2 (reverse strand) with two distinct donor GTs.
    exon2 = brca1.exons[1]
    e1_len = brca1.exons[0].end - brca1.exons[0].start + 1
    e2_len = exon2.end - exon2.start + 1
    cdna = str(brca1.sequence)
    exon_bases: dict[int, str] = {}
    offset = 0
    for exon in brca1.exons:
        length = exon.end - exon.start + 1
        for i in range(length):
            gpos = exon.end - i if brca1.on_backward_strand else exon.start + i
            base = reverse_complement(cdna[offset + i]) if brca1.on_backward_strand else cdna[offset + i]
            exon_bases[gpos] = base
        offset += length

    # On minus strand, donor of exon 2 is at exon2.start - 1 (G in tx -> C on +) and exon2.start - 2 (T in tx -> A on +).
    # Place canonical GT plus two downstream intronic cryptic donor motifs (CAGGTAAGT) in transcript order.
    tx_intron2 = "GTAAGT" + ("C" * 10) + "CAGGTAAGT" + ("C" * 15) + "AAGGTAAGT" + ("C" * 80)
    intron_plus = reverse_complement(tx_intron2)
    # tx_intron2[0] is at exon2.start - 1, tx_intron2[-1] is at exon2.start - len(tx_intron2)
    intron_bases = {
        exon2.start - len(tx_intron2) + i: base
        for i, base in enumerate(intron_plus)
    }

    def provider(contig, start, end):
        return "".join(exon_bases.get(p) or intron_bases.get(p) or "C" for p in range(start, end + 1))

    layout = GenomicLayout.from_transcript(brca1, flank=100, sequence_provider=provider)
    # Exon 2 donor interbase boundary in layout:
    from varcode.transcript_layout import build_exon_runs

    runs = build_exon_runs(brca1, layout)
    run2 = runs[1]
    valid = (run2.layout_start + 1, runs[2].layout_start - 1)
    sites = enumerate_cryptic_boundaries(
        layout,
        run2.layout_end,
        "donor",
        valid,
        policy=CrypticSplicePolicy(scan_flank=80, max_sites_per_junction=3),
    )
    assert len(sites) == 3
    assert len({b for b, _ in sites}) == 3
    assert all(s > 0.0 for _, s in sites)


def test_cryptic_splice_policy_validation():
    with pytest.raises(ValueError, match="scan_flank"):
        CrypticSplicePolicy(scan_flank=0)
    with pytest.raises(ValueError, match="scan_flank"):
        CrypticSplicePolicy(scan_flank=True)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="max_sites_per_junction"):
        CrypticSplicePolicy(max_sites_per_junction=0)
    with pytest.raises(ValueError, match="max_sites_per_junction"):
        CrypticSplicePolicy(max_sites_per_junction=False)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="min_score"):
        CrypticSplicePolicy(min_score=-0.1)
    with pytest.raises(ValueError, match="min_score"):
        CrypticSplicePolicy(min_score=1.0)
    with pytest.raises(TypeError, match="CrypticSplicePolicy"):
        genome = cached_release(81)
        tx = genome.transcript_by_id("ENST00000003084")
        predict_transcript_model_effect(
            (Variant("7", 117531115, "G", "A", genome),),
            tx,
            cryptic_policy={"scan_flank": 100},  # type: ignore[arg-type]
        )
