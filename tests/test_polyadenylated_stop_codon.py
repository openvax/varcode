"""A stop codon completed by the poly(A) tail (#731).

MT-ND1, MT-ND2, MT-CO3, MT-ND3, MT-ND4 and MT-CYB end their reading frame in
``T`` or ``TA``, and polyadenylation adds the rest of a ``TAA`` stop codon
(NCBI NC_012920.1). Ensembl annotates no stop codon for them, so PyEnsembl's
``complete`` is False. Before, every variant in these genes was reported as
``IncompleteTranscript``, including two primary LHON mutations.
"""

import pytest
from pyensembl.transcript import Transcript

from varcode import Variant
from varcode.effects.effect_classes import IncompleteTranscript
from varcode.effects.transcript_helpers import (
    has_complete_coding_sequence, poly_a_completion, three_prime_utr_offset)

from .data import ENSEMBL_GRCH38

pytestmark = pytest.mark.skipif(
    not hasattr(Transcript, "stop_codon_completed_by_polyadenylation"),
    reason="needs PyEnsembl 2.34.1 or later")

POLY_A_GENES = ("MT-ND1", "MT-ND2", "MT-CO3", "MT-ND3", "MT-ND4", "MT-CYB")
ANNOTATORS = ("transcript_model", "fast", "protein_diff")


def transcript_of(gene_name):
    return ENSEMBL_GRCH38.genes_by_name(gene_name)[0].transcripts[0]


def base(transcript, position):
    return transcript.sequence[transcript.spliced_offset(position)]


def effect_on(variant, transcript, annotator):
    effects = [e for e in variant.effects(annotator=annotator)
               if getattr(e, "transcript", None) is not None
               and e.transcript.id == transcript.id]
    assert len(effects) == 1
    return effects[0]


def substitution(transcript, position, alt):
    return Variant("MT", position, base(transcript, position), alt, genome=ENSEMBL_GRCH38)


@pytest.mark.parametrize("annotator", ANNOTATORS)
@pytest.mark.parametrize("gene, position, ref, alt, change", [
    ("MT-ND1", 3460, "G", "A", "p.A52T"),    # LHON
    ("MT-ND4", 11778, "G", "A", "p.R340H"),  # LHON
    ("MT-ND6", 14484, "T", "C", "p.M64V"),   # LHON; MT-ND6 has a stop_codon row
])
def test_lhon_variants_have_their_protein_effect(
        annotator, gene, position, ref, alt, change):
    transcript = transcript_of(gene)
    variant = Variant("MT", position, ref, alt, genome=ENSEMBL_GRCH38)
    effect = effect_on(variant, transcript, annotator)
    assert (type(effect).__name__, effect.short_description) == ("Substitution", change)


@pytest.mark.parametrize("gene", POLY_A_GENES)
def test_every_poly_a_gene_is_complete_and_has_no_three_prime_utr(gene):
    transcript = transcript_of(gene)
    assert not transcript.complete
    assert has_complete_coding_sequence(transcript)
    assert three_prime_utr_offset(transcript) == len(transcript.sequence)
    # ND1 ends in TA and ND4 in T; one or two A complete a TAA.
    assert poly_a_completion(transcript, transcript.coding_sequence) in ("A", "AA")
    assert (transcript.coding_sequence + poly_a_completion(
        transcript, transcript.coding_sequence)).endswith("TAA")


@pytest.mark.parametrize("annotator", ANNOTATORS)
@pytest.mark.parametrize("gene", POLY_A_GENES)
def test_variants_in_poly_a_genes_are_no_longer_incomplete(gene, annotator):
    transcript = transcript_of(gene)
    # The second base of codon 2 (on the plus strand) to a different base.
    position = transcript.start_codon_positions[0] + 4
    ref = base(transcript, position)
    variant = substitution(transcript, position, "C" if ref != "C" else "G")
    effect = effect_on(variant, transcript, annotator)
    assert not isinstance(effect, IncompleteTranscript)
    assert effect.modifies_protein_sequence


@pytest.mark.parametrize("annotator", ANNOTATORS)
@pytest.mark.parametrize("gene, position, alt, expected", [
    # MT-ND1 ends ...ACC TA: the T and A at 3'-end are two thirds of the stop.
    ("MT-ND1", 4261, "C", "p.*319Q (stop-loss)"),  # CA + A = CAA
    ("MT-ND1", 4262, "G", "p.*319W (stop-loss)"),  # TG + A = TGA, Trp in table 2
    ("MT-ND1", 4262, "T", "p.*319L (stop-loss)"),  # TT + A = TTA
    # MT-ND4 ends ...AGT T: one base of its stop codon.
    ("MT-ND4", 12137, "C", "p.*460Q (stop-loss)"),  # C + AA = CAA
    ("MT-ND4", 12137, "A", "p.*460K (stop-loss)"),  # A + AA = AAA
])
def test_a_change_to_the_encoded_part_of_the_stop_is_a_stop_loss(
        annotator, gene, position, alt, expected):
    transcript = transcript_of(gene)
    effect = effect_on(substitution(transcript, position, alt), transcript, annotator)
    assert (type(effect).__name__, effect.short_description) == ("StopLoss", expected)


@pytest.mark.parametrize("annotator", ANNOTATORS)
def test_the_last_whole_codon_before_the_partial_stop_reads_normally(annotator):
    transcript = transcript_of("MT-ND1")
    # ACC (Thr 318) to ACG: still Thr.
    effect = effect_on(substitution(transcript, 4260, "G"), transcript, annotator)
    assert (type(effect).__name__, effect.short_description) == ("Silent", "p.T318=")


@pytest.mark.parametrize("annotator", ANNOTATORS)
def test_a_frameshift_in_a_poly_a_gene_is_a_frameshift(annotator):
    transcript = transcript_of("MT-ND1")
    position = 4255
    anchor, ref = base(transcript, position - 1), base(transcript, position)
    deletion = Variant("MT", position - 1, anchor + ref, anchor, genome=ENSEMBL_GRCH38)
    effect = effect_on(deletion, transcript, annotator)
    assert (type(effect).__name__, effect.short_description) == ("FrameShift", "p.Q317fs")


def test_a_gene_with_an_annotated_stop_codon_is_untouched():
    transcript = transcript_of("MT-CO1")
    assert transcript.complete
    assert poly_a_completion(transcript, transcript.coding_sequence) == ""
    assert three_prime_utr_offset(transcript) == max(transcript.stop_codon_spliced_offsets) + 1


class _Transcript:
    def __init__(self, complete, **attributes):
        self.complete = complete
        self.__dict__.update(attributes)


def test_older_pyensembl_without_the_property_behaves_as_before():
    assert has_complete_coding_sequence(_Transcript(True))
    assert not has_complete_coding_sequence(_Transcript(False))
    assert has_complete_coding_sequence(
        _Transcript(False, stop_codon_completed_by_polyadenylation=True))
    assert poly_a_completion(_Transcript(False), "ATGTA") == ""
