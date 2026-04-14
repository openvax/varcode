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

"""Tests for mitochondrial codon table handling (openvax/varcode#294).

Covers the four codon differences between NCBI table 1 (Standard) and
NCBI table 2 (Vertebrate Mitochondrial):

  * TGA: standard stop -> mt Trp
  * AGA: standard Arg  -> mt stop
  * AGG: standard Arg  -> mt stop
  * ATA: standard Ile  -> mt Met

Also covers the local codon-table data structures (which exist so that
varcode doesn't need BioPython just to enumerate 64 codons — see #293).
"""

from pyensembl import cached_release

from varcode import Variant
from varcode.effects.codon_tables import (
    STANDARD,
    VERTEBRATE_MITOCHONDRIAL,
    codon_table_for_transcript,
    is_mitochondrial_contig,
    translate_sequence,
)


ensembl_grch38 = cached_release(81)

# MT-CO1-201: complete protein-coding transcript on MT (+ strand),
# uses TAA as stop (so pyensembl considers it complete, unlike MT-ND1
# which stops on AGA). Exon 5904-7445, start codon at 5904-5906.
MT_CO1_TRANSCRIPT_ID = "ENST00000361624"


# ====================================================================
# Codon-table data structures
# ====================================================================


def test_standard_and_mt_tables_cover_all_64_codons():
    for table in (STANDARD, VERTEBRATE_MITOCHONDRIAL):
        seen = set(table.forward_table) | set(table.stop_codons)
        assert len(seen) == 64, (
            "Table %s covers %d codons; expected 64" % (table.name, len(seen)))


def test_standard_table_stops():
    assert STANDARD.stop_codons == frozenset({"TAA", "TAG", "TGA"})


def test_mt_table_stops():
    # Key difference: TGA is NOT a stop, AGA/AGG ARE stops.
    assert VERTEBRATE_MITOCHONDRIAL.stop_codons == frozenset({
        "TAA", "TAG", "AGA", "AGG",
    })


def test_tga_codes_trp_in_mt():
    assert STANDARD.forward_table.get("TGA") is None  # it's a stop
    assert VERTEBRATE_MITOCHONDRIAL.forward_table["TGA"] == "W"


def test_aga_agg_are_stops_in_mt():
    assert STANDARD.forward_table["AGA"] == "R"
    assert STANDARD.forward_table["AGG"] == "R"
    assert "AGA" in VERTEBRATE_MITOCHONDRIAL.stop_codons
    assert "AGG" in VERTEBRATE_MITOCHONDRIAL.stop_codons


def test_ata_codes_met_in_mt():
    assert STANDARD.forward_table["ATA"] == "I"
    assert VERTEBRATE_MITOCHONDRIAL.forward_table["ATA"] == "M"


def test_mt_start_codons():
    assert VERTEBRATE_MITOCHONDRIAL.start_codons == frozenset({
        "ATT", "ATC", "ATA", "ATG", "GTG",
    })


def test_standard_start_codons():
    assert STANDARD.start_codons == frozenset({"TTG", "CTG", "ATG"})


# ====================================================================
# Contig-name detection
# ====================================================================


def test_mitochondrial_contig_detection():
    for name in ("MT", "chrMT", "chrM", "M", "mt", "chrmt", "mtDNA", "mito"):
        assert is_mitochondrial_contig(name), (
            "Expected %r to be recognized as mitochondrial" % name)


def test_non_mitochondrial_contigs():
    for name in ("1", "chr1", "X", "chrX", "Y", None, ""):
        assert not is_mitochondrial_contig(name)


def test_codon_table_for_mt_transcript_is_vertebrate_mt():
    mt_t = ensembl_grch38.transcript_by_id(MT_CO1_TRANSCRIPT_ID)
    assert codon_table_for_transcript(mt_t) is VERTEBRATE_MITOCHONDRIAL


def test_codon_table_for_nuclear_transcript_is_standard():
    # CFTR is on chr7
    nuclear_t = ensembl_grch38.transcript_by_id("ENST00000003084")
    assert codon_table_for_transcript(nuclear_t) is STANDARD


# ====================================================================
# Sequence translation
# ====================================================================


def test_translate_sequence_standard_tga_stops():
    # ATG + TGG + TGA stops at TGA under standard
    assert translate_sequence("ATGTGGTGA") == "MW"


def test_translate_sequence_mt_tga_becomes_trp():
    # Same sequence under mt: TGA -> W, no stop here
    assert translate_sequence(
        "ATGTGGTGA", codon_table=VERTEBRATE_MITOCHONDRIAL) == "MWW"


def test_translate_sequence_mt_aga_stops():
    # Under mt, AGA is a stop; standard keeps translating
    assert translate_sequence(
        "ATGAGATTT", codon_table=VERTEBRATE_MITOCHONDRIAL) == "M"
    assert translate_sequence("ATGAGATTT") == "MRF"


def test_translate_sequence_to_stop_false_emits_star():
    assert translate_sequence("ATGTAAGGG", to_stop=False) == "M*G"


def test_translate_sequence_rejects_non_multiple_of_three():
    try:
        translate_sequence("ATGT")
    except ValueError:
        return
    raise AssertionError("Expected ValueError for non-multiple-of-3 input")


# ====================================================================
# End-to-end variant effect prediction on MT
# ====================================================================


def test_mt_tga_creation_is_substitution_not_premature_stop():
    """MT-CO1 aa_pos 278 is encoded by TCA (Ser). A C>G at the middle
    base (genome position 6739) changes it to TGA. Under the standard
    codon table that would be a PrematureStop; under the correct
    vertebrate mitochondrial table it's a Ser->Trp substitution.
    """
    variant = Variant("MT", 6739, "C", "G", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(MT_CO1_TRANSCRIPT_ID)
    effect = variant.effect_on_transcript(transcript)
    assert type(effect).__name__ == "Substitution", (
        "Expected Substitution, got %s (%s). Under mt codon table, "
        "TGA is Trp not stop." % (type(effect).__name__, effect.short_description))
    assert effect.aa_ref == "S"
    assert effect.aa_alt == "W"


def test_mt_aga_creation_is_premature_stop():
    """MT-CO1 aa_pos 37 is encoded by CGA (Arg in both tables). A
    C>A at the first base of that codon (genome 6015) makes AGA.
    Under the standard table, AGA is still Arg (Silent). Under the
    vertebrate mitochondrial table, AGA is a stop (PrematureStop).
    """
    variant = Variant("MT", 6015, "C", "A", ensembl_grch38)
    transcript = ensembl_grch38.transcript_by_id(MT_CO1_TRANSCRIPT_ID)
    effect = variant.effect_on_transcript(transcript)
    assert type(effect).__name__ == "PrematureStop", (
        "Expected PrematureStop, got %s (%s). Under mt codon table, "
        "AGA is a stop." % (type(effect).__name__, effect.short_description))


def test_mt_ata_reference_is_met_not_ile():
    """MT-CO1 aa_pos 64 is encoded by ATA. Under standard that's Ile;
    under mt it's Met (matches the Ensembl-provided protein_sequence,
    which confirms varcode's effect annotation must also use Met).

    Changing ATA to ATT (which is Ile in both tables) is therefore a
    Met->Ile Substitution under mt, but would look Silent under
    standard.
    """
    transcript = ensembl_grch38.transcript_by_id(MT_CO1_TRANSCRIPT_ID)
    # Sanity: the reference protein has Met at aa 64.
    assert transcript.protein_sequence[64] == "M"

    variant = Variant("MT", 6098, "A", "T", ensembl_grch38)
    effect = variant.effect_on_transcript(transcript)
    assert type(effect).__name__ == "Substitution", (
        "Expected Substitution, got %s (%s). Under mt codon table, "
        "ATA is Met, so ATA->ATT is Met->Ile."
        % (type(effect).__name__, effect.short_description))
    assert effect.aa_ref == "M"
    assert effect.aa_alt == "I"


# ====================================================================
# Back-compat: module-level constants still point at standard table
# ====================================================================


def test_legacy_translate_module_constants_are_standard():
    from varcode.effects.translate import (
        DNA_CODON_TABLE,
        START_CODONS,
        STOP_CODONS,
    )
    # Back-compat: callers that imported these directly still see the
    # standard (nuclear) table.
    assert "ATG" in START_CODONS
    assert "TTG" in START_CODONS  # alternate start in standard
    assert "ATA" not in START_CODONS  # only in mt
    assert "TGA" in STOP_CODONS
    assert "AGA" not in STOP_CODONS  # not a stop in standard
    assert DNA_CODON_TABLE["AGA"] == "R"
