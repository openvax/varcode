"""Varcode's NCBI genetic codes and their translation, checked independently.

Expected meanings come from Biotite's own copy of NCBI's tables, never from
Varcode's packaged snapshot or translation code. Biotite's copies (1.6 on
Python 3.11, 2.x later) predate NCBI's current gc.prt (version 4.6) in places,
each corrected below from gc.prt rather than from Varcode:

- both lack code 32, which is code 11 with TAG read as Trp;
- both lack code 33, which is code 24 with TAA read as Tyr;
- both omit GTG from code 3's starts (yeast mitochondrial; gc.prt 4.4);
- Biotite 1.6 reads CTG as Ala in codes 27-30, which gc.prt 4.3 made Leu.
  Codes 27 and 28 are refused, so only 29 and 30 need the correction.
"""

from itertools import product
import random

import pytest
from biotite.sequence import CodonTable as BiotiteCodonTable

from varcode.effects.codon_tables import codon_table_by_id, initiator_as_met, translate_sequence

# NCBI's current genetic codes. 7 and 8 were deleted (merged into 4 and 1);
# 17-20 were never assigned.
NCBI_CODES = (1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16,
              21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33)
# A stop codon in these is also read as sense, depending on context; Varcode refuses them.
CONTEXT_DEPENDENT = (27, 28, 31)
SUPPORTED = tuple(code for code in NCBI_CODES if code not in CONTEXT_DEPENDENT)
# NCBI code -> (Biotite code it starts from, changed meanings, added starts).
FROM_GC_PRT = {32: (11, {"TAG": "W"}, ()), 33: (24, {"TAA": "Y"}, ()), 3: (3, {}, ("GTG",)),
               29: (29, {"CTG": "L"}, ()), 30: (30, {"CTG": "L"}, ())}

CODONS = tuple("".join(bases) for bases in product("TCAG", repeat=3))
# IUPAC nucleotide codes, plus X, which Varcode also accepts as any base.
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT",
         "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
         "X": "ACGT"}
AMBIGUOUS = tuple(codon for codon in ("".join(symbols) for symbols in product(IUPAC, repeat=3))
                  if codon not in CODONS)


def ncbi(code):
    """NCBI code ``code`` as (codon -> residue, ``*`` for a stop; start codons), from Biotite."""
    base, changes, added_starts = FROM_GC_PRT.get(code, (code, {}, ()))
    table = BiotiteCodonTable.load(base)
    return {**table.codon_dict(), **changes}, frozenset(table.start_codons()) | set(added_starts)


def test_every_supported_code_loads_and_retired_codes_are_refused():
    # Refusal of the context-dependent codes is checked in test_species_genetic_codes.
    for code in SUPPORTED:
        assert codon_table_by_id(code).ncbi_id == code
    for code in (7, 8, 17, 18, 19, 20, 34):
        with pytest.raises(ValueError, match="Unsupported"):
            codon_table_by_id(code)


@pytest.mark.parametrize("code", SUPPORTED)
def test_codon_meanings_stops_and_starts_match_ncbi(code):
    residues, starts = ncbi(code)
    table = codon_table_by_id(code)
    assert table.forward_table == {c: aa for c, aa in residues.items() if aa != "*"}
    assert table.stop_codons == {c for c, aa in residues.items() if aa == "*"}
    assert table.start_codons == starts


@pytest.mark.parametrize("code", SUPPORTED)
def test_translation_reads_every_codon_and_stops_at_the_first_stop(code):
    residues, _ = ncbi(code)
    table = codon_table_by_id(code)
    order = list(CODONS)
    random.Random(code).shuffle(order)
    expected = "".join(residues[c] for c in order)
    assert translate_sequence("".join(order), table, to_stop=False) == expected
    assert translate_sequence("".join(order).lower(), table, to_stop=False) == expected
    assert translate_sequence("".join(order), table) == expected.split("*")[0]


@pytest.mark.parametrize("code", SUPPORTED)
def test_ambiguous_codons_read_as_their_one_meaning_or_x(code):
    residues, _ = ncbi(code)
    expected = []
    for codon in AMBIGUOUS:
        meanings = {residues["".join(bases)] for bases in product(*(IUPAC[s] for s in codon))}
        expected.append(meanings.pop() if len(meanings) == 1 else "X")
    actual = translate_sequence("".join(AMBIGUOUS), codon_table_by_id(code), to_stop=False)
    assert len(actual) == len(AMBIGUOUS)
    assert [(c, got, want) for c, got, want in zip(AMBIGUOUS, actual, expected) if got != want] == []


@pytest.mark.parametrize("code", SUPPORTED)
def test_only_ncbi_start_codons_initiate_as_met(code):
    residues, starts = ncbi(code)
    table = codon_table_by_id(code)
    for codon in CODONS:
        if residues[codon] == "*":
            continue
        protein = translate_sequence(codon * 2, table, to_stop=False)
        first = "M" if codon in starts else residues[codon]
        # Only the initiator changes; the same codon later keeps its internal meaning.
        assert initiator_as_met(protein, codon, table) == first + residues[codon], codon
