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

"""Codon translation tables as simple dict/frozenset lookups.

Built by hand from NCBI Translation Tables 1 (Standard) and 2
(Vertebrate Mitochondrial). Keeps codon translation inside varcode so
that downstream callers don't need BioPython just to enumerate 64
codons. See openvax/varcode#293 for the broader move to drop the
biopython dependency.

Refs:
  https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=t#SG1
  https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=t#SG2
"""

from dataclasses import dataclass

_BASES = "TCAG"
_CODONS = tuple(b1 + b2 + b3 for b1 in _BASES for b2 in _BASES for b3 in _BASES)

# AA strings are in NCBI codon order (base1 = T,T,...,G ; base2 cycles
# T,C,A,G within each base1 ; base3 cycles T,C,A,G within each base2).
# So zip(_CODONS, _AA_STRING) gives the forward mapping.

_STANDARD_AAS = \
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

# Vertebrate mitochondrial (table 2) differs from standard at four codons:
#   TGA: * -> W (Trp)
#   ATA: I -> M (Met)
#   AGA: R -> * (stop)
#   AGG: R -> * (stop)
_VERT_MT_AAS = \
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"


def _build_forward_table(aa_string):
    return {codon: aa for codon, aa in zip(_CODONS, aa_string) if aa != "*"}


def _build_stop_codons(aa_string):
    return frozenset(codon for codon, aa in zip(_CODONS, aa_string) if aa == "*")


@dataclass(frozen=True)
class CodonTable:
    """Codon table with amino acid lookup + start/stop sets.

    ``forward_table`` maps sense codons to single-letter AAs; stop
    codons are omitted (consult ``stop_codons`` for those). Matches
    the BioPython ``CodonTable.unambiguous_dna_by_id[n]`` shape so
    callers can be ported in-place.
    """
    name: str
    ncbi_id: int
    forward_table: dict
    start_codons: frozenset
    stop_codons: frozenset


STANDARD = CodonTable(
    name="Standard",
    ncbi_id=1,
    forward_table=_build_forward_table(_STANDARD_AAS),
    start_codons=frozenset({"TTG", "CTG", "ATG"}),
    stop_codons=_build_stop_codons(_STANDARD_AAS),
)

VERTEBRATE_MITOCHONDRIAL = CodonTable(
    name="Vertebrate Mitochondrial",
    ncbi_id=2,
    forward_table=_build_forward_table(_VERT_MT_AAS),
    start_codons=frozenset({"ATT", "ATC", "ATA", "ATG", "GTG"}),
    stop_codons=_build_stop_codons(_VERT_MT_AAS),
)


# Names commonly used for the mitochondrial contig across genome builds
# and VCF dialects: Ensembl "MT", UCSC "chrM", occasionally "chrMT" or
# "M". Lowercased for case-insensitive matching.
_MITOCHONDRIAL_CONTIG_NAMES = frozenset({
    "mt", "chrmt", "m", "chrm", "mtdna", "mito",
})


def is_mitochondrial_contig(contig):
    """True if the given contig name refers to a mitochondrial sequence."""
    if contig is None:
        return False
    return str(contig).lower() in _MITOCHONDRIAL_CONTIG_NAMES


def codon_table_for_transcript(transcript):
    """Select the codon table for a transcript based on its contig.

    Vertebrate mitochondrial genes use NCBI translation table 2;
    every other transcript uses the standard table.
    """
    if transcript is None:
        return STANDARD
    contig = getattr(transcript, "contig", None)
    if is_mitochondrial_contig(contig):
        return VERTEBRATE_MITOCHONDRIAL
    return STANDARD


def translate_sequence(
        nucleotide_sequence,
        codon_table=STANDARD,
        to_stop=True):
    """Translate a cDNA string to amino acids using ``codon_table``.

    Parameters
    ----------
    nucleotide_sequence : str or convertible-to-str
        cDNA. Length must be a multiple of 3.
    codon_table : CodonTable
        Defaults to the standard table.
    to_stop : bool
        If True (default), stop translation at the first stop codon and
        return the protein up to (but not including) the stop. If False,
        include ``'*'`` characters for in-sequence stops.
    """
    seq = str(nucleotide_sequence).upper()
    if len(seq) % 3 != 0:
        raise ValueError(
            "Expected nucleotide sequence length to be a multiple of 3, "
            "got %d (sequence=%r...)" % (len(seq), seq[:30]))
    forward = codon_table.forward_table
    stops = codon_table.stop_codons
    result = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if codon in stops:
            if to_stop:
                break
            result.append("*")
            continue
        aa = forward.get(codon)
        if aa is None:
            raise ValueError(
                "Unknown codon %r at position %d of sequence" % (codon, i))
        result.append(aa)
    return "".join(result)
