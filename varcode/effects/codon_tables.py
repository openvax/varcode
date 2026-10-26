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

"""Offline NCBI genetic codes and species/compartment assignments.

The packaged snapshot records its NCBI source hashes. Refresh it with
``notes/build_genetic_code_data.py``. Context-dependent stop/sense codes
remain unsupported rather than silently losing their recoding semantics.
"""

from dataclasses import dataclass, replace
from importlib.resources import files
import json

_BASES = "TCAG"
_CODONS = tuple(b1 + b2 + b3 for b1 in _BASES for b2 in _BASES for b3 in _BASES)

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


_DATA = json.loads(files("varcode").joinpath("data/genetic_codes.json").read_text())
_TABLES = {
    int(identifier): CodonTable(
        name=row["name"], ncbi_id=int(identifier),
        forward_table=_build_forward_table(row["amino_acids"]),
        start_codons=frozenset(c for c, aa in zip(_CODONS, row["starts_stops"]) if aa == "M"),
        stop_codons=_build_stop_codons(row["starts_stops"]))
    for identifier, row in _DATA["codes"].items()
}
_TAXA_BY_NAME = {alias: taxon for taxon in _DATA["taxa"].values() for alias in taxon["aliases"]}
STANDARD = _TABLES[1]
VERTEBRATE_MITOCHONDRIAL = _TABLES[2]


def codon_table_by_id(code):
    """Return an unambiguous NCBI code without downloading data."""
    if type(code) is not int or code not in _TABLES:
        raise ValueError("Unsupported explicit NCBI translation table: %s" % code)
    table = _TABLES[code]
    if table.stop_codons.intersection(table.forward_table):
        raise ValueError("NCBI translation table %s requires context-dependent stop recoding" % code)
    return table


def _reference_code(transcript):
    assignment = getattr(transcript, "genetic_code", None)
    if getattr(assignment, "source", None) == "reference" and assignment.table_id is None:
        raise ValueError(
            "Explicitly unset %s genetic code for contig %r; "
            "supply a translation table in the annotation/model or reference metadata"
            % (assignment.compartment, getattr(transcript, "contig", None)))
    contig = str(getattr(transcript, "contig", "")).lower()
    compartment = ("mitochondrial" if is_mitochondrial_contig(contig) else
                   "plastid" if contig in {"pt", "cp", "chrc", "chloroplast", "plastid", "apicoplast"} else
                   "nuclear")
    genome = getattr(transcript, "genome", None)
    declared_compartments = getattr(genome, "contig_compartments", {})
    compartment = declared_compartments.get(getattr(transcript, "contig", None), compartment)
    if compartment not in {"nuclear", "mitochondrial", "plastid"}:
        raise ValueError("Unsupported reference compartment: %r" % compartment)
    # Custom/reference providers can supply these without an annotation change.
    explicit = getattr(genome, compartment + "_translation_table", None)
    if explicit is not None:
        return explicit
    species = getattr(genome, "species", None)
    taxon_id = getattr(genome, "taxon_id", None)
    if species is None and taxon_id is None:
        # Source-explicit PyEnsembl Genomes (including with_mirror()) retain
        # the assembly but not the Species object. Resolve exact registered
        # assemblies only; never guess a species from a filename or prefix.
        from pyensembl.species import Species
        species = Species._reference_names_to_species.get(getattr(genome, "reference_name", None))
    name = str(getattr(species, "latin_name", species)).lower().replace(" ", "_")
    if name not in _TAXA_BY_NAME and species is not None:
        from pyensembl.species import find_species_by_name
        try:
            name = find_species_by_name(name).latin_name
        except ValueError:
            pass
    taxon = _DATA["taxa"].get(str(taxon_id)) if taxon_id is not None else _TAXA_BY_NAME.get(name)
    named = _TAXA_BY_NAME.get(name)
    if taxon_id is not None and named is not None and named is not taxon:
        raise ValueError("Reference species and taxon_id disagree when selecting a genetic code")
    if taxon is not None and taxon.get(compartment):
        return taxon[compartment]
    # Preserve the established standard-code default for unscoped nuclear
    # transcripts. Organelle codes cannot be inferred from the contig alone.
    if compartment == "nuclear" and species is None and taxon_id is None:
        return 1
    raise ValueError(
        "No %s genetic-code assignment for species %r (taxon_id=%r); "
        "supply transcript.translation_table or genome.%s_translation_table"
        % (compartment, species, taxon_id, compartment))


# Names commonly used for the mitochondrial contig across genome builds
# and VCF dialects: Ensembl "MT", UCSC "chrM", occasionally "chrMT" or
# "M". Lowercased for case-insensitive matching.
_MITOCHONDRIAL_CONTIG_NAMES = frozenset({
    "mt", "chrmt", "m", "chrm", "mtdna", "mito",
    "mitochondrion_genome", "mitochondrial_genome", "mitochondrion", "mitochondrial",
})


def is_mitochondrial_contig(contig):
    """True if the given contig name refers to a mitochondrial sequence."""
    if contig is None:
        return False
    return str(contig).lower() in _MITOCHONDRIAL_CONTIG_NAMES


def codon_table_for_transcript(transcript):
    """Select this coding model's genetic code and initiation codons.

    Explicit ``translation_table`` and source-declared ``start_codon`` take
    precedence over NCBI species/compartment assignments. A declared start extends only this model's
    initiation set; its internal amino acid and the global tables are unchanged.
    This recognizes a retained source start, not evidence for mutant translation.
    An explicitly unset reference assignment raises ValueError instead of
    falling back to a species default; providers without this metadata retain
    the usual fallback.
    """
    if transcript is None:
        return STANDARD
    # An explicit alternative-ORF model owns its genetic code. Do not infer
    # table 2 from an MT-like contig name for a different species.
    code = getattr(transcript, "translation_table", None)
    table = codon_table_by_id(_reference_code(transcript) if code is None else code)
    declared = getattr(transcript, "start_codon", None)
    if declared is not None:
        if declared not in table.forward_table:
            raise ValueError("Unsupported declared initiation codon: %r" % declared)
        if declared not in table.start_codons:
            return replace(table, start_codons=table.start_codons | {declared})
    return table


def initiator_as_met(protein, first_codon, codon_table):
    """Read a retained, recognized initiation codon as Met.

    NCBI's ``Starts`` assignment is distinct from internal codon decoding.
    This changes only the first residue; callers establish that the codon
    is at the retained CDS start, rather than a downstream candidate start.
    """
    if protein and first_codon.upper() in codon_table.start_codons:
        return "M" + protein[1:]
    return protein


def translate_sequence(
        nucleotide_sequence,
        codon_table=STANDARD,
        to_stop=True,
        selenocysteine=()):
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
    selenocysteine : collection of int
        Offsets into ``nucleotide_sequence`` of TGA codons to read as
        selenocysteine (``'U'``), e.g. annotated Sec codons whose SECIS
        element is retained. See :mod:`varcode.effects.selenocysteine`.
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
        if codon == "TGA" and i in selenocysteine:
            result.append("U")
            continue
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
