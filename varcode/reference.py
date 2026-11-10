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

from __future__ import annotations

import os
from warnings import warn
import re
from typing import TYPE_CHECKING, cast

from pyensembl import EnsemblRelease, Genome
from pyensembl import cached_release as cached_ensembl_release
from pyensembl import genome_for_reference_name as get_genome_for_ensembl_reference_name

import pyensembl.species

from .common import memoize
from .pinned_references import check_supported, find_pinned_name, load_pinned_genome, pinned_name

if TYPE_CHECKING:
    from .genome import Genome as VarcodeGenome

canonical_reference_names = set(
    pyensembl.species.Species._reference_names_to_species.keys())

def normalize_reference_name(name):
    """
    Simplify potentially variable ways to write out a reference name
    by lowercasing it and replacing all dashes and spaces with underscores.

    Parameters
    ----------
    name : str

    Returns
    -------
    str
    """
    return name.strip().lower().replace("-", "_").replace(" ", "_")

normalized_canonical_reference_names = {
    normalize_reference_name(name) for name in canonical_reference_names}


# this dict contains reference name aliases which preserve contig names
ensembl_reference_aliases = {
    # human assemblies
    "NCBI36": ["B36", "GRCh36"],
    "GRCh37": ["B37", "NCBI37"],
    "GRCh38": ["B38", "NCBI38"],

    # Mouse assemblies, as NCBI Assembly's paired GenBank (GCA) and RefSeq
    # (GCF) accessions. Only versioned accessions name one build: the
    # unversioned GCA_000001635 is the series of every mouse build.
    "NCBIM37": [
        "GCA_000001635.1", "GCF_000001635.18",  # MGSCv37
    ],
    "GRCm38": [
        "GCA_000001635.2", "GCF_000001635.20",  # GRCm38
        "GCA_000001635.3", "GCF_000001635.21",  # GRCm38.p1
        "GCA_000001635.4", "GCF_000001635.22",  # GRCm38.p2
        "GCA_000001635.5", "GCF_000001635.23",  # GRCm38.p3
        "GCA_000001635.6", "GCF_000001635.24",  # GRCm38.p4
        "GCA_000001635.7", "GCF_000001635.25",  # GRCm38.p5
        "GCA_000001635.8", "GCF_000001635.26",  # GRCm38.p6
    ],
    "GRCm39": [
        "GCA_000001635.9", "GCF_000001635.27",  # GRCm39
    ],
}

ucsc_to_ensembl_reference_names = {
    # mouse
    "mm9": "NCBIM37",  # UCSC mm9 is NCBI Build 37
    "mm10": "GRCm38",
    "mm39": "GRCm39",

    # human
    "hg18": "NCBI36",
    "hg19": "GRCh37",
    "hg38": "GRCh38",

    # cat
    "felCat9": "Felis_catus_9.0",
    "felCat8": "Felis_catus_8.0",

    # dog. canFam4 to canFam6 are breed assemblies (German Shepherd, Great
    # Dane, Boxer). chrM converts to MT, which their Ensembl builds don't
    # contain. Ensembl's default dog, ROS_Cfam_1.0, has no UCSC name.
    "canFam6": "Dog10K_Boxer_Tasha",
    "canFam5": "UMICH_Zoey_3.1",
    "canFam4": "UU_Cfam_GSD_1.0",
    "canFam3": "CanFam3.1",
    "canFam2": "CanFam2.0",
    "canFam1": "CanFam1.0",

    # rat
    "rn6": "Rnor_6.0",
    "rn5": "Rnor_5.0",

    # chicken
    "galGal6": "Gallus_gallus_6.0",
    "galGal5": "Gallus_gallus_5.0",
    "galGal4": "Gallus_gallus_4.0",

    # TODO: add Zebrafish to PyEnsembl
    #   "danRer11": "GRCz11",
    #   "danRer10": "GRCz10",

    # TODO: add cow to PyEnsembl
    #    "bosTau9": "ARS-UCD1.2",
    #    "bosTau8": "UMD3.1",

    # TODO: add pig to PyEnsembl
    #   "susScr11": "Sscrofa11.1",
    #   "susScr3": "Sscrofa10.2",
    #   "susScr2": "Sscrofa9.2",
}

ucsc_reference_names = set(ucsc_to_ensembl_reference_names.keys())

normalized_ucsc_reference_names = {
    normalize_reference_name(name)
    for name in ucsc_reference_names
}

ensembl_to_ucsc_reference_names = {
    v: k
    for (k, v) in ucsc_to_ensembl_reference_names.items()}

# Widely distributed human FASTA names that don't contain an assembly name,
# to the name inference returns. Builds named with UCSC contigs (chr1, chrM)
# return the UCSC name, as hg19 and hg38 FASTAs do. Matched within the file
# name, like the aliases above.
reference_file_aliases = {
    # 1000 Genomes phase 1 b37 (1..22, X, Y, MT: the rCRS).
    "human_g1k_v37": "GRCh37",
    # 1000 Genomes hs37d5: that build plus EBV and decoys. bwa.kit's hs37 is
    # hs37d5 without the decoys.
    "hs37": "GRCh37",
    "hs37d5": "GRCh37",
    # Broad/GATK b37, also named 1..22, X, Y, MT.
    "Homo_sapiens_assembly19": "GRCh37",
    # Broad/GATK hg38, and bwa.kit's hs38, hs38a and hs38DH (built from NCBI's
    # GRCh38 analysis set with UCSC IDs): chr1..chr22, chrX, chrY, chrM.
    "Homo_sapiens_assembly38": "hg38",
    "hs38": "hg38",
    "hs38a": "hg38",
    "hs38DH": "hg38",
}


@memoize
def is_ucsc_reference_name(name):
    """
    Is the given genome name one of the known UCSC genomes?

    Parameters
    ----------
    name : str

    Returns
    -------
    bool
    """
    return (normalize_reference_name(name) in normalized_ucsc_reference_names)

# merge the UCSC aliases, which are only an approximate correspondence between
# genomes (e.g. hg19 isn't exactly GRCh37) and the more straightforward
# aliases like B37->GRCh37

def _merge_ensembl_aliases_with_ucsc():
    """
    Combine the dictionary of Ensembl reference name aliases with
    the UCSC reference names (which also get included as aliases).

    Returns
    -------
    dict
    """
    result = {}
    for ucsc_name, ensembl_name in ucsc_to_ensembl_reference_names.items():
        result[ensembl_name] = [ucsc_name] + \
            ensembl_reference_aliases.get(ensembl_name, [])
    return result

alias_dict_with_ucsc = _merge_ensembl_aliases_with_ucsc()


def most_recent_assembly_name(assembly_names):
    """
    Given list of (in this case, matched) assemblies, identify the most recent,
    where "recency" is determined by sorting based on the numeric element of
    the assembly name.
    """
    match_recency = []
    for assembly_name in assembly_names:
        match = re.search(r'\d+', assembly_name)
        assert match is not None, "Expected a numbered assembly name: %s" % assembly_name
        match_recency.append(int(match.group()))
    sorted_list_of_names = [
        assembly
        for (number, assembly) in
        sorted(zip(match_recency, assembly_names), reverse=True)]
    most_recent = sorted_list_of_names[0]
    return most_recent


def choose_best_assembly_name(assembly_names):
    """
    Given a list of reference genome names returns the best according to the
    following criteria:
        1) Prefer Ensembl reference names to UCSC
        2) Prefer reference names with higher numbers in them.

    Parameters
    ----------
    assembly_names : list of str

    Returns
    -------
    str
    """
    assembly_names = set(assembly_names)

    if len(assembly_names) == 1:
        return list(assembly_names)[0]

    assembly_names_ucsc = {
        name for name in assembly_names if is_ucsc_reference_name(name)}
    assembly_names_ensembl = assembly_names.difference(assembly_names_ucsc)
    if assembly_names_ensembl:
        # drop the UCSC reference names and pick only between the Ensembl
        # compatible names
        return most_recent_assembly_name(assembly_names_ensembl)
    else:
        return most_recent_assembly_name(assembly_names_ucsc)

def _collect_candidate_matches(reference_name_or_path):
    """
    Generate list of assembly names which occur (or have any alias occur)
    within the input string.

    Parameters
    ----------
    reference_name_or_path : str

    Returns
    -------
    Two lists of strings:
        1) Assembly names which match the filename at the end of a path.
        2) Assembly names which occur anywhere in a path.

    If the input is a reference name and not a filename/path then the second
    list will be empty.
    """
    normalized_reference_name_or_path = \
        normalize_reference_name(reference_name_or_path)

    # if reference_name_or_path is an actual genome name, like GRCh37 then
    # reference_filename_lower will just be the lowercase of the genome name
    # ("grch37) but if it was a full path to a reference file, such as
    # '/path/to/GRCh37.fasta' then reference_filename_lower will be
    # 'grch37.fasta'
    normalized_reference_filename = os.path.basename(normalized_reference_name_or_path)

    file_name_matches = []
    full_path_matches = []
    for assembly_name in canonical_reference_names:
        candidate_list = [assembly_name] + alias_dict_with_ucsc.get(assembly_name, [])
        for candidate in candidate_list:
            normalized_candidate = normalize_reference_name(candidate)
            name_to_add = (
                candidate
                if is_ucsc_reference_name(candidate)
                else assembly_name)
            if normalized_candidate in normalized_reference_filename:
                file_name_matches.append(name_to_add)

            elif normalized_candidate in normalized_reference_name_or_path:
                full_path_matches.append(name_to_add)

    for file_alias, name in reference_file_aliases.items():
        normalized_alias = normalize_reference_name(file_alias)
        if normalized_alias in normalized_reference_filename:
            file_name_matches.append(name)
        elif normalized_alias in normalized_reference_name_or_path:
            full_path_matches.append(name)

    # Names of the genomes ORFdata pins, such as T2T-CHM13v2.0's hs1 and
    # chm13v2.0_maskedY_rCRS: the longest, as a whole word.
    pinned = find_pinned_name(os.path.basename(reference_name_or_path))
    if pinned is not None:
        file_name_matches.append(pinned)
    else:
        pinned = find_pinned_name(reference_name_or_path)
        if pinned is not None:
            full_path_matches.append(pinned)

    # remove duplicate matches (happens due to overlapping aliases)
    file_name_matches = list(set(file_name_matches))
    full_path_matches = list(set(full_path_matches))
    return file_name_matches, full_path_matches

@memoize
def infer_reference_name(reference_name_or_path):
    """
    Given a string containing a reference name (such as a path to
    that reference's FASTA file), return its canonical name
    as used by Ensembl.

    A name of a genome ORFdata pins is returned as ORFdata spells it, so
    ``chm13v2.0_maskedY_rCRS.fa`` gives ``chm13v2.0_maskedY_rCRS``: which
    T2T file it is says which mitochondrion its ``chrM`` is. The longest such
    name wins, each only as a whole word (``hs1.fa``, never ``ths1``).
    T2T-CHM13v1.0 and older drafts raise ``ValueError``, since their
    coordinates differ from T2T-CHM13v2.0's.
    """
    file_name_matches, full_path_matches = _collect_candidate_matches(
        reference_name_or_path)
    # A name of a genome ORFdata pins (T2T-CHM13v2.0's hs1, chm13v2.0*, ...)
    # is a whole file, so it is never one of several candidates.
    matches = file_name_matches or full_path_matches
    pinned = [name for name in matches if pinned_name(name) is not None]
    if pinned:
        if len(matches) > 1:
            raise ValueError(
                "%s names more than one reference (%s); pass the genome explicitly"
                % (reference_name_or_path, ", ".join(sorted(matches))))
        # T2T-CHM13v1.0 and older drafts have other coordinates.
        check_supported(pinned[0])
        if not file_name_matches:
            warn((
                'Reference could not be matched against filename ({}); '
                'using best match against full path ({}).').format(
                    reference_name_or_path, pinned[0]))
        return pinned[0]
    # given set of existing matches, choose one to return
    # (first select based on file_name, then full path. If multiples, use most recent)
    if len(file_name_matches) == 1:
        match = file_name_matches[0]
    elif len(file_name_matches) > 1:
        # separate logic for >1 vs 1 to give informative warning
        match = choose_best_assembly_name(file_name_matches)
        # A UCSC name (canFam3 inside CanFam3.1) is the same assembly, so
        # only matches naming different assemblies are ambiguous.
        assemblies = {
            ucsc_to_ensembl_reference_names.get(name, name)
            for name in file_name_matches}
        if len(assemblies) > 1:
            warn(
                ('More than one reference ({}) matches path in header ({}); '
                 'the most recent one ({}) was used.').format(
                    ','.join(file_name_matches), reference_name_or_path, match))
    elif len(full_path_matches) >= 1:
        # combine full-path logic since warning is the same
        match = choose_best_assembly_name(full_path_matches)
        warn((
            'Reference could not be matched against filename ({}); '
            'using best match against full path ({}).').format(
                reference_name_or_path, match))
    else:
        raise ValueError(
            "Failed to infer genome assembly name for %s" % reference_name_or_path)
    return match


# Numbered releases and annotation dates. Paths never match.
_RELEASE_SUFFIX = re.compile(r"^\s*([^:/\\\s]+?)\s*([:._\-\s])\s*(\d{4}_\d{2}|\d+)\s*$")


def _exact_reference_names():
    """Each known assembly name and alias, normalized, to the name inference returns."""
    names = {}
    for assembly_name in canonical_reference_names:
        for candidate in [assembly_name] + alias_dict_with_ucsc.get(assembly_name, []):
            names[normalize_reference_name(candidate)] = (
                candidate if is_ucsc_reference_name(candidate) else assembly_name)
    names.update((normalize_reference_name(alias), name)
                 for alias, name in reference_file_aliases.items())
    return names


_exact_names = _exact_reference_names()


def _split_release(reference_name):
    """
    The assembly name and requested Ensembl release in ``reference_name``,
    or ``(reference_name, None)`` when it names no release.

    A colon always introduces a release ("GRCh38:2026_04"). Other separators
    do only after an exact assembly name or alias, so names that end in
    digits, such as "Felis_catus_9.0", keep their meaning.
    """
    if normalize_reference_name(reference_name) in _exact_names:
        return reference_name, None
    match = _RELEASE_SUFFIX.match(reference_name)
    if match is None:
        if ":" in reference_name and not any(sep in reference_name for sep in ("/", "\\")):
            raise ValueError(
                "Explicit Ensembl release must be an integer or YYYY_MM: %r" % reference_name)
        return reference_name, None
    name, separator, release = match.groups()
    if separator != ":" and normalize_reference_name(name) not in _exact_names:
        return reference_name, None
    return name, release if "_" in release else int(release)


@memoize
def infer_genome_for_reference_name(reference_name):
    """
    First infer a canonical Ensembl or UCSC reference name
    (e.g. hg19 or GRCh37) and if it's a UCSC genome then map it
    to the equivalent in Ensembl.

    A name may end in a release to choose that Ensembl release, e.g.
    ``"GRCh38:93"``, ``"GRCh38:2026_04"`` or ``"GRCh37:2013_09"``.

    Names of the genomes ORFdata pins open ORFdata's installed annotation:
    ``T2T-CHM13v2.0``, ``hs1``, ``GCA_009914755.4``, ``GCF_009914755.1``, the
    T2T Consortium's ``chm13v2.0*`` analysis sets and T2T-CHM13v1.1 all open
    ``orfdata.load_genome("t2t-chm13")`` (``orfdata.reference_names()`` lists
    them). Its only release is the pinned one, ``"T2T-CHM13v2.0:2022_07"``.
    T2T-CHM13v1.0 and older drafts have other coordinates and raise.
    Selection is offline and validates the requested assembly. Otherwise the
    most recent locally installed numbered release is used, or the newest
    numbered release if none is installed.

    Parameters
    ----------
    reference_name : str

    Returns
    -------
    pyensembl.EnsemblRelease and a boolean flag indicating whether
    UCSC genome was converted to Ensembl.

    Raises
    ------
    ValueError
        If the name is not recognized, or the requested release does not
        provide that assembly.
    """
    converted_ucsc_to_ensembl = False
    reference_name, release = _split_release(reference_name)
    source = reference_name
    reference_name = infer_reference_name(reference_name)
    if pinned_name(reference_name) is not None:
        # A genome ORFdata pins, such as T2T-CHM13v2.0. Its files spell
        # contigs as UCSC or NCBI do, so they are converted, as for hg19.
        return load_pinned_genome(reference_name, release, source=source), True
    if is_ucsc_reference_name(reference_name):
        if reference_name not in ucsc_to_ensembl_reference_names:
            raise ValueError(
                "Unrecognized UCSC reference name '%s'" % reference_name)
        reference_name = ucsc_to_ensembl_reference_names[reference_name]
        converted_ucsc_to_ensembl = True
    if release is None:
        genome = get_genome_for_ensembl_reference_name(reference_name)
    else:
        species = pyensembl.species.Species._reference_names_to_species[reference_name]
        # PyEnsembl resolves historical dated assemblies as well as numbered
        # releases, validates the selection and preserves its cache identity.
        genome = EnsemblRelease.cached(release, species=species, assembly=reference_name)
    return genome, converted_ucsc_to_ensembl

def infer_genome(
        genome_object_string_or_int: Genome | VarcodeGenome | str | int,
) -> tuple[Genome | VarcodeGenome, bool]:
    """
    If given an integer, get the human EnsemblRelease object for that
    Ensembl version.

    If given a string, return the most recent installed EnsemblRelease
    which has an equivalent reference, or the release the string names
    (e.g. "GRCh38:93" or "GRCh38:2026_04"). If the given name is a UCSC genome (e.g. hg19) then
    convert it to the equivalent Ensembl reference (e.g. GRCh37). A name of a
    genome ORFdata pins (e.g. "hs1" or "chm13v2.0_maskedY_rCRS.fa" for
    T2T-CHM13v2.0) returns ORFdata's installed annotation; ORFdata's
    ``NotInstalledError`` says how to install it.

    If given a PyEnsembl or Varcode Genome, return that same object. Do not
    memoize object pass-through: native genomes with different attached DNA
    can compare equal. Integer/string resolvers below already cache results.

    Returns a pair of (Genome, bool) where the bool corresponds to whether
    the input requested a UCSC genome (e.g. "hg19") and an Ensembl (e.g. GRCh37)
    was returned as a substitute, or named a file of a pinned genome (e.g.
    "hs1"), whose contig names are converted the same way.
    """
    # Import at call time to avoid the Genome wrapper's import of this module.
    from .genome import Genome as VarcodeGenome

    converted_ucsc_to_ensembl = False
    # Known genome types must bypass duck-typing: probing .fasta can open
    # and index native DNA before a caller has asked for any sequence.
    if isinstance(genome_object_string_or_int, (Genome, VarcodeGenome)):
        genome = genome_object_string_or_int
    elif (hasattr(genome_object_string_or_int, "transcripts_at_locus")
            and hasattr(genome_object_string_or_int, "fasta")):
        # Preserve legacy duck-typed providers. The concrete public input
        # union cannot express this runtime-validated pass-through (#628).
        genome = cast("Genome | VarcodeGenome", genome_object_string_or_int)
    elif isinstance(genome_object_string_or_int, int):
        genome = cached_ensembl_release(genome_object_string_or_int)
    elif isinstance(genome_object_string_or_int, str):
        genome, converted_ucsc_to_ensembl = \
            infer_genome_for_reference_name(genome_object_string_or_int)
    else:
        raise TypeError(
            ("Expected genome to be an int, string, or pyensembl.Genome "
                "instance, got %s : %s") % (
                str(genome_object_string_or_int),
                type(genome_object_string_or_int)))
    return genome, converted_ucsc_to_ensembl
