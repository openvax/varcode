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

from __future__ import print_function, division, absolute_import

import os
from warnings import warn
import re

from pyensembl import Genome
from pyensembl import cached_release as cached_ensembl_release
from pyensembl import genome_for_reference_name as get_genome_for_ensembl_reference_name

import pyensembl.species
from typechecks import is_string, is_integer

from .common import memoize

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

    # mouse assemblies
    "GRCm38": [
        "GCF_000001635.24",  # GRCm38.p4
        "GCF_000001635.23",  # GRCm38.p3
        "GCF_000001635.22",  # GRCm38.p2
        "GCF_000001635.21",  # GRCm38.p1
        "GCF_000001635.20",  # GRCm38
        "GCA_000001635", # TODO: what's the difference between GCA_ and GCF_?
    ]
}

ucsc_to_ensembl_reference_names = {
    # mouse
    "mm9": "GRCm37",
    "mm10": "GRCm38",

    # human
    "hg18": "NCBI36",
    "hg19": "GRCh37",
    "hg38": "GRCh38",

    # cat
    "felCat9": "Felis_catus_9.0",
    "felCat8": "Felis_catus_8.0",

    # dog
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
    match_recency = [
        int(re.search(r'\d+', assembly_name).group())
        for assembly_name in assembly_names
    ]
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
    """
    file_name_matches, full_path_matches = _collect_candidate_matches(
        reference_name_or_path)
    # given set of existing matches, choose one to return
    # (first select based on file_name, then full path. If multiples, use most recent)
    if len(file_name_matches) == 1:
        match = file_name_matches[0]
    elif len(file_name_matches) > 1:
        # separate logic for >1 vs 1 to give informative warning
        match = choose_best_assembly_name(file_name_matches)
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


# An assembly name followed by an Ensembl release: "GRCh38:93", or
# "GRCh38.93" as in Ensembl's file names. Paths (which contain "/") never match.
_RELEASE_SUFFIX = re.compile(r"^\s*([^:/\\\s]+?)\s*([:._\-\s])\s*(\d+)\s*$")


def _exact_reference_names():
    """Each known assembly name and alias, normalized, to the name inference returns."""
    names = {}
    for assembly_name in canonical_reference_names:
        for candidate in [assembly_name] + alias_dict_with_ucsc.get(assembly_name, []):
            names[normalize_reference_name(candidate)] = (
                candidate if is_ucsc_reference_name(candidate) else assembly_name)
    return names


_exact_names = _exact_reference_names()


def _split_release(reference_name):
    """
    The assembly name and requested Ensembl release in ``reference_name``,
    or ``(reference_name, None)`` when it names no release.

    A colon always introduces a release ("GRCh38.p13:93"). Other separators
    do only after an exact assembly name or alias, so names that end in
    digits, such as "Felis_catus_9.0", keep their meaning.
    """
    if normalize_reference_name(reference_name) in _exact_names:
        return reference_name, None
    match = _RELEASE_SUFFIX.match(reference_name)
    if match is None:
        return reference_name, None
    name, separator, release = match.groups()
    if separator != ":" and normalize_reference_name(name) not in _exact_names:
        return reference_name, None
    return name, int(release)


@memoize
def infer_genome_for_reference_name(reference_name):
    """
    First infer a canonical Ensembl or UCSC reference name
    (e.g. hg19 or GRCh37) and if it's a UCSC genome then map it
    to the equivalent in Ensembl.

    A name may end in a release to choose that Ensembl release, e.g.
    ``"GRCh38:93"``, ``"GRCh38.93"`` or ``"hg19:75"``. Otherwise the most
    recent locally installed release of the assembly is used, or the newest
    release if none is installed.

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
    reference_name = infer_reference_name(reference_name)
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
        # Raises ValueError when the species has no genome in that release.
        provided = species.which_reference(release)
        if provided != reference_name:
            raise ValueError(
                "Ensembl release %d of %s provides %s, not %s" % (
                    release, species.latin_name, provided, reference_name))
        genome = cached_ensembl_release(release, species=species)
    return genome, converted_ucsc_to_ensembl

@memoize
def infer_genome(genome_object_string_or_int):
    """
    If given an integer, get the human EnsemblRelease object for that
    Ensembl version.

    If given a string, return the most recent installed EnsemblRelease
    which has an equivalent reference, or the release the string names
    (e.g. "GRCh38:93"). If the given name is a UCSC genome (e.g. hg19) then
    convert it to the equivalent Ensembl reference (e.g. GRCh37).

    If given a PyEnsembl Genome, simply use it.

    Returns a pair of (Genome, bool) where the bool corresponds to whether
    the input requested a UCSC genome (e.g. "hg19") and an Ensembl (e.g. GRCh37)
    was returned as a substitute.
    """
    converted_ucsc_to_ensembl = False
    # varcode.Genome wraps a pyensembl Genome with optional FASTA.
    # Recognize it via duck-typing (has ``transcripts_at_locus`` and a
    # ``fasta`` attribute) so we don't have to import it here and
    # introduce a cycle.
    if (hasattr(genome_object_string_or_int, "transcripts_at_locus")
            and hasattr(genome_object_string_or_int, "fasta")):
        genome = genome_object_string_or_int
    elif isinstance(genome_object_string_or_int, Genome):
        genome = genome_object_string_or_int
    elif is_integer(genome_object_string_or_int):
        genome = cached_ensembl_release(genome_object_string_or_int)
    elif is_string(genome_object_string_or_int):
        genome, converted_ucsc_to_ensembl = \
            infer_genome_for_reference_name(genome_object_string_or_int)
    else:
        raise TypeError(
            ("Expected genome to be an int, string, or pyensembl.Genome "
                "instance, got %s : %s") % (
                str(genome_object_string_or_int),
                type(genome_object_string_or_int)))
    return genome, converted_ucsc_to_ensembl