# Copyright (c) 2016. Mount Sinai School of Medicine
#
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
from pyensembl import (
    Genome,
    cached_release,
    genome_for_reference_name,
)
from typechecks import is_string, is_integer
import os 
from warnings import warn
import re

# NCBI builds and hg releases aren't identical
# but the differences are all on chrM and unplaced contigs
reference_alias_dict = {
    # human assemblies
    "NCBI36": ["hg18", "B36", "NCBI36"],
    "GRCh37": ["hg19", "B37", "NCBI37"],
    "GRCh38": ["hg38", "B38", "NCBI38"],
    # mouse assemblies
    "GRCm37": ["mm9"],
    "GRCm38": [
        "mm10",
        "GCF_000001635.24",  # GRCm38.p4
        "GCF_000001635.23",  # GRCm38.p3
        "GCF_000001635.22",  # GRCm38.p2
        "GCF_000001635.21",  # GRCm38.p1
        "GCF_000001635.20",  # GRCm38
    ]
}

def infer_reference_name(reference_name_or_path):
    """
    Given a string containing a reference name (such as a path to
    that reference's FASTA file), return its canonical name
    as used by Ensembl.
    """
    # identify all cases where reference name or path matches candidate aliases
    reference_file_name = os.path.basename(reference_name_or_path)
    matches = {'file_name': list(), 'full_path': list()}
    for assembly_name in sorted(reference_alias_dict.keys(), reverse=True):
        candidate_list = [assembly_name] + reference_alias_dict[assembly_name]
        for candidate in candidate_list:
            if candidate.lower() in reference_file_name.lower():
                matches['file_name'].append(assembly_name)
            elif candidate.lower() in reference_name_or_path.lower():
                matches['full_path'].append(assembly_name)
    if len(matches['file_name']) == 1:
        match = matches['file_name'][0]
    elif len(matches['file_name']) >=1:
        # identify most recent matching reference
        match_recency = [int(re.search('\d+', assembly_name).group()) for assembly_name in matches['file_name']]
        most_recent = [x for (y,x) in sorted(zip(match_recency, matches['file_name']), reverse=True)][0] 
        warn('More than one reference matches path in header ({path});' +
            ' the most recent one ({selected}) was used.'.format(
                path=reference_name_or_path, selected=most_recent))
        match = most_recent
    elif len(matches['full_path']) >= 1:
        match_recency = [int(re.search('\d+', assembly_name).group()) for assembly_name in matches['file_name']]
        most_recent = [x for (y,x) in sorted(zip(match_recency, matches['file_name']), reverse=True)][0] 
        warn('Reference could not be matched against filename ({path});' +
            ' using the most recent match ({selected}) against the full path.'.format(
                path=reference_name_or_path, selected=most_recent))
        match = most_recent
    else:
        raise ValueError(
            "Failed to infer genome assembly name for %s" % reference_name_or_path)
    return match


def infer_genome(genome_object_string_or_int):
    """
    If given an integer, return associated human EnsemblRelease for that
    Ensembl version.

    If given a string, return latest EnsemblRelease which has a reference
    of the same name.

    If given a PyEnsembl Genome, simply return it.
    """
    if isinstance(genome_object_string_or_int, Genome):
        return genome_object_string_or_int
    if is_integer(genome_object_string_or_int):
        return cached_release(genome_object_string_or_int)
    elif is_string(genome_object_string_or_int):
        # first infer the canonical reference name, e.g. mapping hg19 -> GRCh37
        # and then get the associated PyEnsembl Genome object
        reference_name = infer_reference_name(genome_object_string_or_int)
        return genome_for_reference_name(reference_name)
    else:
        raise TypeError(
            ("Expected genome to be an int, string, or pyensembl.Genome "
                "instance, got %s : %s") % (
                str(genome_object_string_or_int),
                type(genome_object_string_or_int)))
