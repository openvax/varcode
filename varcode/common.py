from __future__ import print_function, division, absolute_import

from collections import defaultdict
import pyfaidx

def reverse_complement(x):
    """
    Reverse complement of a nucleotide string.

    Parameters
    ----------

    x : str
        Original nucleotide string
    """
    return pyfaidx.Sequence(seq=x).reverse.complement.seq


def trim_shared_prefix(ref, alt):
    """
    Sometimes mutations are given with a shared prefix between the reference
    and alternate strings. Examples: C>CT (nucleotides) or GYFP>G (amino acids).

    This function trims the common prefix and returns the disjoint ref
    and alt strings, along with the shared prefix.
    """
    n_ref = len(ref)
    n_alt = len(alt)
    n_min = min(n_ref, n_alt)
    i = 0
    while i < n_min and ref[i] == alt[i]:
        i += 1

    # guaranteed that ref and alt agree on all the characters
    # up to i'th position, so it doesn't matter which one we pull
    # the prefix out of
    prefix = ref[:i]
    ref_suffix = ref[i:]
    alt_suffix = alt[i:]
    return ref_suffix, alt_suffix, prefix

def trim_shared_suffix(ref, alt):
    """
    Reuse the `trim_shared_prefix` function above to implement similar
    functionality for string suffixes.

    Given ref='ABC' and alt='BC', we first revese both strings:
        reverse_ref = 'CBA'
        reverse_alt = 'CB'
    and then the result of calling trim_shared_prefix will be:
        ('A', '', 'CB')
    We then reverse all three of the result strings to get back
    the shared suffix and both prefixes leading up to it:
        ('A', '', 'BC')
    """
    reverse_ref = ref[::-1]
    reverse_alt = alt[::-1]
    results = trim_shared_prefix(reverse_ref, reverse_alt)
    return tuple(s[::-1] for s in results)

def trim_shared_flanking_strings(ref, alt):
    """
    Given two nucleotide or amino acid strings, identify
    if they have a common prefix, a common suffix, and return
    their unique components along with the prefix and suffix.

    For example, if the input ref = "SYFFQGR" and alt = "SYMLLFIFQGR"
    then the result will be:
        ("F", "MLLFI", "SY", "FQGR")
    """
    ref, alt, prefix = trim_shared_prefix(ref, alt)
    ref, alt, suffix = trim_shared_suffix(ref, alt)
    return ref, alt, prefix, suffix

def group_by(records, field_name):
    # create an empty list for every new key
    groups = defaultdict(list)
    for record in records:
        value = getattr(record, field_name)
        if value is not None:
            groups[value].append(record)
    return dict(groups)
