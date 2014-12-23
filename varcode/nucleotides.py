
# include all pseudonucleotides encoding repeats and uncertain bases
VALID_NUCLEOTIDES = {'A', 'C', 'T', 'G'}

EXTENDED_NUCLEOTIDES = {
    'A', 'C', 'T', 'G',
    'Y', # Pyrimidine (C or T)
    'R', # Purine (A or G)
    'W', # weak (A or T)
    'S', # strong (G or C)
    'K', # keto (T or G)
    'M', # amino (C or A)
    'D', # A, G, T (not C)
    'V', # A, C, G (not T)
    'H', # A, C, T (not G)
    'B', # C, G, T (not A)
    'X', # any base
    'N', # any base
}

def normalize_nucleotide_string(nucleotides, allow_extended_nucleotides=False):
    """
    Normalizes a nucleotide string by converting various ways of encoding empty
    strings into "", making all letters upper case, and checking to make sure
    all letters in the string are actually nucleotides.

    Parameters
    ----------
    nucleotides : str
        Sequence of nucleotides, e.g. "ACCTG"

    extended_nucleotides : bool
        Allow non-canonical nucleotide characters like 'X' for unknown base
    """
    # some MAF files represent deletions/insertions with NaN ref/alt values
    if isinstance(nucleotides, float) and np.isnan(nucleotides):
        return ""
    # VCFs sometimes have '.' ref/alt for insertions and deletions
    elif nucleotides == '.':
        return ""

    if not isinstance(nucleotides, (str, unicode)):
        raise TypeError(
                "Expected nucleotide string, got %s : %s" % (
                    nucleotides, type(nucleotides)))

    nucleotides = nucleotides.upper()

    for letter in set(nucleotides):
        if allow_extended_nucleotides:
            valid_nucleotides = EXTENDED_NUCLEOTIDES
        else:
            valid_nucleotides = VALID_NUCLEOTIDES
        if letter not in valid_nucleotides:
            raise ValueError(
                "Invalid character in nucleotide string: %s" % letter)

    return nucleotides
