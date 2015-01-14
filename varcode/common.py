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