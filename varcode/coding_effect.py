from collections import namedtuple

CodingEffect = namedtuple(
    "CodingEffect",
    [
        "OriginalProteinSequence",

        "MutantProteinSequence",

        # unless a frameshift results in a stop codon within the same
        # exon as the mutation, we leave the sequence incomplete
        # which is marked here
        "MutantSequenceComplete",

        "CodingVariantDescription",
    ]
)
