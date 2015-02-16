
from varcode import VariantAnnotator, Variant, CodingSequenceMutation

annot = VariantAnnotator(ensembl_release=75)


"""
def test_peptide_from_transcript_variant_CASP9():
    transcript_id = 'ENST00000333868'
    peptide, start, stop, annot = \
        transcript_variant.peptide_from_transcript_variant(
            transcript_id,
            15820483,
            ref='C',
            alt='G',
            padding = None,
            max_length = None)
    assert peptide is not None
    assert len(peptide) == 416
    # TODO: actually look up what this variant ought to be

def test_peptide_from_transcript_variant_CASP9_ref():
    transcript_id = 'ENST00000424908'
    peptide, start, stop, annot = \
        transcript_variant.peptide_from_transcript_variant(
            transcript_id,
            15820483,
            ref='C',
            alt='G',
            padding = None,
            max_length = None)
    assert peptide is not None
    assert len(peptide) == 198
    assert peptide[0] == 'X'

def test_peptide_from_transcript_variant_VWA3A():
    transcript_id = 'ENST00000389397'
    peptide, start, stop, annot = \
        transcript_variant.peptide_from_transcript_variant(
            transcript_id,
            22128096,
            ref='G',
            alt='A',
            padding = None,
            max_length = None)
    # expected mutation in the UTR to yield None
    assert peptide is None

def test_peptide_from_transcript_variant_RET():
    transcript_id = 'ENST00000340058'
    peptide, start, stop, annot = \
        transcript_variant.peptide_from_transcript_variant(
            transcript_id,
            43617416,
            ref='T',
            alt='C',
            padding = None,
            max_length = None)
    assert peptide is not None

def test_peptide_from_transcript_variant_PAR():

    transcript_id = 'ENST00000371279'
    peptide, start, stop, annot = \
        transcript_variant.peptide_from_transcript_variant(
            transcript_id,
            55224569,
            ref='T',
            alt='G',
            padding = None,
            max_length = None)
    assert peptide is not None
"""