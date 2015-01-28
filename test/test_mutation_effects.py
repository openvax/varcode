from varcode import (
    VariantAnnotator,
    Variant,
    infer_transcript_effect,
    NoncodingTranscript,
    IncompleteTranscript,
    FivePrimeUTR,
    ThreePrimeUTR,
    Intronic,
    Silent,
    Insertion,
    Deletion,
    Substitution,
    StopLoss,
    StartLoss,
    PrematureStop,
    FrameShift,
    # TODO: SpliceDonor, SpliceReceptor
)

annot = VariantAnnotator(ensembl_release=77)

def test_incomplete():
    # transcript EGFR-009 (ENST00000450046 in Ensembl 77)
    # has an incomplete 3' end
    # chrom. 7 starting at 55,109,723
    # first exon begins: ATCATTCCTTTGGGCCTAGGA

    # change the first nucleotide of the 5' UTR A>T
    variant = Variant("7", 55109723, "A", "T")

    transcript = annot.ensembl.transcript_by_id("ENST00000450046")
    effect = infer_transcript_effect(variant, transcript)
    assert isinstance(effect, IncompleteTranscript), \
        "Expected %s on %s to be IncompleteTranscript, got %s" % (
            variant, transcript, effect)

def test_start_loss():
    # transcript EGFR-005 (ENST00000420316 in Ensembl 77)
    # location: chrom 7 @ 55,019,034-55,156,951 forward strand

    # CDS starts at position 244 of the first exon,
    # which is 55,019,034 + 244 of chr7 = 55019278
    # change the first nucleotide of the 5' UTR A>T
    # making what used to be a start codon into TTG (Leucine)
    variant = Variant("7", 55019278, "A", "T")
    transcript = annot.ensembl.transcript_by_id("ENST00000420316")
    effect = infer_transcript_effect(variant, transcript)
    assert isinstance(effect, StartLoss), \
        "Expected StartLoss, got %s for %s on %s" % (
            effect, variant, transcript, )