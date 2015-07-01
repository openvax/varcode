Varcode
=======

Varcode is a library for working with genomic variant data in Python and predicting the impact of those variants on protein sequences.



Effect Types
------------
  `Intergenic`: Occurs outside of any annotated gene.

  `Intragenic` Within the annotated boundaries of a gene but not in a region that's transcribed into pre-mRNA.

  - *IncompleteTranscript*: Can't determine effect since transcript annotation is incomplete (often missing either the start or stop codon).
  - *NoncodingTranscript*: Transcript doesn't code for a protein.
  - *SpliceDonor*: Mutation in the first two nucleotides of an intron, likely to affect splicing.
  - *SpliceAcceptor*: Mutation in the last two nucleotides of an intron,
  likely to affect splicing.
  - *IntronicSpliceSite*: Mutation near the beginning or end of an intron but less likely to affect splicing than donor/acceptor mutations.
  - *ExonicSpliceSite*: Mutation at the beginning or end of an exon, may affect splicing.
  - *Intronic*: Variant occurs between exons and is unlikely to affect splicing.
  - *ThreePrimeUTR*: Variant affects 3' untranslated region after stop codon of mRNA.
  - *FivePrimeUTR*: Variant affects 5' untranslated region before start codon.
  - *Silent*: Mutation in coding sequence which does not change the amino acid sequence of the translated protein.
  - *Substitution*: Coding mutation which causes simple substitution of one amino acid for another.
  - *Insertion*: Coding mutation which causes insertion of amino acid(s).
  - *Deletion*: Coding mutation which causes deletion of amino acid(s).
  - *ComplexSubstitution*: Insertion and deletion of multiple amino acids.
  - *StartLoss*: Mutation causes loss of start codon, likely result is that an alternate start codon will be used down-stream (possibly in a different frame).
  - *AlternateStartCodon*: Replace annotated start codon with alternative start codon (e.g. ATG>CAG).
  - *StopLoss*: Loss of stop codon, causes extension of protein by translation of nucleotides from 3' UTR.
  - *PrematureStop*: Insertion of stop codon, truncates protein.
  - *FrameShift*: Out-of-frame insertion or deletion of nucleotides, causes novel protein sequence and often premature stop codon.
  - *FrameShiftTruncation*: A frameshift which leads immediately to a stop codon (no novel amino acids created).

  - FrameShift
  - ExonLoss


Coordinate System
-----------------
Varcode currently uses a "base counted, one start" genomic coordinate system, to match the Ensembl annotation database. We are planning to switch over to "space counted, zero start" (interbase) coordinates, since that system allows for more uniform logic (no special cases for insertions). To learn more about genomic coordinate systems, read this [blog post](http://alternateallele.blogspot.com/2012/03/genome-coordinate-conventions.html).




