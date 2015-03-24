Varcode
=======

Varcode is a Python library for predicting the impact of genomic variants on protein sequences. 

Effect Types 
------------
  - Intergenic
  - Intragenic
  - IncompleteTranscript
  - NoncodingTranscript
  - Intronic
  - ThreePrimeUTR
  - FivePrimeUTR
  - Silent
  - Substitution
  - Insertion
  - Deletion
  - ComplexSubstitution
  - AlternateStartCodon
  - IntronicSpliceSite
  - ExonicSpliceSite
  - StopLoss
  - SpliceDonor
  - SpliceAcceptor
  - PrematureStop
  - FrameShiftTruncation
  - StartLoss
  - FrameShift
  - ExonLoss
   
  
Coordinate System
-----------------
Varcode currently uses a "base counted, one start" genomic coordinate system, to match the Ensembl annotation database. We are planning to switch over to "space counted, zero start" (interbase) coordinates, since that system allows for more uniform logic (no special cases for insertions). To learn more about genomic coordinate systems, read this [blog post](http://alternateallele.blogspot.com/2012/03/genome-coordinate-conventions.html). 




