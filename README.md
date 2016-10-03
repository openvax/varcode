<a href="https://travis-ci.org/hammerlab/varcode">
    <img src="https://travis-ci.org/hammerlab/varcode.svg?branch=master" alt="Build Status" />
</a>
<a href="https://coveralls.io/github/hammerlab/varcode?branch=master">
    <img src="https://coveralls.io/repos/hammerlab/varcode/badge.svg?branch=master&service=github" alt="Coverage Status" />
</a>
<a href="https://zenodo.org/badge/latestdoi/18834/hammerlab/varcode">
    <img src="https://zenodo.org/badge/18834/hammerlab/varcode.svg" alt="DOI" />
</a>
<a href="https://pypi.python.org/pypi/varcode/">
    <img src="https://img.shields.io/pypi/v/varcode.svg?maxAge=1000" alt="PyPI" />
</a>

Varcode
=======

Varcode is a library for working with genomic variant data in Python and predicting the impact of those variants on protein sequences.

Installation
------------

You can install varcode using [pip](https://pip.pypa.io/en/latest/quickstart.html):

```bash
pip install varcode
```

Optionally, you can pre-populate metadata caches through [PyEnsembl](https://github.com/hammerlab/pyensembl) as follows:

```bash
# Downloads and installs the Ensembl releases (75 and 76)
pyensembl install --release 75 76
```

This will eliminate a potential delay of several minutes required to install the relevant data
when using the Varcode for the first time.


Example
-------


```python
import varcode

# Load TCGA MAF containing variants from their
variants = varcode.load_maf("tcga-ovarian-cancer-variants.maf")

print(variants)
### <VariantCollection from 'tcga-ovarian-cancer-variants.maf' with 6428 elements>
###  -- Variant(contig=1, start=69538, ref=G, alt=A, genome=GRCh37)
###  -- Variant(contig=1, start=881892, ref=T, alt=G, genome=GRCh37)
###  -- Variant(contig=1, start=3389714, ref=G, alt=A, genome=GRCh37)
###  -- Variant(contig=1, start=3624325, ref=G, alt=T, genome=GRCh37)
###  ...

# you can index into a VariantCollection and get back a Variant object
variant = variants[0]

# groupby_gene_name returns a dictionary whose keys are gene names
# and whose values are themselves VariantCollections
gene_groups = variants.groupby_gene_name()

# get variants which affect the TP53 gene
TP53_variants = gene_groups["TP53"]

# predict protein coding effect of every TP53 variant on
# each transcript of the TP53 gene
TP53_effects = TP53_variants.effects()

print(TP53_effects)
### <EffectCollection with 789 elements>
### -- PrematureStop(variant=chr17 g.7574003G>A, transcript_name=TP53-001, transcript_id=ENST00000269305, effect_description=p.R342*)
### -- ThreePrimeUTR(variant=chr17 g.7574003G>A, transcript_name=TP53-005, transcript_id=ENST00000420246)
### -- PrematureStop(variant=chr17 g.7574003G>A, transcript_name=TP53-002, transcript_id=ENST00000445888, effect_description=p.R342*)
### -- FrameShift(variant=chr17 g.7574030_7574030delG, transcript_name=TP53-001, transcript_id=ENST00000269305, effect_description=p.R333fs)
### ...

premature_stop_effect = TP53_effects[0]

print(str(premature_stop_effect.mutant_protein_sequence))
### 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMF'

print(premature_stop_effect.aa_mutation_start_offset)
### 341

print(premature_stop_effect.transcript)
### Transcript(id=ENST00000269305, name=TP53-001, gene_name=TP53, biotype=protein_coding, location=17:7571720-7590856)

print(premature_stop_effect.gene.name)
### 'TP53'
```

If you are looking for a quick start guide, you can check out [this iPython book](./examples/varcode-quick_start.ipynb) that demonstrates simple use cases of Varcode

Effect Types
------------

Effect type  | Description
-----------: | :-----------
*AlternateStartCodon* | Replace annotated start codon with alternative  start codon (*e.g.* "ATG>CAG").
*ComplexSubstitution* | Insertion and deletion of multiple amino acids.
*Deletion* | Coding mutation which causes deletion of amino acid(s).
*ExonLoss* | Deletion of entire exon, significantly disrupts protein.
*ExonicSpliceSite* | Mutation at the beginning or end of an exon, may affect splicing.
*FivePrimeUTR* | Variant affects 5' untranslated region before start codon.
*FrameShiftTruncation* | A frameshift which leads immediately to a stop codon (no novel amino acids created).
*FrameShift* | Out-of-frame insertion or deletion of nucleotides, causes novel protein sequence and often premature stop codon.
*IncompleteTranscript* | Can't determine effect since transcript annotation is incomplete (often missing either the start or stop codon).
*Insertion* | Coding mutation which causes insertion of amino acid(s).
*Intergenic* | Occurs outside of any annotated gene.
*Intragenic* |Within the annotated boundaries of a gene but not in a region that's transcribed into pre-mRNA.
*IntronicSpliceSite* | Mutation near the beginning or end of an intron but less likely to affect splicing than donor/acceptor mutations.
*Intronic* | Variant occurs between exons and is unlikely to affect splicing.
*NoncodingTranscript* | Transcript doesn't code for a protein.
*PrematureStop* | Insertion of stop codon, truncates protein.
*Silent* | Mutation in coding sequence which does not change the amino acid sequence of the translated protein.
*SpliceAcceptor* | Mutation in the last two nucleotides of an intron, likely to affect splicing.
*SpliceDonor* | Mutation in the first two nucleotides of an intron, likely to affect splicing.
*StartLoss* | Mutation causes loss of start codon, likely result is that an alternate start codon will be used down-stream (possibly in a different frame).
*StopLoss* | Loss of stop codon, causes extension of protein by translation of nucleotides from 3' UTR.
*Substitution* | Coding mutation which causes simple substitution of one amino acid for another.
*ThreePrimeUTR* | Variant affects 3' untranslated region after stop codon of mRNA.


Coordinate System
-----------------
Varcode currently uses a "base counted, one start" genomic coordinate system, to match the Ensembl annotation database. We are planning to switch over to "space counted, zero start" (interbase) coordinates, since that system allows for more uniform logic (no special cases for insertions). To learn more about genomic coordinate systems, read this [blog post](http://alternateallele.blogspot.com/2012/03/genome-coordinate-conventions.html).
