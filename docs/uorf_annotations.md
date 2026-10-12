# Reference ORF consequences

Install `varcode[orf]` for the optional [orfdata](https://github.com/openvax/orfdata)
reference layer. PyEnsembl supplies reference transcripts, orfdata supplies ORF
models and evidence, and Varcode compares those objects with altered products.
This includes upstream and overlapping ORFs, internal/downstream frames, and
ORFs on transcripts annotated as noncoding. The bundled catalogue is human
GRCh38; explicit species-specific catalogues use the same API.

```python
from pyensembl import cached_release
from orfdata import load_catalog
from varcode import Variant, annotate_orf_variants

genome = cached_release(95)  # explicitly provision this reference beforehand
transcript = genome.transcript_by_id("ENST00000304842")
reference, = load_catalog().for_transcript(transcript)
variant = Variant("7", 66205486, "C", "G", genome=genome)
annotation = annotate_orf_variants([variant], reference)
assert annotation.reference is reference
result = annotation.to_dict()
```

The TPST1 +4 substitution is P2A in `c7riboseqorf56` and changes an initiation
context base. Canonical `variant.effects()` still reports its ordinary UTR
consequence. `variant.orf_effects()` discovers overlapping reference models and
reports unavailable/mismatched models as unresolved.

`ReferenceORF.model` holds the ORF identity, species, assembly, coding blocks and
protein; `.evidence` holds source observations, and `.transcript` is the original
PyEnsembl object. Varcode's `ORFAnnotation.reference` retains that base object;
`.effect` and exported `outcomes` describe modifications, including alternative
phase/splice products. A private coding view adapts this ORF to the existing
transcript-model engine. It does not relabel a noncoding PyEnsembl transcript.

Binding checks the entire spliced path, protein, stop, assembly and available
species identity. Literal genomic REF mismatches raise, including VCF padding.
Sequence unavailable for intronic edits remains unresolved. The initial binder
supports complete, stopped ORFs with NCBI translation table 1; other tables are
explicitly unsupported, never inferred from an MT-like contig name.

For RNA products, `annotate_orf_product(reference, product)` accepts a
`RealizedTranscriptProduct` with one `(contig, position, strand)` or `None` in
`evidence['base_origins']` per cDNA base. It validates the retained reference
start and compares the conditional translation with the base ORF. Isovar owns
read support and RNA reconstruction; Varcode owns the sequence consequence.

Reference peptides are projected through codon ancestry, retaining every copied
occurrence and identifying lost/altered peptides. Retention remains reference
MS evidence, never mutant translation or MHC presentation evidence. All Sid
TPST1–CRCP alternatives lack the native KIILSKISV peptide; each keeps its own RNA
support and unresolved polarity/initiation status.

The compatibility module `varcode.uorf` preserves the original upstream-specific
call shapes and `varcode[uorf]` remains an installation alias. New code should
use `varcode.orf` and the bound reference API above. The model catalogue is not
exhaustive ORF discovery: new starts and downstream CDS regulation are unassessed.

## Reference identity

The ORF extras require ORFdata 0.2. A bound `ReferenceORF` retains its species,
strain, assembly, annotation version and sequence fingerprints in both native
and modified-product exports. RNA product metadata must agree with that context.
For custom Genomes without taxonomy, use `record.bind(transcript, context=...)`
before calling the generic annotation API. The compatibility uORF wrappers also
accept `reference_context=...`. See
[ORFdata's reference guide](https://github.com/openvax/orfdata/blob/main/docs/species.md).
