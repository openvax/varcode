# Reference ORF consequences

Install `varcode[orf]` for the optional [orfdata](https://github.com/openvax/orfdata)
reference layer. PyEnsembl supplies reference transcripts, orfdata supplies ORF
models and evidence, and Varcode compares those objects with altered products.
This includes upstream and overlapping ORFs, internal/downstream frames, and
ORFs on transcripts annotated as noncoding. The bundled catalogue is human
GRCh38; explicit species-specific catalogues use the same API.

## Annotate with ORFdata's ORFs for the variant's genome

Pass `extra_orfs="all"`, an evidence name or an ORFdata rule, and ORFdata
chooses the catalogues for the variant's genome and places each ORF on every
transcript that carries its exact coding path and protein:

```python
import orfdata
from pyensembl import cached_release
from varcode import Variant

variant = Variant("7", 66205486, "C", "G", genome=cached_release(115))
effects = variant.effects(extra_orfs="all")       # every catalogue ORF
effects = variant.effects(extra_orfs="ms")        # only MS-supported ones
effects = variant.effects(extra_orfs=orfdata.RiboSeq() & orfdata.MS())
```

| `extra_orfs=` | ORFs used |
| --- | --- |
| `None` (default) | none |
| `"all"` | every ORF ORFdata has for the genome |
| `"ribosome"` or `"riboseq"` | translated in ribosome profiling |
| `"ms"` | any supporting MS peptide |
| `"proteomics"` | peptides from digested cellular proteins |
| `"immunopeptidomics"` | peptides eluted from MHC molecules (HLA in human) |
| a list of names, e.g. `["proteomics", "immunopeptidomics"]` | ORFs with any of them |
| an ORFdata rule, e.g. `orfdata.RiboSeq() & orfdata.MS()` | ORFs meeting every part (`&`) or any part (`\|`) |
| an ORFdata catalogue or bound ORFs | exactly those |

`True` raises rather than guessing which of these you mean.

The ORFs are placed once per genome and reused by later calls and by
`VariantCollection.effects`. The same call works on other assemblies and
species once ORFdata has what they need, and fails with the command that
installs anything missing:

- human GRCh38: the bundled catalogue, nothing to install;
- human GRCh37, NCBI36 or T2T-CHM13: `orfdata setup --species human` for the
  chain files (and, for CHM13, `orfdata genomes install t2t-chm13` for the
  annotation);
- mouse: `pip install "orfdata[data,build-data]"`, then `orfdata setup
  --species mouse`, which builds the three mouse catalogues from their
  sources (add `--accept-terms` for UCSC's GRCm38 and NCBIM37 chains).

ORFs ORFdata could not place on the genome, for example because an older
annotation lacks their transcript, are left out rather than failing the
call; `varcode.extra_orfs.prepare_extra_orfs("all").excluded(genome)` lists
them with the reason.

## Annotate ordinary and extra ORFs together

Supply a catalogue to the usual effects API. Varcode finds the overlapping
records and resolves their transcripts using the variant's genome:

```python
from orfdata import load_catalog
from pyensembl import cached_release
from varcode import Variant

genome = cached_release(95)  # provision the reference before annotation
variant = Variant("7", 66205486, "C", "G", genome=genome)
effects = variant.effects(extra_orfs=load_catalog())
for annotation in effects.orf_annotations:
    print(annotation.orf_id, annotation.effect.short_description)
    for outcome in annotation.outcomes:
        if outcome.mutant is not None:
            print(outcome.mutant.protein_sequence)
```

This returns the ordinary UTR consequences alongside the `p.P2A` consequence in
`c7riboseqorf56`. Omitting `extra_orfs` retains ordinary annotation without loading
ORFdata. `extra_orfs` also accepts a bound reference ORF or an iterable of bound
ORFs, including ORFdata collections; `VariantCollection.effects(extra_orfs=...)`
reuses an iterable across all variants. Several ORFs on one RNA remain separate.
No MS threshold is imposed by Varcode.

Added entries are ordinary Varcode effect classes and participate in collection
sorting and filtering. Their `effect.orf` is the bound reference;
`effect.reference_transcript` is its original PyEnsembl RNA. The existing engine
uses a private coding view in `effect.transcript`; its ID distinguishes the ORF
for candidate grouping. Public `effect.transcript_id` and collection grouping /
expression lookups use the original RNA identity. `effect.orf_annotation` holds native products, source evidence,
and initiation-context changes. `effects.orf_annotations` selects these detailed
annotations in collection order, including failures. Export them with
`annotation.to_dict()` when JSON is needed. Generic effect/collection JSON
serialization rejects these entries because its schema cannot preserve their
coding models. Dataframe/CSV summaries add `orf_id`, `orf_source` and
`orf_unresolved_reason`; `from_csv` rejects ORF rows rather than recomputing them
against the ordinary CDS. Full generic round trips are tracked in #593.

Bound inputs must match the variant genome's reference and annotation context.
Catalogue records resolve against that genome and retain their source identities.
Reference-scope mismatches raise even with `raise_on_error=False`; per-record
binding/calculation errors raise by default and become `Failure` entries with
`orf_annotation.unresolved_reason` when that flag is false. ORFdata collection
exclusions overlapping the variant are reported in the same way.

The `annotator` argument selects ordinary transcript annotation. Extra ORFs reuse
Varcode's existing transcript-model engine, including patient germline hypotheses
and known-cis collection effects. DNA phase queries use the original RNA identity.
Observed RNA outcomes for the ordinary CDS cannot be directly reused for an extra
ORF: such observations raise or produce a `Failure` until explicitly mapped with
`annotate_orf_rna_product`, described below. This interface does not discover new
ORFs or infer IRES activity, and unsupported translation models remain explicit.

## Work with a selected reference and its products

```python
from pyensembl import cached_release
from orfdata import load_catalog
from varcode import Variant, annotate_orf_variants

genome = cached_release(95)  # explicitly provision this reference beforehand
transcript = genome.transcript_by_id("ENST00000304842")
catalog = load_catalog()
reference = catalog.get("c7riboseqorf56").bind(transcript)
variant = Variant("7", 66205486, "C", "G", genome=genome)
annotation = annotate_orf_variants([variant], reference)
assert annotation.reference is reference
for outcome in annotation.outcomes:
    print(outcome.effect)
    if outcome.mutant is not None:
        print(outcome.mutant.protein_sequence)
        print(outcome.mutant.initiation_context.start_codon)
```

Select a known ORF explicitly by ID as above. To annotate every ORF on a
transcript, iterate over `catalog.for_transcript(transcript)`; a transcript
can have zero, one or several ORFs.

```python
for reference in catalog.for_transcript(transcript):
    annotation = annotate_orf_variants([variant], reference)
    print(annotation.orf_id, annotation.effect)
```

The TPST1 +4 substitution is P2A in `c7riboseqorf56` and changes an initiation
context base. Canonical `variant.effects()` still reports its ordinary UTR
consequence. `variant.orf_effects()` discovers overlapping reference models and
reports unavailable/mismatched models as unresolved.

`ReferenceORF.model` holds the ORF identity, species, assembly, coding blocks and
protein; `.evidence` holds source observations, and `.transcript` is the original
PyEnsembl object. Varcode's `ORFAnnotation.reference` retains that base object;
`.effect` and `.outcomes` describe modifications, including alternative
phase/splice products. Each `ORFOutcome` carries native `.effect`, `.hypothesis`,
`.baseline` and `.mutant` objects. Products expose sequences, termination and
an orfdata `InitiationContext` directly. Unavailable structural products remain
`None` with an `.unresolved_reason`; failed reference binding is reported by
`annotation.unresolved_reason`. An empty outcome tuple does not prove no effect.
Use `annotation.to_dict()` only when exporting JSON. A private coding view
adapts this ORF to the existing transcript-model engine. It does not relabel
a noncoding PyEnsembl transcript.

Binding checks the entire spliced path, protein, stop, assembly and available
species identity. Literal genomic REF mismatches raise, including VCF padding.
Sequence unavailable for intronic edits remains unresolved. The initial binder
supports complete, stopped ORFs with NCBI translation table 1; other tables are
explicitly unsupported, never inferred from an MT-like contig name.

For RNA products, `annotate_orf_rna_product(reference, product)` accepts a
`RealizedTranscriptProduct` with one `(contig, position, strand)` or `None` in
`evidence['base_origins']` per cDNA base. It validates the retained reference
start and compares the conditional translation with the base ORF. Its
`ORFProductAnnotation` exposes `.reference`, `.effect`, `.product`, `.status`
and `.reason`. An unresolved result retains the supplied `.cdna_sequence`
and `.evidence` but has no translated `.product`. Isovar owns
read support and RNA reconstruction; Varcode owns the sequence consequence.

Reference peptides are projected through codon ancestry, retaining every copied
occurrence and identifying lost/altered peptides. Retention remains reference
MS evidence, never mutant translation or MHC presentation evidence. All Sid
TPST1–CRCP alternatives lack the native KIILSKISV peptide; each keeps its own RNA
support and unresolved polarity/initiation status.

The compatibility module `varcode.uorf` preserves the original upstream-specific
call shapes and `varcode[uorf]` remains an installation alias. New code should
use `varcode.orf` and the bound reference API above. `annotate_orf_product`
retains its original dictionary return for existing RNA exporters; new Python
code should use `annotate_orf_rna_product`. `ORFAnnotation.data` remains an
export alias; construct annotations with native fields rather than a `data`
dictionary. The model catalogue is not exhaustive ORF discovery: new starts
and downstream CDS regulation are unassessed.

## Reference identity

The ORF extras require ORFdata 0.2. A bound `ReferenceORF` retains its species,
strain, assembly, annotation version and sequence fingerprints in both native
and modified-product exports. RNA product metadata must agree with that context.
For custom Genomes without taxonomy, use `record.bind(transcript, context=...)`
before calling the generic annotation API. The compatibility uORF wrappers also
accept `reference_context=...`. See
[ORFdata's reference guide](https://github.com/openvax/orfdata/blob/main/docs/species.md).
