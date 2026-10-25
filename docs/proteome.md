# Protein collections and FASTA

`to_proteome()` converts Varcode's native results into ORFdata's `Proteome`.
Use it to combine predicted proteins with a reference proteome, export FASTA,
group equal sequences, or save the predictions and their provenance together.
Install the optional integration with `pip install 'varcode[orf]'`.

```python
import varcode
from orfdata import Proteome

variants = varcode.load_vcf("tumor.vcf", genome="GRCh38:95")
effects = variants.effects()
predicted = varcode.to_proteome(effects, sample_id="tumor-1")

predicted.write_fasta("predicted.fasta")
predicted.save("predicted.proteome.json.gz")
restored = Proteome.load("predicted.proteome.json.gz")  # offline
```

`sample_id` and `run_id` are optional. Conversion reads the supplied native
results and their reference transcript properties; it does not rerun variant
annotation. Loading the saved collection requires neither a genome database
nor the original ORF catalogue, and never reconstructs Python objects from
the embedded native archive.

## Alternatives and missing products

The adapter accepts an `EffectCollection` or an iterable of effects,
`EffectCandidate`, `RealizedEffectCandidate`, or `ClassifiedOutcome` objects.
Every realized outcome is visited, including full proteins that differ outside
the local windows used to group candidates. Identical mutant sequences with
different baselines, splice assumptions, phase, or evidence remain separate.

For standalone ORF annotations, use the explicit input below. It retains failed
bindings whose `effect` is `None`, as well as the successful predictions.

```python
predicted = varcode.to_proteome(orf_annotations=annotations)
```

The collection's `records` contain `varcode.proteome.v1` events and comparisons.
Its `links` connect events to comparisons and comparisons to protein entries.
Comparison payloads retain baseline and mutant products, focal variants,
hypotheses, affected intervals, coordinate conventions, and candidate context.
Baseline proteins need not be members of the search collection. Native graphs
retain reference models, edits, ordered fusion segments and evidence.

Available sequences become entries. Missing sequences, known empty products,
annotation failures, failed ORF bindings, unresolved results and hypothesis
limits have distinct statuses and omission reasons. Explicit partial boundaries
and unknown boundaries remain distinct. The manifest describes the supplied
results; it does not claim exhaustive hypothesis enumeration. Ordinal rank,
probability, and phase mass remain separate fields. RNA or reference-ORF
evidence is not relabeled as experimental support for a mutant protein.

## Sequence grouping and composition

```python
sources_by_sequence = predicted.by_sequence()
sources_by_accession = predicted.by_accession()

# The default FASTA keeps one entry per distinct source prediction.
# Opt into exact sequence grouping and retain its token-to-source mapping:
predicted.write_fasta("unique.fasta", deduplicate=True)
token_to_sources = predicted.fasta_index(deduplicate=True)

reference = Proteome.from_genome(variants[0].genome)
combined = reference.merge(predicted)
combined.save("combined.proteome.json.gz")
```

The same `merge()` accepts an ORFdata `Proteome` containing Isovar producer
records. An Isovar adapter is tracked separately in
[Isovar #542](https://github.com/openvax/isovar/issues/542); this API does not
require that unfinished integration. The Varcode tests compose reference,
Varcode, and Isovar-shaped records through the shared interface.

Known conflicting taxa, assemblies, annotation contexts and source identities
raise errors. A collection with missing or mixed source references remains
unscoped. Merging it with a reference requires an explicit shared policy, such
as `reference.merge(predicted, compatibility="allow_unknown")`; this does not
override known conflicts or assign the host's taxon to unknown sources.

FASTA uses ORFdata's readable tokens and headers, with available origin,
producer, gene, transcript/ORF identifiers, boundaries and change attribution.
Supplied amino acids are exported unchanged, including any terminal `*`.
Keep the snapshot alongside FASTA: a header cannot represent every baseline
comparison, alternative hypothesis or evidence record. Unsupported nonportable
extension metadata raises an error instead of disappearing during conversion.

## API

::: varcode.to_proteome
