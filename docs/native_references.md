# Native mouse, dog and T2T-CHM13 references

Annotation needs a matched assembly, transcript annotation and sequence source.
An assembly name alone is not an annotation release. Changing `genome=` never
performs liftover; general variant conversion remains tracked in
[#243](https://github.com/openvax/varcode/issues/243).

## Install and annotate CHM13

Use Python 3.11+ in an isolated environment. ORFdata supplies the pinned Ensembl
rapid **2022_07** annotation for **T2T-CHM13v2.0, GCA_009914755.4**. PyEnsembl
2.32.3+ supplies chromosome DNA and saved taxonomy/translation metadata.

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade 'varcode[orf]' 'orfdata>=0.43.0' 'pyensembl>=2.32.3'
orfdata genomes install t2t-chm13 --genome-fasta \
  https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/genome/Homo_sapiens-GCA_009914755.4-unmasked.fa.gz
```

Installation downloads/indexes GTF, transcript/protein FASTA and the requested
chromosome FASTA. Allow several GB of disk space. Without `--genome-fasta`,
arbitrary intronic sequence is unavailable. Loading is offline.

```python
from orfdata import load_genome
from varcode import Variant, iter_candidate_effects

dna = (
    "https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/"
    "GCA_009914755.4/ensembl/genome/"
    "Homo_sapiens-GCA_009914755.4-unmasked.fa.gz"
)
genome = load_genome("t2t-chm13", genome_fasta_path_or_url=dna)
assert genome.sequence("22", 21798486, 21798486) == "T"
variant = Variant("22", 21798486, "T", "C", genome=genome)
effects = variant.effects()
for root in effects:
    for effect in iter_candidate_effects(root):
        if effect is None:
            continue  # Unresolved traversal, not an unchanged prediction.
        print(effect.transcript_id, effect.short_description,
              effect.mutant_protein_sequence)
effects.to_csv("chm13-effects.csv")
with open("chm13-effects.json", "w") as handle:
    handle.write(effects.to_json())
```

This native 1KGP call overlaps the newly resolved sequence track. Distinct
RIMBP3B and RIMBP3C placements must retain their identities even when their
proteins coincide. Predicted products are not evidence of mutant translation.

**Contig names must match.** Ensembl DNA uses `1`…`22`, `X`, `Y`, matching its
annotation. The consortium analysis FASTA uses `chr1`… and requires explicit,
unambiguous name remapping. Renaming does not convert coordinates. Do not combine
CHM13 variants with GRCh38 transcripts to fill newly resolved regions.

For local DNA, pass its path at installation and load time. Python's `cache_root=`
selects a cache directory. The genome works with splice, structural and supplied-ORF
APIs and with `varcode.Genome`. Missing sequence stays unresolved.

## Native mouse and dog

```sh
pyensembl install --release 115 --species mouse
pyensembl install --release 2020_10 --species dog
```

```python
from pyensembl import cached_release
from varcode import Variant
mouse = cached_release(115, species="mouse")  # GRCm39
dog = cached_release("2020_10", species="dog")  # ROS_Cfam_1.0
print(Variant("X", 105135710, "G", "A", mouse).effects())
print(Variant("12", 60469829, "C", "T", dog).effects())
```

Chromosome FASTA is a separate installation choice here too. Use matching
assembly and contig names. Nuclear/mitochondrial translation tables come from
saved reference metadata, not a universal standard-code assumption.

## Public validation and speed

The [native fixture manifest](https://github.com/openvax/varcode/blob/main/tests/data/native_references/README.md)
describes 100 variant/transcript cases, independent sequences and source limits.
Tests cover both strands, indels, actual intronic DNA, mouse mitochondrial
translation, the two cis changes of the mouse brindled allele, and CHM13
duplicated/newly resolved loci. A controlled CHM13 deletion checks structural
products; it is not a published patient rearrangement. Six further cases cover
public mouse/dog deletions and dog mitochondrial variants on the separately
pinned NCBI NC_002008.4 reference. ROS_Cfam 2020_10 has no mitochondrial
transcripts; the fixture preserves that distinction and source discrepancies.

The dog tumour audit retains all 101 published rows: 74 reproduce the published
amino-acid change with complete native models; four supply splice predictions.
Others retain model/position/format limitations. RefSeq and Ensembl consequences
can differ. The observed Ap2s1 cryptic junction remains outside the DNA-only
heuristic's coverage ([#675](https://github.com/openvax/varcode/issues/675)); RNA
evidence preserves its product and the limits of reference completion.

From a checkout with full references installed:

```sh
python examples/benchmark_native_references.py --genome chm13 \
  --genome-fasta /path/to/Homo_sapiens-GCA_009914755.4-unmasked.fa.gz
python examples/benchmark_native_references.py --genome mouse
python examples/benchmark_native_references.py --genome dog
```

The JSON report includes software/reference versions, distinct variants,
first-pass/warm batch times and exact-protein matches. Annotation includes all
overlapping transcripts; download/index time is excluded. This curated set is
not a whole-genome throughput estimate. `fast` remains compatibility behavior;
the default composes baseline, phase and splice products. See
[annotation choices](experimental_annotators.md).

## Export identity

Tables retain `reference_name`, `assembly_accession`, `species`, `annotation_name`
and `annotation_version` per row. Headers retain shared fields. Mixed files do
not borrow the first row's annotation. Exact ordinary Ensembl releases can be
reconstructed. Custom annotations need their matching installed genome:

```python
from varcode import EffectCollection
reloaded = EffectCollection.from_csv("chm13-effects.csv", genome=genome)
```

Conflicting row/header/supplied metadata raises. JSON/pickle retains the resolved
annotation instead of reinterpreting an assembly alias against a later default.
Filesystem-backed `varcode.Genome` wrappers retain their explicit DNA path;
arbitrary in-memory readers cannot be archived as portable installations. Saved
proteomes retain identity and product sequences independently of local files.
