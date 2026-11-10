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
python -m pip install --upgrade varcode 'orfdata>=0.52' 'pyensembl>=2.32.3'
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

### Reading T2T VCFs by name

`load_vcf` reads the assembly from a VCF's `##reference` line, and `genome=`,
`Variant(..., genome=)` and `varcode.Genome(...)` take the same names. Any name
ORFdata lists for T2T-CHM13v2.0 (`orfdata genomes names`) opens its installed
`t2t-chm13` annotation, so the two packages accept the same names:

| Name, found in a file name or path | Files | `chrM` is |
|---|---|---|
| `T2T-CHM13v2.0`, `GCA_009914755.4` | NCBI GenBank FASTA; Ensembl's DNA | CHM13's own, `MT` |
| `GCF_009914755.1` | NCBI RefSeq FASTA | none |
| `hs1` | UCSC | CHM13's own, `MT` |
| `chm13v2.0`, `chm13v2.0_maskedY`, `chm13v2.0_noY` | T2T Consortium analysis sets | CHM13's own, `MT` |
| `chm13v2.0_maskedY_rCRS`, `chm13v2.0_maskedY.rCRS` | the analysis set the T2T README recommends for short-read variant calling | the rCRS, `NC_012920.1` |
| `T2T-CHM13v1.1`, `GCA_009914755.3`, `chm13.draft_v1.1` | v2.0 without chrY | CHM13's own, `MT` |

The longest name in a file name wins, so `chm13v2.0_maskedY_rCRS.fa` is the
rCRS set, not `chm13v2.0`. `hs1` counts only as a whole word (`hs1.fa`,
`/refs/hs1/`). Bare `chm13` names every version and selects none.
T2T-CHM13v1.0 (`chm13.draft_v1.0`, `GCA_009914755.2`) and older drafts raise:
their coordinates differ. The only release is the pinned one,
`T2T-CHM13v2.0:2022_07`. If the annotation isn't installed, the error says which
name it found and gives the `orfdata genomes install t2t-chm13` command.

```python
from varcode import load_vcf
variants = load_vcf("sample.vcf.gz")  # ##reference=.../chm13v2.0_maskedY_rCRS.fa
print(variants.effects())
```

**Contigs.** Contig names convert by ORFdata's table: `chr1`, `CP068277.2` and
`NC_060925.1` are all `1`. If the header's `##contig` lengths disagree with
T2T-CHM13v2.0's, `load_vcf` raises rather than put another assembly's positions
on T2T transcripts. A GRCh38 VCF is the usual case: chr1 is 248,956,422 bp in
GRCh38 and 248,387,328 bp in T2T.

**Mitochondrion.** A T2T `chrM` is one of two molecules, both 16,569 bp, and
only the file name tells which:

- **CHM13's own** (CP068254.1) is in every file but the `_rCRS` set. Its
  position 1 is what the rCRS numbers 577. Ensembl's annotation has no
  mitochondrial genes, so its variants, on contig `MT`, are kept as failures,
  with one warning that says why.
- **The rCRS** (NC_012920.1) is the `_rCRS` set's `chrM`, and it is GRCh38's
  `MT`. Its variants go on contig `NC_012920.1`, which ORFdata annotates with
  Ensembl GRCh38's MT genes and the vertebrate mitochondrial code. They get the
  same consequences as GRCh38's `MT`: 8993T>G is MT-ATP6 p.L156R.
  `Variant.chrm_is_rcrs` is set; pass `chrm_is_rcrs=True` yourself when a
  `chrM` is the rCRS but no file name says so.

The rCRS is the one exception to keeping CHM13 and GRCh38 annotation apart,
because it is the same sequence in the same coordinates. Don't combine CHM13
variants with GRCh38 transcripts anywhere else, such as to fill newly resolved
regions; renaming does not convert coordinates.

**Chromosome DNA.** Ensembl's DNA uses `1`…`22`, `X`, `Y`, matching its
annotation. PyEnsembl's own `genome_fasta_path_or_url`, used at installation
above, reads names exactly, so give it Ensembl's DNA. To use another T2T FASTA,
attach it with `varcode.Genome(genome, fasta="chm13v2.0.fa")`. That reads each
annotation contig by the FASTA's spelling (`22` from `chr22` or `CP068256.2`).
A FASTA's `chrM` is paired with `MT` or `NC_012920.1` only when its file name
says which mitochondrion it is.

For local DNA, pass its path at installation and load time. Python's `cache_root=`
selects a cache directory. The genome works with splice, structural and supplied-ORF
APIs and with `varcode.Genome`. Missing sequence stays unresolved.

### What is converted, and what isn't

Three different things move between GRCh38 and T2T. Varcode itself does none of them:

- **Annotation.** Ensembl made this annotation by projecting GRCh38 gene
  models onto CHM13. Each GTF record names its GRCh38 parent; SLC5A9's
  `ENST05220073379` carries `projection_parent_transcript "ENST00000533824.5"`.
  A projected model can differ from its parent, so a consequence on CHM13 can
  differ from GRCh38's for the same gene. The annotation names genes but not
  transcripts, and it has no mitochondrial genes.
- **ORF placement.** With `extra_orfs`, catalogues made on GRCh38 reach T2T
  through ORFdata's pinned GRCh38-to-T2T-CHM13v2.0 chain file
  (`orfdata chains download`). An ORF is lifted only when each of its coding
  blocks lies inside one ungapped alignment block. Placement then checks
  the exact coding path and protein again, and leaves out any ORF whose
  sequence changed.
- **Variants.** Varcode converts no variant coordinates
  ([#243](https://github.com/openvax/varcode/issues/243)). Renaming contigs
  (`chr1` or `CP068277.2` to `1`), choosing `genome=` and attaching a FASTA
  change names, not positions. Lift GRCh38 calls to T2T before annotating
  them, and check their REF against T2T.

These are refused, not annotated:

- A VCF whose `##reference` names another assembly, read with `genome=` a
  T2T genome, such as `GRCh38`, `Homo_sapiens_assembly38.fasta`, `hg38.fa` or
  T2T-CHM13v1.0. A T2T VCF read with `genome=` GRCh38 is refused too. If the
  header is wrong, pass `reference_vcf_key=None` to ignore it.
- A VCF whose `##contig` lengths aren't T2T's.
- An allele whose REF disagrees with T2T's transcript or DNA raises
  `ReferenceMismatchError` when it is annotated.
- GRCh38's unplaced scaffolds and EBV (`chrUn_KI270302v1`, `KI270302.1`,
  `chrEBV`), which T2T doesn't have. `Variant.effects()` raises
  `InvalidContigError`, and a collection keeps them as failures, with a warning.

A GRCh38 position whose REF happens to match T2T's base can't be told apart by
its allele. Only the header can say which assembly the VCF is on.

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

UCSC dog names select their Ensembl assembly and convert `chr` contig names,
as `mm10` and `mm39` do for mouse: `canFam3` is CanFam3.1, `canFam4` is
UU_Cfam_GSD_1.0 (German Shepherd), `canFam5` is UMICH_Zoey_3.1 (Great Dane)
and `canFam6` is Dog10K_Boxer_Tasha (Boxer). `chrM` becomes `MT`, which
Ensembl's builds of the last three don't contain. ROS_Cfam_1.0, Ensembl's
default dog since release 105, has no UCSC name.

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
