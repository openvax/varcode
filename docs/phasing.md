# Phasing

Two nearby variants, for example in the same codon, can have a different
combined effect than either one alone, but only if they sit on the same copy
of the chromosome (*in cis*). Variants on different copies (*in trans*) act
separately. Phasing evidence tells Varcode which is the case.

Pass `phase_resolver=` when you have that evidence, from a phased VCF or RNA
reads. Varcode then predicts combined effects for variants in cis. To classify
somatic variants against the patient's own germline sequence, also pass
[germline context](germline.md).

## Phased VCF

Load the jointly phased calls and select the sample whose phase should be used:

```python
from varcode import load_vcf, VCFPhaseResolver

phased_variants = load_vcf("merged_phased.vcf", genome=81)
phaser = VCFPhaseResolver(phased_variants, sample="TUMOR")
effects = phased_variants.effects(phase_resolver=phaser)
```

The resolver uses phased `GT` and `PS` fields from the loaded collection.
Use one resolver per sample. Phase-set tags are local to their source VCF:
matching tags in independently produced tumor and normal files do not establish
relative phase. A joint call set or molecular evidence is needed.

Individual effects remain alongside joint results. The selected annotator owns
the joint prediction: the default returns `HaplotypeEffect` for supported point
edits, while `annotator="transcript_model"` returns classified phase/splice
candidates. Unsupported groups are retained as `Unresolved` with their
`variants` and `phase_source`; they are never silently dropped. The default
declines joint germline composition rather than ignoring the patient baseline.
See the [annotator contract](annotator_contract.md#joint-haplotypes).

For germline-aware somatic annotation, the resolver's collection should contain
both sets of alleles, but call `somatic_variants.effects(germline=...,
phase_resolver=phaser)` on the somatic variants you want to classify.

## RNA evidence

For direct read/fragment co-occurrence in an RNA-seq BAM:

```python
from varcode import MolecularPhaseResolver, RNAReadPhasingSource

source = RNAReadPhasingSource("tumor.rna.bam")
phaser = MolecularPhaseResolver(source)
effects = variants.effects(phase_resolver=phaser)
```

Here `variants` is the loaded collection to annotate. Short- and long-read
data can both leave phase unknown; coverage and linked alleles determine whether
a pair is resolved. Trans needs fragments that carry one variant's alt allele
and the other's reference allele. A source that only reports co-observed
partners establishes cis, and leaves every other pair unknown; sources that see
reference alleles too, such as `RNAReadPhasingSource`, report trans through
their own `in_cis`. Raw BAM phasing does not provide an assembled
`MutantTranscript`. Assembly-backed sources can provide one; see the
[source protocols](api_phasing.md#phasing).
`ReadPhaseResolver` remains a compatibility name for `MolecularPhaseResolver`.

<a id="rna-phase-from-isovar"></a>

### RNA phase from Isovar

Isovar's results can serve as the phasing source. Give `run_isovar` the matched
germline variants so they are recognized in the assembled RNA:

```python
from isovar import IsovarReadPhasing, run_isovar
from varcode import MolecularPhaseResolver

results = run_isovar(
    variants=somatic_variants,
    alignment_file="tumor.rna.bam",
    germline_variants=germline_ctx.variants,
)
phaser = MolecularPhaseResolver(IsovarReadPhasing(results))
effects = somatic_variants.effects(germline=germline_ctx, phase_resolver=phaser)
```

With Isovar 1.36 or later, `IsovarReadPhasing.in_cis` answers from fragments
that cover both variants: cis when they carry both alt alleles, trans when they
carry one alt allele with the other's reference allele. A call needs at least
`min_shared_fragments_for_phasing` fragments and a majority.

For a matched germline variant, Isovar 1.38 or later reads the germline site
in fragments that carry the somatic alt allele: the germline alt allele there
means cis, its reference allele trans. Fragments with the somatic reference
allele are not counted, since the germline alt allele also comes from normal
cells and homozygous sites. When those fragments do not decide, a germline edit
in the somatic variant's assembled RNA is cis; when they say trans but the
assembly has the edit, the answer is unknown. Isovar 1.36 and 1.37 only give
the assembly-based cis.

Any pair Isovar cannot call stays unknown, and Varcode keeps every phase
hypothesis for it.

## Unknown phase

When a somatic and germline variant share a codon and relative phase is unknown,
the result can be a `PhaseCandidateSet`, with one classified effect per
haplotype hypothesis. Read `.candidates` to retain the alternatives. Resolved
phase can reduce the set to one effect, as in the example below.

## Two variants in one codon

CFTR has a somatic at GRCh38 7:117531100 `T→A`. A neighbouring
germline at 7:117531101 `T→C` lands in the same codon.

```python
from pyensembl import cached_release
from varcode import Variant, VariantCollection, GermlineContext

g = cached_release(81)
cftr = g.transcript_by_id("ENST00000003084")

somatic = Variant("7", 117_531_100, "T", "A", genome=g)
germline = Variant("7", 117_531_101, "T", "C", genome=g)
ctx = GermlineContext.from_variants([germline], reference_name="GRCh38")
```

### Without germline context

```python
eff = somatic.effect_on_transcript(cftr)
print(type(eff).__name__, eff.short_description)
# Substitution p.L159M
```

### With germline context, phase unknown

```python
eff = somatic.effect_on_transcript(cftr, germline=ctx)
print(type(eff).__name__, eff.short_description)
# PhaseCandidateSet ?p.L159M

for c in eff.candidates:
    ev = c.evidence
    print(f"  haplotype={ev['haplotype']:<2} "
          f"germline_in_cis={[v.short_description for v in ev['germline_variants']]} "
          f"=> {c.effect.short_description}")
# haplotype=B  germline_in_cis=[]                            => p.L159M
# haplotype=A  germline_in_cis=['chr7 g.117531101T>C']       => p.S159T
```

The `?` prefix on `?p.L159M` flags the description as the
most-likely candidate of a `PhaseCandidateSet`. In application code,
check `isinstance(eff, PhaseCandidateSet)` and iterate
`eff.candidates` (a tuple of `EffectCandidate` objects carrying
per-hypothesis evidence keys); use `eff.effects` if only the
inner classified `MutationEffect`s are needed.

### With known phase

A real pipeline gets the cis/trans answer from a phased VCF or an
RNA assembly. For this demo, two stub resolvers force each answer so you can
see the candidate set collapse to a single effect:

```python
class ForceCis:
    source = "demo"
    def in_cis(self, v1, v2, transcript=None): return True

class ForceTrans:
    source = "demo"
    def in_cis(self, v1, v2, transcript=None): return False

eff_cis = somatic.effect_on_transcript(
    cftr, germline=ctx, phase_resolver=ForceCis())
print("cis   ->", type(eff_cis).__name__, eff_cis.short_description)
# cis   -> Substitution p.S159T

eff_trans = somatic.effect_on_transcript(
    cftr, germline=ctx, phase_resolver=ForceTrans())
print("trans ->", type(eff_trans).__name__, eff_trans.short_description)
# trans -> Substitution p.L159M
```

For real evidence, use the [phased VCF](#phased-vcf) or [RNA](#rna-evidence)
resolvers above. The stubs implement only `in_cis(...)`, which is all the
codon-collapse path consults. Richer resolvers implement more of the
duck-typed interface, including `mutant_transcript` and `phased_partners`; the
[phasing API](api_phasing.md#phasing) documents the built-in resolvers and
source protocols. No public `PhaseResolver` base class is required.

## Combining evidence

Given a somatic `VariantCollection`, patient context, phase resolver, and RNA
resolver:

```python
effects = somatic_variants.effects(
    germline=germline_ctx,
    phase_resolver=phaser,
    rna_resolver=rna,
)
```

Order: germline modifies the transcript first, phase collapses the
candidate set next, and RNA refines multi-candidate effects. Splice
mechanism sets are reconciled to observed mechanisms; other
multi-outcome effects append observed-only candidates. Cross-axis key is
`EffectCandidate.evidence["haplotype"]`, so an
RNA observation tagged with the same haplotype tag aligns with the
right germline-aware outcome.

## Known deletion haplotypes in RNA alignments

Splice-aware aligners can encode a known deletion's RNA sequence with `N`
rather than `D`, or split the gap among mismatches and smaller deletions.
Opt in to a local **sequence hypothesis** instead of interpreting the gap as
proof of a DNA deletion:

```python
from varcode import MolecularPhaseResolver, RNAReadPhasingSource

source = RNAReadPhasingSource("tumor.rna.bam", min_alt_reads=1)
source.register_haplotype([known_deletion, adjacent_snv_1, adjacent_snv_2])
resolver = MolecularPhaseResolver(source)
effects = variants.effects(phase_resolver=resolver)
```

The variants must share one genome dataset and contig. The genome must provide
reference sequence across the interval, through a FASTA or annotated transcript.
Every retained base between five-base reference flanks must match one alignment
and pass the configured quality/read-edge filters. No mates or separate partial
haplotypes are stitched together. Registration invalidates cached counts.

Registration does not assert cis phase: it supplies a candidate to test. For
registered variants, only full-context matches count as alternate support;
nonmatches remain unknown, not trans evidence. The current local mode supports
nonoverlapping substitutions and deletions. Register competing combinations
separately. Unregistered variants keep the ordinary CIGAR-based behavior.

For the audited GRCh38 MAP2 example, the combination of `209694769 C>A`,
`209694770 T>G`, and the 28-base deletion at `209694773` gives the same local
sequence under `28D`, `28N`, and `22D2M6D` alignments. An unrelated 19,449-base
exon skip lacks the required local anchors and is not supporting evidence.

## Implementation limits

See [germline limitations](germline.md#limitations) for the phase-hypothesis
cap and normalization requirements, and the [phasing API](api_phasing.md) for
custom sources.
