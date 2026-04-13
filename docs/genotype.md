# Genotypes and sample-aware queries

*New in varcode 2.3.0 ([#267](https://github.com/openvax/varcode/issues/267)).*

When you load a multi-sample VCF, varcode captures each sample's genotype
(GT, AD, DP, GQ, PS) automatically. As of 2.3.0 this data is available
as structured `Genotype` objects and via sample-aware filtering helpers
on `VariantCollection`.

## The basics

```python
from varcode import load_vcf

vc = load_vcf("tumor_normal.vcf", genome="GRCh38")

vc.samples
# ['normal', 'tumor']

vc.has_sample_data()
# True
```

### One variant, one sample

```python
variant = vc[0]
gt = vc.genotype(variant, "tumor")
# Genotype(raw_gt='0/1', alleles=(0, 1), phased=False, phase_set=None,
#          allele_depths=(10, 5), total_depth=15, genotype_quality=99)

gt.is_called              # True
gt.alleles                # (0, 1)
gt.phased                 # False
gt.allele_depths          # (10, 5)  —  (ref_depth, alt_depth)
gt.total_depth            # 15
gt.depth_for_alt(1)       # 5
```

### Zygosity is relative to *this variant's* alt

```python
from varcode import Zygosity

vc.zygosity(variant, "tumor")
# <Zygosity.HETEROZYGOUS: 'het'>

vc.zygosity(variant, "normal")
# <Zygosity.ABSENT: 'absent'>  — normal is ref/ref at this site
```

Four states, relative to the alt of the `Variant` you query:

| State | Meaning |
|---|---|
| `HETEROZYGOUS` | Sample has at least one copy of this alt, and not all copies are this alt |
| `HOMOZYGOUS` | All called copies are this alt |
| `ABSENT` | Sample was called, but doesn't carry this alt (ref/ref, or a different alt at a multi-allelic site) |
| `MISSING` | Sample's call was `./.` or the sample wasn't called |

## Filtering

These return filtered `VariantCollection`s:

```python
# Variants where this sample carries the alt (het or hom).
vc.for_sample("tumor")

# Finer distinctions:
vc.heterozygous_in("tumor")
vc.homozygous_alt_in("tumor")
```

Typoed sample names fail fast:

```python
vc.for_sample("typo")
# SampleNotFoundError: Sample 'typo' not found. Available samples: ['normal', 'tumor']
```

## Cross-sample queries

Tumor/normal, trio, and cohort queries fall out of set operations on
the primitives:

```python
# Somatic candidates: in tumor, not in normal.
tumor = set(vc.for_sample("tumor"))
normal = set(vc.for_sample("normal"))
somatic = tumor - normal

# De novo candidates in a trio.
de_novo = (set(vc.for_sample("child"))
           - set(vc.for_sample("mom"))
           - set(vc.for_sample("dad")))
```

## Multi-allelic sites

VCF rows like `REF=A ALT=T,G` are split into two `Variant` objects
(one per alt). `zygosity` is computed relative to the specific alt of
the variant you query, so a sample with `GT=1/2` is reported as
heterozygous for both A→T and A→G (one copy of each, with the other
copy being a *different* alt). A sample with `GT=0/1` (one ref, one T)
is heterozygous for A→T and `ABSENT` for A→G.

This means `vc.heterozygous_in("tumor")` on a row where tumor is `1/2`
yields two entries, one per alt — the correct biological
interpretation.

## What this doesn't do yet

Three capabilities are deliberately out of scope for 2.3.0:

1. **Effect prediction that uses zygosity.** `variant.effects()` still
   returns the same output regardless of zygosity. Compound-het
   reasoning and germline-aware effect prediction depend on the
   `Genotype` API landing first and will arrive in the personal-genome
   work tracked under [#270](https://github.com/openvax/varcode/issues/270).

2. **Phased effects.** Phased GTs (`0|1`) and `PS` phase-set tags are
   preserved on `Genotype` but not yet used to predict joint effects
   of nearby variants on the same haplotype. Tracked in
   [#269](https://github.com/openvax/varcode/issues/269).

3. **Dedicated joint-analysis helpers** like `vc.de_novo_in(...)` or
   `vc.somatic(...)`. The set-operation pattern above covers them;
   built-in shortcuts will be added if a clear use case emerges.
