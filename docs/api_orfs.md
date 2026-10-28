# ORFs API

For installation, reference binding, genetic codes, and evidence limitations,
see [ORF annotation](uorf_annotations.md). Ordinary transcript annotation uses
`variant.effects()`; add `extra_orfs=` to that workflow when you want an
EffectCollection containing additional ORF consequences.

## Choose an entry point

| Task | Entry point | Result |
|---|---|---|
| Load reference records | `load_orf_catalog(...)` | ORFdata catalogue |
| Annotate a validated reference ORF | `annotate_orf_variants(variants, reference)` | ORFAnnotation |
| Annotate overlapping catalogue records | `predict_orf_effects(variants, ...)` | Tuple of ORFAnnotation objects |
| Compare a mapped observed RNA product | `annotate_orf_rna_product(reference, product)` | ORFProductAnnotation |
| Export the observed-product comparison directly | `annotate_orf_product(reference, product)` | Dictionary for JSON export |

The tuple from `predict_orf_effects` is not an EffectCollection. Empty results
mean no catalogue model was selected; they do not establish that no ORF can be
affected. Keep each annotation's unresolved reason and alternative outcomes.
The supplied reference and observed product must retain their coordinate and
reference context; a protein string alone cannot supply that mapping.

## Annotation functions

::: varcode.load_orf_catalog

::: varcode.annotate_orf_variants

::: varcode.predict_orf_effects

::: varcode.annotate_orf_rna_product

::: varcode.annotate_orf_product

## Result models

::: varcode.ORFAnnotation

::: varcode.ORFOutcome

::: varcode.ORFProduct

::: varcode.ORFProductAnnotation

## Upstream-ORF compatibility API

These names preserve the earlier record-plus-transcript interface. Prefer
the bound-reference functions above for new integrations. The upstream-only
catalogue default is retained; the wrappers are not all interchangeable aliases
for the general catalogue API.

`UORFAnnotation` is an alias of ORFAnnotation.

::: varcode.UORFAnnotation
    options:
      members: false
      skip_local_inventory: true

::: varcode.ORFTranscript

::: varcode.load_uorf_catalog

::: varcode.annotate_uorf_variants

::: varcode.annotate_uorf_product

::: varcode.predict_uorf_effects
