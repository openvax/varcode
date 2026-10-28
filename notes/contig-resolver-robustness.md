# External reference vocabulary robustness (#650)

Keep one public spelling normalization path for canonical chromosome case
aliases; preserve the explicit normalization/conversion flags on Variant.
External reference lookup additionally accepts the valid name `0`, independent
of annotation's zero-chromosome sentinel. Exact names win; normalized collisions
remain ambiguous and iteration-order independent. No viral/accession synonym
catalogue, accession-version guessing, alternate-to-primary projection or
coordinate equivalence is introduced.

Validate canonical aliases, conflicts, mixed custom/viral/decoy/alternative
vocabularies in an indexed BAM, unchanged annotation/Variant validation and
existing UCSC hg19/hg38 cases. Run lint/mypy, full tests and CI; publish from
clean main. Then Isovar #608 can adopt the released resolver without duplicating
these spelling rules. Reference: SAM v1.6 sections1.2.1/1.3,
https://samtools.github.io/hts-specs/SAMv1.pdf.
