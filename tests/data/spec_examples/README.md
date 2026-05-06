# Public VCF spec example fixtures

These VCF files are reproduced from the official VCF specifications published
by the samtools/hts-specs project (MIT licensed):

  https://github.com/samtools/hts-specs

They serve as canonical, widely-known examples that exercise a broad range of
header / INFO / FORMAT shapes. The accompanying tests in
`tests/test_vcf_spec_examples.py` pin our parser's output against
hand-verified expected values, so a regression here surfaces as an obvious
diff against the published spec.

| File | Source | What it exercises |
|---|---|---|
| `vcf42_spec_trio.vcf` | VCFv4.2 §1.1 | Three samples, multiallelic ALTs, microsatellite indel, `Number=A`, `Number=2`, phased and unphased GTs, missing fields, quoted contig metadata, `ALT=.` row |
| `vcf43_spec_sv.vcf` | VCFv4.3 §5.4 | Symbolic ALT alleles (`<DEL>`, `<DUP>`, `<INS>`, `<INV>`), `END=`, confidence intervals, `IMPRECISE` flag, structural-variant INFO fields |

The oracle tests in `tests/test_vcf_parsing.py` automatically pick up any
`*.vcf` placed under `tests/data/**` (it walks subdirectories), so adding more
spec fixtures here extends both the oracle comparison and the dedicated
spec-example assertions.
