# Bundled fixture inventory

`SHA256SUMS` records every file in this directory tree, including the provenance
notes, except the checksum file itself. The fixtures total about 5 MB; the
largest is the 3.6 MiB compressed historical Osteosarc metadata snapshot.
No BAM/CRAM files or shared RNA bundles belong here. See
[the test-data guide](../README.md) for reference installation and external
`openvax-v2` assets.

## Sources and recipes

| Files | Recorded source and reproduction information |
|---|---|
| `osteosarc_variants.json`, `osteosarc_snapshot_2026-09-18t.zip` | [Snapshot identity, receipts, archive hash and offline export recipe](../README.md#regenerate-from-the-verified-snapshot). `../collect_osteosarc_variants.py` exports the reviewed snapshot with the pinned Osteosarc version. |
| `osteosarc_observed_junctions.json` | Each record retains its source page, ONT read name, reference identity and source-sequence SHA256. [Shared read verification](../README.md#shared-test-data-openvax-v2) uses the external bundle; only the bounded junction windows are stored here. |
| `osteosarc_esvee_somatic.vcf`, `osteosarc_esvee_shared_svs.vcf` | Published Sid esvee records with FORMAT/sample columns removed. The subset selection and reference releases are documented in [the fusion tests](../test_osteosarc_fusions.py) and [shared SV tests](../test_osteosarc_shared_svs.py). |
| `real_callers/*.vcf` | [Pinned GATK source and synthetic Strelka2/VEP construction notes](real_callers/README.md). |
| `spec_examples/*.vcf` | [VCF specification sections and source](spec_examples/README.md). |
| `documentation/*.vcf` | Small inputs constructed for [executable documentation tests](../test_documentation_examples.py), introduced in [14789b1](https://github.com/openvax/varcode/commit/14789b1) and [8108acd](https://github.com/openvax/varcode/commit/8108acd). |
| `dbnsfp_validation_set.csv`, `somatic_hg19_14muts*.vcf*` | Legacy validation inputs introduced in [5d0d268](https://github.com/openvax/varcode/commit/5d0d268); the compressed VCF was added in [85f232c](https://github.com/openvax/varcode/commit/85f232c). The dbNSFP comparisons use Ensembl 75. |
| `tcga_ov.head*.maf`, `ov.wustle.subset5.maf` | Legacy MAF subsets introduced in [679167a](https://github.com/openvax/varcode/commit/679167a) and [30c9082](https://github.com/openvax/varcode/commit/30c9082); record headers retain their source centers and assemblies. |
| `mouse_vcf_dbsnp_chr1_partial.vcf` | Legacy mouse dbSNP subset introduced in [863167b](https://github.com/openvax/varcode/commit/863167b). |
| `mutect-example*.vcf`, `strelka-example.vcf` | Legacy caller outputs introduced in [302ef49](https://github.com/openvax/varcode/commit/302ef49), retaining available caller/reference header metadata. |
| `simple.*.vcf`, `duplicate-id.*.vcf`, `different-samples.*.vcf`, `same-samples.*.vcf`, `multiallelic.vcf`, `duplicates.vcf`, `duplicates.maf` | Parser/collection regression inputs introduced in [f82edd5](https://github.com/openvax/varcode/commit/f82edd5), [6b53cc7](https://github.com/openvax/varcode/commit/6b53cc7), and [22dcb0f](https://github.com/openvax/varcode/commit/22dcb0f). |

The original external acquisition recipes for the legacy inputs were not
recorded. Their repository history is the available provenance; no new upstream
identity or scientific validation is inferred here. All fixture payloads in
this release are byte-for-byte those in
[commit 6f29cf8](https://github.com/openvax/varcode/tree/6f29cf8c5201189ae2803ad319aacca8053f7ddf/tests/data).
To recover those exact payloads, use the raw files under that commit's
`tests/data/` tree and verify them against `SHA256SUMS`. The new inventory notes
and checksums are part of the source distribution itself.

## Verify and update

From the repository or unpacked source-archive root:

```sh
python -m pytest -q tests/test_source_distribution.py
```

After deliberately changing a fixture or its provenance notes, review the
source, size and test expectations, then regenerate the sorted inventory:

```python
from hashlib import sha256
from pathlib import Path

root = Path("tests/data")
files = sorted(p for p in root.rglob("*")
               if p.is_file() and p.name != "SHA256SUMS")
(root / "SHA256SUMS").write_text("".join(
    f"{sha256(p.read_bytes()).hexdigest()}  {p.relative_to(root).as_posix()}\n"
    for p in files))
```

Hashes identify reviewed bytes; they do not establish the scientific correctness
of an expected annotation. Keep independent biological expectations in the
corresponding tests and retain source receipts when adding captured data.
