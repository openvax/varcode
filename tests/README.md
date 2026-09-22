# Osteosarc test variants

`tests/data/osteosarc_variants.json` contains 177 ready site variants and five
unresolved entries collected through `osteosarc==0.1.0` from the pinned public
snapshot below. The ready variants comprise 158 SNVs, 13 deletions, three
insertions, and three complex alleles. Original alleles, assemblies, source
IDs, correction IDs, source receipts, and hashes of complete osteosarc entries
are retained. The fixture is approximately 115 KiB; large count and peptide
tables remain in the source snapshot.

Use the collection in tests without installing osteosarc or opening its cache:

```python
from tests.osteosarc_variants import load_variants, read_fixture

variants = load_variants()  # native VariantCollection, GRCh38 / Ensembl 81
for variant in variants:
    source_entries = variants.metadata[variant]["entries"]
unresolved = [entry for entry in read_fixture()["entries"]
              if entry["status"] != "ready"]
```

The ordinary suite annotates every ready allele with both `fast` and
`protein_diff` and checks for errors or unresolved annotator results. These
are real-input coverage checks; they do not certify source somatic status or
use source protein labels as expected consequences. Dedicated regressions
still provide independent expected effects.

The mitochondrial allele keeps its original `chrM` spelling; Varcode's
explicit UCSC conversion maps it to `MT`. Osteosarc 0.1.4's native
`to_varcode` adapter uses the same conversion and preserves the source name.
Corrected MAP2 and its separately published split representation retain distinct source IDs
and correction notes; they are not assumed to be independent events.

## Regenerate from the verified snapshot

The historical fixture records Osteosarc 0.1.0. To reproduce it byte for byte,
install that version in a separate environment on Python 3.10+, then collect offline:

```sh
python -m pip install -e . 'osteosarc==0.1.0'
python -m tests.collect_osteosarc_variants \
  --cache /path/to/shared/cache --snapshot 2026-09-18t
```

The exporter checks the package version and snapshot identity before writing
the fixture. `--output` selects a different destination. To collect a new
snapshot deliberately, supply its identity with `--expected-snapshot-id`,
review the changed data, and update the test pins. Source acquisition remains
separate from export and tests.

## Optional snapshot integration checks

The targeted GPX4 and BRCA1 regressions use Ensembl 81. Additional offline
integration checks use the public osteosarc dataset through the published
`osteosarc==0.1.4` adapter from the optional `test-data` extra (Python 3.10+):

```sh
python -m pip install -e '.[test-data]'
OSTEOSARC_TEST_CACHE=/path/to/shared/cache \
OSTEOSARC_TEST_SNAPSHOT=2026-09-18t \
pytest -q tests/test_osteosarc_dataset.py
```

The cache must already contain snapshot
`9b34ea0e13f9c1c35c3c88b7e646c0e608b86a143dee0e668bf3f74b909f815c`.
The snapshot name is local; its content identity is checked independently.
Acquisition is explicit and separate from testing: use osteosarc's
`Dataset.sync` for new snapshots, or reuse the verified shared cache for this
snapshot. Creating a new snapshot from current remote sources does not
guarantee this pinned identity. Tests never download data or refresh sources.
Without `OSTEOSARC_TEST_SNAPSHOT`, these optional checks are skipped; when
it is set, a missing package, cache object, or mismatched snapshot fails.

These additional checks preserve the historical fixture's 177 ready alleles
and verify native conversion of all 179 entries now ready in the same snapshot
(182 total). The two newly resolved entries are not silently added to the
checked-in scientific corpus. Checks preserve source provenance and reference
identity, exercise mitochondrial annotation and the corrected MAP2 allele,
reopen the shared cache offline, and reject missing or tampered objects using
a private copy. Package updates do not change the pinned snapshot identity.
They do not treat upstream peptide annotations or vaccine membership as a
protein oracle, or change the existing RNA-read fixtures. Broader RNA fixture
adoption is tracked in [#464](https://github.com/openvax/varcode/issues/464).

Dataset: [osteosarc.com](https://osteosarc.com/data/), snapshot acquired
2026-09-18, collected 2026-09-21; public data listed as CC0-1.0 by the
[AWS Open Data Registry](https://registry.opendata.aws/sid-osteosarc/).
