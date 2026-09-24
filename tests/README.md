# Test data and reproducibility

Keep small synthetic fixture constructors in the test module or a shared
helper under `tests/`. Document the expected behavior and why the inputs
exercise it. For captured public data, check in the acquisition/subsetting
recipe, source identifiers, reference release, and integrity hashes with the
fixture. Do not make an untracked local cache or sibling checkout a test
requirement. RNA fixture harmonization is tracked in
[#464](https://github.com/openvax/varcode/issues/464).

See [CONTRIBUTING.md](../CONTRIBUTING.md) for environment and reference setup.
The reference-cache probe has a separate download issue
([#493](https://github.com/openvax/varcode/issues/493)); the offline guarantees
below apply to the bundled Osteosarc snapshot checks.

## Osteosarc test variants

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
install that version in a separate environment on Python 3.10+, unpack the
bundled snapshot into a new cache directory, then collect offline:

```sh
python -m pip install -e . 'osteosarc==0.1.0'
python -m zipfile -e tests/data/osteosarc_snapshot_2026-09-18t.zip /path/to/new/cache
python -m tests.collect_osteosarc_variants \
  --cache /path/to/new/cache --snapshot 2026-09-18t
```

The exporter checks the package version and snapshot identity before writing
the fixture. `--output` selects a different destination. To collect a new
snapshot deliberately, supply its identity with `--expected-snapshot-id`,
review the changed data, and update the test pins. Source acquisition remains
separate from export and tests.

## Snapshot integration checks

The targeted GPX4 and BRCA1 regressions use Ensembl 81. Offline integration
checks use the public osteosarc dataset through the published
`osteosarc==0.2.3` adapter from the optional `test-data` extra (Python 3.9+):

```sh
python -m pip install -e '.[test-data]'
python -m pytest -q tests/test_osteosarc_dataset.py
```

These five tests run in the ordinary suite and in CI on Python 3.9/3.10/3.11.
They unpack `tests/data/osteosarc_snapshot_2026-09-18t.zip` into a pytest
temporary directory, check the archive SHA256, and open it offline through
Osteosarc, which verifies all source objects against their receipts. No
sibling checkout, pre-existing Osteosarc cache, or environment variables are
required. As elsewhere in the suite, Ensembl 81 reference data must already
be installed; CI provisions this in its existing reference-data step.

The archive contains the unchanged snapshot manifest and its 20 public
metadata objects in the standard OpenVax cache layout. It is approximately
3.6 MiB compressed (53 MiB unpacked), with SHA256
`cbb688ca5cbe775fc4e6a826124f861d6d605c65e34473582c5d86853a6018fa`.
The manifest retains original URLs, acquisition timestamps, sizes and hashes.
It contains no BAMs, generated reads, local cache bookkeeping, or absolute
paths. Its snapshot identity remains
`9b34ea0e13f9c1c35c3c88b7e646c0e608b86a143dee0e668bf3f74b909f815c`.

To test an existing shared cache explicitly, set both
`OSTEOSARC_TEST_CACHE=/path/to/shared/cache` and
`OSTEOSARC_TEST_SNAPSHOT=2026-09-18t`. The name may differ, but the content
identity must match. Missing packages, objects, or a mismatched snapshot then
fail; the tests never fall back to the bundled snapshot. Without an explicit
snapshot, these checks skip only when the optional Osteosarc dependency is
absent. Osteosarc 0.2.3 supports every Python version in the CI matrix.

These snapshot integration tests never download data or refresh sources. `Dataset.sync` acquires new
snapshots separately; current remote sources do not reproduce the historical
identity. Change this fixture only by deliberately reviewing a new snapshot
and updating its source provenance, identity, archive hash, and expectations.

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
