# Optional osteosarc corpus checks

The targeted GPX4 and BRCA1 regressions use Ensembl 81. Additional offline
integration checks use the public osteosarc dataset through the published
`osteosarc==0.1.0` adapter (Python 3.10+):

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

The checks preserve all 177 ready site entries and their source provenance
through native Varcode conversion, using an explicit GRCh38 Ensembl release.
They also distinguish the corrected MAP2 complex allele from the old deletion.
They do not treat upstream peptide annotations or vaccine membership as a
protein oracle, or change the existing RNA-read fixtures. Broader RNA fixture
adoption is tracked in [#464](https://github.com/openvax/varcode/issues/464).

Dataset: [osteosarc.com](https://osteosarc.com/data/), snapshot acquired
2026-09-18, checked 2026-09-19; public data listed as CC0-1.0 by the
[AWS Open Data Registry](https://registry.opendata.aws/sid-osteosarc/).
