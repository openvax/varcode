# Contributing to Varcode

Varcode supports Python 3.9 and later. Read [AGENTS.md](AGENTS.md) for the
repository workflow and [RELEASING.md](RELEASING.md) before preparing a release.
The [2026-09-24 audit](docs/development_audit.md) records known quality gaps and
their issue links.

## Report a problem

Search the [open issues](https://github.com/openvax/varcode/issues) first.
Include a minimal reproduction, the Varcode, PyEnsembl and Python versions,
the annotation release, and the variant's coordinates, REF/ALT and assembly.
For file parsing problems, include a small synthetic input with the relevant
header and sample fields. Report substantial new work as an issue before
implementation, and link related issues from the PR.

## Set up a development environment

From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[rna,test-data]' pytest pytest-cov pytest-xdist ruff
python -m pip install mkdocs mkdocs-material 'mkdocstrings[python]'
```

`rna` installs the optional BAM reader; `test-data` installs the pinned
Osteosarc adapter. `pytest-xdist` is optional: `test.sh` also runs serially.
Activate this environment before invoking the scripts so Python and its tools
come from the same installation. `./develop.sh` installs just the editable
library into the active environment.

Install the reference datasets used by the full suite (these are downloads):

```bash
pyensembl install --release 75 --species human --custom-mirror https://github.com/openvax/ensembl-data/releases/download/GRCh37.75/
pyensembl install --release 81 --species human --custom-mirror https://github.com/openvax/ensembl-data/releases/download/GRCh38.81/
pyensembl install --release 95 --species human --custom-mirror https://github.com/openvax/ensembl-data/releases/download/GRCh38.95/
pyensembl install --release 95 --species mouse --custom-mirror https://github.com/openvax/ensembl-data/releases/download/GRCm38.95/
```

The test reference resolver currently probes uninstalled releases in a way
that can trigger downloads; see [#493](https://github.com/openvax/varcode/issues/493).
The bundled Osteosarc snapshot checks themselves open their data offline.
See [tests/README.md](tests/README.md) for fixture provenance and regeneration.

## Make and verify a change

Create a feature branch before editing; land changes through a reviewed PR.
Use [PEP 8](https://peps.python.org/pep-0008/) and NumPy-style docstrings.
Include a focused regression for a bug fix and keep the fixture construction
or reproducible source acquisition alongside the tests. Scientific expected
values should have an independent basis, rather than comparing a function
with itself. Preserve public APIs where practical: Isovar, Topiary and Vaxrank
consume Varcode results.

```bash
./lint.sh
./test.sh
mkdocs build --strict
```

For a focused test run, use `python -m pytest tests/test_vcf_limits.py`.
`test.sh` accepts pytest arguments and chooses a memory-aware worker count;
set `TEST_SH_MAX=2` to cap it when running other suites concurrently.
The CI parity check is `./test.sh --annotator=protein_diff -m parity`.

Every PR includes a version bump in `varcode/version.py`, relevant changelog
notes, and issue links. Wait for green GitHub checks before merging; follow
[the release procedure](RELEASING.md) to publish and verify the package.

## Licensing

Varcode and contributions are licensed under the Apache License 2.0.

## Writing transforms

Every transform owes three things, documented in its docstring:

| Field | Meaning |
|---|---|
| **Cardinality** | `preserves`, `reduces`, or `increases`. |
| **Provenance** | Every output variant carries `source_variants: tuple[Variant, ...]`. Empty tuple for pass-through; one element for derived-from-one; two or more for combined. Not part of hash/equality. |
| **Metadata behavior** | Explicit rule for how `source_to_metadata_dict` entries flow through (which fields are inherited from which source, which require agreement, what happens on disagreement). |

Transforms are **idempotent on inputs they don't recognize**. Running
`pair_breakends` twice produces the same VC; the second pass finds no
unpaired BNDs to combine because every combined row's `source_variants`
is already populated.


The proposed `combine_cis_snvs` transform is tracked in
[#368](https://github.com/openvax/varcode/issues/368); it is not a shipped API.
User-facing examples belong in the [transforms guide](https://openvax.github.io/varcode/transforms/).
