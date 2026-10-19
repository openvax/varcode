# Python 3.11 support floor (#607)

## Specification

Require Python 3.11 or newer without a future upper bound. Keep the existing
pandas 2.2/3 compatibility policy and test both on Python 3.11; add Python
3.13 and 3.14 pandas 3 jobs while retaining the Python 3.12 job. Update
current support claims and release notes, then bump the patch version to
10.17.2. Preserve variant annotation behavior and dependency bounds.

## Plan

- [x] Read issue, package metadata, CI, release guidance, and current docs.
- [x] Check scope before implementation: change interpreter support policy
  without altering pandas or genome behavior.
- [x] Update metadata, CI, support docs, changelog, and version.
- [x] Install the source with optional test dependencies on Python 3.11 and
  3.14; verify dependency resolution.
- [x] Run `./lint.sh` and `./test.sh` with writable Ensembl cache copies and
  review source-distribution contents.
- [ ] Open PR, pass CI, merge, deploy from clean main, verify PyPI and tag.

## Review

`./lint.sh` passed. The full suite passed 3,446 tests with one skip on Python
3.11.12 and Python 3.14.6. The first 3.11 run found that this isolated venv
lacked the repo's documented `setuptools>=64` test prerequisite; after adding
it, the source-distribution test and full suite passed. Both environments
resolve the optional Osteosarc and ORF dependencies, and dependency checks pass.
After upgrading to released PyEnsembl 2.28.1, the documented CFTR example
returns `p.L159M` on both interpreters. Installed metadata reports varcode
10.17.2 and `Requires-Python: >=3.11`. The source-distribution integrity test
passed in the full suite.
