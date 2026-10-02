# Releasing Varcode

A release is complete after the reviewed PR is merged, CI passes, and the
new version is available on PyPI. Use a feature branch for all code and
version changes; never commit directly to `main`.

## Prepare the PR

1. Choose a [semantic version](https://semver.org/), update `__version__` in
   `varcode/version.py`, and add release notes to `CHANGELOG.md`. Every PR
   needs at least a patch bump. Recheck the version against `origin/main`
   before merging if another release has shipped.
2. Run `./lint.sh`, `./test.sh`, and `mkdocs build --strict` for documentation
   changes. Review the diff and wait for green GitHub CI on the final commit.
3. Merge the PR through GitHub.

## Publish from clean main

Activate the development environment described in [CONTRIBUTING.md](CONTRIBUTING.md),
with build/upload credentials configured for PyPI. Use a clean checkout:

```bash
git switch main
git pull --ff-only origin main
git status --porcelain
python -c 'import runpy; print(runpy.run_path("varcode/version.py")["__version__"])'
```

The script requires a clean `main` or `master` whose HEAD matches the live
branch on `origin`. The release tag must be absent both locally and remotely.
Then run:

```bash
./deploy.sh
```

Use `./deploy.sh X.Y.Z` to also assert the already-merged version.

`deploy.sh` runs lint and tests, builds a wheel and source distribution in a
new directory under `dist/`, and runs `twine check --strict` on that exact pair.
It rechecks the checkout and remote refs before uploading, then tags the
original merged commit with `v<version>` and pushes only that tag. Version
changes belong in the PR; a different version argument fails without editing
or committing files. Activate the development environment first; the script
uses `python` and the installed `build` and `twine` tools.

`./deploy.sh --dry-run [version]` performs the same guards, lint, tests, build,
and distribution checks, then stops before uploading or changing Git history.
It needs read access to `origin` and leaves its build artifacts under `dist/`.

After a successful upload, verify the exact version on
[PyPI](https://pypi.org/project/varcode/), compare its wheel/sdist SHA256 hashes
with the files in the printed artifact directory, and smoke-test the published
package in an isolated environment. Verify that the pushed tag targets the
merged commit. Check the main-branch CI and documentation deployment, then
review open issues for the next dependency or correctness blocker.

If publication fails, keep the printed artifact directory and inspect PyPI
before retrying. Upload only missing files from that original pair after
verifying the published hashes; rebuilding may produce different bytes.
Finish tagging the original commit and pushing the tag once both files are
verified. If only the tag push failed, verify the local tag's target and retry
`git push origin refs/tags/vX.Y.Z`. Never move an existing release tag.
