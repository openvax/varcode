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
python -c 'from varcode.version import __version__; print(__version__)'
```

Confirm that the working tree is clean, the branch is `main`, and the imported
version is the merged release. Then run:

```bash
./deploy.sh
```

The current script runs lint and tests, installs/upgrades `build` and `twine`,
recreates `dist/`, builds distributions, and uploads them to PyPI. It does
**not** check the branch or working tree, accept a version argument, bump the
version, create a tag, or push commits. Those missing safeguards are tracked
in [#414](https://github.com/openvax/varcode/issues/414); enforce the clean-main
requirement yourself until the script implements it. If a step fails, fix its
root cause and verify the publication state before retrying an upload.

After a successful upload, verify the exact version on
[PyPI](https://pypi.org/project/varcode/), compare its wheel/sdist SHA256 hashes
with the files in `dist/`, and smoke-test the published package in an isolated
environment. Tag the same merged commit with `v<version>` and push that tag:

```bash
# Replace X.Y.Z with the verified published version.
git tag vX.Y.Z
git push origin vX.Y.Z
```

If the tag already exists, verify its target instead of moving it. Check the
main-branch CI and documentation deployment, then review open issues for the
next dependency or correctness blocker.
