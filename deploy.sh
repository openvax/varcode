#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: ./deploy.sh [--dry-run] [version]"
    echo "Publish the merged version; an optional version must match varcode/version.py."
}

fail() { echo "deploy.sh: $*" >&2; exit 1; }

dry_run=0
expected_version=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) dry_run=1 ;;
        -h|--help) usage; exit 0 ;;
        -*) fail "Unknown option: $1" ;;
        *)
            [[ -z "$expected_version" ]] || fail "Expected at most one version."
            expected_version="${1#v}"
            [[ -n "$expected_version" ]] || fail "Expected a release version."
            ;;
    esac
    shift
done

cd "$(dirname "$0")"
branch="$(git symbolic-ref --quiet --short HEAD)" || fail "Deploy from main or master, not a detached HEAD."
[[ "$branch" == main || "$branch" == master ]] || fail "Deploys require main or master (current: $branch)."
commit="$(git rev-parse HEAD)"
version="$(python -c 'import runpy; print(runpy.run_path("varcode/version.py")["__version__"])')"
[[ "$version" =~ ^(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)$ ]] || fail "Invalid release version: $version"
[[ -z "$expected_version" || "$expected_version" == "$version" ]] || fail "Requested $expected_version, but the merged version is $version. Bump versions in a PR."
tag="v$version"

# Check live origin refs, not a possibly stale remote-tracking branch.
check_release() {
    local status remote_branch remote_tag
    status="$(git status --porcelain --untracked-files=all)"
    [[ -z "$status" ]] || fail "Working tree is not clean."
    [[ "$(git symbolic-ref --quiet --short HEAD)" == "$branch" && "$(git rev-parse HEAD)" == "$commit" ]] || fail "Checkout changed during release."
    remote_branch="$(git ls-remote --exit-code origin "refs/heads/$branch")" || fail "Cannot read origin/$branch."
    [[ "${remote_branch%%$'\t'*}" == "$commit" ]] || fail "HEAD must match live origin/$branch; pull the merged release first."
    if git show-ref --verify --quiet "refs/tags/$tag"; then
        fail "Tag $tag already exists locally. Verify the publication state before retrying."
    fi
    remote_tag="$(git ls-remote origin "refs/tags/$tag")" || fail "Cannot check release tags on origin."
    [[ -z "$remote_tag" ]] || fail "Tag $tag already exists on origin."
}

check_release
./lint.sh
./test.sh

# Keep each attempt's original bytes, including after a failed/partial upload.
mkdir -p dist
artifact_dir="$(mktemp -d "dist/$tag.XXXXXX")"
echo "Building $tag at $commit; artifacts: $artifact_dir"
python -m build --outdir "$artifact_dir"
distributions=("$artifact_dir/varcode-$version-py3-none-any.whl" "$artifact_dir/varcode-$version.tar.gz")
python -m twine check --strict "${distributions[@]}"
check_release

if [[ "$dry_run" -eq 1 ]]; then
    echo "Dry run passed for $tag; no upload, commit, tag, or push. Artifacts: $artifact_dir"
    exit 0
fi

python -m twine upload --non-interactive --disable-progress-bar "${distributions[@]}"
git tag "$tag" "$commit"
git push origin "refs/tags/$tag"
echo "Deployed varcode $version: https://pypi.org/project/varcode/$version/"
