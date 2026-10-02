"""Exercise deployment with real local Git repos and offline build/upload stubs."""

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
VERSION = "1.2.3"
TAG = "v" + VERSION


def git(repo, *args):
    return subprocess.run(
        ["git", "-C", str(repo), *args], check=True, capture_output=True,
        text=True).stdout.strip()


@pytest.fixture
def release(tmp_path, monkeypatch):
    repo = tmp_path / "release repo"
    repo.mkdir()
    git(repo, "init", "-b", "main")
    git(repo, "config", "user.name", "Release Test")
    git(repo, "config", "user.email", "release@example.invalid")
    git(repo, "config", "commit.gpgsign", "false")
    git(repo, "config", "tag.gpgsign", "false")
    shutil.copy(ROOT / "deploy.sh", repo / "deploy.sh")
    (repo / "varcode").mkdir()
    (repo / "varcode/version.py").write_text('__version__ = "1.2.3"\n')
    (repo / ".gitignore").write_text("dist/\n__pycache__/\n")
    for step in ("lint", "test"):
        script = repo / (step + ".sh")
        script.write_text(
            '#!/bin/sh\n'
            'echo %s >> "$RELEASE_LOG"\n'
            '[ "${FAIL_STEP:-}" != %s ]\n' % (step, step))
        script.chmod(0o755)
    git(repo, "add", ".")
    git(repo, "commit", "-m", "Release fixture")
    origin = tmp_path / "origin.git"
    git(tmp_path, "init", "--bare", str(origin))
    git(repo, "remote", "add", "origin", str(origin))
    git(repo, "push", "-u", "origin", "main")
    hook = origin / "hooks/pre-receive"
    hook.write_text(
        '#!/bin/sh\ncat >> "$PUSH_LOG"\n[ "${FAIL_STEP:-}" != push ]\n')
    hook.chmod(0o755)

    # Only Python package build/upload commands are faked. Version reading and
    # all Git operations use the real executables; no network service is used.
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    python = bin_dir / "python"
    python.write_text("#!" + sys.executable + "\n" + '''
import json
import os
from pathlib import Path
import runpy
import subprocess
import sys

args = sys.argv[1:]
if args[:1] == ["-c"]:
    os.execv(sys.executable, [sys.executable, *args])
if args[:2] == ["-m", "build"]:
    step = "build"
else:
    assert args[:2] == ["-m", "twine"], args
    step = args[2]
with open(os.environ["RELEASE_LOG"], "a") as log:
    log.write(step + "\\n")
with open(os.environ["RELEASE_ARGS"], "a") as log:
    log.write(json.dumps(args) + "\\n")
if os.environ.get("FAIL_STEP") == step:
    sys.exit(17)
if step == "build":
    output = Path(args[args.index("--outdir") + 1])
    version = runpy.run_path("varcode/version.py")["__version__"]
    for suffix in ("-py3-none-any.whl", ".tar.gz"):
        (output / ("varcode-" + version + suffix)).write_bytes(b"original bytes")
    if os.environ.get("BUILD_HOOK"):
        subprocess.run(["bash", os.environ["BUILD_HOOK"]], check=True)
else:
    paths = [Path(arg) for arg in args[3:] if not arg.startswith("--")]
    assert len(paths) == 2 and all(path.is_file() for path in paths), args
    if step == "upload":
        assert not subprocess.check_output(["git", "tag", "--list", "v1.2.3"])
''')
    python.chmod(0o755)
    (bin_dir / "python3").symlink_to(python)
    monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ["PATH"])
    monkeypatch.setenv("RELEASE_LOG", str(tmp_path / "commands"))
    monkeypatch.setenv("RELEASE_ARGS", str(tmp_path / "arguments"))
    monkeypatch.setenv("PUSH_LOG", str(tmp_path / "pushes"))
    monkeypatch.delenv("FAIL_STEP", raising=False)
    monkeypatch.delenv("BUILD_HOOK", raising=False)
    return repo


def deploy(repo, *args):
    # Running from outside the repository must still release the script's repo.
    return subprocess.run(
        ["bash", str(repo / "deploy.sh"), *args], cwd=repo.parent,
        capture_output=True, text=True)


def steps(repo):
    log = repo.parent / "commands"
    return log.read_text().splitlines() if log.exists() else []


@pytest.mark.parametrize("branch", ["main", "master"])
@pytest.mark.parametrize("args", [[], [VERSION], [TAG]])
def test_publish_tags_only_the_merged_commit(release, branch, args):
    if branch == "master":
        git(release, "branch", "-m", "master")
        git(release, "push", "-u", "origin", "master")
        (release.parent / "pushes").unlink()
    before = git(release, "rev-parse", "HEAD")
    result = deploy(release, *args)
    assert result.returncode == 0, result.stderr
    assert steps(release) == ["lint", "test", "build", "check", "upload"]
    assert git(release, "rev-parse", "HEAD", TAG).splitlines() == [before, before]
    assert git(release, "ls-remote", "origin", "refs/tags/" + TAG).startswith(before)
    assert (release.parent / "pushes").read_text().split()[2:] == ["refs/tags/" + TAG]
    assert not git(release, "status", "--porcelain")
    files = list((release / "dist").glob("*/*"))
    assert {path.name for path in files} == {
        "varcode-1.2.3-py3-none-any.whl", "varcode-1.2.3.tar.gz"}
    commands = [json.loads(line) for line in
                (release.parent / "arguments").read_text().splitlines()]
    assert commands[1][:4] == ["-m", "twine", "check", "--strict"]
    assert set(commands[2][5:]) == {str(path.relative_to(release)) for path in files}


def test_dry_run_builds_without_publishing_or_changing_refs(release):
    before = git(release, "show-ref")
    result = deploy(release, "--dry-run", VERSION)
    assert result.returncode == 0, result.stderr
    assert steps(release) == ["lint", "test", "build", "check"]
    assert git(release, "show-ref") == before
    assert not (release.parent / "pushes").exists()
    assert not git(release, "status", "--porcelain")


@pytest.mark.parametrize("args", [
    ["1.2.4"], ["--dry-run", "1.2.4"], ["--unknown"],
    [VERSION, "1.2.4"], [""], ["v"], ["01.2.3"],
])
def test_bad_arguments_stop_before_gates_and_never_bump(release, args):
    before = git(release, "rev-parse", "HEAD")
    result = deploy(release, *args)
    assert result.returncode != 0
    assert not steps(release)
    assert git(release, "rev-parse", "HEAD") == before
    assert not git(release, "status", "--porcelain")
    assert not (release.parent / "pushes").exists()


@pytest.mark.parametrize("state", ["topic", "detached", "dirty", "staged", "untracked"])
def test_invalid_checkout_stops_before_gates(release, state):
    if state == "topic":
        git(release, "checkout", "-b", "topic")
    elif state == "detached":
        git(release, "checkout", "--detach")
    elif state == "untracked":
        (release / "new-file").touch()
    else:
        (release / ".gitignore").write_text("# edited\n")
        if state == "staged":
            git(release, "add", ".gitignore")
    result = deploy(release)
    assert result.returncode != 0
    assert ("main or master" if state in ("topic", "detached") else "not clean") in result.stderr
    assert not steps(release)


@pytest.mark.parametrize("state", ["ahead", "behind", "unreachable", "missing"])
def test_requires_live_origin_branch(release, state):
    if state in ("ahead", "behind"):
        git(release, "commit", "--allow-empty", "-m", "Unreleased change")
        if state == "behind":
            git(release, "push", "origin", "main")
            git(release, "reset", "--hard", "HEAD~")
            # Even a stale remote-tracking ref matching HEAD must not pass.
            git(release, "update-ref", "refs/remotes/origin/main", "HEAD")
    elif state == "unreachable":
        git(release, "remote", "set-url", "origin", str(release / "missing"))
    else:
        git(release.parent / "origin.git", "update-ref", "-d", "refs/heads/main")
    result = deploy(release)
    assert result.returncode != 0
    assert "origin/main" in result.stderr
    assert not steps(release)


@pytest.mark.parametrize("remote", [False, True])
def test_existing_tags_are_never_overwritten(release, remote):
    git(release, "tag", TAG)
    if remote:
        git(release, "push", "origin", "refs/tags/" + TAG)
        git(release, "tag", "-d", TAG)
    result = deploy(release)
    assert result.returncode != 0
    assert "already exists" in result.stderr
    assert not steps(release)


@pytest.mark.parametrize("step", ["lint", "test", "build", "check", "upload", "push"])
def test_failed_step_stops_release_and_preserves_artifacts(release, monkeypatch, step):
    monkeypatch.setenv("FAIL_STEP", step)
    result = deploy(release)
    assert result.returncode != 0
    pipeline = ["lint", "test", "build", "check", "upload", "push"]
    assert steps(release) == pipeline[:min(pipeline.index(step) + 1, 5)]
    assert bool(git(release, "tag", "--list", TAG)) == (step == "push")
    assert not git(release, "ls-remote", "origin", "refs/tags/" + TAG)
    if step in ("check", "upload", "push"):
        files = list((release / "dist").glob("*/*"))
        assert len(files) == 2
        assert all(path.read_bytes() == b"original bytes" for path in files)


@pytest.mark.parametrize("command", [
    "touch untracked",
    "git checkout -b topic",
    "git commit --allow-empty -m changed",
    "git tag v1.2.3",
    "git push origin HEAD:refs/tags/v1.2.3",
    'git --git-dir="$PWD/../origin.git" update-ref -d refs/heads/main',
])
def test_changes_during_build_block_upload(release, monkeypatch, command):
    hook = release.parent / "during-build.sh"
    hook.write_text(command + "\n")
    monkeypatch.setenv("BUILD_HOOK", str(hook))
    result = deploy(release)
    assert result.returncode != 0
    assert steps(release) == ["lint", "test", "build", "check"]


def test_each_attempt_preserves_earlier_artifacts(release):
    for _ in range(2):
        result = deploy(release, "--dry-run")
        assert result.returncode == 0, result.stderr
    assert len(list((release / "dist").iterdir())) == 2
    assert len(list((release / "dist").glob("*/*"))) == 4
