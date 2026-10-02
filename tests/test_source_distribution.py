"""The source archive must carry the inputs needed to reproduce its tests."""

import hashlib
from pathlib import Path
import subprocess
import sys
import tarfile


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests/data"


def test_bundled_fixture_checksums():
    lines = (DATA / "SHA256SUMS").read_text().splitlines()
    expected = {name: digest for digest, name in
                (line.split("  ", 1) for line in lines)}
    assert len(expected) == len(lines), "Duplicate fixture checksum entries"
    actual = {path.relative_to(DATA).as_posix(): path for path in DATA.rglob("*")
              if path.is_file() and path.name != "SHA256SUMS"}
    assert actual.keys() == expected.keys(), "Update the fixture inventory after reviewing changed data"
    for name, path in actual.items():
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected[name], name


def test_source_distribution_contains_test_support(tmp_path):
    # Run the configured backend without downloads. This also works from an
    # extracted sdist, which has no Git metadata or surrounding checkout.
    result = subprocess.run(
        [sys.executable, "-c",
         "import sys; from setuptools.build_meta import build_sdist; build_sdist(sys.argv[1])",
         str(tmp_path)], cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    archives = list(tmp_path.glob("*.tar.gz"))
    assert len(archives) == 1

    required = {ROOT / name for name in (
        "README.md", "LICENSE", "pyproject.toml", "requirements.txt", "MANIFEST.in",
        "AGENTS.md", "CHANGELOG.md", "CONTRIBUTING.md", "RELEASING.md",
        "develop.sh", "lint.sh", "test.sh", "deploy.sh", "mkdocs.yml", ".coveragerc",
    )}
    for pattern in ("tests/**/*.py", "tests/**/*.md", "docs/**/*.md", "examples/**/*.py"):
        required.update(ROOT.glob(pattern))
    required.update(path for path in DATA.rglob("*") if path.is_file())

    with tarfile.open(archives[0]) as archive:
        members = {Path(member.name).relative_to(archives[0].name[:-7]).as_posix(): member
                   for member in archive.getmembers() if member.isfile()}
        missing = {path.relative_to(ROOT).as_posix() for path in required} - members.keys()
        assert not missing, "Missing source-distribution test support: %s" % sorted(missing)
        for path in required:
            name = path.relative_to(ROOT).as_posix()
            assert archive.extractfile(members[name]).read() == path.read_bytes(), name
        assert not any("__pycache__" in Path(name).parts or name.endswith(".pyc")
                       for name in members)
