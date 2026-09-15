"""CLI imports must not depend on setuptools' removed pkg_resources API."""

from pathlib import Path
import subprocess
import sys
import textwrap

import pytest


@pytest.mark.parametrize("module", ["effects_script", "genes_script"])
def test_cli_help_without_pkg_resources(module):
    # Run in a fresh interpreter so previous test imports cannot hide the
    # dependency. Blocking the import also covers older test environments
    # where setuptools still happens to provide pkg_resources.
    code = textwrap.dedent("""
        import importlib
        import importlib.abc
        import sys

        class NoPkgResources(importlib.abc.MetaPathFinder):
            def find_spec(self, fullname, path=None, target=None):
                if fullname == "pkg_resources" or fullname.startswith("pkg_resources."):
                    raise ModuleNotFoundError("No module named 'pkg_resources'")

        sys.meta_path.insert(0, NoPkgResources())
        assert "pkg_resources" not in sys.modules
        module = importlib.import_module(sys.argv[1])
        module.main(["--help"])
    """)
    result = subprocess.run(
        [sys.executable, "-c", code, "varcode.cli." + module],
        cwd=Path(__file__).resolve().parents[1],
        capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "usage:" in result.stdout.lower()
    assert "--output-csv" in result.stdout
    assert "pkg_resources" not in result.stderr
