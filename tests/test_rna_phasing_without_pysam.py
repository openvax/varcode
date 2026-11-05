"""Without the optional pysam (varcode[rna]), Varcode still imports and RNA phasing explains why it can't run."""

import importlib.abc
from pathlib import Path
import subprocess
import sys

import pytest

from varcode import RNAReadPhasingSource


def test_missing_pysam_gives_the_install_command_and_platforms(monkeypatch, tmp_path):
    monkeypatch.setitem(sys.modules, "pysam", None)  # makes `import pysam` raise ModuleNotFoundError
    with pytest.raises(ImportError, match=r"pip install 'varcode\[rna\]'\. pysam does not support Windows"):
        RNAReadPhasingSource(str(tmp_path / "reads.bam"))


def test_a_broken_pysam_install_reports_its_own_error(monkeypatch, tmp_path):
    class BrokenPysam(importlib.abc.MetaPathFinder):
        def find_spec(self, name, path=None, target=None):
            if name == "pysam":
                raise ImportError("libhts.so.3: cannot open shared object file", name="pysam")
            return None

    monkeypatch.delitem(sys.modules, "pysam", raising=False)
    monkeypatch.setattr(sys, "meta_path", [BrokenPysam(), *sys.meta_path])
    with pytest.raises(ImportError, match="libhts") as caught:
        RNAReadPhasingSource(str(tmp_path / "reads.bam"))
    assert "pip install" not in str(caught.value)


def test_varcode_and_every_module_import_without_pysam():
    # A fresh interpreter, so pysam is blocked before anything in Varcode is imported.
    code = (
        "import importlib, pkgutil, sys\n"
        "sys.modules['pysam'] = None\n"
        "import varcode\n"
        "names = [m.name for m in pkgutil.walk_packages(varcode.__path__, 'varcode.')]\n"
        "for name in names:\n"
        "    importlib.import_module(name)\n"
        "print(len(names))\n")
    result = subprocess.run([sys.executable, "-c", code], cwd=Path(__file__).resolve().parents[1],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert int(result.stdout) > 50
