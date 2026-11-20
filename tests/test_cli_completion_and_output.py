"""CLI output streams, unified CSV columns, lazy startup, argcomplete, and help styling (#432, #697)."""

from __future__ import annotations

import io
import os
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest

from varcode.cli import effects_script, genes_script, main as cli_main
from varcode.cli.help import HelpArgumentParser, color_enabled, style_help, style_table_header

REPO_ROOT = Path(__file__).resolve().parents[1]
HEAVY_MODULES = ("pandas", "numpy", "pyarrow", "pyensembl", "orfdata", "pysam")


def _clean_env(**extra):
    env = os.environ.copy()
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = (
        f"{REPO_ROOT}{os.pathsep}{existing}" if existing else str(REPO_ROOT)
    )
    env.update(extra)
    return env


@pytest.mark.parametrize(
    "module_name",
    [
        "varcode",
        "varcode.cli.main",
        "varcode.cli.effects_script",
        "varcode.cli.genes_script",
        "varcode.cli.sample_checks",
    ],
)
def test_cli_and_package_import_avoid_heavy_dependencies(module_name):
    code = (
        "import importlib, sys\n"
        f"importlib.import_module({module_name!r})\n"
        f"loaded = [m for m in {HEAVY_MODULES!r} if m in sys.modules]\n"
        "assert not loaded, f'Heavy modules loaded on import: {loaded}'\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=str(REPO_ROOT),
        env=_clean_env(),
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize(
    "entry_module,comp_line,expected_items",
    [
        ("varcode.cli.main", "varcode ", {"effects", "genes", "check-samples"}),
        ("varcode.cli.main", "varcode effects --out", {"--output-csv", "--output"}),
        ("varcode.cli.main", "varcode genes --v", {"--vcf", "--variant", "--version"}),
        ("varcode.cli.main", "varcode check-samples --m", {"--manifest", "--mixed"}),
        ("varcode.cli.genes_script", "varcode-genes --out", {"--output-csv", "--output"}),
    ],
)
def test_argcomplete_returns_expected_completions(
    tmp_path, entry_module, comp_line, expected_items
):
    out_file = tmp_path / "argcomplete_out.txt"
    env = _clean_env(
        _ARGCOMPLETE="1",
        _ARGCOMPLETE_IFS="\n",
        _ARGCOMPLETE_STDOUT_FILENAME=str(out_file),
        COMP_LINE=comp_line,
        COMP_POINT=str(len(comp_line)),
    )
    result = subprocess.run(
        [sys.executable, "-c", f"import {entry_module} as m; m.main()"],
        cwd=str(REPO_ROOT),
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    completions = {
        line.strip()
        for line in out_file.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    assert expected_items <= completions, (
        f"Expected {expected_items} in completions {completions}"
    )


def test_entry_modules_have_argcomplete_marker_within_first_1024_bytes():
    for rel in (
        "varcode/cli/main.py",
        "varcode/cli/effects_script.py",
        "varcode/cli/genes_script.py",
        "varcode/cli/sample_checks.py",
    ):
        head = (REPO_ROOT / rel).read_bytes()[:1024]
        assert b"PYTHON_ARGCOMPLETE_OK" in head, f"Missing marker in {rel}"


def test_help_styling_plain_in_format_help_and_cyan_only_on_tty(monkeypatch):
    class _TTYBuffer(io.StringIO):
        def isatty(self):
            return True

    parser = HelpArgumentParser(prog="varcode-test", description="Test CLI")
    parser.add_argument("--vcf", help="Input ``VCF`` path")
    plain = parser.format_help()
    assert "\033[" not in plain
    assert "``" not in plain

    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")

    non_tty = io.StringIO()
    assert not color_enabled(non_tty)
    parser.print_help(file=non_tty)
    assert "\033[" not in non_tty.getvalue()

    tty = _TTYBuffer()
    assert color_enabled(tty)
    parser.print_help(file=tty)
    styled = tty.getvalue()
    assert "\033[1musage:\033[0m" in styled
    assert "\033[36m--vcf\033[0m" in styled
    assert "\033[1;36m" in styled
    assert style_table_header("contig  start\n7       100\n", tty).startswith("\033[1;36m")

    monkeypatch.setenv("NO_COLOR", "1")
    assert not color_enabled(tty)
    assert "\033[" not in style_help(plain, tty)

    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "dumb")
    assert not color_enabled(tty)


def test_help_and_version_do_not_print_legacy_banner(capsys):
    for mod in (effects_script, genes_script):
        with pytest.raises(SystemExit) as exc:
            mod.main(["--help"])
        assert exc.value.code == 0
        out = capsys.readouterr().out
        assert "Path:" not in out
        assert "usage:" in out.lower()

        with pytest.raises(SystemExit) as exc:
            mod.main(["--version"])
        assert exc.value.code == 0
        out = capsys.readouterr().out
        assert out.strip().startswith("varcode ")
        assert "Path:" not in out

    assert cli_main.main(["--help"]) == 0
    out = capsys.readouterr().out
    assert "effects" in out and "genes" in out and "check-samples" in out
    assert "Path:" not in out


def test_effects_and_genes_stdout_piped_csv_and_unified_columns(tmp_path, capsys):
    # Insertion using VCF-padded alleles (C -> CA) to verify trimmed normalization parity.
    vcf_path = tmp_path / "ins.vcf"
    vcf_path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "17\t7577548\t.\tC\tCA\t.\tPASS\t.\n",
        encoding="utf-8",
    )

    # Without --output-csv, piped (non-TTY) stdout emits clean parseable CSV with no logger lines.
    effects_script.main(["--genome", "GRCh37:75", "--vcf", str(vcf_path), "--one-per-variant"])
    eff_captured = capsys.readouterr()
    assert "INFO -" not in eff_captured.out
    assert "<EffectCollection" not in eff_captured.out
    eff_df = pd.read_csv(io.StringIO(eff_captured.out), dtype={"contig": str}, keep_default_na=False)

    genes_script.main(["--genome", "GRCh37:75", "--vcf", str(vcf_path)])
    gene_captured = capsys.readouterr()
    assert "INFO -" not in gene_captured.out
    gene_df = pd.read_csv(io.StringIO(gene_captured.out), dtype={"contig": str}, keep_default_na=False)

    # Shared variant columns & identical trimmed allele normalization.
    for col in ("contig", "start", "ref", "alt", "reference_name", "annotation_name", "annotation_version"):
        assert col in eff_df.columns
        assert col in gene_df.columns
    assert "chr" not in gene_df.columns
    assert tuple(eff_df.iloc[0][["contig", "start", "ref", "alt"]]) == (
        tuple(gene_df.iloc[0][["contig", "start", "ref", "alt"]])
    )
    assert tuple(gene_df.iloc[0][["contig", "start", "ref", "alt"]]) == ("17", 7577548, "", "A")
    assert "TP53" in gene_df.iloc[0]["gene_name"]
